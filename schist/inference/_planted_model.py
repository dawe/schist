from typing import Optional, Tuple, Sequence, Type, Union, Dict

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from scanpy import logging as logg
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency

from .._utils import get_graph_tool_from_adjacency, prune_groups

try:
    import graph_tool.all as gt
except ImportError:
    raise ImportError(
        """Please install the graph-tool>=2.33 library either visiting

        https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

        or by conda: `conda install -c conda-forge graph-tool`
        """
    )

from ._utils import get_cell_loglikelihood

from ._helpers import *

def planted_model(
    adata: AnnData,
    n_sweep: int = 1000,
    beta: float = np.inf, 
    max_iterations: int = 1000000,
    epsilon: float = 0,
    equilibrate: bool = False,
    wait: int = 1000,
    nbreaks: int = 2,
    collect_marginals: bool = False,
    niter_collect: int = 10000,
    deg_corr: bool = True,
    n_init: int = 1,
    beta_range: Tuple[float] = (1., 100.),
    steps_anneal: int = 5,
    resume: bool = False,
    *,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    random_seed: Optional[int] = None,
    key_added: str = 'ppbm',
    adjacency: Optional[sparse.spmatrix] = None,
    neighbors_key: Optional[str] = 'neighbors',
    directed: bool = False,
    use_weights: bool = False,
    copy: bool = False,
    minimize_args: Optional[Dict] = {},
    equilibrate_args: Optional[Dict] = {},    
) -> Optional[AnnData]:
    """\
    Cluster cells into subgroups [Peixoto14]_.

    Cluster cells using the  Stochastic Block Model [Peixoto14]_, performing
    Bayesian inference on node groups. This function, in particular, uses
    the Planted Block Model, which is particularly suitable in case of
    assortative graphs and it returns the optimal number of communities

    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    n_sweep
        Number of MCMC sweeps to get the initial guess
    beta
        Inverse temperature for the initial MCMC sweep            
    max_iterations
        Maximal number of iterations to be performed by the equilibrate step.
    epsilon
        Relative changes in entropy smaller than epsilon will
        not be considered as record-breaking.
    equilibrate
        Whether or not perform the mcmc_equilibrate step.
        Equilibration should always be performed. Note, also, that without
        equilibration it won't be possible to collect marginals.
    collect_marginals
        Whether or not collect node probability of belonging
        to a specific partition.
    niter_collect
        Number of iterations to force when collecting marginals. This will
        increase the precision when calculating probabilites
    wait
        Number of iterations to wait for a record-breaking event.
        Higher values result in longer computations. Set it to small values
        when performing quick tests.
    nbreaks
        Number of iteration intervals (of size `wait`) without
        record-breaking events necessary to stop the algorithm.
    deg_corr
        Whether to use degree correction in the minimization step. In many
        real world networks this is the case, although this doesn't seem
        the case for KNN graphs used in scanpy.
    n_init
        Number of initial minimizations to be performed. The one with smaller
        entropy is chosen
    beta_range
        Inverse temperature at the beginning and the end of the equilibration
    steps_anneal
        Number of steps in which the simulated annealing is performed
    resume
        Start from a previously created model, if any, without initializing a novel
        model    
    key_added
        `adata.obs` key under which to add the cluster labels.
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']` in case of scanpy<=1.4.6 or
        `adata.obsp[neighbors_key][connectivity_key]` for scanpy>1.4.6
    neighbors_key
        The key passed to `sc.pp.neighbors`
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges). Note that this
        increases computation times
    copy
        Whether to copy `adata` or modify it inplace.
    random_seed
        Random number to be used as seed for graph-tool

    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of samples) that stores the subgroup id
        (`'0'`, `'1'`, ...) for each cell.
    `adata.uns['sbm']['params']`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    `adata.uns['sbm']['stats']`
        A dict with the values returned by mcmc_sweep
    `adata.uns['sbm']['cell_affinity']`
        A `np.ndarray` with cell probability of belonging to a specific group
    `adata.uns['sbm']['state']`
        The BlockModel state object
    """

    # first things first
    check_gt_version()
    
    if resume: 
        equilibrate=True
        
    if resume and (key_added not in adata.uns or 'state' not in adata.uns[key_added]):
        # let the model proceed as default
        logg.warning('Resuming has been specified but a state was not found\n'
                     'Will continue with default minimization step')

        resume=False

    if random_seed:
        np.random.seed(random_seed)
        gt.seed_rng(random_seed)

    if collect_marginals:
        logg.warning('Collecting marginals has a large impact on running time')
        if not equilibrate:
            raise ValueError(
                "You can't collect marginals without MCMC equilibrate "
                "step. Either set `equlibrate` to `True` or "
                "`collect_marginals` to `False`"
            )

    start = logg.info('minimizing the Stochastic Block Model')
    adata = adata.copy() if copy else adata
    # are we clustering a user-provided graph or the default AnnData one?
    if adjacency is None:
        if neighbors_key not in adata.uns:
            raise ValueError(
                'You need to run `pp.neighbors` first '
                'to compute a neighborhood graph.'
            )
        elif 'connectivities_key' in adata.uns[neighbors_key]:
            # scanpy>1.4.6 has matrix in another slot
            conn_key = adata.uns[neighbors_key]['connectivities_key']
            adjacency = adata.obsp[conn_key]
        else:
            # scanpy<=1.4.6 has sparse matrix here
            adjacency = adata.uns[neighbors_key]['connectivities']
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories,
            adjacency,
        )
    # convert it to igraph
    g = get_graph_tool_from_adjacency(adjacency, directed=directed)

    recs=[]
    rec_types=[]
    if use_weights:
        # this is not ideal to me, possibly we may need to transform
        # weights. More tests needed.
        recs=[g.ep.weight]
        rec_types=['real-normal']

    if resume:
        # create the state and make sure sampling is performed
        state = adata.uns[key_added]['state'].copy()
        g = state.g
    else:
        if n_init < 1:
            n_init = 1
        
        # initialize  the block states
        states = [gt.PPBlockState(g) for x in range(n_init)]
        
        # perform a mcmc sweep on each 
        # no list comprehension as I need to collect stats
        
        _dS = np.zeros(n_init)
        _nattempts = np.zeros(n_init)
        _nmoves = np.zeros(n_init)
        for x in range(n_init):
            _dS[x], _nattempts[x], _nmoves[x] = states[x].multiflip_mcmc_sweep(beta=beta, niter=n_sweep)

        _amin = np.argmin([s.entropy() for s in states])            
        state = states[_amin]
        dS = _dS[_amin]
        nattempts = _nattempts[_amin]
        nmoves = _nmoves[_amin]
        

        logg.info('    done', time=start)
    
    # equilibrate the Markov chain
    if equilibrate:
        logg.info('running MCMC equilibration step')
        equilibrate_args['wait'] = wait
        equilibrate_args['nbreaks'] = nbreaks
        equilibrate_args['max_niter'] = max_iterations
        equilibrate_args['mcmc_args'] = {'niter':10}
        
        dS, nattempts, nmoves = gt.mcmc_anneal(state, 
                                               mcmc_equilibrate_args=equilibrate_args,
                                               niter=steps_anneal,
                                               beta_range=beta_range)

    if collect_marginals and equilibrate:
        # we here only retain level_0 counts, until I can't figure out
        # how to propagate correctly counts to higher levels
        # I wonder if this should be placed after group definition or not
        logg.info('    collecting marginals')
        group_marginals = np.zeros(g.num_vertices() + 1)
        def _collect_marginals(s):
            group_marginals[s.get_B()] += 1

        gt.mcmc_equilibrate(state, wait=wait, nbreaks=nbreaks, epsilon=epsilon,
                            max_niter=max_iterations, multiflip=True,
                            force_niter=niter_collect, mcmc_args=dict(niter=10),
                            callback=_collect_marginals)
        logg.info('    done', time=start)

    # everything is in place, we need to fill all slots
    # first build an array with
    groups = pd.Series(state.get_blocks().get_array()).astype('category')
    new_cat_names = dict([(cx, u'%s' % cn) for cn, cx in enumerate(groups.cat.categories)])
    groups.cat.rename_categories(new_cat_names, inplace=True)

    if restrict_to is not None:
        groups.index = adata.obs[restrict_key].index
    else:
        groups.index = adata.obs_names

    # add column names
    adata.obs.loc[:, key_added] = groups

    # add some unstructured info

    adata.uns[key_added] = {}
    adata.uns[key_added]['stats'] = dict(
    dS=dS,
    nattempts=nattempts,
    nmoves=nmoves,
    modularity=gt.modularity(g, state.get_blocks())
    )
    adata.uns[key_added]['state'] = state

    # now add marginal probabilities.

    if collect_marginals:
        # cell marginals will be a list of arrays with probabilities
        # of belonging to a specific group
        adata.uns[key_added]['group_marginals'] = group_marginals

    # calculate log-likelihood of cell moves over the remaining levels
    
    # adata.uns[key_added]['cell_affinity'] = {'1':get_cell_loglikelihood(state, as_prob=True, rescale=True)}
    
    # last step is recording some parameters used in this analysis
    adata.uns[key_added]['params'] = dict(
        epsilon=epsilon,
        wait=wait,
        nbreaks=nbreaks,
        equilibrate=equilibrate,
        collect_marginals=collect_marginals,
        random_seed=random_seed
    )


    logg.info(
        '    finished',
        time=start,
        deep=(
            f'found {state.get_B()} clusters and added\n'
            f'    {key_added!r}, the cluster labels (adata.obs, categorical)'
        ),
    )
    return adata if copy else None

