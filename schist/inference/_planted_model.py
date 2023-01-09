from typing import Optional, Tuple, Sequence, Type, Union, Dict

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from joblib import delayed, Parallel
from natsort import natsorted
from scanpy import logging as logg
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency
from scanpy._utils import get_igraph_from_adjacency

import graph_tool.all as gt


def planted_model(
    adata: AnnData,
    n_sweep: int = 10,
    beta: float = np.inf, 
    tolerance = 1e-6,
    collect_marginals: bool = True,
    deg_corr: bool = True,
    n_init: int = 10,
    n_jobs: int = -1,
    refine_model: bool = True,
    refine_iter: int = 100,
    *,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    random_seed: Optional[int] = None,
    key_added: str = 'ppbm',
    adjacency: Optional[sparse.spmatrix] = None,
    neighbors_key: Optional[str] = 'neighbors',
    directed: bool = False,
    use_weights: bool = False,
    copy: bool = False,
    save_model: Union[str, None] = None,
#    minimize_args: Optional[Dict] = {},
    dispatch_backend: Optional[str] = 'processes',
) -> Optional[AnnData]:
    """\
    Cluster cells into subgroups [Peixoto14]_.

    Cluster cells using the  Planted Partition Block Model [Peixoto14]_, performing
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
    tolerance
        Difference in description length to stop MCMC sweep iterations        
    collect_marginals
        Whether or not collect node probability of belonging
        to a specific partition.
    deg_corr
        Whether to use degree correction in the minimization step. In many
        real world networks this is the case, although this doesn't seem
        the case for KNN graphs used in scanpy.
    n_init
        Number of initial minimizations to be performed. This influences also the 
        precision for marginals
    refine_model
    	Wether to perform a further mcmc step to refine the model
    refine_iter
    	Number of refinement iterations.
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
    save_model
        If provided, this will be the filename for the PartitionModeState to 
        be saved    
    random_seed
        Random number to be used as seed for graph-tool
    n_jobs
        Number of parallel computations used during model initialization

    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of cells) that stores the subgroup id
        (`'0'`, `'1'`, ...) for each cell.
    `adata.uns['schist']['params']`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    `adata.uns['schist']['stats']`
        A dict with the values returned by mcmc_sweep
    `adata.obsm['CM_ppbm']`
        A `np.ndarray` with cell probability of belonging to a specific group
    `adata.uns['schist']['state']`
        The BlockModel state object
    """

    if random_seed:
        np.random.seed(random_seed)
    
    seeds = np.random.choice(range(n_init**2), size=n_init, replace=False)

    if collect_marginals and not refine_model:
        if n_init < 100:
            logg.warning('Collecting marginals without refinement requires sufficient number of n_init\n'
                     f'It is now set to {n_init} and should be at least 100\n')
    elif refine_model and refine_iter < 100:                     
        logg.warning('Collecting marginals with refinement requires sufficient number of iterations\n'
                     f'It is now set to {refine_iter} and should be at least 100\n')

    start = logg.info('minimizing the Planted Partition Block Model')
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
    # convert it to igraph and graph-tool
    g = get_igraph_from_adjacency(adjacency, directed=directed)
    g = g.to_graph_tool()
    gt.remove_parallel_edges(g)

    recs=[]
    rec_types=[]
    if use_weights:
        # this is not ideal to me, possibly we may need to transform
        # weights. More tests needed.
        recs=[g.ep.weight]
        rec_types=['real-normal']

    if n_init < 1:
        n_init = 1
        
    # initialize  the block states
    def fast_min(state, beta, n_sweep, fast_tol, seed=None):
        if seed:
            gt.seed_rng(seed)
        dS = 1
        while np.abs(dS) > fast_tol:
            dS, _, _ = state.multiflip_mcmc_sweep(beta=beta, niter=n_sweep)
        return state

    states = [gt.PPBlockState(g) for x in range(n_init)]
        
    # perform a mcmc sweep on each 
    # no list comprehension as I need to collect stats
        
    states = Parallel(n_jobs=n_jobs, prefer=dispatch_backend)(
             delayed(fast_min)(states[x], beta, n_sweep, tolerance, seeds[x]) for x in range(n_init)
             )
    logg.info('        minimization step done', time=start)        
    pmode = gt.PartitionModeState([x.get_blocks().a for x in states], converge=True)
        
    bs = pmode.get_max(g)
    logg.info('        consensus step done', time=start)
    state = gt.PPBlockState(g, b=bs)

    if refine_model:
        # we here reuse pmode variable, so that it is consistent
        logg.info('        Refining model')
        bs = []
        def collect_partitions(s):
            bs.append(s.get_blocks().a)
        gt.mcmc_equilibrate(state, force_niter=refine_iter, 
                            multiflip=True, 
                            mcmc_args=dict(niter=n_sweep, beta=beta),
                            callback=collect_partitions)
        pmode = gt.PartitionModeState(bs, converge=True)
        bs = pmode.get_max(g)
        state = gt.PPBlockState(g, b=bs)
        logg.info('        refinement complete', time=start)
        
    if save_model:
        import pickle
        fname = save_model
        if not fname.endswith('pkl'):
            fname = f'{fname}.pkl'
        logg.info(f'Saving model into {fname}')    
        with open(fname, 'wb') as fout:
            pickle.dump(pmode, fout, 2)
    
    
    logg.info('    done', time=start)

    groups = np.array(bs.get_array())
    u_groups = np.unique(groups)
    n_groups = len(u_groups)
    last_group = np.max(u_groups) + 1
    if collect_marginals:
        pv_array = pmode.get_marginal(g).get_2d_array(range(last_group)).T[:, u_groups] / n_init

    rosetta = dict(zip(u_groups, range(len(u_groups))))
    groups = np.array([rosetta[x] for x in groups])

    if restrict_to is not None:
        if key_added == 'ppbm':
            key_added += '_R'
        groups = rename_groups(
            adata,
            key_added,
            restrict_key,
            restrict_categories,
            restrict_indices,
            groups,
        )

    # add column names
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(map(str, np.unique(groups))),
    )

    # now add marginal probabilities.

    if collect_marginals:
        # cell marginals will be a list of arrays with probabilities
        # of belonging to a specific group
        adata.obsm[f"CM_{key_added}"] = pv_array

    # add some unstructured info
    if not 'schist' in adata.uns:
        adata.uns['schist'] = {}

    adata.uns['schist'][f'{key_added}'] = {}
    adata.uns['schist'][f'{key_added}']['stats'] = dict(
    entropy=state.entropy(),
    modularity=gt.modularity(g, state.get_blocks())
    )

    # record state as list of blocks
    # for compatibility with nested model, use a dictionary with a single key here
    # although a np.array would be ok
    adata.uns['schist'][f'{key_added}']['blocks'] = {'0':np.array(state.get_blocks().a)}

    # last step is recording some parameters used in this analysis
    adata.uns['schist'][f'{key_added}']['params'] = dict(
        model='planted',
        use_weights=use_weights,
        neighbors_key=neighbors_key,
        key_added=key_added,
        n_init=n_init,
        refine_model=refine_model,
        refine_iter=refine_iter,
        collect_marginals=collect_marginals,
        random_seed=random_seed,
        deg_corr=deg_corr,
        recs=recs,
        rec_types=rec_types
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

