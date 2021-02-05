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
        """Please install the graph-tool library either visiting

        https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

        or by conda: `conda install -c conda-forge graph-tool`
        """
    )

from ._utils import get_cell_loglikelihood

def nested_model(
    adata: AnnData,
#    max_iterations: int = 1000000,
#    epsilon: float = 0,
#    equilibrate: bool = False,
#    wait: int = 1000,
#    nbreaks: int = 2,
#    collect_marginals: bool = False,
#    niter_collect: int = 10000,
    hierarchy_length: int = 10,
    deg_corr: bool = True,
#    multiflip: bool = True,
    fast_model: bool = False,
    fast_tol: float = 1e-6,
    n_sweep: int = 10,
    beta: float = np.inf,
    n_init: int = 1,
#    beta_range: Tuple[float] = (1., 1000.),
#    steps_anneal: int = 3,
    resume: bool = False,
    *,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    random_seed: Optional[int] = None,
    key_added: str = 'nsbm',
    adjacency: Optional[sparse.spmatrix] = None,
    neighbors_key: Optional[str] = 'neighbors',
    directed: bool = False,
    use_weights: bool = False,
    prune: bool = False,
    return_low: bool = False,
    calculate_affinity: bool = True,
    copy: bool = False,
    minimize_args: Optional[Dict] = {},
#    equilibrate_args: Optional[Dict] = {},
) -> Optional[AnnData]:
    """\
    Cluster cells into subgroups [Peixoto14]_.

    Cluster cells using the nested Stochastic Block Model [Peixoto14]_,
    a hierarchical version of Stochastic Block Model [Holland83]_, performing
    Bayesian inference on node groups. NSBM should circumvent classical
    limitations of SBM in detecting small groups in large graphs
    replacing the noninformative priors used by a hierarchy of priors
    and hyperpriors.

    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    hierarchy_length
        Initial length of the hierarchy. When large values are
        passed, the top-most levels will be uninformative as they
        will likely contain the very same groups. Increase this valus
        if a very large number of cells is analyzed (>100.000).
    deg_corr
        Whether to use degree correction in the minimization step. In many
        real world networks this is the case, although this doesn't seem
        the case for KNN graphs used in scanpy.
    fast_model
        Whether to skip initial minization step and let the MCMC find a solution. 
        This approach tend to be faster and consume less memory, but may be
        less accurate.
    fast_tol
        Tolerance for fast model convergence.
    n_sweep 
        Number of iterations to be performed in the fast model MCMC greedy approach
    beta
        Inverse temperature for MCMC greedy approach    
    n_init
        Number of initial minimizations to be performed. The one with smaller
        entropy is chosen
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
    prune
        Some high levels in hierarchy may contain the same information in terms of 
        cell assignments, even if they apparently have different group names. When this
        option is set to `True`, the function only returns informative levels.
        Note, however, that cell affinities are still reported for all levels. Pruning
        does not rename group levels
    return_low
        Whether or not return nsbm_level_0 in adata.obs. This level usually contains
        so many groups that it cannot be plot anyway, but it may be useful for particular
        analysis. By default it is not returned
    calculate_affinity
        Whether or not calculate cell affinities. For large datasets it may save some time
    copy
        Whether to copy `adata` or modify it inplace.
    random_seed
        Random number to be used as seed for graph-tool

    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of samples) that stores the subgroup id
        (`'0'`, `'1'`, ...) for each cell. 
    `adata.uns['schist']['params']`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    `adata.uns['schist']['stats']`
        A dict with the values returned by mcmc_sweep
    `adata.obsm['CA_nsbm_level_{n}']`
        A `np.ndarray` with cell probability of belonging to a specific group
    `adata.uns['schist']['state']`
        The NestedBlockModel state object
    """

#    if resume: 
#        # if the fast_model is chosen perform equilibration anyway
#        # also if a model has previously created
#        equilibrate=True
        
    if resume and ('schist' not in adata.uns or 'state' not in adata.uns['schist']):
        # let the model proceed as default
        logg.warning('Resuming has been specified but a state was not found\n'
                     'Will continue with default minimization step')

        resume=False
        
    if random_seed:
        np.random.seed(random_seed)
        gt.seed_rng(random_seed)

#    if collect_marginals:
#        logg.warning('Collecting marginals has a large impact on running time')
#        if not equilibrate:
#            raise ValueError(
#                "You can't collect marginals without MCMC equilibrate "
#                "step. Either set `equlibrate` to `True` or "
#                "`collect_marginals` to `False`"
#            )

    start = logg.info('minimizing the nested Stochastic Block Model')
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

    if n_init < 1:
        n_init = 1

    if fast_model:
        # do not minimize, start with a dummy state and perform only a sweep
        # also, use consensus partitio from all models

        states = [gt.NestedBlockState(g=g,
                                    state_args=dict(deg_corr=deg_corr,
                                    recs=recs,
                                    rec_types=rec_types
                                    )) for n in range(n_init)]
        for x in range(n_init):
            dS = 1
            while np.abs(dS) > fast_tol:
                # perform sweep until a tolerance is reached
                dS, _, _ = states[x].multiflip_mcmc_sweep(beta=beta, niter=n_sweep, c=0.5)

        logg.info(f'Getting consensus over {len(bs)} models', time=start)

        state = states[np.argmin([s.entropy() for s in states])]
        bs = state.get_bs()
        
        logg.info('    done', time=start)
        
    elif resume:
        # create the state and make sure sampling is performed
        state = adata.uns['schist']['state'].copy(sampling=True)
        bs = state.get_bs()
        # get the graph from state
        g = state.g
    else:
        states = [gt.minimize_nested_blockmodel_dl(g, deg_corr=deg_corr, 
                  state_args=dict(recs=recs,  rec_types=rec_types), 
                  **minimize_args) for n in range(n_init)]
                  
        state = states[np.argmin([s.entropy() for s in states])]    
#        state = gt.minimize_nested_blockmodel_dl(g, deg_corr=deg_corr, 
#                                                 state_args=dict(recs=recs,
#                                                 rec_types=rec_types), 
#                                                 **minimize_args)
        logg.info('    done', time=start)
        bs = state.get_bs()
        if len(bs) <= hierarchy_length:
            # increase hierarchy length up to the specified value
            # according to Tiago Peixoto 10 is reasonably large as number of
            # groups decays exponentially
            bs += [np.zeros(1)] * (hierarchy_length - len(bs))
        else:
            logg.warning(f'A hierarchy length of {hierarchy_length} has been specified\n'
                         f'but the minimized model contains {len(bs)} levels')
            pass    
        # create a new state with inferred blocks   
        state = gt.NestedBlockState(g, bs, state_args=dict(recs=recs,
                                    rec_types=rec_types), sampling=True)
    
    # equilibrate the Markov chain
#    if equilibrate:
#        logg.info('running MCMC equilibration step')
#        # equlibration done by simulated annealing
#        
#        equilibrate_args['wait'] = wait
#        equilibrate_args['nbreaks'] = nbreaks
#        equilibrate_args['max_niter'] = max_iterations
#        equilibrate_args['multiflip'] = multiflip
#        equilibrate_args['mcmc_args'] = {'niter':10}
#        
#        dS, nattempts, nmoves = gt.mcmc_anneal(state, 
#                                               mcmc_equilibrate_args=equilibrate_args,
#                                               niter=steps_anneal,
#                                               beta_range=beta_range)
#    if collect_marginals and equilibrate:
#        # we here only retain level_0 counts, until I can't figure out
#        # how to propagate correctly counts to higher levels
#        # I wonder if this should be placed after group definition or not
#        logg.info('    collecting marginals')
#        group_marginals = [np.zeros(g.num_vertices() + 1) for s in state.get_levels()]
#        def _collect_marginals(s):
#            levels = s.get_levels()
#            for l, sl in enumerate(levels):
#                group_marginals[l][sl.get_nonempty_B()] += 1
#
#        gt.mcmc_equilibrate(state, wait=wait, nbreaks=nbreaks, epsilon=epsilon,
#                            max_niter=max_iterations, multiflip=True,
#                            force_niter=niter_collect, mcmc_args=dict(niter=10),
#                            callback=_collect_marginals)
#        logg.info('    done', time=start)

    # everything is in place, we need to fill all slots
    # first build an array with
    groups = np.zeros((g.num_vertices(), len(bs)), dtype=int)

    for x in range(len(bs)):
        # for each level, project labels to the vertex level
        # so that every cell has a name. Note that at this level
        # the labels are not necessarily consecutive
        groups[:, x] = state.project_partition(x, 0).get_array()

    groups = pd.DataFrame(groups).astype('category')

    # rename categories from 0 to n
    for c in groups.columns:
        ncat = len(groups[c].cat.categories)
        new_cat = [u'%s' % x for x in range(ncat)]
        groups[c].cat.rename_categories(new_cat, inplace=True)

    levels = groups.columns
    
    if restrict_to is not None:
        groups.index = adata.obs[restrict_key].index
    else:
        groups.index = adata.obs_names

    # add column names
    groups.columns = [f"{key_added}_level_{level}" for level in range(len(bs))]

    # remove any column with the same key
    keep_columns = [x for x in adata.obs.columns if not x.startswith('%s_level_' % key_added)]
    adata.obs = adata.obs[keep_columns]
    # concatenate obs with new data, skipping level_0 which is usually
    # crap. In the future it may be useful to reintegrate it
    # we need it in this function anyway, to match groups with node marginals
    if return_low:
        adata.obs = pd.concat([adata.obs, groups], axis=1)
    else:
        adata.obs = pd.concat([adata.obs, groups.iloc[:, 1:]], axis=1)

    # add some unstructured info

    adata.uns['schist'] = {}
    adata.uns['schist']['stats'] = dict(
    level_entropy=np.array([state.level_entropy(x) for x in range(len(state.levels))]),
    modularity=np.array([gt.modularity(g, state.project_partition(x, 0))
                         for x in range(len((state.levels)))])
    )
#    if equilibrate:
#        adata.uns['schist']['stats']['dS'] = dS
#        adata.uns['schist']['stats']['nattempts'] = nattempts
#        adata.uns['schist']['stats']['nmoves'] = nmoves


    adata.uns['schist']['state'] = state

    # now add marginal probabilities.

#    if collect_marginals:
#        # refrain group marginals. We collected data in vector as long as
#        # the number of cells, cut them into appropriate length data
#        adata.uns['schist']['group_marginals'] = {}
#        for nl, level_marginals in enumerate(group_marginals):
#            idx = np.where(level_marginals > 0)[0] + 1
#            adata.uns['schist']['group_marginals'][nl] = np.array(level_marginals[:np.max(idx)])

    # prune uninformative levels, if any
    if prune:
        to_remove = prune_groups(groups)
        logg.info(
            f'    Removing levels f{to_remove}'
        )
        adata.obs.drop(to_remove, axis='columns', inplace=True)
        
    # calculate log-likelihood of cell moves over the remaining levels
    # we have to calculate events at level 0 and propagate to upper levels
    # Should I use virtual vertex moves for all levels, instead?
    if calculate_affinity:
        logg.info('    calculating cell affinity to groups')
#        levels = [int(x.split('_')[-1]) for x in adata.obs.columns if x.startswith(f'{key_added}_level')]    
        p0 = get_cell_loglikelihood(state, level=0, as_prob=True)
    
        adata.obsm[f'CA_{key_added}_level_0'] = p0
        l0 = "%s_level_0" % key_added
        for nl, level in enumerate(groups.columns[1:]):
            cross_tab = pd.crosstab(groups[l0], groups[level])
            cl = np.zeros((p0.shape[0], cross_tab.shape[1]), dtype=p0.dtype)
            for x in range(cl.shape[1]):
                # sum counts of level_0 groups corresponding to
                # this group at current level
                cl[:, x] = p0[:, np.where(cross_tab.iloc[:, x] > 0)[0]].sum(axis=1)
            adata.obsm[f'CA_{key_added}_level_{nl + 1}'] = cl / np.sum(cl, axis=1)[:, None]
    
    # last step is recording some parameters used in this analysis
    adata.uns['schist']['params'] = dict(
        model='nested',
#        epsilon=epsilon,
#        wait=wait,
#        nbreaks=nbreaks,
#        equilibrate=equilibrate,
        fast_model=fast_model,
#        collect_marginals=collect_marginals,
        hierarchy_length=hierarchy_length,
        random_seed=random_seed,
        calculate_affinity=calculate_affinity,
        prune=prune,
    )


    logg.info(
        '    finished',
        time=start,
        deep=(
            f'found {state.get_levels()[1].get_nonempty_B()} clusters at level_1, and added\n'
            f'    {key_added!r}, the cluster labels (adata.obs, categorical)'
        ),
    )
    return adata if copy else None

