from typing import Optional, Tuple, Sequence, List, Type, Union, Dict

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from joblib import delayed, Parallel

from scanpy import logging as logg
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency
from scanpy._utils import get_igraph_from_adjacency


try:
    import graph_tool.all as gt
except ImportError:
    raise ImportError(
        """Please install the graph-tool library either visiting

        https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

        or by conda: `conda install -c conda-forge graph-tool`
        """
    )

def nested_model_multi(
    adatas: List[AnnData],
    deg_corr: bool = True,
    tolerance: float = 1e-6,
    n_sweep: int = 10,
    beta: float = np.inf,
    n_init: int = 10,
    collect_marginals: bool = True,
    n_jobs: int = -1,
    refine_model: bool = False,
    refine_iter: int = 1000,
    *,
    random_seed: Optional[int] = None,
    key_added: str = 'multi_nsbm',
    adjacency: Optional[List[sparse.spmatrix]] = None,
    neighbors_key: Optional[List[str]] = ['neighbors'],
    directed: bool = False,
    use_weights: bool = False,
    save_model: Union[str, None] = None,
    copy: bool = False,
#    minimize_args: Optional[Dict] = {},
    dispatch_backend: Optional[str] = 'processes',
#    equilibrate_args: Optional[Dict] = {},
) -> Optional[List[AnnData]]:
    """\
    Cluster cells into subgroups using multiple modalities.

    Cluster cells using the nested Stochastic Block Model [Peixoto14]_,
    performing Bayesian inference on node groups. This function takes multiple
    experiments, possibly across different modalities, and perform joint
    clustering.
    

    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first. It also requires cells having the same
    names if coming from paired experiments

    Parameters
    ----------
    adatas
        A list of processed AnnData. Neighbors must have been already
        calculated.
    deg_corr
        Whether to use degree correction in the minimization step. In many
        real world networks this is the case, although this doesn't seem
        the case for KNN graphs used in scanpy.
    tolerance
        Tolerance for fast model convergence.
    n_sweep 
        Number of iterations to be performed in the fast model MCMC greedy approach
    beta
        Inverse temperature for MCMC greedy approach    
    n_init
        Number of initial minimizations to be performed. The one with smaller
        entropy is chosen
    refine_model
    	Wether to perform a further mcmc step to refine the model
    refine_iter
    	Number of refinement iterations.
    n_jobs
        Number of parallel computations used during model initialization
    key_added
        `adata.obs` key under which to add the cluster labels.
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']` in case of scanpy<=1.4.6 or
        `adata.obsp[neighbors_key][connectivity_key]` for scanpy>1.4.6
    neighbors_key
        The key passed to `sc.pp.neighbors`. If all AnnData share the same key, one
        only has to be specified, otherwise the full tuple of all keys must 
        be provided
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges). Note that this
        increases computation times
    save_model
        If provided, this will be the filename for the PartitionModeState to 
        be saved    
    copy
        Whether to copy `adata` or modify it inplace.
    random_seed
        Random number to be used as seed for graph-tool

    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of cells) that stores the subgroup id
        (`'0'`, `'1'`, ...) for each cell. 
    `adata.uns['schist']['multi_level_params']`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    `adata.uns['schist']['multi_level_stats']`
        A dict with the values returned by mcmc_sweep
    `adata.obsm['CA_multi_nsbm_level_{n}']`
        A `np.ndarray` with cell probability of belonging to a specific group
    `adata.uns['schist']['multi_level_state']`
        The NestedBlockModel state object
    """

    if random_seed:
        np.random.seed(random_seed)
    
    seeds = np.random.choice(range(n_init**2), size=n_init, replace=False)
        

    if collect_marginals and n_init < 100:
        logg.warning('Collecting marginals requires sufficient number of n_init\n'
                     f'It is now set to {n_init} and should be at least 100')
        

    start = logg.info('minimizing the nested Stochastic Block Model')
    
    if copy:
        adatas = [x.copy() for x in adatas]

    n_keys = len(neighbors_key)
    n_data = len(adatas)
    # are we clustering a user-provided graph or the default AnnData one?
    if adjacency is None:
        adjacency = []
        if n_keys > 1 and n_keys < n_data:
            raise ValueError(
                'The number of neighbors keys does not match'
                'the number of data matrices. Either fix this'
                'or pass a neighbor key that is shared across all modalities'
            )
        if n_keys == 1:
            neighbors_key = [neighbors_key[0] for x in range(n_data)]    
        for x in range(n_data):
            logg.info(f'getting adjacency for data {x}', time=start)
            if neighbors_key[x] not in adatas[x].uns:
                raise ValueError(
                    'You need to run `pp.neighbors` first '
                    'to compute a neighborhood graph. for'
                    f'data entry {x}'
                )
            elif 'connectivities_key' in adatas[x].uns[neighbors_key[x]]:
                # scanpy>1.4.6 has matrix in another slot
                conn_key = adatas[x].uns[neighbors_key[x]]['connectivities_key']
                adjacency.append(adatas[x].obsp[conn_key])
            else:
                # scanpy<=1.4.6 has sparse matrix here
                adjacency.append(adatas[x].uns[neighbors_key[x]]['connectivities'])


    # convert it to igraph and graph-tool
    
    graph_list = []
    for x in range(n_data):
        g = get_igraph_from_adjacency(adjacency[x], directed=directed)
        g = g.to_graph_tool()
        gt.remove_parallel_edges(g)
        # add cell names to graph, this will be used to create
        # layered graph 
        g_names = g.new_vertex_property('string') 
        d_names = adatas[x].obs_names
        for xn in range(len(d_names)):
            g_names[xn] = d_names[xn]
        g.vp['cell'] = g_names
        graph_list.append(g)
    
# skip weights for now    
#    recs=[]
#    rec_types=[]
#    if use_weights:
        # this is not ideal to me, possibly we may need to transform
        # weights. More tests needed.
#        recs=[g.ep.weight]
#        rec_types=['real-normal']
    
    # get a non-redundant list of all cell names across all modalities
    all_names = set(adatas[0].obs_names)
    [all_names.update(adatas[x].obs_names) for x in range(1, n_data)]
    all_names = list(all_names)
    # create the shared graph
    union_g = gt.Graph(directed=False)
    union_g.add_vertex(len(all_names))
    u_names = union_g.new_vertex_property('string')
    for xn in range(len(all_names)):
        u_names[xn] = all_names[xn]
    union_g.vp['cell'] = u_names
    
    # now handle in a non elegant way the index mapping across all 
    # modalities and the unified Graph
    
    u_cell_index = dict([(union_g.vp['cell'][x], x) for x in range(union_g.num_vertices())])
    # now create layers
    layer = union_g.new_edge_property('int')
    for ng in range(n_data):
        for e in graph_list[ng].edges():
            S, T = e.source(), e.target()
            Sn = graph_list[ng].vp['cell'][S]
            Tn = graph_list[ng].vp['cell'][T]
            Sidx = u_cell_index[Sn]
            Tidx = u_cell_index[Tn]
            ne = union_g.add_edge(Sidx, Tidx)
            layer[ne] = ng + 1 # this is the layer label

    union_g.ep['layer'] = layer
    # DONE! now proceed with standard minimization, ish
    
    if n_init < 1:
        n_init = 1

    states = [gt.NestedBlockState(g=union_g,
                                  base_type=gt.LayeredBlockState,
                                  state_args=dict(deg_corr=deg_corr,
                                  ec=union_g.ep.layer,
                                  layers=True
                                  )) for n in range(n_init)]

    def fast_min(state, beta, n_sweep, fast_tol, seed=None):
        if seed:
            gt.seed_rng(seed)
        dS = 1
        while np.abs(dS) > fast_tol:
            dS, _, _ = state.multiflip_mcmc_sweep(beta=beta, niter=n_sweep, c=0.5)
        return state                            
            
    states = Parallel(n_jobs=n_jobs, prefer=dispatch_backend)(
        delayed(fast_min)(states[x], beta, n_sweep, tolerance, seeds[x]) for x in range(n_init)
    )
    logg.info('        minimization step done', time=start)
    pmode = gt.PartitionModeState([x.get_bs() for x in states], converge=True, nested=True)
    bs = pmode.get_max_nested()
    logg.info('        consensus step done', time=start)
        
    # prune redundant levels at the top
    bs = [x for x in bs if len(np.unique(x)) > 1]
    bs.append(np.array([0], dtype=np.int32)) #in case of type changes, check this
    state = gt.NestedBlockState(union_g, bs=bs,
                                  base_type=gt.LayeredBlockState,
                                  state_args=dict(deg_corr=deg_corr,
                                  ec=union_g.ep.layer,
                                  layers=True
                                  ))
    
    if refine_model:
        # we here reuse pmode variable, so that it is consistent
        logg.info('        Refining model')
        bs = []
        def collect_partitions(s):
            bs.append(s.get_bs())
        gt.mcmc_equilibrate(state, force_niter=refine_iter, 
                            multiflip=True, 
                            mcmc_args=dict(niter=n_sweep, beta=beta),
                            callback=collect_partitions)
        pmode = gt.PartitionModeState(bs, nested=True, converge=True)
        bs = [x for x in pmode.get_max_nested() if len(np.unique(x)) > 1]
        bs.append(np.array([0], dtype=np.int32)) #in case of type changes, check this
        state = gt.NestedBlockState(union_g, bs=bs,
                                  base_type=gt.LayeredBlockState,
                                  state_args=dict(deg_corr=deg_corr,
                                  ec=union_g.ep.layer,
                                  layers=True
                                  ))
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
    u_groups = np.unique(bs[0])
    n_groups = len(u_groups)
    last_group = np.max(u_groups) + 1

    if collect_marginals:
        # note that the size of this will be equal to the number of the groups in Mode
        # but some entries won't sum to 1 as in the collection there may be differently
        # sized partitions
        pv_array = pmode.get_marginal(union_g).get_2d_array(range(last_group)).T[:, u_groups] / n_init	
         
    groups = np.zeros((union_g.num_vertices(), len(bs)), dtype=int)

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
        groups[c] = groups[c].cat.rename_categories(new_cat)

    levels = groups.columns
    
    # recode block names to have consistency with group names
    i_groups = groups.astype(int)
    bs = [i_groups.iloc[:, 0].values]
    for x in range(1, groups.shape[1]):
        bs.append(np.where(pd.crosstab(i_groups.iloc[:, x - 1], i_groups.iloc[:, x])> 0)[1])
    state = gt.NestedBlockState(union_g, bs)
    del(i_groups)

    groups.index = all_names

    # add column names
    groups.columns = [f"{key_added}_level_{level}" for level in range(len(bs))]

    # remove any column with the same key
    for xn in range(n_data):
        drop_columns = groups.columns.intersection(adatas[xn].obs.columns)
        adatas[xn].obs.drop(drop_columns, 'columns', inplace=True)
        adatas[xn].obs = pd.concat([adatas[xn].obs, groups.loc[adatas[xn].obs_names]], axis=1)

        # now add marginal probabilities.

        if collect_marginals:
            # add marginals for level 0, the sum up according to the hierarchy
            _groups = groups.loc[adatas[xn].obs_names]
            _pv_array = pd.DataFrame(pv_array, index=all_names).loc[adatas[xn].obs_names].values
            adatas[xn].obsm[f"CM_{key_added}_level_0"] = _pv_array
            for group in groups.columns[1:]:
                ct = pd.crosstab(_groups[_groups.columns[0]], _groups[group], 
                                 normalize='index', dropna=False)
                adatas[xn].obsm[f'CM_{group}'] = _pv_array @ ct.values

        # add some unstructured info
        if not 'schist' in adatas[xn].uns:
            adatas[xn].uns['schist'] = {}

        adatas[xn].uns['schist'][f'{key_added}'] = {}
        adatas[xn].uns['schist'][f'{key_added}']['stats'] = dict(
        level_entropy=np.array([state.level_entropy(x) for x in range(len(state.levels))]),
        modularity=np.array([gt.modularity(union_g, state.project_partition(x, 0))
                             for x in range(len((state.levels)))])
        )

        bl_d = {}
        levels = state.get_levels()
        for nl in range(len(levels)):
            bl_d[str(nl)] = np.array(levels[nl].get_blocks().a)
        adatas[xn].uns['schist'][f'{key_added}']['blocks'] = bl_d

        # last step is recording some parameters used in this analysis
        adatas[xn].uns['schist'][f'{key_added}']['params'] = dict(
            model='multiome_nested',
            use_weights=use_weights,
            neighbors_key=neighbors_key[xn],
            key_added=key_added,
            n_init=n_init,
            collect_marginals=collect_marginals,
            random_seed=random_seed,
            deg_corr=deg_corr,
            refine_model=refine_model,
            refine_iter=refine_iter
#            recs=recs,
#            rec_types=rec_types
        )


    logg.info(
        '    finished',
        time=start,
        deep=(
            f'and added\n'
            f'    {key_added!r}, the cluster labels (adata.obs, categorical)'
        ),
    )
    return adatas if copy else None

    