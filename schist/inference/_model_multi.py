from typing import Optional, Tuple, Sequence, Type, Union, Dict, List, Literal

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData
from scipy import sparse
from natsort import natsorted
from joblib import delayed, Parallel, parallel_config
from tqdm import tqdm
from scanpy import logging as logg
from .._utils import get_graph_tool_from_adjacency


try:
    import graph_tool.all as gt
except ImportError:
    raise ImportError(
        """Please install the graph-tool library either visiting

        https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

        or by conda: `conda install -c conda-forge graph-tool`
        """
    )

def fast_min(state, beta=np.inf, n_sweep=10, fast_tol=1e-4, max_iter=1000, seed=None):
    if seed:
        gt.seed_rng(seed)
    dS = 1e9
    n = 0
    while (np.abs(dS) > fast_tol) and (n < max_iter):
        dS, _, _ = state.multiflip_mcmc_sweep(beta=beta, niter=n_sweep, c=0.5)
        n += 1
    return state                            


def fit_model_multi(
    mdata: Union[List[AnnData], MuData],
    deg_corr: bool = True,
    tolerance: float = 1e-4,
    n_sweep: int = 10,
    beta: float = np.inf,
    n_init: int = 100,
    model: Literal["nsbm", "sbm"] = "nsbm",
    max_iter: int = 1000,
    collect_marginals: bool = True,
    refine_model: bool = False,
    refine_iter: int = 100,
    n_jobs: int = -1,
    overlap: bool = False,
    key_added: str | None = None,
    adjacency: Optional[List[sparse.spmatrix]] = None,
    neighbors_key: Optional[List[str]] = ['neighbors'],
    directed: bool = False,
    use_weights: bool = False,
    save_model: Union[str, None] = None,
    copy: bool = False,
    dispatch_backend: Optional[str] = 'loky',
    random_seed: Optional[int] = None,
) -> [Union[List[AnnData]], MuData, None]:
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
    mdata
        A list of processed AnnData. Neighbors must have been already
        calculated. If a MuData object is passed, a model on the layered graph
        will be fitted. If you want to fit a model on the shared graph representation, 
        e.g. WNN graph or a graph built on MOFA latent factors, you still can use
        the standard ``scs.inference.model()`` function.
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
    max_iter
        Maximum number of iterations during minimization, set to infinite to stop 
        minimization only on tolerance
    overlap
        Whether the different layers are dependent (overlap=True) or not (overlap=False)
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

    gt_model = {'nsbm':gt.NestedBlockState, 
                'sbm':gt.LayeredBlockState,
    }[model]
        
    # if key is not set, use the model name
    key_added = f'multi_{model}' if key_added is None else key_added

    if random_seed:
        np.random.seed(random_seed)
    
    if n_init < 1:
        n_init = 1
    seeds = np.random.choice(range(n_init**2), size=n_init, replace=False)
        
    # the following lines are for compatibility
    if dispatch_backend == 'threads':
        dispatch_backend = 'threading'
    elif dispatch_backend == 'processes':
        dispatch_backend = 'loky'

    if dispatch_backend == 'threading' and float(gt.__version__.split()[0]) > 2.55:
        logg.warning('Threading backend does not work with this version of graph-tool\n'
                     'Switching to loky backend')
        dispatch_backend = 'loky'

    if collect_marginals and not refine_model:
        if n_init < 100:
            logg.warning('Collecting marginals without refinement requires sufficient number of n_init\n'
                     f'It is now set to {n_init} and should be at least 100\n')
    elif refine_model and refine_iter < 100:                     
        logg.warning('Collecting marginals with refinement requires sufficient number of iterations\n'
                     f'It is now set to {refine_iter} and should be at least 100\n')
        

    start = logg.info('minimizing the Block Model')

    is_mudata = False
    adata_list = []
    if type(mdata) == MuData:
        is_mudata = True
        # treat MuData as a list of anndatas
        # just keep it aside as we need to reconstruct it later
        if copy:
            mdata = mdata.copy()
        # adata_list will keep reference to mdata, whether it's a copy or not
        # remember, in case, to put things back    
        adata_list = list(mdata.mod.values())
    else:
        if copy:
            adata_list = [x.copy() for x in mdata]
        else:
            adata_list = mdata

    n_keys = len(neighbors_key)
    n_data = len(adata_list)
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
            if neighbors_key[x] not in adata_list[x].uns:
                raise ValueError(
                    'You need to run `pp.neighbors` first '
                    'to compute a neighborhood graph. for'
                    f'data entry {x}'
                )
            elif 'connectivities_key' in adata_list[x].uns[neighbors_key[x]]:
                # scanpy>1.4.6 has matrix in another slot
                conn_key = adata_list[x].uns[neighbors_key[x]]['connectivities_key']
                adjacency.append(adata_list[x].obsp[conn_key])
            else:
                # scanpy<=1.4.6 has sparse matrix here
                adjacency.append(adata_list[x].uns[neighbors_key[x]]['connectivities'])


    # create a union graph with layers
        
    graph_list = []
    for x in range(n_data):
        g = get_graph_tool_from_adjacency(adjacency[x], directed=directed, use_weights=use_weights)
        # add cell names to graph, this will be used to create
        # layered graph 
        g_names = g.new_vertex_property('string') 
        d_names = adata_list[x].obs_names
        for xn in range(len(d_names)):
            g_names[xn] = d_names[xn]
        g.vp['cell'] = g_names
        graph_list.append(g)
       
    # get a non-redundant list of all cell names across all modalities
    all_names = set(adata_list[0].obs_names)
    [all_names.update(adata_list[x].obs_names) for x in range(1, n_data)]
    all_names = list(all_names)
    # create the shared graph
    union_g = gt.Graph(directed=False)
    union_g.add_vertex(len(all_names))
    u_names = union_g.new_vertex_property('string')
    for xn in range(len(all_names)):
        u_names[xn] = all_names[xn]
    union_g.vp['cell'] = u_names
    
    # check that there are overlapping nodes, otherwise exit
    if union_g.num_vertices() == sum(adata_list[xn].shape[0] for xn in range(n_data)):
        raise ValueError(
                'The number of nodes in the merged graph is the same as'
                'the total number of cells across all datasets, it seems'
                'there are no shared cells across modalities.'
                'Check if this is the case and change cell names so that'
                'shared cells have the same name across modalities.'
            )
        
    
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
    
    if model == "nsbm":
        states = [gt_model(g=union_g,
                          base_type=gt.LayeredBlockState,
                          state_args=dict(deg_corr=deg_corr,
                          ec=union_g.ep.layer,
                          layers=True,
                          overlap=overlap
                          )) for n in range(n_init)]
    
    else:
        states = [gt_model(g=union_g,
                          deg_corr=deg_corr,
                          ec=union_g.ep.layer,
                          layers=True,
                          overlap=overlap
                          ) for n in range(n_init)]

    with parallel_config(backend=dispatch_backend,
                         max_nbytes=None,
                         n_jobs=n_jobs):
        states = list(tqdm(Parallel(return_as="generator")(
            delayed(fast_min)(states[x], beta, n_sweep, tolerance, max_iter, seeds[x]) for x in range(n_init)
                      ),
                      total=n_init))

    if model == "nsbm":
        pmode = gt.PartitionModeState([x.get_bs() for x in states], converge=True, nested=True)
        bs = pmode.get_max_nested()
        bs = [x for x in bs if len(np.unique(x)) > 1]
        bs.append(np.array([0], dtype=np.int32)) #in case of type changes, check this
        state = gt_model(union_g, bs=bs,
                         base_type=gt.LayeredBlockState,
                         state_args=dict(deg_corr=deg_corr,
                         ec=union_g.ep.layer,
                         layers=True,
                         overlap=overlap
                         ))


    else:
        pmode = gt.PartitionModeState([x.get_blocks().a for x in states], converge=True, nested=False)
        bs = pmode.get_max(union_g)
        state = gt_model(union_g, b=bs,
                         deg_corr=deg_corr,
                         ec=union_g.ep.layer,
                         layers=True,
                         overlap=overlap
                         )


    logg.info('        consensus step done', time=start)
        
    # prune redundant levels at the top
    
    if refine_model:
        # we here reuse pmode variable, so that it is consistent
        logg.info('        Refining model')
        bs = []
        if model == "nsbm":
            def collect_partitions(s):
                bs.append(s.get_bs())
        else:
            def collect_partitions(s):
                bs.append(s.get_blocks().a)

        gt.mcmc_equilibrate(state, force_niter=refine_iter, 
                            multiflip=True, 
                            mcmc_args=dict(niter=n_sweep, beta=beta),
                            callback=collect_partitions)

        if model == "nsbm":
            pmode = gt.PartitionModeState(bs, nested=True, converge=True)
            bs = [x for x in pmode.get_max_nested() if len(np.unique(x)) > 1]
            bs.append(np.array([0], dtype=np.int32)) #in case of type changes, check this
            state = gt_model(union_g, bs=bs,
                             base_type=gt.LayeredBlockState,
                             state_args=dict(deg_corr=deg_corr,
                             ec=union_g.ep.layer,
                             layers=True,
                             overlap=overlap
                             ))
        else:
            pmode = gt_model(bs, converge=True)
            bs = pmode.get_max(union_g)
            state = gt_model(union_g, b=bs,
                             deg_corr=deg_corr,
                             ec=union_g.ep.layer,
                             layers=True,
                             overlap=overlap
                             )

        logg.info('        refinement complete', time=start)
    
    
    if save_model:
        import pickle
        fname = save_model
        if not fname.endswith('pkl'):
            fname = f'{fname}.pkl'
        logg.info(f'Saving model into {fname}')    
        with open(fname, 'wb') as fout:
            dump = {'PartitionModeState':pmode,
                    'Graph':g
                    }
            pickle.dump(dump, fout, 2)

    logg.info('    done', time=start)
    # reorganize things so that groups are ordered literals
    if model == "nsbm":
        groups = np.zeros((union_g.num_vertices(), len(bs)), dtype=int)
        u_groups = np.unique(bs[0])
    else:
        groups = np.array(bs.get_array())
        u_groups = np.unique(groups)

    last_group = np.max(u_groups) + 1
    n_groups = len(u_groups)
    
    if model == "nsbm":
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
        state = gt_model(union_g, bs)
        del(i_groups)

        groups.index = all_names

        # add column names
        groups.columns = [f"{key_added}_level_{level}" for level in range(len(bs))]

        # remove any column with the same key
        for xn in range(n_data):
            drop_columns = groups.columns.intersection(adata_list[xn].obs.columns)
            if len(drop_columns) > 0:
                adata_list[xn].obs.drop(drop_columns, axis='columns', inplace=True)
            adata_list[xn].obs = pd.concat([adata_list[xn].obs, groups.loc[adata_list[xn].obs_names]], axis=1)
    else:
        # for ppbm and sbm is simpler
        rosetta = dict(zip(u_groups, range(len(u_groups))))
        groups = np.array([rosetta[x] for x in groups])
        groups = groups.astype('U')
        groups = pd.Series(groups, index=all_names)
        for xn in range(n_data):
            adata_list[xn].obs[key_added] = pd.Categorical(groups.loc[adata_list[xn].obs_names], 
                                                       categories=natsorted(np.unique(groups)),
                                                       )

    # now add marginal probabilities.

    if collect_marginals:
        # note that the size of this will be equal to the number of the groups in Mode
        # but some entries won't sum to 1 as in the collection there may be differently
        # sized partitions
        pv_array = pmode.get_marginal(union_g).get_2d_array(range(last_group)).T[:, u_groups] / n_init    
        for xn in range(n_data):
            # add marginals for level 0, the sum up according to the hierarchy
            _groups = groups.loc[adata_list[xn].obs_names]
            _pv_array = pd.DataFrame(pv_array, index=all_names).loc[adata_list[xn].obs_names].values
            if model == "nsbm":
                adata_list[xn].obsm[f"CM_{key_added}_level_0"] = _pv_array
                for group in groups.columns[1:]:
                    ct = pd.crosstab(_groups[_groups.columns[0]], _groups[group], 
                                     normalize='index', dropna=False)
                    adata_list[xn].obsm[f'CM_{group}'] = _pv_array @ ct.values
            else:
                adata_list[xn].obsm[f"CM_{key_added}"] = _pv_array

    # add some unstructured info
    if model == "nsbm":
        modularity=np.array([gt.modularity(union_g, state.project_partition(x, 0))
                         for x in range(len((state.levels)))])
    else:
        modularity=np.array([gt.modularity(union_g, state.get_blocks())])
    
    for xn in range(n_data):
        if not 'schist' in adata_list[xn].uns:
            adata_list[xn].uns['schist'] = {}

        adata_list[xn].uns['schist'][key_added] = {}
        adata_list[xn].uns['schist'][key_added]['stats'] = dict(
              entropy=state.entropy(),
              modularity=modularity
              )
        if model == "nsbm":
            adata_list[xn].uns['schist'][key_added]['stats']['level_entropy']=np.array([state.level_entropy(x) for x in range(len(state.levels))])

        if model == "nsbm":
            # record state as list of blocks
            # unfortunately this cannot be a list of lists but needs to be a dictionary
            bl_d = {}
            levels = state.get_levels()
            for nl in range(len(levels)):
                bl_d[str(nl)] = np.array(levels[nl].get_blocks().a)
        else:
            bl_d = {'0':np.array(state.get_blocks().a)}        
    
        adata_list[xn].uns['schist'][key_added]['blocks'] = bl_d

        # last step is recording some parameters used in this analysis
        adata_list[xn].uns['schist'][key_added]['params'] = dict(
            model=model,
            use_weights=use_weights,
            neighbors_key=neighbors_key[xn],
            key_added=key_added,
            n_init=n_init,
            collect_marginals=collect_marginals,
            random_seed=random_seed,
            deg_corr=deg_corr,
            refine_model=refine_model,
            refine_iter=refine_iter,
            overlap=overlap,
            directed=directed
        )

    logg.info(
        '    finished',
        time=start,
        deep=(
            f'and added\n'
            f'    {key_added!r}, the cluster labels (adata.obs, categorical)'
        ),
    )
    if copy:
        if is_mudata:
            return mdata
        return adata_list
    return None

    