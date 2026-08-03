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
        np.random.seed(seed)
        gt.seed_rng(seed)
    dS = 1e9
    n = 0
    while (np.abs(dS) > fast_tol) and (n < max_iter):
        dS, _, _ = state.multiflip_mcmc_sweep(beta=beta, niter=n_sweep, c=0.5)
        n += 1
    return state                            


def fit_model_multi(
    mdata: Union[List[AnnData], MuData],
    nested: bool = True,
    collect_marginals: bool = True,
    n_samples: int = 100,
    key_added: str | None = None,
    adjacency: Optional[List[sparse.spmatrix]] = None,
    neighbors_key: Optional[List[str]] = ['neighbors'],
    constraint_key: Optional[str] = None,
    deg_corr: bool = True,
    directed: bool = False,
    use_weights: bool = False,
    bisection: bool = True,
    simple_init: bool = False,
    overlap: bool = False,
    n_jobs: int = -1,
    n_iter: int = 10,
    beta: float = 1.0,
    save_model: Union[str, None] = None,
    copy: bool = False,
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
    nested
        Wether to use the hierarchical version of SBM (default) or a simple SBM
    collect_marginals
        Collect marginal distribution of cells, that is the probability
        to belong to any cluster. Note that enabling this option requires sampling
        from the posterior and modifying the community structure.
    n_samples
        If marginals are collected, this number of samples is taken from the posterior.
        These are then used to compute the consensus and the marginals.
    key_added
        `adata.obs` key under which to add the cluster labels.
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']` in case of scanpy<=1.4.6 or
        `adata.obsp[neighbors_key][connectivity_key]` for scanpy>1.4.6
    neighbors_key
        The key passed to `sc.pp.neighbors`
    constraint_key
        Use this annotation in adata.obs as constraint when minimizing the model, 
        that is, avoid grouping cells with different values in this field. 
    deg_corr
        Whether to use degree correction in the minimization step. In many
        real world networks this is the case, although this doesn't seem
        the case for KNN graphs used in scanpy.
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges). Note that this
        increases computation times
    bisection
        Bisection search is enabled by default to determine the appropriate number of
        groups. Disable this to obtain a faster computation by agglomerative search.
    simple_init
        This parameter estimates an upper bound on the number of groups at each level, 
        the space of solutions is reduced and so is the time needed to converge.
    overlap
        Whether the different layers are dependent (overlap=True) or not (overlap=False)
    n_jobs
        Number of OMP threads used during model minimization
    n_iter
        When sampling from the posterior, set this number of iterations for the 
        MCMC sweep
    beta
        When sampling from the posterior, set this as the inverse of the temperature.
        Higher values make computation faster but may be less effective
    save_model
        If provided, this will be the filename for the PartitionModeState to 
        be saved. The PartitionModeState contains all the models minimized during 
        inference.
    copy
        Whether to copy data or modify them inplace.
    random_seed
        Random number to be used as seed for graph-tool. Note that setting this 
        to 0 is equivalent to not setting it at all.

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
    `adata.obsm['CM_multi_nsbm_level_{n}']`
        A `np.ndarray` with cell probability of belonging to a specific group
    `adata.uns['schist']['multi_level_state']`
        The NestedBlockModel state object
    """

    # define the base model.
    base_state=gt.LayeredBlockState
    if overlap:
        base_state=gt.LayeredOverlapBlockState
    if use_weights:
        base_state=gt.LayeredWeightedBlockState
        if overlap:
            base_state=gt.LayeredWeightedOverlapBlockState
                
    if nested:
        key_added = 'multi_nsbm' if key_added is None else key_added
    else:
        key_added = 'multi_sbm' if key_added is None else key_added

    # define the function that will be used to minimize    
    f_minimize = gt.minimize_blockmodel_dl
    if nested:
        f_minimize = gt.minimize_nested_blockmodel_dl
    
    if random_seed:
        logg.warning('Setting random seed disables threading during minimization\n'
                     f'due to non deterministic behaviour of thread allocation\n')
    
    if n_samples < 1:
        n_samples = 1
    if collect_marginals and n_samples < 100:
        logg.warning('Collecting marginals requires sufficient number of n_samples\n'
                     f'It is now set to {n_samples} and should be at least 100\n')
        
    start = logg.info('minimizing the Block Model')
    # here we have to build the graph
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
        logg.info(f'getting adjacency for data {x}', time=start)
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
    if union_g.num_vertices() == sum([adata_list[xn].shape[0] for xn in range(n_data)]):
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
        logg.info(f'getting cells for data {ng}', time=start)
        for e in graph_list[ng].edges():
            S, T = e.source(), e.target()
            Sn = graph_list[ng].vp['cell'][S]
            Tn = graph_list[ng].vp['cell'][T]
            Sidx = u_cell_index[Sn]
            Tidx = u_cell_index[Tn]
            ne = union_g.add_edge(Sidx, Tidx)
            layer[ne] = ng + 1 # this is the layer label, must be >0

    union_g.ep['layer'] = layer
    # DONE! now proceed with standard minimization, ish

    if union_g.num_vertices() > 1e5:
        if not save_model:
            logg.warning('When working with large networks it is a good idea\n'
                         'to save the model and load it separately for further analysis.\n'
                        )
        if collect_marginals:
            logg.warning('When working with large networks it may be appropriate to\n'
                         'skip sampling and do not collect marginals.\n'
                         'You can save the model and perform sampling in a separate experiment.\n'
                        )
        if not simple_init or bisection:
            logg.warning('When working with large networks it may be appropriate to\n'
                         'enable fast minimization setting the ```bisection``` and\n'
                         '```simple_init``` accordingly.\n'
                        )

    state_args={}
    base_state_args={}
    f_args={}

    state_args['deg_corr']=deg_corr
    base_state_args['ec']=union_g.ep.layer
    base_state_args['overlap']=overlap
    if use_weights:
#### TODO ####
#### ADD weights to union_g
        base_state_args['rec']=[union_g.ep.weight] 
        base_state_args['rec_types']=['real-exponential']
    
    f_args['state_args']=state_args
    if nested:
        f_args['base_state']=base_state
        f_args['base_state_args']=base_state_args
    else:
        f_args['state']=base_state
        f_args['state_args'].update(base_state_args) #merge...
    
    if simple_init and nested:
        bisection=False
    f_args['multilevel_mcmc_args'] = {'bisection':bisection}
    if nested:
        f_args['simple_init']=simple_init


    if constraint_key:
        # the key must be present in all adata. Moreover, we require the same
        # categories to be in all adata. Collect all info in a single dict
        # then use it to create the pclabel
        d_constraint = dict.from_keys(all_names)
        for data in adata_list:
            if not constraint_key in data.obs.columns:
                raise NameError(f"{constraint_key} was not found all your datasets")
            if data.obs[constraint_key].dtype.name != 'category':
                raise AttributeError(f"{constraint_key} must be categorical in all datasets")
            if data.obs[constraint_key].cat.categories != adata_list[0].obs[constraint_key].cat.categories:
                raise AttributeError(f"{constraint_key} must contain the same cateogires in all datasets")
            _s = data.obs[constraint_key].cat.codes
            for c in _s:
                d_constraint[c] = _s[c]
        pclabel = g.new_vp('int64_t')
        pclabel.a = np.array([d_constraint[c]] for c in all_names)
        if nested:
            raise NotImplementedError("Constraints do not (yet) work with NSBM")
#            base_state_args['pclabel'] = pclabel
#            base_state_args['clabel'] = pclabel
        elif assortative:
            raise NotImplementedError("Constraints do not (yet) work with PPBM")
#            base_state_args['pclabel'] = pclabel
#            base_state_args['clabel'] = pclabel
        else:
            state_args['pclabel'] = pclabel
            state_args['clabel'] = pclabel

    current_omp_threads = gt.openmp_get_num_threads()
    if n_jobs <= 0:
        # this because context does not understand negative values
        n_threads_set=current_omp_threads
        n_jobs=current_omp_threads
    
    if random_seed:
        gt.seed_rng(random_seed)
        n_threads_set=1 #disable threading to have deterministic behaviour

    # do the actual minimization
    with gt.openmp_context(nthreads=n_threads_set, schedule='guided'):
        state=f_minimize(g, **f_args)
        if nested:
            # this is needed in case marginals are not collected
            bs=[np.array(x) for x in state.get_bs() if len(np.unique(x)) > 1]
            bs.append(np.array([0], dtype=np.int32)) 
        else:
            bs=np.array(state.b.a)

    logg.info('        done', time=start)


    if collect_marginals:
        logg.info('Sampling posterior and getting cell marginals')
        bs = []
        if random_seed:
            gt.seed_rng(random_seed)
        # create seeds for each MCMC sweep
        with gt.openmp_context(nthreads=n_jobs, schedule='guided'):
            for n in tqdm(range(n_samples)):
                state.multiflip_mcmc_sweep(niter=n_iter, beta=beta)
                if nested:
                    bs.append(state.get_bs())
                else:
                    bs.append(state.b.a.copy())

        logg.info('        done', time=start)
        logg.info('Computing the consesus')
        pmode=gt.PartitionModeState(bs, converge=True, nested=nested)
        if nested:
            bs=pmode.get_max_nested()
            # prune redundant levels at the top
            bs = [x for x in bs if len(np.unique(x)) > 1]
            bs.append(np.array([0], dtype=np.int32)) 
            
        else:
            bs=pmode.get_max(union_g)
        pv=pmode.get_marginal(union_g)

    logg.info('        done', time=start)

    if nested:
        state=gt.NestedBlockState(union_g, bs=bs,
                                  base_state=base_state,
                                  base_state_args=base_state_args,

        )
    else:
        if use_weights:
            state=gt.WeightedBlockState(g=union_g, b=bs, 
                                        **state_args,
                                        **base_state_args)
        else:
            state=gt.BlockState(g=union_g, b=bs, **state_args)
        
    # prune redundant levels at the top
    
    if save_model:
        import pickle
        fname = save_model
        if not fname.endswith('pkl'):
            fname = f'{fname}.pkl'
        logg.info(f'Saving model into {fname}')    
        dump = {'State':state,
                'Graph':union_g
                }
        if collect_marginals:
            dump['PartitionModeState']=pmode
            
        with open(fname, 'wb') as fout:
            pickle.dump(dump, fout, 2)

    logg.info('    done', time=start)
    # reorganize things so that groups are ordered literals
    if nested:
        groups = np.zeros((union_g.num_vertices(), len(bs)), dtype=int)
        u_groups = np.unique(bs[0])
    else:
        groups = np.array(bs.get_array())
        u_groups = np.unique(groups)

    last_group = np.max(u_groups) + 1
    n_groups = len(u_groups)
    
    if nested:
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
        state = gt.NestedBlockState(union_g, bs=bs,
                             base_state=base_state,
                             base_state_args=base_state_args
                             )

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
        # for sbm is simpler
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
        pv_array = pv.get_2d_array(range(last_group)).T[:, u_groups] / n_samples    
        for xn in range(n_data):
            # add marginals for level 0, the sum up according to the hierarchy
            _groups = groups.loc[adata_list[xn].obs_names]
            _pv_array = pd.DataFrame(pv_array, index=all_names).loc[adata_list[xn].obs_names].values
            if nested:
                adata_list[xn].obsm[f"CM_{key_added}_level_0"] = _pv_array
                for group in groups.columns[1:]:
                    ct = pd.crosstab(_groups[_groups.columns[0]], _groups[group], 
                                     normalize='index', dropna=False)
                    adata_list[xn].obsm[f'CM_{group}'] = _pv_array @ ct.values
            else:
                adata_list[xn].obsm[f"CM_{key_added}"] = _pv_array

    # add some unstructured info
    if nested:
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
        if nested:
            adata_list[xn].uns['schist'][key_added]['stats']['level_entropy']=np.array([state.level_entropy(x) for x in range(len(state.levels))])

        if nested:
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
            nested=nested,
            neighbors_key=neighbors_key[xn],
            key_added=key_added,
            use_weights=use_weights,
            n_init=n_init,
            collect_marginals=collect_marginals,
            random_seed=random_seed,
            deg_corr=deg_corr,
            overlap=overlap,
            directed=directed,
            n_iter=n_iter,
            beta=beta
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

    