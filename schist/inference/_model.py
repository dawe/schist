from typing import Optional, Tuple, Sequence, Type, Union, Dict, Literal

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from natsort import natsorted
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
    
def fit_model(
    adata: AnnData,
    nested: bool = True,
    assortative: bool = False,
    collect_marginals: bool = True,
    n_samples: int = 100,
    key_added: str | None = None,
    adjacency: Optional[sparse.spmatrix] = None,
    neighbors_key: Optional[str] = 'neighbors',
    constraint_key: Optional[str] = None,
    deg_corr: bool = True,
    directed: bool = False,
    use_weights: bool = False,
    bisection: bool = True,
    simple_init: bool = False,
    n_jobs: int = -1,
    n_iter: int = 10,
    beta: float = 1.0,
    save_model: Union[str, None] = None,
    copy: bool = False,
    random_seed: Optional[int] = None,
) -> Optional[AnnData]:
    """\
    Cluster cells using the nested Stochastic Block Model [Peixoto14]_,
    performing Bayesian inference on node groups. 
    
    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    nested
        Wether to use the hierarchical version of SBM (default) or a simple SBM
    assortative
        Wether to use the planted partition model, using only a prior on assorativity.
        Note that setting this option disables the ability to use hierarchical models.
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
        Whether to copy `adata` or modify it inplace.
    random_seed
        Random number to be used as seed for graph-tool. Note that setting this 
        to 0 is equivalent to not setting it at all.
    
    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of cells) that stores the subgroup id
        (`'0'`, `'1'`, ...) for each cell. 
    `adata.uns['schist'][model]['stats']`
        A dict with entropy and modularity values
    `adata.uns['schist'][model]['params']`
        A dict with the values for the parameters used
    `adata.obsm['CM_nsbm_level_{n}']` or `adata.obsm['CM_model']`
        A `np.ndarray` with cell probability of belonging to a specific group if
        marginals are computed
    `adata.uns['schist'][model]['state']`
        The block model, to be used in case a gt state should be initialized
    """

    # define the block model class to be used
    if assortative:
        base_state = gt.PPBlockState
        key_added = 'ppbm' if key_added is None else key_added      
        nested=False #force this  
    else:
        base_state = gt.BlockState
        if nested:
            key_added = 'nsbm' if key_added is None else key_added
        else:
            key_added = 'sbm' if key_added is None else key_added

        
    if use_weights:
        # if weighted, it cannot be assortative, overrides any previous choice
        base_state = gt.WeightedBlockState
        if simple_init==False:
            logg.warning('Working with weighted graphs usually requires\n'
                     f'considerably more time, consider enabling speed options\n'
                     f'as ```simple_init``` or disabling bisection')

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
                     
    start = logg.info('Minimizing the Model')
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

    # convert it to igraph and graph-tool
    g = get_graph_tool_from_adjacency(adjacency, directed=directed, use_weights=use_weights)
    
    if g.num_vertices() > 1e5:
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

    if not assortative:
        state_args['deg_corr']=deg_corr          

    if use_weights:
        base_state_args['rec']=[g.ep.weight]
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
        if not constraint_key in adata.obs.columns:
            raise NameError(f"{constraint_key} was not found in your dataset")
        if adata.obs[constraint_key].dtype.name != 'category':
            raise AttributeError(f"{constraint_key} must be categorical")
        pclabel = g.new_vp('int64_t')
        pclabel.a = np.array(adata.obs[constraint_key].cat.codes)
        if nested:
            base_state_args['pclabel'] = pclabel
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
            bs=pmode.get_max(g)
        pv=pmode.get_marginal(g)

    logg.info('        done', time=start)

    if nested:
        state=gt.NestedBlockState(g, bs=bs,
                                  base_state=base_state,
                                  base_state_args=base_state_args,

        )
    elif assortative:
        state=gt.PPBlockState(g, b=bs, **state_args)
    else:
        if use_weights:
            state=gt.WeightedBlockState(g=g, b=bs, 
                                        **state_args,
                                        **base_state_args)
        else:
            state=gt.BlockState(g=g, b=bs, **state_args)

    if save_model:
        import pickle
        fname = save_model
        if not fname.endswith('pkl'):
            fname = f'{fname}.pkl'
        logg.info(f'Saving model into {fname}')    
        dump = {'State':state,
                'Graph':g
                }
        if collect_marginals:
            dump['PartitionModeState']=pmode
            
        with open(fname, 'wb') as fout:
            pickle.dump(dump, fout, 2)
    
    # reorganize things so that groups are ordered literals
    if nested:
        groups = np.zeros((g.num_vertices(), len(bs)), dtype=int)
        u_groups = np.unique(bs[0])
    else:
        groups = np.array(bs)
        u_groups = np.unique(groups)

    last_group = np.max(u_groups) + 1
    n_groups = len(u_groups)

    if nested:
        # this is harder as we have to consider and refactor the entire hierarchy
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

        state=gt.NestedBlockState(g, bs=bs,
                                  base_state=base_state,
                                  base_state_args=base_state_args
        )

        del(i_groups)
    
        groups.index = adata.obs_names
    
        # add column names
        groups.columns = [f"{key_added}_level_{level}" for level in range(len(bs))]
    
        # remove any column with the same key
        keep_columns = [x for x in adata.obs.columns if not x.startswith('%s_level_' % key_added)]
        adata.obs = adata.obs[keep_columns]
        adata.obs = pd.concat([adata.obs, groups], axis=1)

    else:
        # for ppbm and sbm is simpler
        rosetta = dict(zip(u_groups, range(len(u_groups))))
        groups = np.array([rosetta[x] for x in groups])
        groups = groups.astype('U')
        adata.obs[key_added] = pd.Categorical(values=groups,
                                              categories=natsorted(np.unique(groups)),
                                              )

    # now add marginal probabilities.
    if collect_marginals:
        # note that the size of this will be equal to the number of the groups in Mode
        # but some entries won't sum to 1 as in the collection there may be differently
        # sized partitions
        pv_array = pv.get_2d_array(range(last_group)).T[:, u_groups] / n_samples
        if nested:
            # add marginals for level 0, the sum up according to the hierarchy
            adata.obsm[f"CM_{key_added}_level_0"] = pv_array
            for group in groups.columns[1:]:
                ct = pd.crosstab(groups[groups.columns[0]], groups[group], normalize='index')
                adata.obsm[f'CM_{group}'] = pv_array @ ct.values
        else:
            adata.obsm[f"CM_{key_added}"] = pv_array

    # add some unstructured info
    if nested:
        modularity=np.array([gt.modularity(g, state.project_partition(x, 0))
                         for x in range(len((state.levels)))])
    else:
        modularity=np.array([gt.modularity(g, state.get_blocks())])

    if not 'schist' in adata.uns:
        adata.uns['schist'] = {}
    adata.uns['schist'][key_added] = {}
    adata.uns['schist'][key_added]['stats'] = dict(
              entropy=state.entropy(),
              modularity=modularity
              )
    if nested:
        adata.uns['schist'][key_added]['stats']['level_entropy']=np.array([state.level_entropy(x) for x in range(len(state.levels))])

    if nested:
        # record state as list of blocks
        # unfortunately this cannot be a list of lists but needs to be a dictionary
        bl_d = {}
        levels = state.get_levels()
        for nl in range(len(levels)):
            bl_d[str(nl)] = np.array(levels[nl].get_blocks().a)
    else:
        bl_d = {'0':np.array(state.get_blocks().a)}        
    
    adata.uns['schist'][key_added]['blocks'] = bl_d

    # last step is recording some parameters used in this analysis
    adata.uns['schist'][key_added]['params'] = dict(
        nested=nested,
        assortative=assortative,
        neighbors_key=neighbors_key,
        use_weights=use_weights,
        key_added=key_added,
        n_samples=n_samples,
        collect_marginals=collect_marginals,
        random_seed=random_seed,
        deg_corr=deg_corr,
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
    del g
    return adata if copy else None

