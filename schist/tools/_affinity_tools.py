from typing import Optional, Tuple, Sequence, Type, Union, Dict

import numpy as np
from anndata import AnnData
import scipy.stats
from scipy import sparse

from scanpy import logging as logg
import graph_tool.all as gt
import pandas as pd
from .._utils import get_cell_loglikelihood, get_cell_back_p, state_from_blocks
from scanpy._utils import get_igraph_from_adjacency, _choose_graph


def calculate_affinity(
    adata: AnnData,
    level: int = 1,
    block_key: Optional[str] = 'nsbm',
    group_by: Optional[str] = None,
    state: Optional = None,
    neighbors_key: Optional[str] = 'neighbors',
    adjacency: Optional[sparse.spmatrix] = None,
    directed: bool = False,
    use_weights: bool = False,
    obsp: Optional[str] = None,
    back_prob: bool = False,
    copy: bool = False
    
) -> Optional[AnnData]:
    
    """\
    Calculate cell affinity given a partition scheme. It can be used for 
    partitions calculated using schist or for any partition scheme, given
    for example by cell annotations.
    
    Parameters
    ----------
    adata:
        The AnnData object. Should have been already processed with schist
    level:
        The level to calculate affinity. This parameter is effective
        only for Nested partitions
    block_key:
        The prefix for partitions. This parameter is ignored if the state
        is not gt.NestedBlockState
    group_by:
        The key for group names used for calculations. Setting this will override
        level and block_key. This is effective only for NestedBlockState partitions
    state:
        Optionally calculate affinities on this state.
    neighbors_key
        Use neighbors connectivities as adjacency.
        If not specified, leiden looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, leiden looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    adjacency
        Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    copy:
        Return a new object or do everything in place
        

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with affinity values 
    in adata.obsm[f'CA_{block_key}_level_{level}']
        
"""    

    matrix_key = f'CA_{block_key}_level_{level}' # the default name of the matrix
    if group_by:
        logg.info(f'Calculating cell affinity to {group_by}')
    else:
        logg.info(f'Calculating cell affinity to level {level}')
        
    if not state:
        # if no state is provided, use the default to retrieve graph
        if 'schist' in adata.uns and 'blocks' in adata.uns['schist'][f'{block_key}']:
            params = adata.uns['schist'][f'{block_key}']['params']
            if 'neighbors_key' in params:
                neighbors_key=params['neighbors_key']
            if 'use_weights' in params:
                use_weights=params['use_weights']
            if 'deg_corr' in params:
                deg_corr=params['deg_corr']
            state = state_from_blocks(adata, 
                                  state_key=block_key,
                                  neighbors_key=neighbors_key,
                                  adjacency=adjacency,
                                  directed=directed,
                                  use_weights=use_weights,
                                  deg_corr=deg_corr
                                  )
            g = state.g                                  
        elif not neighbors_key:
            # no state and no adjacency provided, raise an error
            raise ValueError("A state or an adjacency matrix should be given"
                             "Otherwise a graph cannot be computed")
        else:
            # get the graph from the adjacency    
            adjacency = _choose_graph(adata, obsp, neighbors_key)
            g = get_igraph_from_adjacency(adjacency, directed=directed)
            g = g.to_graph_tool()
            gt.remove_parallel_edges(g)
            state = gt.BlockState(g)
    else:
        g = state.g        
    
    if group_by:
        matrix_key = f'CA_{group_by}'
        # if groups are given, we generate a new BlockState and work on that
        if group_by in adata.obs.columns and adata.obs[group_by].dtype.name == 'category':
            partitions = adata.obs[group_by].cat.codes.values
            state = gt.BlockState(g, b=partitions)
            if back_prob:
                ca_matrix = get_cell_back_p(state)
            else:
                ca_matrix = get_cell_loglikelihood(state, as_prob=True)
        else:
            raise ValueError(f"{group_by} should be a categorical entry in adata.obs")    
    else:        
        # use precomputed blocks and states
        if type(state) == gt.NestedBlockState:
            if back_prob:
                p0 = get_cell_back_p(state, level=0)
            else:
                p0 = get_cell_loglikelihood(state, level=0, as_prob=True)
            group_col = None
            if group_by and group_by in adata.obs.columns:
                group_col = group_by
            else:
                g_name = f'{block_key}_level_{level}'
                if g_name in adata.obs.columns:
                    group_col = g_name
            if not group_col:
                raise ValueError("The provided groups or level/blocks do not exist")
            
            g0 = pd.Categorical(state.project_partition(0, 0).a)
            cross_tab = pd.crosstab(g0, adata.obs[group_col], normalize='index')
            ca_matrix = (p0 @ cross_tab).values

        elif type(state) == gt.PPBlockState:
            if back_prob:
                ca_matrix = get_cell_back_p(state)
            else:
	            ca_matrix = get_cell_loglikelihood(state, as_prob=True)
            matrix_key = 'CA_ppbm'
    
    adata.obsm[matrix_key] = ca_matrix 
    
    return adata if copy else None


def cluster_consistency(
    adata: AnnData,
    groups: str = None,
    key_added: Optional[str] = 'cluster_consistency',
    use_marginals: Optional[bool] = False,
    copy: bool = False
) -> Optional[AnnData]:
    """\
    Calculate cluster consistency at a given level
    
    Parameters
    ----------
    adata
        Annotated data matrix. 
    groups
        The key for clusters in adata.obs
    key_added
        The name of obs values that will be added to the adata
    use_marginals
        By default it uses cell affinities for the analysis, but if group marginals
        are available from the inference, those can be used here.
    copy
        Return a copy instead of writing to adata.


    Returns
    -------
    Depending on `copy`, returns or updates `adata` with consistency values 
    in adata.uns['cluster_consistency'] and adata.obs['cluster_consistency']
"""    

    matrix_prefix = 'CA'
    if use_marginals:
        matrix_prefix = 'CM'

    if not groups or not groups in adata.obs.columns:
        raise ValueError("Valid groups should be specified")
    else:
         ca_key = f'{matrix_prefix}_{groups}'
         if not ca_key in adata.obsm.keys():
             msg = "Affinities for the provided group were not calculated"
             if use_marginals:
                 msg = "Marginals for the provided group were not calculated"
             raise ValueError(msg)

    affinity = adata.obsm[ca_key]
    entropy = scipy.stats.entropy(affinity, axis=0) / np.log(adata.shape[0]) #normalized entropy

    adata.uns['cluster_consistency'] = entropy

    # now assign consistency to each cell, according to their group
    e_dict = dict(zip(adata.obs[groups].cat.categories, entropy))
    g = adata.obs[groups].values
    adata.obs['cluster_consistency'] = [e_dict[g[x]] for x in range(adata.shape[0])]
    
    return adata if copy else None

def cell_stability(
    adata: AnnData,
    block_key: Optional[str] = 'nsbm', # dummy default
    key_added: Optional[str] = 'cell_stability',
    use_marginals: Optional[bool] = False,
    neighbors_key: Optional[str] = 'neighbors',
    adjacency: Optional[sparse.spmatrix] = None,
    directed: bool = False,
    use_weights: bool = False,
    obsp: Optional[str] = None,    
    state: Optional = None,
    back_prob: bool = False,
    copy: bool = False
) -> Optional[AnnData]:
    """\
    Calculate cell stability given cell affinity.
    
    Parameters
    ----------
    adata
        Annotated data matrix. 
    key
        The prefix of CA matrices in adata.obsm to evaluate.
    copy
        Return a copy instead of writing to adata.
    use_marginals
        Whether to use marginals in place of affinities

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with stability values 
    in adata.obs['cell_stability']
"""    

    matrix_prefix = 'CA'
    if use_marginals:
        matrix_prefix = 'CM'

    if not state:
        if not adata.uns['schist'][f'{block_key}']['blocks']:
            raise ValueError("No state detected")
        else:
            params = adata.uns['schist'][f'{block_key}']['params']
            if 'neighbors_key' in params:
                neighbors_key=params['neighbors_key']
            if 'use_weights' in params:
                use_weights=params['use_weights']
            if 'deg_corr' in params:
                deg_corr=params['deg_corr']
            state = state_from_blocks(adata, 
                                  state_key=block_key,
                                  neighbors_key=neighbors_key,
                                  adjacency=adjacency,
                                  directed=directed,
                                  use_weights=use_weights,
                                  deg_corr=deg_corr
                                  )

    # check if we have levels we want to prune
    n_effective_levels = sum([x.get_nonempty_B() > 1 for x in state.get_levels()])
    n_effective_levels = min(n_effective_levels, len(state.get_levels()))
    obsm_names = [x for x in adata.obsm if x.startswith(f"{matrix_prefix}_{block_key}_level")]    
    if len(obsm_names) < n_effective_levels:
        logg.warning("Your dataset doesn't contain all the required matrices\n")
        if use_marginals:
            logg.warning("Marginals cannot be recomputed from current data, switching to affinities")
            matrix_prefix='CA'
        logg.warning("They will be recalculated from scratch")
        if back_prob:
            adata.obsm[f'{matrix_prefix}_{block_key}_level_0'] = get_cell_back_p(state, level=0)
        else:
            adata.obsm[f'{matrix_prefix}_{block_key}_level_0'] = get_cell_loglikelihood(state, level=0, 
                                                                      as_prob=True)
        obsm_names = [f'{matrix_prefix}_{block_key}_level_0']
        for n in range(n_effective_levels):
            calculate_affinity(adata, level = n+1, block_key=block_key, state=state)
            obsm_names.append(f'{matrix_prefix}_{block_key}_level_{n}')

    # take only matrices with at least 2 groups
    obsm_names = [x for x in obsm_names if adata.obsm[x].shape[1] > 1] 
    
    # take the max value for each matrix
    _M = np.array([np.max(adata.obsm[x], axis=1) for x in obsm_names]).T 
    
    # set a threshold given by a uniform distribution 
    # this is questionable, may be improved
    thr = np.array([1 - 1/adata.obsm[x].shape[1] for x in obsm_names])
    
    # use the fraction of levels that are over the level specific threshold
    _S = np.sum(_M > thr, axis=1) / _M.shape[1]
    adata.obs[f'{key_added}'] = _S
#    _S = np.array([scipy.stats.entropy(adata.obsm[x], axis=1) /np.log(adata.obsm[x].shape[1]) for x in obsm_names]).T
#    adata.obs[f'{key_added}'] = 1-np.nanmax(_S, axis=1) #/ np.nanmean(EE, axis=1)


    return adata if copy else None


def cell_similarity(
    adata: AnnData,
    key_added: Optional[str] = 'cell_similarity',
    sim_type: Optional[str] = 'hub-promoted',
    use_weights: Optional[bool] = True,
    copy: bool = False,
    **neighbors_kwds
) -> Optional[AnnData]:
    """\
    Calculate cell similarity score based on the kNN graph. Higher scores
    are associated to cells mostly close to similar cells.
    
    Parameters
    ----------
    adata
        Annotated data matrix. 
    key_added
        The name of the entry in adata.obs with calculated values.
    copy
        Return a copy instead of writing to adata.
    sim_type:
    	Similarity function. Can be one in 'dice', 'salton', 'hub-promoted','hub-suppressed', 'jaccard', 'inv-log-weight', 'resource-allocation','leight-holme-newman'. For more information check here https://graph-tool.skewed.de/static/doc/topology.html?highlight=distance#graph_tool.topology.vertex_similarity
    state
        A separate block state object

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with stability values 
    in adata.obs['cell_stability']
"""    
    from .._utils import get_graph_tool_from_adata
    logg.info("Adding cell similarity scores")
    g = get_graph_tool_from_adata(adata, use_weights=use_weights, **neighbors_kwds)
    n_cells = g.num_vertices()
    S = gt.vertex_similarity(g, sim_type=sim_type).get_2d_array(range(n_cells))
    D = np.dot(S, S)
    D = np.diag(D / np.max(D)) # take the scaled diagonal 
    adata.obs[f'{key_added}'] = D
    return adata if copy else None


def label_transfer(
    adata: AnnData,
    adata_ref: Optional[AnnData] = None,
    obs: Optional[str] = None,
    label_unk: Optional[str] = 'unknown',
    use_best: Optional[bool] = False,
    neighbors_key: Optional[str] = 'neighbors',
    adjacency: Optional[sparse.spmatrix] = None,
    directed: bool = False,
    use_weights: bool = False,
    pca_args: Optional[dict] = {},
    hamony_args: Optional[dict] = {},
    copy: bool = False
    
) -> Optional[AnnData]:
    
    """\
    Transfer annotation from one dataset to another using cell affinities
    on the kNN graph. If two datasets are given, it uses harmony to perform
    integration and then the kNN graph. If only no reference is given, it is assumed
    that the only adata already contains the proper kNN graph and that
    labels to be reassigned have a specified value.
    
    Parameters
    ----------
    adata:
        The AnnData object. 
    adata_ref
        The optional reference dataset. If None, then all the needed information
        should be included in `adata` (i.e. the kNN graph and the labels)
    obs
        The label that needs to be transfered. Should be in `adata_ref.obs` or in 
        `adata.obs` if no `adata_ref` is given
    label_unk
        The label for unassigned cells. If no `adata_ref` is given, this label 
        identifies cells to be assigned in `adata`. If `adata_ref` is given, this
        label will be given to all cells that cannot be assigned.
    use_best
        When assigning labels, some cells may have not enough evidence and, therefore, 
        left `unknown`. If this parameter is set to `True`, all cells will be assigned
        to the best possible, even if it may not be optimal
    neighbors_key
        Use neighbors connectivities as adjacency.
        If not specified, leiden looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, leiden looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    adjacency
        Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    pca_args
        Parameters to be passed to `sc.tl.pca` before harmony is issued
    harmony_args
    	Parameters to be passed to `sc.external.pp.harmony_integrate`
    copy:
        Return a new object or do everything in place
        

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with added labels 
    in adata.obs[f'{label_ref}']
        
"""    

    if adata_ref:
        # we have to create a merged dataset and integrate
        # before that check that the labels are not in the recipient, in case drop
        
        if obs in adata.obs_keys():
            logg.warning(f'{obs} was found in dataset 1, it will be wiped')
            adata.obs.drop(obs, inplace=True, axis='columns')

        if not obs in adata_ref.obs_keys():
            raise ValueError(
                f'Annotation {obs} is not present in reference dataset.'
            )         

        # now do the merge, so that the empty category is now created
        adata_merge = adata.concatenate(adata_ref, batch_categories=['_unk', '_ref'],
                                        batch_key='_label_transfer')
        # 
        adata_merge.obs[obs] = adata_merge.obs[obs].cat.add_categories(label_unk).fillna(label_unk)
        
        # perform integration using harmony
        sc.tl.pca(adata_merge, **pca_args)
		sc.external.pp.harmony_integrate(adata_merge, 
		                                 key='_label_transfer', 
		                                 **harmony_args)
        # now calculate the kNN graph		                                 
        n_neighbors = int(np.sqrt(adata_merge.shape[0])/2)
        sc.pp.neighbors(adata_merge, use_rep='X_pca_harmony', 
                        n_neighbors=n_neighbors, key_added=neighbors_key) 
    else:
        adata_merge = adata#.copy()
        if not obs in adata_merge.obs_keys():
            raise ValueError(
                f'Annotation {obs} is not present in dataset.'
            )         
        if not label_unk in adata_merge.obs[obs].cat.categories:
            raise ValueError(
                f'Label {label_unk} is not present in {obs}.'
            ) 

    # calculate affinity
    
    calculate_affinity(mdata, group_by=obs, neighbors_key=neighbors_key)