from typing import Optional#, Tuple, Sequence, Type, Union, Dict

import numpy as np
from anndata import AnnData
import scipy.stats
from scipy import sparse

from scanpy import logging as logg
import graph_tool.all as gt
import pandas as pd
from ._utils import get_cell_loglikelihood
from scanpy._utils import get_igraph_from_adjacency, _choose_graph


def calculate_affinity(
    adata: AnnData,
    level: int = 1,
    block_key: Optional[str] = 'nsbm',
    group_by: Optional[str] = None,
    state: Optional = None,
    neighbors_key: Optional[str] = None,
    adjacency: Optional[sparse.spmatrix] = None,
    directed: bool = True,
    use_weights: bool = True,
    obsp: Optional[str] = None,
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
        if 'schist' in adata.uns and 'state' in adata.uns['schist']:
            state = adata.uns['schist']['state']
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
            ca_matrix = get_cell_loglikelihood(state, as_prob=True)
        else:
            raise ValueError(f"{group_by} should be a categorical entry in adata.obs")    
    else:        
        # use precomputed blocks and states
        if type(state) == gt.NestedBlockState:
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
    state: Optional = None,
    copy: bool = False
) -> Optional[AnnData]:
    """\
    Calculate cell stability given cell affinity
    Parameters
    ----------
    adata
        Annotated data matrix. 
    key
        The prefix of CA matrices in adata.obsm to evaluate
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with stability values 
    in adata.obs['cell_stability']
"""    

    matrix_prefix = 'CA'
    if use_marginals:
        matrix_prefix = 'CM'

    if not state:
        if not adata.uns['schist']['state']:
            raise ValueError("No state detected")
        else:
            state = adata.uns['schist']['state']

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
    are associated to cells mostly close to similar cells
    Parameters
    ----------
    adata
        Annotated data matrix. 
    key_added
        The name of the entry in adata.obs with calculated values
    copy
        Return a copy instead of writing to adata.
    sim_type:
        Similarity function. Can be one in 'dice', 'salton', 'hub-promoted', 
        'hub-suppressed', 'jaccard', 'inv-log-weight', 'resource-allocation',
        'leight-holme-newman'. For more information check here
        https://graph-tool.skewed.de/static/doc/topology.html?highlight=distance#graph_tool.topology.vertex_similarity
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


