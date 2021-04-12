from typing import Optional#, Tuple, Sequence, Type, Union, Dict

import numpy as np
from anndata import AnnData
import scipy.stats

from scanpy import logging as logg
import graph_tool.all as gt
import pandas as pd
from ._utils import get_cell_loglikelihood


def calculate_affinity(
    adata: AnnData,
    level: int = 1,
    block_key: Optional[str] = 'nsbm',
    group_by: Optional[str] = None,
    state: Optional = None,
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
    copy:
        Return a new object or do everything in place
        
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with affinity values 
    in adata.obsm[f'CA_{block_key}_level_{level}']
        
"""    

    if groups:
        logg.info(f'Calculating cell affinity to {groups}')
    else:
        logg.info(f'Calculating cell affinity to level {level}')
        
    if group_by:
        if group_by in adata.obs.columns and adata.obs[group_by].dtype.name == 'category':
            partitions = adata.obs[group_by].cat.codes
            

    if not state:
        if not adata.uns['schist']['state']:
            raise ValueError("No state detected")
        else:
            state = adata.uns['schist']['state']
    
    if type(state) == gt.NestedBlockState:
        p0 = get_cell_loglikelihood(state, level=0, as_prob=True)
        group_col = None
        if groups and groups in adata.obs.columns:
            group_col = groups
        else:
            g_name = f'{block_key}_level_{level}'
            if g_name in adata.obs.columns:
                group_col = g_name
        if not group_col:
            raise ValueError("The provided groups or level/blocks do not exist")
            
        g0 = pd.Categorical(state.project_partition(0, 0).a)
        cross_tab = pd.crosstab(g0, adata.obs[group_col], normalize='index')
        ca_matrix = p0 @ cross_tab

    elif type(state) == gt.PPBlockState:
        ca_matrix = get_cell_loglikelihood(state, as_prob=True)
        level = 1
        block_key = 'ppbm'
    
    adata.obsm[f'CA_{block_key}_level_{level}'] = ca_matrix 
    
    return adata if copy else None


def cluster_consistency(
    adata: AnnData,
    level: int = 1,
    group: Optional[str] = None,
    key: Optional[str] = 'nsbm',
    copy: bool = False
) -> Optional[AnnData]:
    """\
    Calculate cluster consistency at a given level
    Parameters
    ----------
    adata
        Annotated data matrix. 
    level
        The NSBM level, as an alternative of full group name
    group
        The name of the NSBM level for which consistency should be calculated        
    key
        The key used to store NSBM groupings
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with consistency values 
    in adata.uns['cluster_consistency'] and adata.obs['cluster_consistency']
"""    

    if group:
        level = group.split('_')[-1]
    else:
        group = f'{key}_level_{level}'

    if not group and not level:
        raise ValueError("You should specify at least one of group or level")

    if not f'CA_{key}_level_{level}' in adata.obsm_keys():
        raise ValueError(
            f"Affinitity for the specfified level {level} do not exist"
        )
        

    affinity = adata.obsm[f'CA_{key}_level_{level}']
    entropy = scipy.stats.entropy(affinity, axis=0) / np.log(adata.shape[0]) #normalized entropy

    adata.uns['cluster_consistency'] = entropy

    # now assign consistency to each cell, according to their group
    e_dict = dict(zip(adata.obs[group].cat.categories, entropy))
    g = adata.obs[group].values
    adata.obs['cluster_consistency'] = [e_dict[g[x]] for x in range(adata.shape[0])]
    
    return adata if copy else None


def cell_stability(
    adata: AnnData,
    block_key: Optional[str] = 'nsbm', # dummy default
    key_added: Optional[str] = 'cell_stability',
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

    if not state:
        if not adata.uns['schist']['state']:
            raise ValueError("No state detected")
        else:
            state = adata.uns['schist']['state']

    n_effective_levels = sum([x.get_nonempty_B() > 1 for x in state.get_levels()])
    n_effective_levels = min(n_effective_levels, len(state.get_levels()))
    obsm_names = [x for x in adata.obsm if x.startswith(f"CA_{block_key}_level")]    
    if len(obsm_names) < n_effective_levels:
        logg.warning("Your dataset doesn't contain all the required affinities\n"
                     "They will be recalculated from scratch")
        adata.obsm[f'CA_{block_key}_level_0'] = get_cell_loglikelihood(state, level=0, 
                                                                      as_prob=True)
        obsm_names = [f'CA_{block_key}_level_0']
        for n in range(n_effective_levels):
            calculate_affinity(adata, level = n+1, block_key=block_key, state=state)
            obsm_names.append(f'CA_{block_key}_level_{n}')

    _S = np.array([scipy.stats.entropy(adata.obsm[x], axis=1) /np.log(adata.obsm[x].shape[1]) for x in obsm_names]).T
    adata.obs[f'{key_added}'] = 1-np.nanmax(_S, axis=1) #/ np.nanmean(EE, axis=1)

    return adata if copy else None

def max_marginal(
    adata: AnnData,
    key: Optional[str] = 'nsbm', # dummy default
    key_added: Optional[str] = 'max_marginal',
    copy: bool = False,
    n_iter: int = 100,
    state: Optional = None,
    level: Optional[int] = 0
) -> Optional[AnnData]:
    """\
    Perform a MCMC sweep and calculate the maximal marginal probability
    for a cell. For Nested Model returns the max probability at the lowest level.
    Parameters
    ----------
    adata
        Annotated data matrix. 
    key
        The prefix of CA matrices in adata.obsm to evaluate
    key_added
        The name of the entry in adata.obs with calculated values
    copy
        Return a copy instead of writing to adata.
    n_iter:
        Number of iterations to collect. The higher this number, the higher the 
        precision
    state
        A separate block state object

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with stability values 
    in adata.obs['cell_stability']
"""    
    if not adata.uns['schist']['state'] and not state:
        raise ValueError(
            "A BlockState should be passed to this function"
        )
    if not state:
        state = adata.uns['schist']['state']

    logg.info(f'Collecting marginals for {n_iter} iterations')
    
    nested = False
    if type(state) == gt.NestedBlockState:
        nested = True
    bs = []
    def collect_partitions(s):
        if type(s) == gt.NestedBlockState:
            bs.append(s.get_bs())
        if type(s) == gt.PPBlockState:
            bs.append(s.get_blocks().a.copy())
    
    gt.mcmc_equilibrate(state, force_niter=n_iter, 
                        mcmc_args=dict(niter=10),
                        callback=collect_partitions)
                    
    # Disambiguate partitions and obtain marginals
    pmode = gt.PartitionModeState(bs, converge=True, nested=nested)
    pv = pmode.get_marginal(state.g)
    
    if nested:
        n_groups = state.get_levels()[0].get_nonempty_B()
    else:   
        n_groups = state.get_nonempty_B()
    pv_array = pv.get_2d_array(np.arange(n_groups)) / (n_iter - 1)

    if nested and level > 0:
        p0 = pd.Categorical(state.project_partition(0, 0))
        pL = pd.Categorical(state.project_partition(level, 0))
        ct = pd.crosstab(p0, pL, normalize='index')
        pv_array = (pv_array.T @ ct.values)
        pv_array = pv_array/np.sum(pv_array, axis=1)[:, None]
        pv_array = pv_array.T

    
    adata.obs[f'{key_added}'] = np.max(pv_array, axis=0)
    adata.obsm[f"CM_{key_added}"] = pv_array.T
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


