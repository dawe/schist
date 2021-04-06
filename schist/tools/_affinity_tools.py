from typing import Optional#, Tuple, Sequence, Type, Union, Dict

import numpy as np
from anndata import AnnData
import scipy.stats

from scanpy import logging as logg
import graph_tool.all as gt


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
    key: Optional[str] = 'nsbm', # dummy default
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

    obsm_names = [x for x in adata.obsm_keys() if x.startswith(f'CA_{key}_level')]
    if len(obsm_names) == 0:
        raise KeyError(
            f"Your dataset does not contain cell affinities, did you run nSBM?"
        )

    _S = np.array([scipy.stats.entropy(adata.obsm[x], axis=1) /np.log(adata.obsm[x].shape[1]) for x in obsm_names]).T
    adata.obs['cell_stability'] = 1-np.nanmax(_S, axis=1) #/ np.nanmean(EE, axis=1)

    return adata if copy else None

def max_marginal(
    adata: AnnData,
    key: Optional[str] = 'nsbm', # dummy default
    key_added: Optional[str] = 'max_marginal',
    copy: bool = False,
    n_iter: int = 100,
    state: Optional = None
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
    adata.obs[f'{key_added}'] = np.max(pv_array, axis=0)
    return adata if copy else None


def cell_similarity(
    adata: AnnData,
    key_added: Optional[str] = 'cell_similarity',
    sim_type: Optional[str] = 'salton',
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


