from typing import Optional#, Tuple, Sequence, Type, Union, Dict

import numpy as np
from anndata import AnnData
import scipy.stats

from scanpy import logging as logg


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

    _S = np.array([scipy.stats.entropy(np.log(adata.obsm[x]), axis=1) /np.log(adata.obsm[x].shape[1]) for x in obsm_names]).T
    adata.obs['cell_stability'] = np.nanmax(_S, axis=1) #/ np.nanmean(EE, axis=1)

    return adata if copy else None
