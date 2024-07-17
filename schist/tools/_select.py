from typing import Optional

import numpy as np
from anndata import AnnData
import pandas as pd
from scanpy import logging as logg
from scanpy._compat import Literal

import graph_tool.all as gt

def select_affinity(
    adata: AnnData,
    level: str = '1',
    threshold: float = 0.9999,
    inverse: bool = False,
    key: Optional[str] = 'nsbm',
    update_state: Optional[bool] = False,
    filter: Optional[bool] = True,
    copy: bool = False
):
    """\
    Selects cells based on the affinity values at a specified level.
    
    Parameters
    ----------
    adata
        Annotated data matrix. A NestedBlockState object needs to be saved
    level
        The level to be used for selection
    threshold
        The maximal affinity to be used. Cells with affinities lower than the
        threshold will be discarded
    inverse
        Whether to return cells with affinity lower than the threshold
    key
        key of the groupings used to evaluate the model
    update_state
        Whether to update the state removing unselected cells
    filter
        If False, cells are not filtered and only marked in `adata.obs['selected']`
    copy
        Whether to perform selection in place or return a subsetted object

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with selected cells.

    """

    #this function is needed as subsetting doesn't work on marginals
    raise ValueError(
                'This function has been deprecated'
            )
    
    adata = adata.copy() if copy else adata

    level = str(level) # ensure it is a string
    
    if f'CA_{key}_level_{level}' not in adata.obsm_keys():
        logg.error(f'Level {level} was not found in your data')
        raise
    
    affinities = adata.obsm[f'CA_{key}_level_{level}']
    max_aff = np.max(affinities, axis=1)
    if inverse:
        mask = max_aff < threshold
    else:
        mask = max_aff >= threshold
    
    adata.obs['selected'] = mask#pd.Categorical(mask)
    
    if filter:
        adata = adata[adata.obs['selected']] #actually filter cells
    
        if update_state and adata.uns['schist'][f'{key}']:
            logg.warning('Removing a vertex from a BlockState may result in inconsistent data')
            v_idx = np.where(np.bitwise_not(mask)) #vertex to be removed
            adata.uns['schist'][key]['state'].remove_vertex(v_idx)
    
    return adata if copy else None
