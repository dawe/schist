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
        key of the `adata.uns` slot storing the affinities
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
    
    adata = adata.copy() if copy else adata

    if level not in adata.uns[key]['cell_affinity']:
        logg.error(f'Level {level} was not found in your data')
        raise
    
    affinities = adata.uns[key]['cell_affinity'][level]
    max_aff = np.max(affinities, axis=1)
    if inverse:
        mask = max_aff < threshold
    else:
        mask = max_aff >= threshold
    
    adata.obs['selected'] = pd.Categorical(mask)
    
    if filter:
        for l in adata.uns[key]['cell_affinity'].keys():
            # filter affinities
            adata.uns[key]['cell_affinity'][l] = adata.uns[key]['cell_affinity'][l][mask]
    
        adata = adata[mask] #actually filter cells
    
        if update_state and adata.uns[key]['state']:
            logg.warning('Removing a vertex from a BlockState may result in inconsistent data')
            v_idx = np.where(np.bitwise_not(mask)) #vertex to be removed
            adata.uns[key]['state'].remove_vertex(v_idx)
    
    return adata if copy else None
