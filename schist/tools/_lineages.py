from typing import Optional#, Tuple, Sequence, Type, Union, Dict

import numpy as np
from anndata import AnnData
from scanpy import logging as logg
#import pandas as pd

def cr_lineages(
    adata: AnnData,
    model_key: Optional[str] = 'nsbm',    
    level: int = 1,
    use_affinity: bool = False,
    
) -> Optional[AnnData]:
    
    """\
    Return marginals/affinities as cellrank Lineage.
    
    Parameters
    ----------
    adata:
        The AnnData object. Should have been already processed with schist
    model_key:
        The partition model used. It matches the "key_added" parameter of
        `scs.inference.model()`. If no key was provided when modeling, it matches the
        name of the model (i.e. 'nsbm', 'sbm' or 'ppbm')
    level:
        The level to consider. This parameter is effective
        only for Nested partitions
    use_weights
        Default mode uses cell marginals (adata.obsm['CM_{model_key}']) to work. 
        Set this to True if affinities should be uses instead.

    Returns
    -------
    A cellrank Lineage class instance, where each block model cluster corresponds
    to a single lineage. 
        
"""    
    try:
        from cellrank._utils._lineage import Lineage
    except ModuleNotFoundError:
        from cellrank.tl import Lineage
    except ImportError:
        raise ImportError(
        """Please install cellrank first. You can do it by pip:
        
        `pip install cellrank`

        or by conda: 
        
        `conda install -c conda-forge cellrank`
        
        Refer to

        https://cellrank.readthedocs.io/en/latest/installation.html
        """
    )
    
    if 'schist' in adata.uns and model_key in adata.uns['schist']:
        params = adata.uns['schist'][model_key]['params']
    else:
        raise KeyError(
        """
        I can't find the parameters for the specified model key
        """
        )
    key_added = params['key_added']
    model = params['model']
    obs_df = adata.obs.filter(like=key_added)
    if obs_df.shape[1] == 0:
        raise ValueError(
        f"""
        It seems that the {key_added}, specified in parameters,
        can't be found in adata.obs
        """
        )

    prefix = 'CA' if use_affinity else 'CM'
    group_key = f'{key_added}_level_{level}' if model == "nsbm" else key_added
    matrix_key = f'{prefix}_{group_key}'

    if matrix_key not in adata.obsm:
        raise KeyError(
        f"""
        I can't find the specified matrix in adata.obsm ({matrix_key})
        """
        )

    lineage_matrix = adata.obsm[matrix_key]
    lineage_names = adata.obs[group_key].cat.categories
    # For multi_model the number of groups and dimensions in the matrix may
    # not match, this because it may happen that some groups are only given
    # to a specific modality. 
    # We prefer to keep the matrices with the same dimensions as they can be used
    # for other puroposes. Here instead we have to match. We can either subset the 
    # matrix or augment the lineages adding back the missing ones. 
    # Choose the second
    if len(lineage_names) != lineage_matrix.shape[1]:
        # this can be done because the groups have been previously sorted
        lineage_names = [str(x) for x in range(lineage_matrix.shape[1])]

    lineages = Lineage(lineage_matrix,  names=lineage_names)
    
    return lineages    
    

