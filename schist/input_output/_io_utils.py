from typing import Union, Optional

import pickle
from anndata import AnnData, read_h5ad
from scanpy import logging as logg

def read(
    prefix: str = 'adata', 
    key: str = 'schist', 
    h5ad_fname: Optional[str] = None, 
    pkl_fname: Optional[str] = None
    ) -> Optional[AnnData]:
    """\
    Read anndata object when a NestedBlockState has been saved separately.
    This function reads the h5ad and the pkl files, then rebuilds the `adata` properly,
    returning it to the user. Note that if pkl is not found, an AnnData object
    is returned anyway

    Parameters
    ----------
    prefix
        The prefix for .h5ad and .pkl files, it is supposed to be the same for 
        both. If this is not, specify file names (see below)
    key
        The slot in `AnnData.uns` in which nsbm information is placed
    h5ad_filename
        If `prefix` is not shared between h5ad and pkl, specify the h5ad file here
    pkl_filename
        If `prefix` is not shared between h5ad and pkl, specify the pkl file here
"""
    if not h5ad_fname:
        h5ad_fname = "%s.h5ad" % prefix
    if not pkl_fname:
        pkl_fname = "%s.pkl" % prefix
    
    # read the anndata
    adata = read_h5ad(h5ad_fname)
    
    if not key in adata.uns:
        raise KeyError(
            f"Your dataset does not contain {key}, did you run schist?"
        )
    
    try:
        with open(pkl_fname, 'rb') as fh:
            state = pickle.load(fh)
            if type(state) == dict:
                if 'state' in state:
                    adata.uns[key]['state'] = state['state']
                if 'multi_level_state' in state:
                    adata.uns[key]['multi_level_state'] = state['multi_level_state']
            else:
                adata.uns[key]['state'] = state
    except IOError:
        logg.warning(
            f'The specified file for state {pkl_fname} does not exist. '
            'Proceeding anyway'
        )
        pass            
    return adata


def write(
    adata: AnnData, 
    prefix: str = 'adata', 
    key: str = 'schist', 
    h5ad_fname: Optional[str] = None, 
    pkl_fname: Optional[str] = None
    ) -> Optional[AnnData]:
    """Save anndata object when a NestedBlockState has been retained during
    inference. The `state` object is stripped out the `adata.uns` and saved as pickle
    separately.

    Parameters
    ----------
    adata
        The AnnData object to be saved
    prefix
        The prefix for .h5ad and .pkl files. Two files (prefix.h5ad, prefix.pkl) 
        will be saved
    key
        The slot in `AnnData.uns` in which nsbm information is placed
    h5ad_filename
        Specify a file name for AnnData
    pkl_filename
        Specify a file name for `state` pickle
"""
    state = {}
    if 'state' in adata.uns[key]:
        state['state'] = adata.uns[key].pop('state')
    if 'multi_level_state' in adata.uns[key]:
        state['multi_level_state'] = adata.uns[key].pop('multi_level_state')
    if not h5ad_fname:
        h5ad_fname = "%s.h5ad" % prefix
    if not pkl_fname:
        pkl_fname = "%s.pkl" % prefix
    # write the anndata
    adata.write(h5ad_fname)
    if state:
        with open(pkl_fname, 'wb') as fh:
            pickle.dump(state, fh, 2)
    # restore state        
    adata.uns[key]['state'] = state            
