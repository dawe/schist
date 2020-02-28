import pickle
import anndata
from scanpy import logging as logg

def read(prefix='adata', key='nsbm', h5ad_fname=None, pkl_fname=None):
    """Read anndata object when a NestedBlockState has been saved separately.
    This function reads the h5ad and the pkl files, then rebuilds the `adata` properly,
    returning it to the user.
    """
    if not h5ad_fname:
        h5ad_fname = "%s.h5ad" % prefix
    if not pkl_fname:
        pkl_fname = "%s.pkl" % prefix
    
    # read the anndata
    adata = anndata.read_h5ad(h5ad_fname)
    
    try:
        with open(pkl_fname, 'rb') as fh:
            state = pickle.load(fh)
            adata.uns[key]['state'] = state
    except IOError:
        logg.warning(
            f'The specified file for state {pkl_fname} does not exist. '
            'Proceeding anyway'
        )
        pass            
    return adata        