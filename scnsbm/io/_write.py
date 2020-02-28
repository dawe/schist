import pickle

def write(adata, prefix='adata', key='nsbm', h5ad_fname=None, pkl_fname=None):
    """Save anndata object when a NestedBlockState has been retained during
    inference. The `state` object is stripped out the `adata.uns` and saved as pickle
    separately.
    """
    state = None
    if 'state' in adata.uns[key]:
        state = adata.uns[key].pop('state')
    if not h5ad_fname:
        h5ad_fname = "%s.h5ad" % prefix
    if not pkl_fname:
        pkl_fname = "%s.pkl" % prefix
    # write the anndata
    adata.write(h5ad_fname)
    if state:
        with open(pkl_fname, 'wb') as fh:
            pickle.dump(state, fh, 2)