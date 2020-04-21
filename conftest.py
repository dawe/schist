import scnsbm
import scanpy as sc                          
sc.settings.verbosity=2                           
adata = sc.datasets.blobs() 
sc.tl.pca(adata)                                                        
sc.pp.neighbors(adata, n_neighbors=3, key_added='foo')
scnsbm.inference.nested_model(adata, fast_model=True, use_weights=False,  wait=10, equilibrate_args=dict(force_niter=2,verbose=True) , neighbors_key='foo')          
scnsbm.io.write(adata, prefix='test')
adata = scnsbm.io.read(prefix='test')
sc.pp.neighbors(adata, n_neighbors=3)
scnsbm.inference.flat_model(adata, fast_model=True, use_weights=False, wait=10, collect_marginals=True)          
scnsbm.inference.nested_model(adata, resume=True, wait=10)
scnsbm.inference.flat_model(adata, resume=True, wait=10)
scnsbm.inference.nested_model(adata, n_init=2, equilibrate=False) 
