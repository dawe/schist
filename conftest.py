import scnsbm
import scanpy as sc                          
sc.settings.verbosity=2                           
adata = sc.datasets.blobs() 
sc.tl.pca(adata)                                                        
sc.pp.neighbors(adata, n_neighbors=3)                                   
scnsbm.inference.nested_model(adata, fast_model=True, use_weights=False, sweep_iterations=10, wait=10, save_state=True, equilibrate_args=dict(force_niter=2,verbose=True), sweep_args=dict(d=0.05))          
scnsbm.io.write(adata, prefix='test')
adata = scnsbm.io.read(prefix='test')
scnsbm.inference.flat_model(adata, fast_model=True, use_weights=False, sweep_iterations=10, wait=10, save_state=True)          
