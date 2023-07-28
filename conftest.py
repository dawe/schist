import schist
import scanpy as sc                          
sc.settings.verbosity=2                           
adata = sc.datasets.blobs() 
sc.tl.pca(adata)                                                        
try:
    sc.pp.neighbors(adata, n_neighbors=3, key_added='foo')
    schist.inference.nested_model(adata,  neighbors_key='foo', n_init=2)
except TypeError:
    sc.pp.neighbors(adata, n_neighbors=3)
    schist.inference.nested_model(adata, n_init=2)

adata.write('foo.h5ad')
schist.inference.leiden(adata, neighbors_key='foo', save_model='lle', n_init=2)
schist.tools.calculate_affinity(adata, neighbors_key='foo', group_by='leiden')
schist.tools.cell_similarity(adata, neighbors_key='foo')
sc.pp.neighbors(adata, n_neighbors=3)
schist.inference.planted_model(adata, use_weights=False, save_model='foo', n_init=2, refine_iter=2)
schist.inference.flat_model(adata, use_weights=False, save_model='foo', n_init=2, refine_iter=2)
#schist.inference.nested_model(adata, save_model='test', dispatch_backend='threads', n_init=2, refine_iter=2)
schist.tools.calculate_affinity(adata, level=0, back_prob=True)

#test label transfer
d1 = sc.datasets.blobs()
d2 = sc.datasets.blobs()
sc.pp.neighbors(d1)
sc.tl.leiden(d1)
adata = ad.concat([d1, d2], join='outer', index_unique='-', keys=['0','1'])
sc.pp.neighbors(d2)

schist.tools.label_transfer(d2, d1, obs='leiden')

sc.pp.neighbors(adata)
adata.obs['leiden'] = adata.obs['leiden'].cat.add_categories('unknown').fillna('unknown')

schist.tools.label_transfer(adata, obs='leiden')
d1 = sc.datasets.blobs(n_observations=200)
d2 = sc.datasets.blobs(n_observations=100)
sc.pp.neighbors(d1)
sc.pp.neighbors(d2)

schist.inference.nested_model_multi([d1, d2], n_init=2, refine_iter=2)
u = schist._utils.get_multi_graph_from_adata([d1, d2])
print(u.num_vertices())
