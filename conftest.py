import numpy as np
import anndata as ad
import schist as scs
import scanpy as sc            
from mudata import MuData           
import pandas as pd   
sc.settings.verbosity=2                           
adata = sc.datasets.blobs() 
sc.tl.pca(adata)                                                        


print("Test standard hierarchical model")
sc.pp.neighbors(adata,  key_added='test', use_rep='X_pca')
scs.inference.fit_model(adata, nested=True , neighbors_key='test', n_samples=100)
print("Test drawing")
adata.obs['rn'] = np.random.rand(adata.shape[0])
scs.pl.draw_tree(adata, color='rn', save='test.png')
adata.write_h5ad('test.h5ad')
lineages = scs.tl.cr_lineages(adata, level=0)
print("Test cell similarity")
scs.tools.cell_similarity(adata, neighbors_key='test')
sc.pp.neighbors(adata, n_neighbors=3)
print("Test assortative model")
scs.inference.fit_model(adata, assortative=True, use_weights=False, save_model='test', n_samples=100)
lineages = scs.tl.cr_lineages(adata, model_key='ppbm', level=0)
print("Test flat model")
scs.inference.fit_model(adata, nested=False, use_weights=False, save_model='test', n_samples=100)
#scs.inference.nested_model(adata, save_model='test',  n_samples=2, refine_iter=2)
print("Test cell affinity calculation")
scs.tools.calculate_affinity(adata, neighbors_key='test', model_key='sbm')#group_by='sbm')
adata.uns.pop('schist')
scs.tools.calculate_affinity(adata, neighbors_key='test', group_by='sbm')
print("Testing constrained modeling")
adata.obs['batch'] = pd.Categorical(np.random.choice([0, 1], adata.shape[0]))
scs.inference.fit_model(adata, nested=False, constraint_key='batch')

#SKIP THIS UNTIL adata 1.13 is released...
#test label transfer
#d1 = sc.datasets.blobs()
#d2 = sc.datasets.blobs()
#sc.pp.neighbors(d1)
#d1.obs['cluster'] = pd.Categorical(np.random.randint(low=0, high=3, size=d1.shape[0]).astype(str))
#adata = ad.concat([d1, d2], join='outer', index_unique='-', keys=['0','1'])
#sc.pp.neighbors(d2)
#
#scs.tools.label_transfer(d2, d1, obs='cluster')
#
#sc.pp.neighbors(adata)
#adata.obs['cluster'] = adata.obs['cluster'].cat.add_categories('unknown').fillna('unknown')
#
#scs.tools.label_transfer(adata, obs='cluster')

print("Test multimodal model")
d1 = sc.datasets.blobs(n_observations=200)
d2 = sc.datasets.blobs(n_observations=100)
sc.pp.neighbors(d1)
sc.pp.neighbors(d2)

scs.inference.fit_model_multi([d1, d2], nested=True, n_samples=100)
u = scs._utils.get_multi_graph_from_adata([d1, d2])
print(u.num_vertices())
mdata = MuData({'A':d1, 'B':d2})
scs.inference.fit_model_multi(mdata, nested=False, n_samples=100)

