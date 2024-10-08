import numpy as np
import anndata as ad
import schist as scs
import scanpy as sc            
from mudata import MuData           
import pandas as pd   
sc.settings.verbosity=2                           
adata = sc.datasets.blobs() 
sc.tl.pca(adata)                                                        


print("NESTED MODEL")
sc.pp.neighbors(adata,  key_added='foo', use_rep='X_pca')
scs.inference.fit_model(adata,  model='nsbm', neighbors_key='foo', n_init=2)
adata.obs['rn'] = np.random.rand(adata.shape[0])
scs.pl.draw_tree(adata, color='rn', save='test.png')
adata.write('foo.h5ad')
lineages = scs.tl.cr_lineages(adata, level=0)
#comment out as it's broken by numpy
scs.tools.cell_similarity(adata, neighbors_key='foo')
sc.pp.neighbors(adata, n_neighbors=3)
print("PLANTED MODEL")
scs.inference.fit_model(adata, model='ppbm', use_weights=False, save_model='foo', n_init=2, refine_iter=2)
lineages = scs.tl.cr_lineages(adata, model_key='ppbm', level=0)
print("SB MODEL")
scs.inference.fit_model(adata, model='sbm', use_weights=False, save_model='foo', n_init=2, refine_iter=2)
#scs.inference.nested_model(adata, save_model='test',  n_init=2, refine_iter=2)
scs.tools.calculate_affinity(adata, neighbors_key='foo', model_key='sbm')#group_by='sbm')
adata.uns.pop('schist')
scs.tools.calculate_affinity(adata, neighbors_key='foo', group_by='sbm')

#test label transfer
d1 = sc.datasets.blobs()
d2 = sc.datasets.blobs()
sc.pp.neighbors(d1)
d1.obs['cluster'] = pd.Categorical(np.random.randint(low=0, high=3, size=d1.shape[0]).astype(str))
adata = ad.concat([d1, d2], join='outer', index_unique='-', keys=['0','1'])
sc.pp.neighbors(d2)

scs.tools.label_transfer(d2, d1, obs='cluster')

sc.pp.neighbors(adata)
adata.obs['cluster'] = adata.obs['cluster'].cat.add_categories('unknown').fillna('unknown')

scs.tools.label_transfer(adata, obs='cluster')
d1 = sc.datasets.blobs(n_observations=200)
d2 = sc.datasets.blobs(n_observations=100)
sc.pp.neighbors(d1)
sc.pp.neighbors(d2)

scs.inference.fit_model_multi([d1, d2], model='nsbm', n_init=2, refine_iter=2)
u = scs._utils.get_multi_graph_from_adata([d1, d2])
print(u.num_vertices())
mdata = MuData({'A':d1, 'B':d2})
scs.inference.fit_model_multi(mdata, model='sbm', n_init=2, refine_iter=2)

