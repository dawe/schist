
=========
Tutorials
=========

.. highlight:: python

In ored to perform cluster analysis on single cell data, count matrices must be previously prepared for the analysis. For a brief and effective example of single cell RNA sequencing (scRNA-seq) data processing, check the `Scanpy tutorial of 3k PBMCs <https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html>`_.

The first step necessary to perform cluste analyses with ``schist`` is to import the library::
    
    import schist as scs

After that, standard analysis with ``scanpy`` can be performed::
    
    import scanpy as sc
    
    adata=sc.read(adata = sc.read_10x_mtx('data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols', cache=True)  
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    mito_genes = adata.var_names.str.startswith('MT-') 
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata = adata[adata.obs['percent_mito'] < 0.05, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)

----------------------
Clustering of 3k PBMCs
----------------------

In this section, the "3k PBMCs" dataset is preprocessed, using ``scanpy``. Afterwards, cluster analysis is performed, using ``schist``. 

Before the cluster analysis with ``schist``, the *k*NN graph must be built, representing connections between cells of the dataset. The *k*NN graph can be created, using the scanpy function ``sc.pp.neighbors``, which allows the selection of the number of neighbors and the number of pca components considered for the analysis::

    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)

Moreover, in order to spatially visualize the outcome of cluster analysis, we can compute the UMAP embedding of our dataset, using the function ``sc.tl.umap``::
   
    sc.tl.umap(adata)

nested_model
^^^^^^^^^^^^

spiego nested


codice nested
outcome nested


alluvial




planted_model
^^^^^^^^^^^^^

spiego planted

codice planted

esito planted


--------------
Label transfer
--------------

spiego label transfer (problema)
spiego come lo faccio (cell aff)
spiego cosa mi serve (harmony)

mostro codice

mostro outcome con dati buoni (quartzseq)

mostro outcome con dati non buoni (iCEll8)
