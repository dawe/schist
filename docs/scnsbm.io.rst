*****************
scnsbm.io package
*****************

scnsbm.io.read
##############

::

    scnsbm.io.read(prefix: 'adata', key: 'nsbm', h5ad_fname: None, pkl_fname: None)

Read anndata object when a NestedBlockState has been saved separately. 

This function reads the h5ad and the pkl files, then rebuilds the ``adata`` properly, returning it to the user. Note that if pkl is not found, an AnnData object is returned anyway.

Parameters
**********
    *   **prefix**: ``str`` (default: ``'adata'`` )
            The prefix for .h5ad and .pkl files, it is supposed to be the same for both. If this is not, specify file names (see below).
    *   **key**: ``str`` (default: ``'nsbm'`` )
            The slot in ``AnnData.uns`` in which nsbm information is placed.
    *   **h5ad_filename**: Optional [ ``str`` ] (default: ``None`` )
            If ``prefix`` is not shared between h5ad and pkl, specify the h5ad file here.
    *   **pkl_filename**: Optional [ ``str`` ] (default: ``None`` )
            If ``prefix`` is not shared between h5ad and pkl, specify the pkl file here.

scnsbm.io.write
###############

::

    scnsbm.io.write(adata: AnnData, prefix: 'adata', key: 'nsbm', h5ad_fname: None, pkl_fname: None)

Save anndata object when a NestedBlockState has been retained during inference. 

The ``state`` object is stripped out the ``adata.uns`` and saved as pickle separately.

Parameters
**********
    *   **adata**: ``AnnData``
            The AnnData object to be saved.
    *   **prefix**: ``str`` (default: ``'adata'`` )
            The prefix for .h5ad and .pkl files. Two files (prefix.h5ad, prefix.pkl) will be saved.
    *   **key**: ``str`` (default: ``'nsbm'`` )
            The slot in ``AnnData.uns`` in which nsbm information is placed.
    *   **h5ad_filename**: Optional [ ``str`` ] (default: ``None`` )
            Specify a file name for AnnData.
    *   **pkl_filename**: Optional [ ``str`` ] (default: ``None`` )
            Specify a file name for ``state`` pickle.

.. automodule:: scnsbm.io
   :members:
   :undoc-members:
   :show-inheritance:
