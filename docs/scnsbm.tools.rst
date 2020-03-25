********************
scnsbm.tools package
********************

scnsbm.tools._draw_graph
########################


::

    scnsbm.tools._draw_graph(adata: AnnData, layout: 'sfdp', use_tree: False, random_seed: None, adjacency: None, key_added_ext: None, key: 'nsbm', copy: False, **kwds)

Extends scanpy.tools.draw_graph function using some layouts available in graph-tool library. 

Three layouts are available here:
    *   SFDP spring-block layout.
    *   ARF spring-block layout.
    *   Fruchterman-Reingold spring-block layout.

Fruchterman-Reingold is already available in scanpy, but here can be used to render the nested model tree. 
    
In order to use these plotting function, the NestedBlockState needs to be saved when building the model, so `save_state=True` needs to be set.

Parameters:
***********
    *   **adata**: ``AnnData``
            Annotated data matrix. A NestedBlockState object needs to be saved.
    *   **layout**: _Layout (default: ``'sfdp'`` )
            A layout among 'sfdp', 'fr' or 'arf'. Other graph-tool layouts haven't been implemented.
    *   **use_tree**: ``bool`` (default: ``False`` )
            When this is set, the tree of the nested model is used to generate layout, otherwise the layout only accounts for the neighborhood graph.    
    *   **random_seed**: Optional [ ``int`` ] (default: ``None`` )
            Random number to be used as seed for graph-tool.
    *   **adjacency**: Optional [ ``spmatrix`` ] (default: ``None`` )
            Sparse adjacency matrix of the graph, defaults to ``adata.uns['neighbors']['connectivities']`` .
    *   **key_added_ext**: Optional [ ``str`` ] (default: ``None`` )
            By default, append ``layout`` .
    *   **key**: Optional [ ``str`` ] (default: ``'nsbm'`` )
            The slot in ``AnnData.uns`` containing the state. Default is 'nsbm'
    *   **copy**: ``bool`` (default: ``False`` )
            Return a copy instead of writing to adata.
    *   ****kwds**
            Parameters of chosen igraph layout. See e.g. ``fruchterman-reingold`` [Fruchterman91]_. One of the most important ones is ``maxiter`` .

.. _fruchterman-reingold: http://igraph.org/python/doc/igraph.Graph-class.html#layout_fruchterman_reingold

Returns:
********
    *   Depending on ``copy`` , returns or updates ``adata`` with the following field.
    *   **X_draw_graph_layout**: ``adata.obsm``
            Coordinates of graph layout. E.g. for layout='fa' (the default), the field is called 'X_draw_graph_fa'.







.. automodule:: scnsbm.tools
   :members:
   :undoc-members:
   :show-inheritance:
