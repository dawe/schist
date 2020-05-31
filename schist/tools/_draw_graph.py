from typing import Union, Optional

import numpy as np
from anndata import AnnData
from scipy.sparse import spmatrix

from .._utils import get_graph_tool_from_adjacency
from scanpy import logging as logg
from scanpy._compat import Literal

import graph_tool.all as gt

_LAYOUTS = ('fr', 'sfdp', 'arf')
_Layout = Literal[_LAYOUTS]


def draw_graph(
    adata: AnnData,
    layout: _Layout = 'sfdp',
#    init_pos: Union[str, bool, None] = None,
#    root: Optional[int] = None,
    use_tree: bool = False,
    random_seed: Optional[int] = None,
    adjacency: Optional[spmatrix] = None,
    key_added_ext: Optional[str] = None,
    key: Optional[str] = 'nsbm',
    copy: bool = False,
    **kwds,
):
    """\
    Extends scanpy.tools.draw_graph function using some layouts available in 
    graph-tool library. Three layouts are available here:
    
    - SFDP spring-block layout.
    - ARF spring-block layout.
    - Fruchterman-Reingold spring-block layout.
    
    Fruchterman-Reingold is already available in scanpy, but here can be used
    to render the nested model tree. 
    
    In order to use these plotting function, the NestedBlockState needs to be
    saved when building the model, so `save_state=True` needs to be set.
    
    Parameters
    ----------
    adata
        Annotated data matrix. A NestedBlockState object needs to be saved
    layout
        A layout among 'sfdp', 'fr' or 'arf'. Other graph-tool layouts haven't been
        implemented.
    use_tree
        When this is set, the tree of the nested model is used to generate layout, 
        otherwise the layout only accounts for the neighborhood graph.    
    random_seed
        Random number to be used as seed for graph-tool
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']`.
    key_added_ext
        By default, append `layout`.
    key
        The slot in `AnnData.uns` containing the state. Default is 'nsbm'
    copy
        Return a copy instead of writing to adata.
    **kwds
        Parameters of chosen igraph layout. See e.g. `fruchterman-reingold`_
        [Fruchterman91]_. One of the most important ones is `maxiter`.

        .. _fruchterman-reingold: http://igraph.org/python/doc/igraph.Graph-class.html#layout_fruchterman_reingold

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following field.

    **X_draw_graph_layout** : `adata.obsm`
        Coordinates of graph layout. E.g. for layout='fa' (the default),
        the field is called 'X_draw_graph_fa'
    """
    if random_seed:
        np.random.seed(random_seed)
        gt.seed_rng(random_seed)

    n_cells = adata.shape[0]
    start = logg.info(f'drawing single-cell graph using layout {layout!r}')
    if layout not in _LAYOUTS:
        raise ValueError(f'Provide a valid layout, one of {_LAYOUTS}.')
    adata = adata.copy() if copy else adata
    if adjacency is None and 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first '
            'to compute a neighborhood graph.'
        )
    if not key in adata.uns:
        raise ValueError(
            'You need to run `nested_model` before trying to run this function '
        )
        
    if use_tree and 'state' not in adata.uns[key]:    
        raise ValueError(
            'When `use_tree` is set to `True`, a state should be saved'
            'running  `nested_model(adata, save_state=True)`.'
        )
    if adjacency is None:
        adjacency = adata.uns['neighbors']['connectivities']
    
    g = get_graph_tool_from_adjacency(adjacency)
    weights=g.ep['weight']
    if use_tree:
        state = adata.uns[key]['state']
        g, _, _ = gt.get_hierarchy_tree(state, empty_branches=False)
        weights=None
    
    # actual drawing
    positions = np.zeros((n_cells, 2))
    if layout == 'fr':
        positions = gt.fruchterman_reingold_layout(g, weight=weights)
        positions = np.array([x for x in positions][:n_cells])
    elif layout == 'sfdp':
        positions = gt.sfdp_layout(g)
        positions = np.array([x for x in positions][:n_cells])
    elif layout == 'arf':
        positions = gt.arf_layout(g)
        positions = np.array([x for x in positions][:n_cells])

    adata.uns['draw_graph'] = {}
    adata.uns['draw_graph']['params'] = dict(
        layout=layout, random_seed=random_seed
    )
    key_added = f'X_draw_graph_{layout}'
    adata.obsm[key_added] = positions
    logg.info(
        '    finished',
        time=start,
        deep=(
            'added\n'
            f'    {key_added!r}, graph_drawing coordinates (adata.obsm)'
        ),
    )
    return adata if copy else None
