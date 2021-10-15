from typing import Optional, Tuple, Sequence, Type, Union, Dict

import graph_tool
import graph_tool.all as gt
import numpy as np
import numba
from anndata import AnnData
from scipy.sparse import spmatrix

from ._helpers import *


def state_from_blocks(
    adata: AnnData,
    state_key: Optional[str] = 'nsbm',
    neighbors_key: Optional[str] = 'neighbors',
    adjacency: Optional[spmatrix] = None,
    directed: bool = False,
    use_weights: bool = False,
    deg_corr: bool = True,
):
    """
    Returns a gt state object given an AnnData

    Parameters
    ----------
    adata
        The annotated data matrix.
    state_key
        The key under which the state has been saved
    neighbors_key
        The key passed to `sc.pp.neighbors`
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']` in case of scanpy<=1.4.6 or
        `adata.obsp[neighbors_key][connectivity_key]` for scanpy>1.4.6
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges). Note that this
        increases computation times
    deg_corr
        Whether to use degree correction in the minimization step. In many
        real world networks this is the case, although this doesn't seem
        the case for KNN graphs used in scanpy.
        
    """
    blocks = adata.uns['schist'][f'{state_key}']['blocks']
    params = adata.uns['schist'][f'{state_key}']['params']
    if 'deg_corr' in params:
        deg_corr=params['deg_corr']

    recs=[]
    rec_types=[]
    if use_weights:
        # this is not ideal to me, possibly we may need to transform
        # weights. More tests needed.
        recs=[g.ep.weight]
        rec_types=['real-normal']
        
    if 'recs' in params:
        recs=params['recs']
    if 'rec_types' in params:
        rec_types=params['rec_types']
            
    if adjacency is None:
        if neighbors_key not in adata.uns:
            raise ValueError(
                'You need to run `pp.neighbors` first '
                'to compute a neighborhood graph.'
            )
        elif 'connectivities_key' in adata.uns[neighbors_key]:
            # scanpy>1.4.6 has matrix in another slot
            conn_key = adata.uns[neighbors_key]['connectivities_key']
            adjacency = adata.obsp[conn_key]
        else:
            # scanpy<=1.4.6 has sparse matrix here
            adjacency = adata.uns[neighbors_key]['connectivities']

    g = get_igraph_from_adjacency(adjacency, directed=directed)
    g = g.to_graph_tool()
    gt.remove_parallel_edges(g)

    if params['model'] == 'flat':
        state = gt.BlockState(g, b=blocks, 
            state_args=dict(deg_corr=deg_corr,
            recs=recs,
            rec_types=rec_types)
            )
    elif params['model'] == 'ppbm':
        state = gt.PPBlockState(g, b=blocks, 
            state_args=dict(deg_corr=deg_corr,
            recs=recs,
            rec_types=rec_types)
            )
    else:
        state = gt.NestedBlockState(g, bs=blocks, 
            state_args=dict(deg_corr=deg_corr,
            recs=recs,
            rec_types=rec_types)
            )
    return state            


#@numba.jit(forceobj=True, parallel=True)
def get_cell_loglikelihood(
    state: Union[graph_tool.inference.nested_blockmodel.NestedBlockState, graph_tool.inference.planted_partition.PPBlockState],
    level: int = 0,
    rescale: bool = False, 
    as_prob: bool = False,
    
):
    """
    Returns the matrix of log-likelihood differences
    when moving a cell into a different block
    
    Parameters
    ----------
    state
        A graphtool BlockState or NestedBlockState objexct
    level
        The level in NestedBlockState to consider
    rescale
        For some models, moving a cell into a different block may result in a 
        negative log-likelihood, indicating that cells may be better assigned 
        to another group. Set this parameter to `True` if you want 
        every cell to have LL=0 for the best group and avoid negative values
    as_prob
        Return values as probabilites

    Returns
    -------
    `M`
        Array of dim (n_cells, n_blocks) that stores the entropy difference
        of moving a cell into a specific group
    """
    
    # get the graph from state
#    g = state.g
    
    try:
        if level < 0 or level > len(state.get_levels()):
            # by now return the lowest level if invalid 
            level = 0
        B = gt.BlockState(state.g, b=state.project_partition(level, 0))
    except AttributeError:
        B = state
    
    
    n_cells = state.g.num_vertices()
    n_blocks = B.get_nonempty_B()
#    M = np.zeros((n_cells, n_blocks))
#    for v in range(n_cells): #one day this will be parallel
#        for s in range(n_blocks):
#            M[v, s] = B.virtual_vertex_move(v, s)
    shape = (n_cells, n_blocks)
    M = np.array([B.virtual_vertex_move(v, s) for v in range(n_cells) for s in range(n_blocks)]).reshape(shape)

   
    if rescale:
        # some cells may be better in other groups, hence their LL
        # is negative when moved. Rescaling sets the minimum LL in the
        # best group
        M = M - np.min(M, axis=1)[:, None]
        
    if as_prob:
        E = np.exp(-M)
        return (E / np.sum(E, axis=1)[:, None])

            
    return M

def get_cell_back_p(
    state: Union[graph_tool.inference.nested_blockmodel.NestedBlockState, graph_tool.inference.planted_partition.PPBlockState],
    level: int = 0,
):
    """
    Returns the matrix of proabilities of moving a cell back to its
    group from a different block
    
    Parameters
    ----------
    state
        A graphtool BlockState or NestedBlockState objexct
    level
        The level in NestedBlockState to consider

    Returns
    -------
    `M`
        Array of dim (n_cells, n_blocks) that stores the entropy difference
        of moving a cell into a specific group
    """
    
    # get the graph from state
#    g = state.g
    
    try:
        if level < 0 or level > len(state.get_levels()):
            # by now return the lowest level if invalid 
            level = 0
        B = gt.BlockState(state.g, b=state.project_partition(level, 0))
    except AttributeError:
        B = state
    
    
    n_cells = state.g.num_vertices()
    n_blocks = B.get_nonempty_B()
    shape = (n_cells, n_blocks)
    M = np.array([B.get_move_prob(v, s, reverse=True) for v in range(n_cells) for s in range(n_blocks)]).reshape(shape)

    return M    
