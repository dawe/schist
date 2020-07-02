import graph_tool.all as gt
import numpy as np

def get_cell_loglikelihood(
    state,
    level: int = 0,
    rescale: bool = False, 
    as_prob: bool = False
    
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
    g = state.g
    
    if isinstance(state, gt.NestedBlockState):
        if level < 0 or level > len(state.get_levels()):
            # by now return the lowest level if invalid 
            level = 0
        B = gt.BlockState(g, b=state.project_partition(level, 0))
    elif isinstance(state, gt.PPBlockState):
        # PP doesn't support virtual moves, so we transform it into
        # a blockstate 
        B = gt.BlockState(g, b=state.get_blocks())
    else:
        B = state
    
    
    n_cells = g.num_vertices()
    n_blocks = B.get_nonempty_B()
    
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