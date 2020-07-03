# helpers functions that will be deprecated as they will included in graph-tool


def pp_virtual_vertex_move(state, v, s):
    """
    A a function that calculates the difference in entropy when a cell
    is moved in another cluster of a PPBlockState
    
    """
    
    blocks = np.array(state.get_blocks().get_array())
    # IDK what's going on here, but apparently entropy of state and of 
    # reinitialied state is different, take this one
    E0 = gt.PPBlockState(g, b=blocks).entropy() 
    blocks[v] = s 
    E1 = gt.PPBlockState(g, b=blocks).entropy()  
    return E1 - E0

def check_gt_version(min_v=2.33):
    # this should be written in a more general way
    # but gt versioning is M.mm only, so...
    import graph_tool.all as gt
    raw_v = gt.__version__[:4]
    if float(raw_v) < min_v:
        raise VersionError(f"""
        You should install graph-tool version {min_v} to use this Feature.
        You have version {raw_v}
        """)
    