# helpers functions that will be deprecated as they will included in graph-tool


def pp_virtual_vertex_move(state, v, s):
    """
    A a function that calculates the difference in entropy when a cell
    is moved in another cluster of a PPBlockState
    
    """
    
    blocks = state.get_blocks()
    E0 = state.entropy()
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
    