# helpers functions that will be deprecated as they will included in graph-tool
import numpy as np
import graph_tool.all as gt


def check_gt_version(min_v=2.37):
    # this should be written in a more general way
    # but gt versioning is M.mm only, so...
    import graph_tool.all as gt
    raw_v = gt.__version__[:4]
    if float(raw_v) < min_v:
        raise ImportError(f"""
        You should install graph-tool version {min_v} to use this Feature.
        You have version {raw_v}
        """)
    
