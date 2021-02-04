from typing import Optional, Tuple, Sequence, Type, Union, Dict

import numpy as np
from anndata import AnnData
from scanpy import logging as logg
import pickle
from scipy import sparse
from sklearn.metrics import adjusted_mutual_info_score as ami

try:
    import graph_tool.all as gt
except ImportError:
    raise ImportError(
        """Please install the graph-tool library either visiting

        https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

        or by conda: `conda install -c conda-forge graph-tool`
        """
    )

def get_graph_tool_from_adjacency(adjacency, directed=None):
    """Get graph-tool graph from adjacency matrix."""
    idx = np.nonzero(np.triu(adjacency.todense(),1))
    weights = adjacency[idx]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = gt.Graph(directed=directed)
    g.add_edge_list(np.transpose(idx))  # add
    try:
        ew = g.new_edge_property("double")
        ew.a = weights
        g.ep['weight'] = ew
    except:
        pass
    if g.num_vertices() != adjacency.shape[0]:
        logg.warning(
            f'The constructed graph has only {g.num_vertices()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g


def prune_groups(groups, inverse=False):
    """
    Returns the index of informative levels after the nested_model has
    been run. It works by looking at level entropy and, moreover, checks if
    two consecutive levels have the same clustering
    """
    
    n_groups = groups.shape[1]
    
    mi_groups = np.array([ami(groups.iloc[:, x - 1], groups.iloc[:, x]) for x in range(1, n_groups)])
    
    if inverse:
        return groups.columns[np.where(mi_groups != 1)]
    
    return groups.columns[np.where(mi_groups == 1)]

def get_graph_tool_from_adata(adata: AnnData,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    random_seed: Optional[int] = None,
    key_added: str = 'nsbm',
    adjacency: Optional[sparse.spmatrix] = None,
    neighbors_key: Optional[str] = 'neighbors',
    directed: bool = False,
    use_weights: bool = False,

):
    """Get graph-tool graph from adata."""
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
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories,
            adjacency,
        )
    # convert it to igraph
    g = get_graph_tool_from_adjacency(adjacency, directed=directed)
    return g
