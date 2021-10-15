from typing import Optional, Tuple, Sequence, Type, Union, Dict

import numpy as np
from anndata import AnnData
from scanpy import logging as logg
import pickle
from scipy import sparse
from sklearn.metrics import adjusted_mutual_info_score as ami
import pandas as pd

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

def plug_state(adata: AnnData,
    state: Union[gt.NestedBlockState, 
                gt.BlockState, 
                gt.PPBlockState],
    nested: bool = True,
    key_added: str = 'nsbm',
    calculate_affinity: bool = False,
    use_weights=False,
    copy: bool = False

) -> Optional[AnnData]:
    """\
    Add a state to a dataset, populate the AnnData.obs consistenly 

    Parameters
    ----------
    adata
        The annotated data matrix.
    state
        The graph_tool state. Supported types are NestedBlockState
        BlockState and PPBlockState
    nested
        If False plug only the lowest level, otherwise the full hierarchy
    key_added
        The prefix for annotations

    """

    adata = adata.copy() if copy else adata
    g = state.g
    
    model_type = 'nested'
    if type(state) == gt.PPBlockState:
        model_type = 'planted'
    elif type(state) == gt.BlockState:
        model_type = 'flat'
    
    if type(state) == gt.NestedBlockState:
        bs = state.get_bs()
        if not nested:
            bs = bs[:1]
        groups = np.zeros((g.num_vertices(), len(bs)), dtype=int)
        for x in range(len(bs)):
            groups[:, x] = state.project_partition(x, 0).get_array()
        groups = pd.DataFrame(groups).astype('category')
        for c in groups.columns:
            ncat = len(groups[c].cat.categories)
            new_cat = [u'%s' % x for x in range(ncat)]
            groups[c].cat.rename_categories(new_cat, inplace=True)
        levels = groups.columns
        groups.columns = [f"{key_added}_level_{level}" for level in range(len(bs))]
        groups.index = adata.obs_names
        # remove any column with the same key
        keep_columns = [x for x in adata.obs.columns if not x.startswith('%s_level_' % key_added)]
        adata.obs = adata.obs[keep_columns]
        adata.obs = pd.concat([adata.obs, groups], axis=1)
        
        if not 'schist' in adata.uns:
            adata.uns['schist'] = {}

        adata.uns['schist'][f'{key_added}'] = {}
        adata.uns['schist'][f'{key_added}']['stats'] = dict(
        level_entropy=np.array([state.level_entropy(x) for x in range(len(bs))]),
        modularity=np.array([gt.modularity(g, state.project_partition(x, 0))
                         for x in range(len(bs))])
        )
        adata.uns['schist'][f'{key_added}']['blocks'] = [np.array(l.get_blocks().a) for l in state.get_levels()]
        adata.uns['schist'][f'{key_added}']['params'] = dict(
            model='nested',
            use_weights=use_weights,
            key_added=key_added,
            collect_marginals=collect_marginals,
            deg_corr=state.deg_corr,
            recs=state.rec,
            rec_types=state.rec_types
        )
        
        if calculate_affinity:
            p0 = get_cell_loglikelihood(state, level=0, as_prob=True)
            adata.obsm[f'CA_{key_added}_level_0'] = p0
            l0 = "%s_level_0" % key_added
            for nl, level in enumerate(groups.columns[1:]):
                cross_tab = pd.crosstab(groups[l0], groups[level])
                cl = np.zeros((p0.shape[0], cross_tab.shape[1]), dtype=p0.dtype)
                for x in range(cl.shape[1]):
                    # sum counts of level_0 groups corresponding to
                    # this group at current level
                    cl[:, x] = p0[:, np.where(cross_tab.iloc[:, x] > 0)[0]].sum(axis=1)
                adata.obsm[f'CA_{key_added}_level_{nl + 1}'] = cl / np.sum(cl, axis=1)[:, None]

    else:
        
        groups = pd.Series(state.get_blocks().get_array()).astype('category')
        ncat = len(groups.cat.categories)
        new_cat = [u'%s' % x for x in range(ncat)]
        groups.cat.rename_categories(new_cat, inplace=True)
        groups.index = adata.obs_names
        adata.obs[key_added] = groups
        if not 'schist' in adata.uns:
            adata.uns['schist'] = {}
        adata.uns['schist'][f'{key_added}'] = {}
        adata.uns['schist'][f'{key_added}']['stats'] = dict(
            modularity=gt.modularity(g, state.get_blocks())
        )
        adata.uns['schist'][f'{key_added}']['blocks'] = np.array(state.get_blocks().a)
        if calculate_affinity:
            adata.obsm[f'CA_{key_added}_level_1'] = get_cell_loglikelihood(state, as_prob=True)
            
        adata.uns['schist'][f'{key_added}']['params'] = dict(
            model=model_type,
            use_weights=use_weights
            calculate_affinity=calculate_affinity,
            deg_corr=state.deg_corr,
            recs=state.rec,
            rec_types=state.rec_types
        )


    return adata if copy else None