#this is from scanpy _leiden.py code, slightly modified

from typing import Optional, Tuple, Sequence, Type, Union

import numpy as np
import pandas as pd
from natsort import natsorted
from anndata import AnnData
from scipy import sparse

from scanpy import _utils
from scanpy import logging as logg
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency

from scanpy._utils import get_igraph_from_adjacency, _choose_graph
from joblib import delayed, Parallel, parallel_config

try:
    from leidenalg.VertexPartition import MutableVertexPartition
except ImportError:

    class MutableVertexPartition:
        pass

    MutableVertexPartition.__module__ = 'leidenalg.VertexPartition'

try:
    import graph_tool.all as gt
except ImportError:
    raise ImportError(
        """Please install the graph-tool>=2.33 library either visiting

        https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

        or by conda: `conda install -c conda-forge graph-tool`
        """
    )

def leiden(
    adata: AnnData,
    resolution: float = 1,
    n_init: int = 100,
    *,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    random_state: _utils.AnyRandom = 0,
    key_added: str = 'leiden',
    adjacency: Optional[sparse.spmatrix] = None,
    directed: bool = True,
    use_weights: bool = True,
    n_iterations: int = -1,
    partition_type: Optional[Type[MutableVertexPartition]] = None,
    neighbors_key: Optional[str] = None,
    obsp: Optional[str] = None,
    collect_marginals: bool = True,
    n_jobs: int = -1,
    copy: bool = False,
    save_model: Union[str, None] = None,
    dispatch_backend: Optional[str] = 'threads',
    **partition_kwargs,
) -> Optional[AnnData]:
    """\
    Cluster cells into subgroups [Traag18]_.

    Cluster cells using the Leiden algorithm [Traag18]_,
    an improved version of the Louvain algorithm [Blondel08]_.
    It has been proposed for single-cell analysis by [Levine15]_.

    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.


    Parameters
    ----------
    adata
        The annotated data matrix.
    resolution
        A parameter value controlling the coarseness of the clustering.
        Higher values lead to more clusters.
        Set to `None` if overriding `partition_type`
        to one that doesn’t accept a `resolution_parameter`.
    n_init
	The number of random initializations to take for consensus        
    random_state
        Change the initialization of the optimization.
    restrict_to
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain `(obs_key, list_of_categories)`.
    key_added
        `adata.obs` key under which to add the cluster labels.
    adjacency
        Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    n_iterations
        How many iterations of the Leiden clustering algorithm to perform.
        Positive values above 2 define the total number of iterations to perform,
        -1 has the algorithm run until it reaches its optimal clustering.
    partition_type
        Type of partition to use.
        Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
        For the available options, consult the documentation for
        :func:`~leidenalg.find_partition`.
    neighbors_key
        Use neighbors connectivities as adjacency.
        If not specified, leiden looks .obsp['connectivities'] for connectivities
        (default storage place for pp.neighbors).
        If specified, leiden looks
        .obsp[.uns[neighbors_key]['connectivities_key']] for connectivities.
    obsp
        Use .obsp[obsp] as adjacency. You can't specify both
        `obsp` and `neighbors_key` at the same time.
    collect_marginals
    	Wheter to retrieve the marginal probability to belong to a group
    n_jobs
        Number of parallel jobs to calculate partitions
    copy
        Whether to copy `adata` or modify it inplace.
    save_model
        If provided, this will be the filename for the PartitionModeState to 
        be saved    
    **partition_kwargs
        Any further arguments to pass to `~leidenalg.find_partition`
        (which in turn passes arguments to the `partition_type`).


    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of cells) that stores the subgroup id
        (`'0'`, `'1'`, ...) for each cell.
    `adata.uns['leiden']['params']`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    """
    try:
        import leidenalg
    except ImportError:
        raise ImportError(
            'Please install the leiden algorithm: `conda install -c conda-forge leidenalg` or `pip3 install leidenalg`.'
        )
    partition_kwargs = dict(partition_kwargs)

    # the following lines are for compatibility
    if dispatch_backend == 'threads':
        dispatch_backend = 'threading'
    elif dispatch_backend == 'processes':
        dispatch_backend = 'loky'

    if dispatch_backend == 'threading' and float(gt.__version__.split()[0]) > 2.55:
        logg.warning('Threading backend does not work with this version of graph-tool\n'
                     'Switching to loky backend')
        dispatch_backend = 'loky'

    start = logg.info('running Leiden clustering')
    adata = adata.copy() if copy else adata
    # are we clustering a user-provided graph or the default AnnData one?
    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)
    if restrict_to is not None:
        restrict_key, restrict_categories = restrict_to
        adjacency, restrict_indices = restrict_adjacency(
            adata,
            restrict_key,
            restrict_categories,
            adjacency,
        )
    # convert it to igraph
    g = get_igraph_from_adjacency(adjacency, directed=directed)
    g_gt = g.to_graph_tool()
    gt.remove_parallel_edges(g_gt)
    # flip to the default partition type if not overriden by the user
    if partition_type is None:
        partition_type = leidenalg.RBConfigurationVertexPartition
    # Prepare find_partition arguments as a dictionary,
    # appending to whatever the user provided. It needs to be this way
    # as this allows for the accounting of a None resolution
    # (in the case of a partition variant that doesn't take it on input)
    if use_weights:
        partition_kwargs['weights'] = np.array(g.es['weight']).astype(np.float64)
    partition_kwargs['n_iterations'] = n_iterations
    np.random.seed(random_state)
    seeds = np.random.choice(range(0, n_init**2), size=n_init, replace=False)
    

    if resolution is not None:
        partition_kwargs['resolution_parameter'] = resolution
    # clustering proper
    def membership(g, partition_type, seed, **partition_kwargs):
        return leidenalg.find_partition(g, partition_type, 
                                        seed=seed, **partition_kwargs).membership

    with parallel_config(backend=dispatch_backend,
                         max_nbytes=None,
                         n_jobs=n_jobs):
        parts = Parallel()(
            delayed(membership)(g, partition_type, seeds[x], **partition_kwargs) for x in range(n_init)
        )


    pmode = gt.PartitionModeState(parts, converge=True) 

    if save_model:
        import pickle
        fname = save_model
        if not fname.endswith('pkl'):
            fname = f'{fname}.pkl'
        logg.info(f'Saving model into {fname}')    
        with open(fname, 'wb') as fout:
            pickle.dump(pmode, fout, 2)

    groups = np.array(pmode.get_max(g_gt).get_array())     
    u_groups = np.unique(groups)
    n_groups = len(u_groups)
    last_group = np.max(u_groups) + 1
    if collect_marginals:
        pv_array = pmode.get_marginal(g_gt).get_2d_array(range(last_group)).T[:, u_groups] / n_init
    # rename groups to ensure they are a continuous range
    rosetta = dict(zip(u_groups, range(len(u_groups))))
    groups = np.array([rosetta[x] for x in groups])

    # store output into adata.obs
        
    if restrict_to is not None:
        if key_added == 'leiden':
            key_added += '_R'
        groups = rename_groups(
            adata,
            key_added,
            restrict_key,
            restrict_categories,
            restrict_indices,
            groups,
        )
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(map(str, np.unique(groups))),
    )
    if collect_marginals:
        adata.obsm[f"CM_{key_added}"] = pv_array
    # store information on the clustering parameters
    adata.uns['leiden'] = {}
    adata.uns['leiden']['params'] = dict(
        resolution=resolution,
        random_state=random_state,
        n_iterations=n_iterations,
        n_init=n_init,
        collect_marginals=collect_marginals
    )
    logg.info(
        '    finished',
        time=start,
        deep=(
            f'found {len(np.unique(groups))} clusters and added\n'
            f'    {key_added!r}, the cluster labels (adata.obs, categorical)'
        ),
    )
    return adata if copy else None
