from typing import Optional, Tuple, Sequence, Type, Union, Dict

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from joblib import delayed, Parallel
from natsort import natsorted
from scanpy import logging as logg
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency
from scanpy._utils import get_igraph_from_adjacency

import graph_tool.all as gt


def planted_model(
    adata: AnnData,
    n_sweep: int = 10,
    beta: float = np.inf, 
    tolerance = 1e-6,
    collect_marginals: bool = True,
    deg_corr: bool = True,
    samples: int = 100,
    n_jobs: int = -1,
    *,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    random_seed: Optional[int] = None,
    key_added: str = 'ppbm',
    adjacency: Optional[sparse.spmatrix] = None,
    neighbors_key: Optional[str] = 'neighbors',
    directed: bool = False,
    use_weights: bool = False,
    copy: bool = False,
    minimize_args: Optional[Dict] = {},
) -> Optional[AnnData]:
    """\
    Cluster cells into subgroups [Peixoto14]_.

    Cluster cells using the  Planted Partition Block Model [Peixoto14]_, performing
    Bayesian inference on node groups. This function, in particular, uses
    the Planted Block Model, which is particularly suitable in case of
    assortative graphs and it returns the optimal number of communities

    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    n_sweep
        Number of MCMC sweeps to get the initial guess
    beta
        Inverse temperature for the initial MCMC sweep        
    tolerance
        Difference in description length to stop MCMC sweep iterations        
    collect_marginals
        Whether or not collect node probability of belonging
        to a specific partition.
    deg_corr
        Whether to use degree correction in the minimization step. In many
        real world networks this is the case, although this doesn't seem
        the case for KNN graphs used in scanpy.
    samples
        Number of initial minimizations to be performed. This influences also the 
        precision for marginals
    key_added
        `adata.obs` key under which to add the cluster labels.
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']` in case of scanpy<=1.4.6 or
        `adata.obsp[neighbors_key][connectivity_key]` for scanpy>1.4.6
    neighbors_key
        The key passed to `sc.pp.neighbors`
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges). Note that this
        increases computation times
    copy
        Whether to copy `adata` or modify it inplace.
    random_seed
        Random number to be used as seed for graph-tool
    n_jobs
        Number of parallel computations used during model initialization

    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of samples) that stores the subgroup id
        (`'0'`, `'1'`, ...) for each cell.
    `adata.uns['schist']['params']`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    `adata.uns['schist']['stats']`
        A dict with the values returned by mcmc_sweep
    `adata.obsm['CM_ppbm']`
        A `np.ndarray` with cell probability of belonging to a specific group
    `adata.uns['schist']['state']`
        The BlockModel state object
    """

    if random_seed:
        np.random.seed(random_seed)
    
    seeds = np.random.choice(range(samples**2), size=samples, replace=False)

    if collect_marginals and samples < 100:
        logg.warning('Collecting marginals requires sufficient number of samples\n'
                     f'It is now set to {samples} and should be at least 100')

    start = logg.info('minimizing the Planted Partition Block Model')
    adata = adata.copy() if copy else adata
    # are we clustering a user-provided graph or the default AnnData one?
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
    # convert it to igraph and graph-tool
    g = get_igraph_from_adjacency(adjacency, directed=directed)
    g = g.to_graph_tool()
    gt.remove_parallel_edges(g)

    recs=[]
    rec_types=[]
    if use_weights:
        # this is not ideal to me, possibly we may need to transform
        # weights. More tests needed.
        recs=[g.ep.weight]
        rec_types=['real-normal']

    if samples < 1:
        samples = 1
        
    # initialize  the block states
    def fast_min(state, beta, n_sweep, fast_tol, seed=None):
        if seed:
            gt.seed_rng(seed)
        dS = 1
        while np.abs(dS) > fast_tol:
            dS, _, _ = state.multiflip_mcmc_sweep(beta=beta, niter=n_sweep)
        return state

    states = [gt.PPBlockState(g) for x in range(samples)]
        
    # perform a mcmc sweep on each 
    # no list comprehension as I need to collect stats
        
    states = Parallel(n_jobs=n_jobs)(
             delayed(fast_min)(states[x], beta, n_sweep, tolerance, seeds[x]) for x in range(samples)
             )
        
    pmode = gt.PartitionModeState([x.get_blocks().a for x in states], converge=True)
        
    bs = pmode.get_max(g)
    state = gt.PPBlockState(g, b=bs)
    logg.info('    done', time=start)

    groups = np.array(bs.get_array())
    u_groups = np.unique(groups)
    n_groups = len(u_groups)
    if collect_marginals:
        pv_array = pmode.get_marginal(g).get_2d_array(range(n_groups)).T / samples

    rosetta = dict(zip(u_groups, range(len(u_groups))))
    groups = np.array([rosetta[x] for x in groups])

    if restrict_to is not None:
        if key_added == 'ppbm':
            key_added += '_R'
        groups = rename_groups(
            adata,
            key_added,
            restrict_key,
            restrict_categories,
            restrict_indices,
            groups,
        )

    # add column names
    adata.obs[key_added] = pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(map(str, np.unique(groups))),
    )

    # add some unstructured info

    adata.uns['schist'] = {}
    adata.uns['schist']['stats'] = dict(
    entropy=state.entropy(),
    modularity=gt.modularity(g, state.get_blocks())
    )
    adata.uns['schist']['state'] = state

    # now add marginal probabilities.

    if collect_marginals:
        # cell marginals will be a list of arrays with probabilities
        # of belonging to a specific group
        adata.obsm[f"CM_{key_added}"] = pv_array

    # last step is recording some parameters used in this analysis
    adata.uns['schist']['params'] = dict(
        model='planted',
        samples=samples,
        collect_marginals=collect_marginals,
        random_seed=random_seed
    )


    logg.info(
        '    finished',
        time=start,
        deep=(
            f'found {state.get_B()} clusters and added\n'
            f'    {key_added!r}, the cluster labels (adata.obs, categorical)'
        ),
    )
    return adata if copy else None

