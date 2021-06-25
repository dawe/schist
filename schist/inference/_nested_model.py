from typing import Optional, Tuple, Sequence, Type, Union, Dict

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from joblib import delayed, Parallel

from scanpy import logging as logg
from scanpy.tools._utils_clustering import rename_groups, restrict_adjacency
from scanpy._utils import get_igraph_from_adjacency


try:
    import graph_tool.all as gt
except ImportError:
    raise ImportError(
        """Please install the graph-tool library either visiting

        https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

        or by conda: `conda install -c conda-forge graph-tool`
        """
    )

def nested_model(
    adata: AnnData,
    deg_corr: bool = True,
    tolerance: float = 1e-6,
    n_sweep: int = 10,
    beta: float = np.inf,
    samples: int = 100,
    collect_marginals: bool = True,
    n_jobs: int = -1,
    *,
    restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
    random_seed: Optional[int] = None,
    key_added: str = 'nsbm',
    adjacency: Optional[sparse.spmatrix] = None,
    neighbors_key: Optional[str] = 'neighbors',
    directed: bool = False,
    use_weights: bool = False,
    save_model: Union[str, None] = None,
    copy: bool = False,
    minimize_args: Optional[Dict] = {},
    dispatch_backend: Optional[str] = 'processes',
#    equilibrate_args: Optional[Dict] = {},
) -> Optional[AnnData]:
    """\
    Cluster cells into subgroups [Peixoto14]_.

    Cluster cells using the nested Stochastic Block Model [Peixoto14]_,
    a hierarchical version of Stochastic Block Model [Holland83]_, performing
    Bayesian inference on node groups. NSBM should circumvent classical
    limitations of SBM in detecting small groups in large graphs
    replacing the noninformative priors used by a hierarchy of priors
    and hyperpriors.

    This requires having ran :func:`~scanpy.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    deg_corr
        Whether to use degree correction in the minimization step. In many
        real world networks this is the case, although this doesn't seem
        the case for KNN graphs used in scanpy.
    tolerance
        Tolerance for fast model convergence.
    n_sweep 
        Number of iterations to be performed in the fast model MCMC greedy approach
    beta
        Inverse temperature for MCMC greedy approach    
    samples
        Number of initial minimizations to be performed. The one with smaller
        entropy is chosen
    n_jobs
        Number of parallel computations used during model initialization
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
    save_model
        If provided, this will be the filename for the PartitionModeState to 
        be saved    
    copy
        Whether to copy `adata` or modify it inplace.
    random_seed
        Random number to be used as seed for graph-tool

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
    `adata.obsm['CA_nsbm_level_{n}']`
        A `np.ndarray` with cell probability of belonging to a specific group
    `adata.uns['schist']['state']`
        The NestedBlockModel state object
    """

    if random_seed:
        np.random.seed(random_seed)
    
    seeds = np.random.choice(range(samples**2), size=samples, replace=False)
        

    if collect_marginals and samples < 100:
        logg.warning('Collecting marginals requires sufficient number of samples\n'
                     f'It is now set to {samples} and should be at least 100')
        

    start = logg.info('minimizing the nested Stochastic Block Model')
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

    states = [gt.NestedBlockState(g=g,
                                  state_args=dict(deg_corr=deg_corr,
                                  recs=recs,
                                  rec_types=rec_types
                                  )) for n in range(samples)]

    def fast_min(state, beta, n_sweep, fast_tol, seed=None):
        if seed:
            gt.seed_rng(seed)
        dS = 1
        while np.abs(dS) > fast_tol:
            dS, _, _ = state.multiflip_mcmc_sweep(beta=beta, niter=n_sweep, c=0.5)
        return state                            
            
    states = Parallel(n_jobs=n_jobs, prefer=dispatch_backend)(
        delayed(fast_min)(states[x], beta, n_sweep, tolerance, seeds[x]) for x in range(samples)
    )
    logg.info('        minimization step done', time=start)
    pmode = gt.PartitionModeState([x.get_bs() for x in states], converge=True, nested=True)
    bs = pmode.get_max_nested()
    logg.info('        consensus step done', time=start)
    
    if save_model:
        import pickle
        fname = save_model
        if not fname.endswith('pkl'):
            fname = f'{fname}.pkl'
        logg.info(f'Saving model into {fname}')    
        with open(fname, 'wb') as fout:
            pickle.dump(pmode, fout, 2)
    
    # prune redundant levels at the top
    bs = [x for x in bs if len(np.unique(x)) > 1]
    bs.append(np.array([0], dtype=np.int32)) #in case of type changes, check this
    state = gt.NestedBlockState(g, bs)
    

    logg.info('    done', time=start)
    u_groups = np.unique(bs[0])
    n_groups = len(u_groups)
    last_group = np.max(u_groups) + 1

    if collect_marginals:
        # note that the size of this will be equal to the number of the groups in Mode
        # but some entries won't sum to 1 as in the collection there may be differently
        # sized partitions
        pv_array = pmode.get_marginal(g).get_2d_array(range(last_group)).T[:, u_groups] / samples	
         
    groups = np.zeros((g.num_vertices(), len(bs)), dtype=int)

    for x in range(len(bs)):
        # for each level, project labels to the vertex level
        # so that every cell has a name. Note that at this level
        # the labels are not necessarily consecutive
        groups[:, x] = state.project_partition(x, 0).get_array()

    groups = pd.DataFrame(groups).astype('category')

    # rename categories from 0 to n
    for c in groups.columns:
        ncat = len(groups[c].cat.categories)
        new_cat = [u'%s' % x for x in range(ncat)]
        groups[c].cat.rename_categories(new_cat, inplace=True)

    levels = groups.columns
    
    # recode block names to have consistency with group names
    i_groups = groups.astype(int)
    bs = [i_groups.iloc[:, 0].values]
    for x in range(1, groups.shape[1]):
        bs.append(np.where(pd.crosstab(i_groups.iloc[:, x - 1], i_groups.iloc[:, x])> 0)[1])
    state = gt.NestedBlockState(g, bs)
    del(i_groups)

    if restrict_to is not None:
        groups.index = adata.obs[restrict_key].index
    else:
        groups.index = adata.obs_names

    # add column names
    groups.columns = [f"{key_added}_level_{level}" for level in range(len(bs))]

    # remove any column with the same key
    keep_columns = [x for x in adata.obs.columns if not x.startswith('%s_level_' % key_added)]
    adata.obs = adata.obs[keep_columns]
    adata.obs = pd.concat([adata.obs, groups], axis=1)

    # add some unstructured info

    adata.uns['schist'] = {}
    adata.uns['schist']['stats'] = dict(
    level_entropy=np.array([state.level_entropy(x) for x in range(len(state.levels))]),
    modularity=np.array([gt.modularity(g, state.project_partition(x, 0))
                         for x in range(len((state.levels)))])
    )

    adata.uns['schist']['state'] = state

    # now add marginal probabilities.

    if collect_marginals:
        # add marginals for level 0, the sum up according to the hierarchy
        adata.obsm[f"CM_{key_added}_level_0"] = pv_array
        for group in groups.columns[1:]:
            ct = pd.crosstab(groups[groups.columns[0]], groups[group], normalize='index')
            adata.obsm[f'CM_{group}'] = pv_array @ ct.values

    # last step is recording some parameters used in this analysis
    adata.uns['schist']['params'] = dict(
        model='nested',
        key_added=key_added,
        samples=samples,
        collect_marginals=collect_marginals,
        random_seed=random_seed,
    )


    logg.info(
        '    finished',
        time=start,
        deep=(
            f'and added\n'
            f'    {key_added!r}, the cluster labels (adata.obs, categorical)'
        ),
    )
    return adata if copy else None

