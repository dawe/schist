from typing import Optional, Tuple, Sequence, Type, Union, Dict, List
from anndata import AnnData
import numpy as np
import scipy.sparse as sparse
from . import fit_model_multi
from scanpy import logging as logg

try:
    import graph_tool.all as gt
except ImportError:
    raise ImportError(
        """Please install the graph-tool library either visiting

        https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

        or by conda: `conda install -c conda-forge graph-tool`
        """
    )

def nested_model_multi(
    adatas: List[AnnData],
    deg_corr: bool = True,
    tolerance: float = 1e-6,
    n_sweep: int = 10,
    beta: float = np.inf,
    n_init: int = 100,
    collect_marginals: bool = True,
    n_jobs: int = -1,
    refine_model: bool = False,
    refine_iter: int = 100,
    overlap: bool = False,
    max_iter: int = 100000,
    *,
    random_seed: Optional[int] = None,
    key_added: str = 'multi_nsbm',
    adjacency: Optional[List[sparse.spmatrix]] = None,
    neighbors_key: Optional[List[str]] = ['neighbors'],
    directed: bool = False,
    use_weights: bool = False,
    save_model: Union[str, None] = None,
    copy: bool = False,
#    minimize_args: Optional[Dict] = {},
    dispatch_backend: Optional[str] = 'threads',
#    equilibrate_args: Optional[Dict] = {},
) -> Optional[List[AnnData]]:
    """\
    This function has been deprecated and it soon will be removed.
    It now wraps ``scs.inference.fit_model()`` function.
    """

    logg.warning('This function has been deprecated, and soon will be removed\n'
                  'use `scs.inference.fit_model_multi(adata, model="nsbm")` instead')

    return model_multi(adatas, model='nsbm',
                 deg_corr = deg_corr,
                 tolerance = tolerance,
                 n_sweep = n_sweep,
                 beta = beta,
                 n_init = n_init,
                 collect_marginals = collect_marginals,
                 n_jobs = n_jobs,
                 refine_model = refine_model,
                 refine_iter = refine_iter,
                 overlap = overlap,
                 max_iter = max_iter,
                 random_seed = random_seed,
                 key_added = key_added,
                 neighbors_key = neighbors_key, 
                 directed = directed,
                 use_weights = use_weights,
                 save_model = save_model,
                 copy = copy,
                 dispatch_backend = dispatch_backend  
    )
