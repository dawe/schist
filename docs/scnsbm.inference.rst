************************
scnsbm.inference package
************************

scnsbm.inference._flat_model
############################

=,
scnsbm.inference.**_flat_model**(*adata:AnnData*,*max_iterations:1000000*,*epsilon:1e-3*,*equilibrate:True*,*wait:1000*,*nbreaks:2*,*collect_marginals:False*,*niter_collect:10000*,*deg_corr:False*,*multiflip:True*,*fast_model:False*,*beta_range:1.0,100.0*,*steps_anneal:int=5*,*restrict_to:Optional=None*,*random_seed:Optional=None*,*key_added:sbm*,*adjacency:Optional=None*,*directed:False*,*use_weights:False*,*save_state:False*,*copy:False*,*minimize_args:Optional*={},*equilibrate_args:Optional*={})
=,

Cluster cells into subgroups,using the Stochastic Block Model [Peixoto14]_, performing Bayesian inference on node groups. 

This requires having ran :func:`~scanpy.pp.neighbors` or :func:`~scanpy.external.pp.bbknn` first.

Parameters:
***********
    *   **adata**:**AnnData**
            The annotated data matrix.
    *   **max_iterations**:**int**(default:1000000)
            Maximal number of iterations to be performed by the equilibrate step.
    *   **epsilon**:''**float**''(default:''1e-3'')
            Relative changes in entropy smaller than epsilon will not be considered as record-breaking.
    *   **equilibrate**:''**bool**''(default:''True)
            Whether or not perform the mcmc_equilibrate step. Equilibration should always be performed. Note, also, that without equilibration it won't be possible to collect marginals.
    *   **collect_marginals**:''**bool**''(default:''False'')
            Whether or not collect node probability of belonging to a specific partition.
    *   **niter_collect**:''**int**''(default:''10000'')
            Number of iterations to force when collecting marginals. This will increase the precision when calculating probabilites.
    *   **wait**:''**int**''(default:''1000'')
            Number of iterations to wait for a record-breaking event. Higher values result in longer computations. Set it to small values when performing quick tests.
    *   **nbreaks**:''**int**''(default:''2'')
            Number of iteration intervals (of size `wait`) without record-breaking events necessary to stop the algorithm.
    *   **deg_corr**:''**bool**''(default:''False'')
            Whether to use degree correction in the minimization step. In many real world networks this is the case, although this doesn't seem the case for KNN graphs used in scanpy.
    *   **multiflip**:''**bool**''(default:''True'')
            Whether to perform MCMC sweep with multiple simultaneous moves to sample network partitions. It may result in slightly longer runtimes, but under the hood it allows for a more efficient space exploration.
    *   **fast_model**:''**bool**''(default:''False'')
            Whether to skip initial minization step and let the MCMC find a solution. This approach tend to be faster and consume less memory, but less accurate.
    *   **beta_range**:''**Tuple**[**float**]''(default:(''1.0,100.0''))
            Inverse temperature at the beginning and the end of the equilibration.
    *   **steps_anneal**:''**int**''(default:''5'')
            Number of steps in which the simulated annealing is performed.
    *   **key_added**:''**str**''(default:'''sbm''')
            `adata.obs` key under which to add the cluster labels.
    *   **adjacency**:''*Optional*[**sparse.spmatrix**]''(default:''None'')
            Sparse adjacency matrix of the graph, defaults to `adata.uns['neighbors']['connectivities']`.
    *   **directed**:''**bool**''(default:'False'')
            Whether to treat the graph as directed or undirected.
    *   **use_weights**:''**bool**''(default:''False'')
            If `True`, edge weights from the graph are used in the computation (placing more emphasis on stronger edges). Note that this increases computation times.
    *   **save_state**:''**bool**''(default:''False'')
            Whether to keep the block model state saved for subsequent custom analysis with graph-tool. Use only for debug session, state is not (yet) supported for `sc.write` function.
    *   **copy**:''**bool**''(default:''False'')
            Whether to copy `adata` or modify it inplace.
    *   **random_seed**:''**Optional**[**int**]''(default:''None'')
            Random number to be used as seed for graph-tool.

Returns:
********
    *   ''`adata.obs[key_added]`''
            Array of dim (number of samples) that stores the subgroup id (`'0'`, `'1'`, ...) for each cell.
    *   ''`adata.uns['sbm']['params']`''
            A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    *   ''`adata.uns['sbm']['stats']`''
            A dict with the values returned by mcmc_sweep
    *   ''`adata.uns['sbm']['cell_marginals']`''
            A `np.ndarray` with cell probability of belonging to a specific group
    *   ''`adata.uns['sbm']['state']`''
            The BlockModel state object

scnsbm.inference._nested_model
############################

-,
scnsbm.inference.**_flat_model**(*adata:AnnData*,*max_iterations:1000000*,*epsilon:1e-3*,*equilibrate:True*,*wait:1000*,*nbreaks:2*,*collect_marginals:False*,*niter_collect:10000*,*hierarchy_length:10*,*deg_corr:False*,*multiflip:True*,*fast_model:False*,*beta_range:1.0,100.0*,*steps_anneal:int=5*,*restrict_to:Optional=None*,*random_seed:Optional=None*,*key_added:nsbm*,*adjacency:Optional=None*,*directed:False*,*use_weights:False*,*save_state:False*,*prune=False*,*return_low=False*,*copy:False*,*minimize_args:Optional*={},*equilibrate_args:Optional*={})
-,

Cluster cells using the nested Stochastic Block Model [Peixoto14]_, a hierarchical version of Stochastic Block Model [Holland83]_, performing Bayesian inference on node groups. NSBM should circumvent classical limitations of SBM in detecting small groups in large graphs replacing the noninformative priors used by a hierarchy of priors and hyperpriors. 

This requires having ran :func:`~scanpy.pp.neighbors` or :func:`~scanpy.external.pp.bbknn` first.

Parameters:
***********
    *   **adata**:''**AnnData**''
            The annotated data matrix.
    *   **max_iterations**:''**int**''(default:''1000000'')
            Maximal number of iterations to be performed by the equilibrate step.
    *   **epsilon**:''**float**''(default:''1e-3'')
            Relative changes in entropy smaller than epsilon will not be considered as record-breaking.
    *   **equilibrate**:''**bool**''(default:''True)
            Whether or not perform the mcmc_equilibrate step. Equilibration should always be performed. Note, also, that without equilibration it won't be possible to collect marginals.
    *   **collect_marginals**:''**bool**''(default:''False'')
            Whether or not collect node probability of belonging to a specific partition.
    *   **niter_collect**:''**int**''(default:''10000'')
            Number of iterations to force when collecting marginals. This will increase the precision when calculating probabilites.
    *   **wait**:''**int**''(default:''1000'')
            Number of iterations to wait for a record-breaking event. Higher values result in longer computations. Set it to small values when performing quick tests.
    *   **nbreaks**:''**int**''(default:''2'')
            Number of iteration intervals (of size `wait`) without record-breaking events necessary to stop the algorithm.
    *   **hierarchy_length**:''**int**''(default:''10'')
            Initial length of the hierarchy. When large values are passed, the top-most levels will be uninformative as they will likely contain the very same groups. Increase this valus if a very large number of cells is analyzed (>100.000).
    *   **deg_corr**:''**bool**''(default:''False'')
            Whether to use degree correction in the minimization step. In many real world networks this is the case, although this doesn't seem the case for KNN graphs used in scanpy.
    *   **multiflip**:''**bool**''(default:''True'')
            Whether to perform MCMC sweep with multiple simultaneous moves to sample network partitions. It may result in slightly longer runtimes, but under the hood it allows for a more efficient space exploration.
    *   **fast_model**:''**bool**''(default:''False'')
            Whether to skip initial minization step and let the MCMC find a solution. This approach tend to be faster and consume less memory, but less accurate.
    *   **beta_range**:''**Tuple**[**float**]''(default:(''1.0,100.0''))
            Inverse temperature at the beginning and the end of the equilibration.
    *   **steps_anneal**:''**int**''(default:''5'')
            Number of steps in which the simulated annealing is performed.
    *   **key_added**:''**str**''(default:'''nsbm''')
            `adata.obs` key under which to add the cluster labels.
    *   **adjacency**:''*Optional*[**sparse.spmatrix**]''(default:''None'')
            Sparse adjacency matrix of the graph, defaults to `adata.uns['neighbors']['connectivities']`.
    *   **directed**:''**bool**''(default:'False'')
            Whether to treat the graph as directed or undirected.
    *   **use_weights**:''**bool**''(default:''False'')
            If `True`, edge weights from the graph are used in the computation (placing more emphasis on stronger edges). Note that this increases computation times.
    *   **save_state**:''**bool**''(default:''False'')
            Whether to keep the block model state saved for subsequent custom analysis with graph-tool. Use only for debug session, state is not (yet) supported for `sc.write` function.
    *   **prune**:''**bool**''(default:''False'')
            Some high levels in hierarchy may contain the same information in terms of cell assignments, even if they apparently have different group names. When this option is set to `True`, the function only returns informative levels. Note, however, that cell_marginals are still reported for all levels. Pruning does not rename group levels.
    *   **return_low**:''**bool**''(default:''False'')
            Whether or not return nsbm_level_0 in adata.obs. This level usually contains so many groups that it cannot be plot anyway, but it may be useful for particular analysis. By default it is not returned
    *   **copy**:''**bool**''(default:''False'')
            Whether to copy `adata` or modify it inplace.
    *   **random_seed**:''**Optional**[**int**]''(default:''None'')
            Random number to be used as seed for graph-tool.

Returns:
********
    *   ''`adata.obs[key_added]`''
            Array of dim (number of samples) that stores the subgroup id (`'0'`, `'1'`, ...) for each cell. 
    *   ''`adata.uns['nsbm']['params']`''
            A dict with the values for the parameters `resolution`, `random_state`, and `n_iterations`.
    *   ''`adata.uns['nsbm']['stats']`''
            A dict with the values returned by mcmc_sweep
    *   ''`adata.uns['nsbm']['cell_marginals']`''
            A `np.ndarray` with cell probability of belonging to a specific group
    *   ''`adata.uns['nsbm']['state']`''
            The NestedBlockModel state object

.. automodule:: scnsbm.inference
   :members:
   :undoc-members:
   :show-inheritance:
