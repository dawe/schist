.. _large_samples:

==========================
Analysis of large datasets
==========================

It is stated in the `original paper <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04489-7>`__
(and also raised by users) that ``schist`` may be slow, especially for
large datasets. While the underlying library is computationally
efficient, the MCMC problem may be difficult to minimize and scales with
the square of the number of the edges. In lay terms, this means that the
number of cells and the number of neighbors used to build the kNN graph
will, in turn, influence the speed of the procedure. Moreover, it is
difficult to predict the actual execution time due to the specific
nature of the data. In other words, there are datasets in which the
algorithms converges pretty fast because the community structure is
“easy to find”.

As for July 2026, the major library under the hood (``graph-tool``) implemented a 
parallel version of the algorithm that is used to identify the most likely block scheme 
in graphs.
Starting from ``schist`` version 0.10.0, we support this major improvement, analysis of 
large networks is no more prohibitive, given that enough computational resources are 
available. Parallelization relies on `OpenMP <https://www.openmp.org>`_, the more 
threads are available the faster is convergence. For SBM and NSBM models, some options 
are also helpful to speed up computation. Note that memory can be an issue.

The following tutorial goes through the analysis of the HCA Lung
datasets, downloaded from the cellxgene data portal, which includes more
than 2 millions cells. The whole processing is made using 32 threads.

.. code:: python

	>>> import numpy as np
	>>> import schist as scs
	>>> import scanpy as sc
	>>> import time
	>>> import tracemalloc as tm
	>>> sc.settings.verbosity=2

The dataset is rather huge (20 Gb), we will perform operations without
loading it into memory. Whenever we will need to modify its attributes,
we will create a lightweight copy, including only the relevant
information.

.. code:: python

	>>> print(f"Loading data ({time.ctime()})")
	>>> sc.settings.verbosity=2
	>>> tm.start()
	>>> adata = sc.read("local.h5ad", backed='r')
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()


.. parsed-literal::

	Loading data (Mon Aug 10 13:54:26 2026)
	Size: 3734.81 MB, peak: 4298.19 MB


---------------------------
The Planted Partition Model
---------------------------

We first test performance of the fastest model (planted partition), which was
shown to be feasible with the old version (it took ~13 hours to converge). We start 
with the default configuration, without collecting cell marginals.

.. code:: python

	>>> scs.inference.fit_model(adata,
	...	collect_marginals=False,
	...	assortative=True,
	...	bisection=True,
	...	adjacency=adata.obsp['connectivities'],
	...	save_model='ppbm_scs_cm0_b1.pkl'
	...	)
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()


.. parsed-literal::

	Minimizing the Model
			done (0:07:54)
			done (0:07:54)
	Saving model into ppbm_scs_cm0_b1.pkl
			done (0:08:02)
		finished (0:08:55)
	Size: 3756.79 MB, peak: 6066.77 MB


Enabling sampling from the posterior to compute cell marginals may be required, for
example when working with cell trajectories (see :ref:`cellrank_integration`). This 
adds extra effort and time to the overall process, which may be extremely demanding 
(see below). 

.. code:: python

	>>> scs.inference.fit_model(adata,
	...	collect_marginals=True,
	...	assortative=True,
	...	bisection=True,
	...	adjacency=adata.obsp['connectivities'],
	...	save_model='ppbm_scs_cm1_b0.pkl'
	...	)
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()


.. parsed-literal::

	Minimizing the Model
			done (0:08:08)
	Sampling the posterior and getting cell marginals
	100%|██████████| 100/100 [2:35:22<00:00, 93.23s/it] 
			done (2:43:31)
	Computing the consesus partition
			done (2:50:20)
	Saving model into ppbm_scs_cm1_b0.pkl
			done (2:50:33)
		finished (2:54:58)
	Size: 22424.49 MB, peak: 31941.33 MB



As expected, resources required to compute cell marginals are large, but still
manageable. The memory required to process the data increased from 6 Gb to 32 Gb
to store and process the 100 samples.


------------------
The default models
------------------

The planted partition model is an effective approach and returns assortative
communities, in the analysis of kNN graphs derived from single cell data
it is reasonable to expect those communities to reflect the population
structure in terms of cell types. Nevertheless, we may be interested in
the Stochastic Block Model or its Nested formulation as well. Those
approach incorporate other priors in the model, possibly identifying
other properties of the dataset. As a consequence, more computing resources
are needed. We first start fitting a SBM model, without cell marginals.


.. code:: python

	>>> scs.inference.fit_model(adata, 
	...	collect_marginals=False, 
	...	nested=False, 
	...	bisection=True,
	...	adjacency=adata.obsp['connectivities'],
	...	save_model='sbm_scs_cm0_b1.pkl'
	... )	
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()



.. parsed-literal::

	Minimizing the Model
			done (2:16:08)
			done (2:16:08)
	Saving model into sbm_scs_cm0_b1.pkl
			done (2:19:31)
		finished (2:21:07)
	Size: 3757.49 MB, peak: 7748.71 MB


It took less than 2.5 hours to converge to a solution. As a comparison, the
previous version needed almost 4 days to achieve the same results. This is in 
perfect agreement with the fact we are using 32 threads. We can push this further
by enabling `specific search options <https://graph-tool.skewed.de/static/docs/stable/demos/inference/inference.html#dealing-with-large-networks>`_, 
disabling the bisection search in the MCMC step.

.. code:: python

	>>> scs.inference.fit_model(adata, 
	...	collect_marginals=False, 
	...	nested=False, 
	...	bisection=False,
	...	adjacency=adata.obsp['connectivities'],
	...	save_model='sbm_scs_cm0_b0.pkl'
	... )
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()


.. parsed-literal::

	Minimizing the Model
			done (1:19:25)
			done (1:19:25)
	Saving model into sbm_scs_cm0_b0.pkl
			done (1:22:53)
		finished (1:24:30)
	Size: 3757.27 MB, peak: 7770.93 MB


The agglomerative search takes half the time. We can compare results in terms of 
partitioning schemes:

.. code:: python

	>>> import numpy as np
	>>> import graph_tool.all as gt
	>>> import pickle
	>>> import sklearn.metrics as skm
	>>> with open('sbm_scs_cm0_b1.pkl', 'rb') as fh:
	... 	b1_state = pickle.load(fh)['State']
	>>> with open('sbm_scs_cm0_b0.pkl', 'rb') as fh:
	... 	b0_state = pickle.load(fh)['State']
	>>> print("Entropy difference: ", b1_state.entropy() - b0_state.entropy())
	
.. parsed-literal::

	Entropy difference:  -1935484.1506269574


The bisection search generates a model that is much more likely than the other, and
yet we can check how different they are:

.. code:: python

	>>> print(f'Model with bisection search has {b1_state.get_B()} blocks')
	>>> print(f'Model with agglomerative search has {b0_state.get_B()} blocks')
	>>> print('Adjusted Rand score: ', skm.adjusted_rand_score(b1_state.b.a, b0_state.b.a))
	>>> print("Homogeneity: ", h)
	>>> print("Completeness: ", c)
	>>> print("V measure: ", v)

.. parsed-literal::	

	Model with bisection search has 2944 blocks
	Model with agglomerative search has 4019 blocks	
	Adjusted Rand score:  0.5113166674314255
	Homogeneity:  0.8779627143141272
	Completeness:  0.8460632844422329
	V measure:  0.8617178838674843


In the end, agglomerative search produces a more fragmented solution, yet the two
solutions are quite overlapping.
Lastly, we may want to test how cell marginals impact on the required resources:

.. code:: python

	>>> scs.inference.fit_model(adata, 
	...	collect_marginals=True, 
	...	nested=False, 
	...	bisection=False,
	...	adjacency=adata.obsp['connectivities'],
	...	save_model='sbm_scs_cm1_b0.pkl'
	... )
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()


.. parsed-literal::

	Minimizing the Model
			done (1:17:49)
	Sampling the posterior and getting cell marginals
	100%|██████████| 100/100 [21:21:37<00:00, 768.97s/it]
			done (22:39:26)
	Computing the consesus partition
			done (23:51:28)
	Saving model into sbm_scs_cm1_b0.pkl
			done (23:55:09)
		finished (1 day, 0:02:34)
	Size: 69285.43 MB, peak: 102232.87 MB


Computation of cell marginals requires substantially more resources for SBM compared to
PPBM, both in term of time (1 day vs 3 hours) and memory (102Gb vs 32Gb), these 
requirements should always been kept in mind.


As for the hierachical version of the SBM model, we can test the performance with
an additional switch (``simple_init``). First run the default model:

.. code:: python

	>>> scs.inference.fit_model(adata, 
	...	collect_marginals=False, 
	...	nested=True, 
	...	bisection=True,
	...	adjacency=adata.obsp['connectivities'],
	...	save_model='nsbm_scs_cm0_b1_si0.pkl'
	... )
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()


.. parsed-literal::

	Minimizing the Model
			done (16:37:48)
			done (16:37:48)
	Saving model into nsbm_scs_cm0_b1_si0.pkl
			done (16:41:18)
		finished (16:42:26)
	Size: 3800.21 MB, peak: 13193.28 MB


It takes slightly less than 17 hours, which is not bad, considered that the old version
failed to converge in reasonable times (more than a month). When bisection search is
disabled, runtimes are greatly reduced:

.. code:: python

	>>> scs.inference.fit_model(adata, 
	...	collect_marginals=False, 
	...	nested=True, 
	...	bisection=False,
	...	adjacency=adata.obsp['connectivities'],
	...	save_model='nsbm_scs_cm0_b0_si0.pkl'
	... )
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()


.. parsed-literal::

	Minimizing the Model
			done (1:52:30)
			done (1:52:30)
	Saving model into nsbm_scs_cm0_b0_si0.pkl
			done (1:56:01)
		finished (1:56:22)
	Size: 3786.25 MB, peak: 7861.88 MB


Convergence happens in less than 2 hours. When working with nested SBM, another
option is available (``simple_init``), which simplifies the search strategy estimating
an upper bound for the number of groups at each level:

.. code:: python

	>>> scs.inference.fit_model(adata, 
	...	collect_marginals=False, 
	...	nested=True, 
	...	bisection=False,
	... simple_init=True,
	...	adjacency=adata.obsp['connectivities'],
	...	save_model='nsbm_scs_cm0_b0_si1.pkl'
	... )
	>>> size, peak = tm.get_traced_memory()
	>>> size /= 1<<20
	>>> peak /= 1<<20
	>>> print(f'Size: {size:.2f} MB, peak: {peak:.2f} MB')
	>>> tm.reset_peak()


.. parsed-literal::

	Minimizing the Model
			done (2:37:01)
			done (2:37:01)
	Saving model into nsbm_scs_cm0_b0_si1.pkl
			done (2:40:36)
		finished (2:42:41)
	Size: 3831.45 MB, peak: 16942.78 MB


Despite the runtime is greatly reduced compared to the default parameters, this
dataset does not benefit from this strategy, the time required to converge is larger.
It is also interesting to notice that the amount of memory is varying among the 
strategies, due to the differences in number of blocks found:


.. code:: python

	>>> with open('nsbm_scs_cm0_b1_si0.pkl', 'rb') as fh:
	... 	b1_state = pickle.load(fh)['State']
	>>> with open('nsbm_scs_cm0_b0_si0.pkl', 'rb') as fh:
	... 	b0_state = pickle.load(fh)['State']
	>>> with open('nsbm_scs_cm0_b0_si1.pkl', 'rb') as fh:
	... 	s1_state = pickle.load(fh)['State']
	>>> print([x.get_B() for x in b1_state.get_levels()])
	>>> print([x.get_B() for x in b0_state.get_levels()])
	>>> print([x.get_B() for x in s1_state.get_levels()])
	
	
.. parsed-literal::

	[87169, 10915, 3044, 790, 243, 91, 35, 10, 4, 2, 1]
	[3756, 1717, 709, 307, 130, 49, 20, 6, 2, 1]
	[74247, 17676, 5254, 1772, 690, 292, 124, 57, 34, 20, 10, 8, 7, 6, 5, 5, 3, 2, 1, 1, 1, 1, 1, 1]



These results actually reveal that the second option (``bisection=False``, 
``simple_init=True``) converges to a solution that seems to lack fine-grained descriptions
at levels 0 and 1. Indeed there is some correspondence between level 2 of the first 
solution and level 0 of the second solution:

.. code:: python

	>>> l2_B1 = b1_state.project_partition(2, 0).a # project level 2 to single cells
	>>> l0_B0 = b0_state.project_partition(0, 0).a
	>>> h, c, v = skm.homogeneity_completeness_v_measure(l2_B1, l0_B0)
	>>> print("Homogeneity: ", h)
	>>> print("Completeness: ", c)
	>>> print("V measure: ", v)
	
.. parsed-literal::

	Homogeneity:  0.779528789164119
	Completeness:  0.830511899563901
	V measure:  0.8042131357126281

