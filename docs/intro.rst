************
Introduction
************

Understanding NSBM
##################

`schist` is based on Nested Stochastic Block Models (NSBM), a generative process based on the notion of group of nodes. Here, we use NSBM to cluster cells in scRNA-seq experiments.

Bayesian approach
#################

Cell to cell relations are commonly represented as a neighborhood graph (typically KNN or  SNN), cell groups are identified as graph communities, this step is usually performed maximising graph modularity. In `schist`, instead, partitioning is performed maximising the  Bayesian posterior probability **P(b|A)**, that is the likehood of generating a network _A_ with a partition _b_ and it is obtained according to

<img src="https://latex.codecogs.com/gif.latex?P(\boldsymbol&space;b&space;|&space;\boldsymbol&space;A)&space;=&space;\frac{P(\boldsymbol&space;A|\boldsymbol\theta,&space;\boldsymbol&space;b)P(\boldsymbol\theta,&space;\boldsymbol&space;b)}{P(\boldsymbol&space;A)}" title="P(\boldsymbol b | \boldsymbol A) = \frac{\sum_{\boldsymbol\theta}P(\boldsymbol A|\boldsymbol\theta, \boldsymbol b)P(\boldsymbol\theta, \boldsymbol b)}{P(\boldsymbol A)}">

Where **P(A|θ,b)** is the probability of obtaining the network _A_ given the partition _b_ and additional parameters _θ_; **P(θ,b)** is the probability of occurrence of the partition _b_ having observed the netwok _A_; **P(A)** is the “model evidence” and it is the same for all possible partitions. Refer to the excellent `graph-tool` documentation `<https://graph-tool.skewed.de/static/doc/demos/inference/inference.html`>_ for more details. Note that maximising this quantity is equivalent to minimizing the entropy

<img src="https://latex.codecogs.com/gif.latex?\Sigma&space;=&space;-\ln&space;P(\boldsymbol&space;A|\boldsymbol\theta,&space;\boldsymbol&space;b)&space;-&space;\ln&space;P(\boldsymbol\theta,&space;\boldsymbol&space;b)" title="\Sigma = -\ln P(\boldsymbol A|\boldsymbol\theta, \boldsymbol b) - \ln P(\boldsymbol\theta, \boldsymbol b)" />

The nested model introduces a hierarchy of priors used to infer the optimal recursive grouping of single cell groups. If you are familiar with Leiden or Louvain methods to find cell groups, you may think at this multilevel approach as a multiresolution one, except that it is not. Here, not only the cell groups at each hierarchy level are found maximising the equation above, but the hierarchy itself (hence the groups of groups) is part of the model.
Since there may be more than one fit with similar probability, `schist` uses the `graph-tool` routines to apply a Markov chain Monte Carlo sampling of the posterior distribution aiming to converge to the best model. 
One of the main limitations of `schist` is that it requires significantly more time than any other state of the art approach to identify cell groups. This cost comes with the benefit that it is possible to choose between different parameters according to the likelihood of a partition set to be found over a network.

Fast model vs standard approach
*******************************

In the standard approach, the model is initialized by minimizing the description length (entropy). This requires extra time but, in general, returns better results. It is possible to achieve good results with lower memory footprint and in shorter times settting `fast_model`:

::
	nested_model(adata, fast_model=True)

This will seed the model with a dummy partition scheme, then a greedy merge-split MCMC will explore solutions until it converges.

Marginals
*********

When using the following invocation 

::

	nested_model(adata, collect_marginals=True, equilibrate=True)

an additional step (with fixed number of iterations) is added to execution. During this step, `schist` collects the probability of having a certain number of groups for each level of the hierarchy. These are stored into `adata.uns['nsbm']['group_marginals']`:

::

	level = 2
	S = adata.uns['nsbm']['group_marginals'][level].sum()
	p = adata.uns['nsbm']['group_marginals'][level] / S
	ng = range(1, len(p) + 1)
	bar(ng, p)
	xticks(ng)
	xlabel('Number of groups')
	ylabel('Probability')
	title(f'Group marginals for level {level}')

.. image:: ./figures/gm.png
    :width: 600px
    :align: center
    :height: 600px
    :alt: alternate text
