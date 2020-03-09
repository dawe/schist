# Advanced usage

This page explains some advanced usage options for scNSBM. 

- [Understanding NSBM](#Understanding-NSBM)
- [Plotting](#Plotting)

## Understanding NSBM

scNSBM is based on Nested Stochastic Block Models (NSBM), a generative process based on the notion of group of nodes. Here, we use NSBM to cluster cells in scRNA-seq experiments.

### Bayesian approach
Cell to cell relations are commonly represented as a neighborhood graph (typically KNN or  SNN), cell groups are identified as graph communities, this step is usually performed maximising graph modularity. In scNSBM, instead, partitioning is performed maximising the  Bayesian posterior probability **P(b|A)**, that is the likehood of generating a network _A_ with a partition _b_ and it is obtained according to 

<img src="https://latex.codecogs.com/gif.latex?P(\boldsymbol&space;b&space;|&space;\boldsymbol&space;A)&space;=&space;\frac{P(\boldsymbol&space;A|\boldsymbol\theta,&space;\boldsymbol&space;b)P(\boldsymbol\theta,&space;\boldsymbol&space;b)}{P(\boldsymbol&space;A)}" title="P(\boldsymbol b | \boldsymbol A) = \frac{\sum_{\boldsymbol\theta}P(\boldsymbol A|\boldsymbol\theta, \boldsymbol b)P(\boldsymbol\theta, \boldsymbol b)}{P(\boldsymbol A)}">

Where **P(A|θ,b)** is the probability of obtaining the network _A_ given the partition _b_ and additional parameters _θ_; **P(θ,b)** is the probability of occurrence of the partition _b_ having observed the netwok _A_; **P(A)** is the “model evidence” and it is the same for all possible partitions. Refer to the excellent [`graph-tool` documentation](https://graph-tool.skewed.de/static/doc/demos/inference/inference.html) for more details. Note that maximising this quantity is equivalent to minimizing the entropy

<img src="https://latex.codecogs.com/gif.latex?\Sigma&space;=&space;-\ln&space;P(\boldsymbol&space;A|\boldsymbol\theta,&space;\boldsymbol&space;b)&space;-&space;\ln&space;P(\boldsymbol\theta,&space;\boldsymbol&space;b)" title="\Sigma = -\ln P(\boldsymbol A|\boldsymbol\theta, \boldsymbol b) - \ln P(\boldsymbol\theta, \boldsymbol b)" />


The nested model introduces a hierarchy of priors used to infer the optimal recursive grouping of single cell groups. If you are familiar with Leiden or Louvain methods to find cell groups, you may think at this multilevel approach as a multiresolution one, except that it is not. Here, not only the cell groups at each hierarchy level are found maximising the equation above, but the hierarchy itself (hence the groups of groups) is part of the model.
Since there may be more than one fit with similar probability, scNSBM uses the `graph-tool` routines to apply a Markow chain Monte Carlo sampling of the posterior distribution aiming to converge to the best model. 
One of the main limitations of scNSBM is that it requires significantly more time than any other state of the art approach to identify cell groups. This cost comes with the benefit that it is possible to choose between different parameters according to the likelihood of a partition set to be found over a network. 

### Fast model vs standard approach

In the standard approach, the model is initialized by minimizing the description length (entropy) before running the MCMC algorithm. This requires extra time but, in general, returns better results. It is possible to skip the initial minimization and the sweep step of MCMC setting 

```python
nested_model(adata, fast_model=True)
```

This will seed the model with a dummy description of the graph where every cell belongs to its own partition. This approach is generally faster and requires less memory, but it could be less precise

### Marginals

When using the following invocation 

```python
nested_model(adata, collect_marginals=True)
```

an additional step (with fixed number of iterations) is added to execution. During this step, scNSBM collects two types of marginals that can be used to understand the configuration of the single cell experiment. Cell marginals, that is the probability of a cell to be assigned to a group, are stored into `adata.uns['nsbm']['cell_marginals']`. Here, a dictionary keyed with NSBM levels counts the times a cell has been successfully moved to a group. These marginals can be efficiently used as covariates when looking for marker genes, this approach will weight the belief that a cell belongs to a group. We have prepared a [notebook](https://github.com/dawe/scNSBM-notebooks/blob/master/Cell_Marginals.ipynb) showing an example. 

In addition to cell marginals, scNSBM collects the probability of having a certain number of groups for each level of the hierarchy. These are stored into `adata.uns['nsbm']['group_marginals']`:

```python
level = 2
S = adata.uns['nsbm']['group_marginals'][level].sum()
p = adata.uns['nsbm']['group_marginals'][level] / S
ng = range(1, len(p) + 1)
bar(ng, p)
xticks(ng)
xlabel('Number of groups')
ylabel('Probability')
title(f'Group marginals for level {level}')
```

<img src="docs/figures/gm.png" width=400>


## Plotting

### Alluvial Plots

scNSBM provides an interface to `graph-tool` to infer Nested Stochastic Block Models from single cell data in `scanpy`. Once models are built, data are partitioned in multiple groups, linked together in hierarchical way. In order to represent a hierarchy, scNSBM implements a simple plot function that represents data using alluvial plots:

```python
adata = scnsbm.io.read('adata')
scnsbm.pl.alluvial(adata)
```
<img src="docs/figures/alluvial_01.png"  width=400>

This function will plot all levels in hierarchy by default. As many level are uninformative, they can be excluded from the plot

```python
adata = scnsbm.io.read('adata')
scnsbm.pl.alluvial(adata, level_end=5)
```
<img src="docs/figures/alluvial_02.png" width=400>

Leaf levels can be also excluded

```python
adata = scnsbm.io.read('adata')
scnsbm.pl.alluvial(adata, level_end=5, level_start=2)               
```
<img src="docs/figures/alluvial_03.png" width=400>

### Extending `sc.tl.draw_graph()`

`graph-tools` has built-in functionalities to plot graphs. Some of these have been implemented into scNSBM using a syntax compatibile with `scanpy`'s functions. Note that Fruchterman-Reingold spring-block layout is already implemented into `scanpy`, and it gives the same output. 

```python
adata = scnsbm.io.read('adata')
scnsbm.tl.draw_graph(adata, layout='fr') 
sc.pl.draw_graph(adata, layout='fr', color='nsbm_level_2', legend_loc='on data')     
```

<img src="docs/figures/fr_01.png" width=400>

However, scNSBM allows to seed the plot using the graph tree.

```python
adata = scnsbm.io.read('adata')
scnsbm.tl.draw_graph(adata, layout='fr', use_tree=True) 
sc.pl.draw_graph(adata, layout='fr', color='nsbm_level_2')
```

<img src="docs/figures/fr_02.png" width=400>

Default layout is SFDP spring-block layout

```python
adata = scnsbm.io.read('adata')
scnsbm.tl.draw_graph(adata)
sc.pl.draw_graph(adata, layout='sfdp', color='nsbm_level_2', legend_loc='on data')     
```

<img src="docs/figures/sfdp_01.png" width=400>

With tree information

```python
adata = scnsbm.io.read('adata')
scnsbm.tl.draw_graph(adata, use_tree=True)
sc.pl.draw_graph(adata, layout='sfdp', color='nsbm_level_2')     
```

<img src="docs/figures/sfdp_02.png" width=400>
