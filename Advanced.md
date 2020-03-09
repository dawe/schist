# Advanced usage

This page explains some advanced usage options for scNSBM

## Understanding NSBM

scNSBM is based on Nested Stochastic Block Models (NSBM), a generative process based on the notion of group of nodes. Here, we use NSBM to cluster cells in scRNA-seq experiments.

### Bayesian approach
Cell to cell relations are commonly represented as a neighborhood graph (typically KNN or  SNN), cell groups are identified as graph communities, this step is usually performed maximising graph modularity. In scNSBM, instead, partitioning is performed maximising the  Bayesian posterior probability **P(b|A)**, that is the likehood of generating a network _A_ with a partition _b_ and it is obtained according to 

<img src="https://latex.codecogs.com/gif.latex?P(\boldsymbol&space;b&space;|&space;\boldsymbol&space;A)&space;=&space;\frac{P(\boldsymbol&space;A|\boldsymbol\theta,&space;\boldsymbol&space;b)P(\boldsymbol\theta,&space;\boldsymbol&space;b)}{P(\boldsymbol&space;A)}" title="P(\boldsymbol b | \boldsymbol A) = \frac{\sum_{\boldsymbol\theta}P(\boldsymbol A|\boldsymbol\theta, \boldsymbol b)P(\boldsymbol\theta, \boldsymbol b)}{P(\boldsymbol A)}">

Where **P(A|θ,b)** is the probability of obtaining the network _A_ given the partition _b_ and additional parameters _θ_; **P(θ,b)** is the probability of occurrence of the partition _b_ having observed the netwok _A_; **P(A)** is the “model evidence” and it is the same for all possible partitions. Refer to the excellent [`graph-tool` documentation](https://graph-tool.skewed.de/static/doc/demos/inference/inference.html) for more details. Note that maximising this quantity is equivalent to minimizing the entropy

<img src="https://latex.codecogs.com/gif.latex?\Sigma&space;=&space;-\ln&space;P(\boldsymbol&space;A|\boldsymbol\theta,&space;\boldsymbol&space;b)&space;-&space;\ln&space;P(\boldsymbol\theta,&space;\boldsymbol&space;b)" title="\Sigma = -\ln P(\boldsymbol A|\boldsymbol\theta, \boldsymbol b) - \ln P(\boldsymbol\theta, \boldsymbol b)" />


The nested model introduces a hierarchy of priors used to infer the optimal recursive grouping of single cell groups. If you are familiar with Leiden or Louvain methods to find cell groups, you may think at this multilevel approach as a multiresolution one, except that it is not. Here, not only the cell groups at each hierarchy level are found maximising the equation above, but the hierarchy itself (hence the groups of groups) is part of the model.
Since there may be more than one fit with similar probability, scNSBM uses the `graph-tool` routines to apply a Markow chain Monte Carlo sampling of the posterior distribution aiming to converge to the best model. 


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

### Extending `sc.tl.draw_grap()`

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
