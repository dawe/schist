# Advanced usage

This page explains some advanced usage options for scNSBM

## Understanding NSBM

scNSBM is based on Nested Stochastic Block Models (NSBM), a generative process based on the notion of group of nodes. Here, we use NSBM tu cluster cells in scRNA-seq experiments.

### Bayesian approach
The first step of a network analysis is the subdivision of nodes (or cells) into communities. In our situation this is achieved via Bayesian posterior probability, where **P(b|A)** is the likehood of generating a network A with a partition b and it is obtained according to 

<img src="https://latex.codecogs.com/gif.latex?P(\boldsymbol&space;b&space;|&space;\boldsymbol&space;A)&space;=&space;\frac{\sum_{\boldsymbol\theta}P(\boldsymbol&space;A|\boldsymbol\theta,&space;\boldsymbol&space;b)P(\boldsymbol\theta,&space;\boldsymbol&space;b)}{P(\boldsymbol&space;A)}" title="P(\boldsymbol b | \boldsymbol A) = \frac{\sum_{\boldsymbol\theta}P(\boldsymbol A|\boldsymbol\theta, \boldsymbol b)P(\boldsymbol\theta, \boldsymbol b)}{P(\boldsymbol A)}">



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
