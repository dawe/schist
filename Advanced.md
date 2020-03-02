# Advanced usage

This page explains some advanced usage options for scNSBM

## Understanding NSBM

Foo bar baz

## Plotting

### Alluvial Plots

scNSBM provides an interface to `graph-tool` to infer Nested Stochastic Block Models from single cell data in `scanpy`. Once models are built, data are partitioned in multiple groups, linked together in hierarchical way. In order to represent a hierarchy, scNSBM implements a simple plot function that represents data using alluvial plots:

```python
adata = scnsbm.io.read('adata')
scnsbm.pl.alluvial(adata)
```
<img src="figures/alluvial_01.png" width=150>

This function will plot all levels in hierarchy by default. As many level are uninformative, they can be excluded from the plot

```python
adata = scnsbm.io.read('adata')
scnsbm.pl.alluvial(adata, level_end=5)
```
![alluvial_02](figures/alluvial_02.png | width=150)

Leaf levels can be also excluded

```python
adata = scnsbm.io.read('adata')
scnsbm.pl.alluvial(adata, level_end=5, level_start=2)               
```
![alluvial_03](figures/alluvial_03.png | width=150)

### Extending `sc.tl.draw_grap()`

`graph-tools` has built-in functionalities to plot graphs. Some of these have been implemented into scNSBM using a syntax compatibile with `scanpy`'s functions. Note that Fruchterman-Reingold spring-block layout is already implemented into `scanpy`, and it gives the same output. 

```python
adata = scnsbm.io.read('adata')
scnsbm.tl.draw_graph(adata, layout='fr') 
sc.pl.draw_graph(adata, layout='fr', color='nsbm_level_2', legend_loc='on data')     
```

![fr_01](figures/fr_01.png | width=150)

However, scNSBM allows to seed the plot using the graph tree.

```python
adata = scnsbm.io.read('adata')
scnsbm.tl.draw_graph(adata, layout='fr', use_tree=True) 
sc.pl.draw_graph(adata, layout='fr', color='nsbm_level_2')
```

![fr_02](figures/fr_02.png | width=150)

Default layout is SFDP spring-block layout

```python
adata = scnsbm.io.read('adata')
scnsbm.tl.draw_graph(adata)
sc.pl.draw_graph(adata, layout='sfdp', color='nsbm_level_2', legend_loc='on data')     
```

![sfdp_01](figures/sfdp_01.png | width=150)
