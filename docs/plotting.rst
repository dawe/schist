********
Advance Plotting
********

Alluvial Plots
##############

`schist` provides an interface to `graph-tool` to infer Nested Stochastic Block Models from single cell data in `scanpy`. Once models are built, data are partitioned in multiple groups, linked together in hierarchical way. In order to represent a hierarchy, `schist` implements a simple plot function that represents data using alluvial plots:

::

	adata = schist.io.read('adata')
	schist.pl.alluvial(adata)

.. image:: ./figures/alluvial_01.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text


This function will plot all levels in hierarchy by default. As many level are uninformative, they can be excluded from the plot

::

	adata = schist.io.read('adata')
	schist.pl.alluvial(adata, level_end=5)

.. image:: ./figures/alluvial_02.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text

Leaf levels can be also excluded

::

	adata = schist.io.read('adata')
	schist.pl.alluvial(adata, level_end=5, level_start=2)               

.. image:: ./figures/alluvial_03.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text

Extending `sc.tl.draw_graph()`
##############################

`graph-tools` has built-in functionalities to plot graphs. Some of these have been implemented into `schist` using a syntax compatibile with `scanpy`'s functions. Note that Fruchterman-Reingold spring-block layout is already implemented into `scanpy`, and it gives the same output. 

::

	adata = schist.io.read('adata')
	schist.tl.draw_graph(adata, layout='fr') 
	sc.pl.draw_graph(adata, layout='fr', color='nsbm_level_2', legend_loc='on data')     

.. image:: ./figures/fr_01.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text

However, `schist` allows to seed the plot using the graph tree.

::

	adata = schist.io.read('adata')
	schist.tl.draw_graph(adata, layout='fr', use_tree=True) 
	sc.pl.draw_graph(adata, layout='fr', color='nsbm_level_2')

.. image:: ./figures/fr_02.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text

Default layout is SFDP spring-block layout

::

	adata = schist.io.read('adata')
	schist.tl.draw_graph(adata)
	sc.pl.draw_graph(adata, layout='sfdp', color='nsbm_level_2', legend_loc='on data')     

.. image:: ./figures/sfdp_01.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text

With tree information

::

	adata = schist.io.read('adata')
	schist.tl.draw_graph(adata, use_tree=True)
	sc.pl.draw_graph(adata, layout='sfdp', color='nsbm_level_2')     

.. image:: ./figures/sfdp_02.png
    :width: 500px
    :align: center
    :height: 400px
    :alt: alternate text
