<img src='garnet.png' alt='logo' width="200" height="200">

<img src="https://travis-ci.org/dawe/scNSBM.svg?branch=master" title="Build status"> 



# SCHIST
An interface for Nested Stochastic Block Model for single cell analysis.

## Status
- Version 0.7.0 Some changes in the way nested model is calculated. Also, advanced functions have been eliminated, one should use `graph_tool` for that
- Version 0.6.0 PlantedModel fully supported. Also, some changes in the data structure (affinities in obsm)
- Version 0.5.2 Fast model for NSBM (and PP).
- Version 0.5.0 Introduces Planted Partition model
- Version 0.4.2 Introduces functions to calculate cluster consistency.
- Version 0.4.0 Some bugfixes and name change to `schist`
- Version 0.3.1 It is possible to resume operatons on a previously modeled state (i.e. you can do two step analysis)
- Version 0.3.0 Major changes in the way MCMC is performed. SBM is also available.
- Version 0.2.5 extends `scanpy.tools.draw_graph` layouts to the ones available in `graph-tool` library
- Version 0.2 The code for the interface is mature. It is also possible to read/write objects
- Version 0.1 implements a working interface to `graph-tool` library for `scanpy`


## How to use
Once `schist` has been installed, it can be used out of the box on `scanpy` objects:

```python
from schist.inference import nested_model

nested_model(adata)
```

Once the MCMC has converged, the `adata.obs` object will contain additional columns for multiple levels, named `nsbm_level_1`, `nsbm_level_2` and so on (by default up to `nsbm_level_10`). 
More details can be found at [this page](Advanced.md)

## Installation
The key component (`graph-tool`) is not available through pip and requires extra compilation by the user, refer to its [installation page](https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions). Note, however, that a conda package is available from conda-forge, that means you may install it (and `schist` dependencies) just issuing

```
conda create -n schist -c conda-forge -c bioconda numpy scipy anndata pandas 'graph-tool>=2.35' scanpy
conda activate schist
```

After that, `schist` can be installed from source:

```
git clone https://github.com/dawe/schist.git
cd schist
pip install .
```


## Known issues
### Saving objects
`schist` allows to return the `NestedBlockState` object in `adata.uns['nsbm']['state']` slot. Unfortunately, this object cannot be dumped into `.h5ad` files by the `sc.write()` function. If you returned the state, e.g. for debugging, you should pop it out from your dataset before writing:

```python
state = adata.uns['nsbm'].pop('state')
adata.write('myfile.h5ad')

# save the state separately
import pickle
with open('state.pkl', 'wb') as pkl_state:
    pickle.dump(state, pkl_state, 2)
```

Since version 0.2 it is possible to save `AnnData` objects like above simply issuing

```python
schist.io.write(adata, prefix='myfile')
```

This will create two files: `myfile.h5ad`, containing the actual `AnnData`, and 
`myfile.pkl` containing the pickled state. With the same logic, assuming the two files
are in the same place, issuing

```python
adata = schist.io.read('myfile')
```

will read the `.h5ad` and the `.pkl` files and create the proper `AnnData` object

## Cite
We are preparing the manuscript. In the meantime, if you use `schist` you may cite the preprint:

```
@article{morelli_2020,
title = {Nested stochastic block models applied to the analysis of single cell data},
author = {Morelli, Leonardo and Giansanti, Valentina and Cittaro, Davide},
url = {http://biorxiv.org/lookup/doi/10.1101/2020.06.28.176180},
year = {2020},
month = {jun},
day = {29},
urldate = {2020-07-02},
journal = {BioRxiv},
doi = {10.1101/2020.06.28.176180},
}
```


## Name
`schist` is a [type of rock](https://en.wikipedia.org/wiki/Schist). Previous name for this project was `scNSBM`, which was hard to pronounce and caused typos when writing (`scnbsm` or `scbsnm` and so on…). We looked for a name which should have "single cell" in it (sc), something about the stochastic model (st) and something about the hierarchy (hi). That's were `schist` comes from. 
