<img src='garnet.png' alt='logo' width="200" height="200">

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/schist/badges/version.svg)](https://anaconda.org/conda-forge/schist)

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/schist/badges/downloads.svg)](https://anaconda.org/conda-forge/schist)



# SCHIST
An interface for Nested Stochastic Block Model for single cell analysis. The idea behind `schist` is to create a `scanpy`-compatible interface to `graph-tool` methods.

## Installation
As `schist` is available from conda-forge, installation is just easy as

```
conda install -c conda-forge schist
```

This will install `graph-tool-base` package as dependency, which does not include complete graphical capabilities. If you need access to the full `graph-tool` infrastracture, install it with

```
conda install -c conda-forge graph-tool
```

If you want to install `schist` from source, clone this repository and install in the usual way

```
git clone https://github.com/dawe/schist.git
cd schist
pip install .
```

## How to use
Once `schist` has been installed, it can be used out of the box on `scanpy` objects:

```python
import schist as scs

scs.inference.nested_model(adata)
```

Once the MCMC has converged, the `adata.obs` object will contain additional columns for multiple levels, named `nsbm_level_0`, `nsbm_level_1`, `nsbm_level_2` and so on (by default up to `nsbm_level_10`). 
More details can be found at [this page](Advanced.md)


## Known issues
### Saving objects
`schist` allows to return the `NestedBlockState` object in `adata.uns['nsbm']['state']` slot. Unfortunately, this object cannot be dumped into `.h5ad` files by the `sc.write()` function. If you returned the state, e.g. for debugging, you should pop it out from your dataset before writing:

```python
state = adata.uns['schist'].pop('state')
adata.write('myfile.h5ad')

# save the state separately
import pickle
with open('state.pkl', 'wb') as pkl_state:
    pickle.dump(state, pkl_state, 2)
```

Since version 0.2 it is possible to save `AnnData` objects like above simply issuing

```python
scs.io.write(adata, prefix='myfile')
```

This will create two files: `myfile.h5ad`, containing the actual `AnnData`, and 
`myfile.pkl` containing the pickled state. With the same logic, assuming the two files
are in the same place, issuing

```python
adata = scs.io.read('myfile')
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
