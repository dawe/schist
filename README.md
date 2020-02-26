# scNSBM
An interface for Nested Stochastic Block Model for single cell analysis.

## Status
- Version 0.2 The code for the interface is mature. It is also possible to read/write objects
- Version 0.1 implements a working interface to `graph-tool` library for `scanpy`


## How to use
Once scNSBM has been installed, it can be used out of the box on `scanpy` objects:

```python
from scnsbm.inference import nested_model

nested_model(adata)
```

Once the MCMC has converged, the `adata.obs` object will contain additional columns for multiple levels, named `nsbm_level_1`, `nsbm_level_2` and so on (by default up to `nsbm_level_10`). 
More details can be found at [this page](Advanced.md)

## Installation
scNSBM is not (yet) available on PyPI. You have to install from source as:

```
git clone https://github.com/dawe/scNSBM.git
cd scNSBM
pip install .
```

The following dependencies will be managed by pip itself

- `numpy`
- `scipy`
- `anndata`
- `pandas`
- `scanpy`

Alas, the key component (`graph-tool`) is not available through pip and requires extra compilation by the user, refer to its [installation page](https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions). Note, however, that a conda package is available from conda-forge, that means you may install it just issuing

```
conda install -c conda-forge graph-tool
```

## Known issues
### Cairo interface
`graph-tool` requilres `Gtk` to plot graphs. We do not plan to use those capabilities natively. This means that you may safely disregard the following warning:

```python
graph_tool/draw/cairo_draw.py:1494: RuntimeWarning: Error importing Gtk module: No module named 'gi'; GTK+ drawing will not work.
  warnings.warn(msg, RuntimeWarning)
```

### Saving objects
scNSBM allows to return the `NestedBlockState` object in `adata.uns['nsbm']['state']` slot. Unfortunately, this object cannot be dumped into `.h5ad` files by the `sc.write()` function. If you returned the state, e.g. for debugging, you should pop it out from your dataset before writing:

```python
state = adata.uns['nsbm'].pop('state')
adata.write('myfile.h5ad')

# save the state separately
import pickle
with open('state.pkl', 'wb') as pkl_state:
    pickle.dump(state, pkl_state, 2)
```
