=============
Release notes
=============

.. highlight:: python

-------------
version 0.9.0
-------------
In addition to various bugfixes, this version unifies the inteference functions into one (`fit_model`), so that it's easier to use and maintain. Similarly, `fit_model_multi` works for multimodal data, now supporting `sbm` in addition to `nsbm`. Multimodal analysis now supports `MuData` objects.

-------------
version 0.8.4
-------------
This is a minor update over version 0.8.3. In fixes some `scanpy` compatibility issues and limits the version of `graph-tool` to 2.59. We have noticed that `gt>=2.60` crashes on `osx-arm64` platform when working with multimodal data. 