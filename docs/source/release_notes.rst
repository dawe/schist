=============
Release notes
=============

.. highlight:: python

-------------
version 0.10.0
-------------
This update introduces support for ``graph-tool`` 3, which is now the required version. Due to recent updates to ``graph-tool`` it was impossible to make ``schist`` compatible with elder versions. This update introduces some modifications:

    - If marginals are computed, ``n_samples`` must be specified in place of ``n_init``
    - Due to threading policy in ``graph-tool``, setting a ``random_seed`` causes single thread execution of the minimization step. 
    - It is now possible to set a constraint to the model, so that cells that belong to different categories (``constraint_key``) are never grouped together.
    - Different models are not accessed with ``model`` parameter, it is now necessary to set ``nested`` or ``assortative`` flags.

-------------
version 0.9.4
-------------
This is a minor update. It fixes setting random seed for reproducibile cell marginals.

-------------
version 0.9.3
-------------
This is a minor update. It fixes setting random seed for reproducibile modeling.

-------------
version 0.9.2
-------------
This is a minor update. It fixes dumping of state objects and also removes a legacy dependency for numpy <2.

-------------
version 0.9.1
-------------
This version addresses a vicious bug that prevented convergence during minimization.

-------------
version 0.9.0
-------------
In addition to various bugfixes, this version unifies the inteference functions into one (`fit_model`), so that it's easier to use and maintain. Similarly, `fit_model_multi` works for multimodal data, now supporting `sbm` in addition to `nsbm`. Multimodal analysis now supports `MuData` objects.

-------------
version 0.8.4
-------------
This is a minor update over version 0.8.3. In fixes some `scanpy` compatibility issues and limits the version of `graph-tool` to 2.59. We have noticed that `gt>=2.60` crashes on `osx-arm64` platform when working with multimodal data. 