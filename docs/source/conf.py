
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

#-----------------------------------------------------------------------------
# Trying to save autodoc generation on readthedocs, when C packages are imported

#if 'READTHEDOCS' not in os.environ:
#    import cython_generated_ext

import mock
#MOCK_MODULES = ['numpy', 'numpy.dtype','pandas', 'anndata','joblib','scanpy','natsort','scanpy.tools._utils_clustering','scanpy._utils','scanpy._compat','scipy','scipy.sparse','scipy.stats', 'matplotlib', 'matplotlib.patches', 'scipy.interpolate','schist','matplotlib.axes','matplotlib.pyplot','sklearn.preprocessing','graph_tool','graph_tool.all','numba','graph-tool','sklearn.utils.murmurhash','sklearn']
MOCK_MODULES = ['graph-tool','numpy']
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock.Mock()
#------------------------------------------------------------------------------
#sys.path.insert(0, os.path.abspath('..'))
#sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'schist'
copyright = '2021, Morelli Leonardo, Giansanti Valentina, Cittaro Davide'
author = 'Morelli Leonardo, Giansanti Valentina, Cittaro Davide'

# The full version, including alpha/beta/rc tags
release = '2020'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc','sphinx.ext.intersphinx', 'sphinx.ext.coverage', 'sphinx.ext.napoleon', 'sphinx.ext.duration',
              'sphinx.ext.doctest','sphinx.ext.autosummary','sphinx.ext.extlinks','sphinx.ext.viewcode']


numpydoc_show_class_members = False


# generate autosummary even if no references
autosummary_generate = True
autosummary_imported_members = True

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'python': ('https://docs.python.org/3', None),
                       'numpy': ('https://docs.scipy.org/doc/numpy', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
                       'matplotlib': ('https://matplotlib.org', None),
                       'cairo': ('https://www.cairographics.org/documentation/pycairo/3', None),
                       'ipython': ('https://ipython.org/ipython-doc/stable/', None),
                       'pandas': ('https://pandas.pydata.org/pandas-docs/stable/', None)}


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

#html_theme = 'classic'
#html_theme_options = {
#    "rightsidebar": "false",
#    "stickysidebar": "true",
#    "collapsiblesidebar":"false",
#    "externalrefs":"true"
#}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []
