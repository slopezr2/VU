# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
#sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'LEKF'
copyright = '2020, Arjo Segers'
author = 'Arjo Segers'

# The full version, including alpha/beta/rc tags
release = 'v3.0.004'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# name of the master doc:
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# ~~ options for autodoc

# Options for autodoc entries, see:
#   https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html

# List of entities that should be produced automatically:
# * 'show-inheritance'   : name of base class;
# * 'members'            : members of modules (classes) and classes (methods);
# * 'inherited-members'  : inherited members:
autodoc_default_options = {
    'members'           :   True,
    'show-inheritance'  :   True
}
# do not sort methods alphabetically, keep source order:
autodoc_member_order = 'bysource'

# ~~ options for intersphinx

# Options for references to python documenations, see:
#   https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html

# add links to existing online documentation if present:
docs_url = 'https://docs.python.org/%i.%i' % (sys.version_info.major,sys.version_info.minor)
intersphinx_mapping = { 'python' : (docs_url,None) }


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'bizstyle'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


