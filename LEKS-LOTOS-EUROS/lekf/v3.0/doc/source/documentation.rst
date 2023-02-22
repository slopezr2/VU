.. Documentation description.

.. Label between '.. _' and ':' ; use :ref:`text <label>` for reference
.. _documentation:


*************
Documentation
*************


Generation
==========

The documentation of LEKF that you are reading right now is generated
from LEKF source files.
To (re)generate it, run in the main directory::

  make docu

The documentation files are formatted for use with 'Sphinx'; see:

* `Sphinx homepage <http://sphinx-doc.org/index.html>`_

The formatting of documentation is done using 'reStructuredText', see:

* `reStructuredText Primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_

The 'Sphink' python package should be available to be able to (re)generate the documentation.
If not available on your system yet, see below for installation instructions.


Source files
============

The main files of the documentation are located in::

  ./doc/source/

Parts of the documentation are included in the various Python modules and scripts.
The documentation is configured in such a way that these are automatically incorporated.

The source tree of the documentation files is:

  * `index.rst <../../source/index.rst>`_
  * `history.rst <../../source/history.rst>`_
  * `documentation.rst <../../source/documentation.rst>`_


Installation of Sphinx
======================

To check if Sphinx is already installed, try to run the quick-start script::

  sphinx-quickstart --version
  
  sphinx-quickstart 2.4.0

If not available yet, try to install using one of the following methods.

Install Sphinx
--------------

Try if it is possible to download and install Sphinx using 
one of the standard installation commands.

* For an `Anaconda <https://www.anaconda.com/>`_ distribution, try::

    conda install -c anaconda sphinx
    
* For other distributions, try::

    pip install sphinx

Build Spinx from source 
-----------------------

Download the latest version from the
`Sphinx homepage <http://sphinx-doc.org/index.html>`_,
for example::

  Sphinx-3.0.3.tar.gz

To let the package be installed at a local destination, 
define the following environment variable::

  export PYTHON_PREFIX="${HOME}/opt/Python-3.7"

Build and install in the user-defined location using::

  cd Sphinx-3.0.3
  python setup.py install --user

Extend the search path with::

  export PATH="${PYTHON_PREFIX}/bin:${PATH}"


Configuration
=============

The source of the documentation is located in (subdirectories of):

  ./doc/

In this directory, the Sphinx quick start command has been called::

    sphinx-quickstart \
      --sep \
      --project='LEKF' \
      --author='Arjo Segers' \
      --ext-autodoc \
      --ext-intersphinx \
      --ext-mathjax \
      --no-batchfile

The following entities have been created now:

* ``source`` directory to hold the documenation source files;
  initially the following files are created:
  
  * `conf.py <../../source/conf.py>`_ : configuration file;
  * `index.rst <../../source/index.rst>`_ : source of the main page of the documentation;

* ``build`` directory to hold the generated documenation;

* ``Makefile`` : make commands.

In the ``./doc/source/conf.py`` file created in this way,
the following changes were made manually:

* the location of the python modules was added to the search path;
* added options for `autodoc` entries;
* added options for `intersphinx` entries;
* the html theme was changed.

The initial documentation could be created using:

    (cd doc; make html)


