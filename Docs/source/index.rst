.. PeleLM documentation master file, created by
   sphinx-quickstart on Sat Oct 20 12:17:48 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the documentation for `PeleLM`
=========================================

`PeleLM` is an adaptive-mesh low Mach number hydrodynamics code for reacting flows.  `PeleLM` has an
official project `homepage <https://amrex-combustion.github.io/PeleLM/>`_, and can be obtained via
`GitHub <https://github.com/AMReX-Combustion/PeleLM>`_.  If you need help or have questions, please join the users `forum <https://groups.google.com/forum/#!forum/pelelmusers>`_. The documentation pages appearing here
are distributed with the code in the ``Docs`` folder as "restructured text" files.  The html is built
automatically with certain pushes to the `PeleLM` GibHub repository and are maintained online by 
`ReadTheDocs <https://pelelm.readthedocs.io/en/latest>`_.  A local version can also be built as follows ::

    cd ${PELELM_DIR}/Docs
    make html

where ``PELELM_DIR`` is the location of your clone of the `PeleLM` repository.  To view the local pages,
point your web browser at the file ``${PELELM_DIR}/Docs/build/html/index.html``.

**Current docs build status on ReadTheDocs:**

.. image:: https://readthedocs.org/projects/pelelm/badge


.. toctree::
   :maxdepth: 2
   :caption: Documentation contents:

   GettingStarted.rst
   Model.rst
   ProblemSetup.rst
   GNUmakeSystem.rst
   RunningPeleLM.rst
   Visualization.rst
   Tutorials.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
