Installation
============

Requirements
------------

* AmberTools: ``tleap`` and ``sander``, plus ``$AMBERHOME/dat/leap/parm/gaff*.dat``
* Python 3.9+ with `ParmEd <https://parmed.github.io/ParmEd/>`_, `RDKit <https://www.rdkit.org/>`_,
  NumPy and SciPy

Install
-------

The easiest way to get everything is conda-forge, whose ``ambertools`` package includes ParmEd:

.. code-block:: bash

   conda create -n unleap -c conda-forge ambertools rdkit
   conda activate unleap
   pip install unleap

If you already have an AmberTools install, source it (``source $AMBERHOME/amber.sh``) and
``pip install unleap`` into that Python.

From a checkout, for development:

.. code-block:: bash

   git clone https://github.com/piskuliche/unleap.git
   cd unleap
   pip install -e ".[test,docs]"

``AMBERHOME`` must be set. If it isn't, unleap stops with ``AMBERHOME is not set``.

Building these docs
-------------------

.. code-block:: bash

   pip install -r docs/requirements.txt
   make -C docs html          # output in docs/_build/html
