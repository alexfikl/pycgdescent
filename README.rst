.. image:: https://github.com/alexfikl/pycgdescent/workflows/CI/badge.svg
    :alt: Build Status
    :target: https://github.com/alexfikl/pycgdescent/actions?query=branch%3Amain+workflow%3ACI

.. image:: https://readthedocs.org/projects/pycgdescent/badge/?version=latest
    :alt: Documentation
    :target: https://pycgdescent.readthedocs.io/en/latest/?badge=latest

pycgdescent
===========

Python wrapper for `CG_DESCENT <https://users.clas.ufl.edu/hager/papers/Software/>`__.
A previous wrapper can be found `here <https://github.com/martiniani-lab/PyCG_DESCENT>`__.
Some differences:

* This one only depends on `pybind11 <https://github.com/pybind/pybind11>`__.
* Tries to emulate the interface from `scipy <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`__
  (still needs work).

Interesting links:

* `Documentation <https://pycgdescent.readthedocs.io/en/latest/>`__.
* `Code <https://github.com/alexfikl/pycgdescent>`__.
