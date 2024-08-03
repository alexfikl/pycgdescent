.. |badge-ci| image:: https://github.com/alexfikl/pycgdescent/workflows/CI/badge.svg
    :alt: Build Status
    :target: https://github.com/alexfikl/pycgdescent/actions?query=branch%3Amain+workflow%3ACI

.. |badge-rtd| image:: https://readthedocs.org/projects/pycgdescent/badge/?version=latest
    :alt: Documentation
    :target: https://pycgdescent.readthedocs.io/en/latest/?badge=latest

.. |badge-pypi| image:: https://badge.fury.io/py/pycgdescent.svg
    :alt: PyPI
    :target: https://pypi.org/project/pycgdescent/

.. |badge-reuse| image:: https://api.reuse.software/badge/github.com/alexfikl/pycgdescent
    :alt: REUSE
    :target: https://api.reuse.software/info/github.com/alexfikl/pycgdescent

|badge-ci| |badge-rtd| |badge-pypi| |badge-reuse|

pycgdescent
===========

Python wrapper for the `CG_DESCENT <https://people.clas.ufl.edu/hager/software/>`__
algorithm by Hager and Zang (see `DOI <https://doi.org/10.1145/1132973.1132979>`__).
A previous wrapper can be found `here <https://github.com/martiniani-lab/PyCG_DESCENT>`__.
Some differences:

* This one only depends on `pybind11 <https://github.com/pybind/pybind11>`__.
* Tries to emulate the interface from `scipy <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`__
  (still needs work).

Interesting links:

* `Documentation <https://pycgdescent.readthedocs.io/en/latest/>`__.
* `Code <https://github.com/alexfikl/pycgdescent>`__.
