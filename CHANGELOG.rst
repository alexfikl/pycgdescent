pycgdescent 0.10.3 (October 12, 2024)
-------------------------------------

**Dependencies**

* Added official support for Python 3.13.

pycgdescent 0.10.2 (August 28, 2024)
------------------------------------

**Fixes**

* Do not build and publish wheels for Python 3.9.

pycgdescent 0.10.1 (August 28, 2024)
------------------------------------

**Fixes**

* Update CI Python version matrix to exclude 3.9.

pycgdescent 0.10.0 (August 28, 2024)
------------------------------------

**Dependencies**

* Bumped minimum Python version to 3.10 to match the newly released numpy 2.1.

**Features**

* Various maintenance updates (bumping CI dependencies, linting, fixing docs).

pycgdescent 0.9.0 (August 3, 2024)
----------------------------------

**Features**

* Publishing wheels to PyPI.
* Tested and working on numpy 2.0.
* Use ``pybind11-stubgen`` to generate the type annotations for the internal
  ``_cg_descent`` module.
* Various maintenance updates (bumping CI dependencies, linting, fixing docs).

**Fixes**

* Fix ``sphinx-autoapi`` doc build.
* Fix the gradient in a test case.

pycgdescent 0.8.1 (February 7, 2024)
------------------------------------

**Fixes**

* Fix the meson build without BLAS and other related issues.

pycgdescent 0.8.0 (February 7, 2024)
------------------------------------

**Dependencies**

* Requires ``meson-python`` for the build-system.

**Features**

* Export the low level wrapper as used in some examples.

pycgdescent 0.7.0 (January 16, 2024)
------------------------------------

**Dependencies**

* Requires ``scikit-build-core`` for the build-system.

**Features**

* Use the ``src`` layout for the code.
* Use ``ruff format`` for all our formatting needs.
* Switch to ``scikit-build-core`` for the build-system. This effectively leaves
  ``setup.py`` as an empty shim now.

pycgdescent 0.6.0 (October 4, 2023)
-----------------------------------

**Dependencies**

* Support Python 3.12.

**Fixes**

* Update README with more links.
* Update `.readthedocs.yml` config: remove deprecated ``system_packages``.
* Fix latest linting and typing errors.

pycgdescent 0.5.0 (June 20, 2023)
---------------------------------

**Dependencies**

* Bump minimum version to Python 3.9.

**Features**

* Switch to ``ruff`` for all static linting (``flake8`` and ``isort``).
* Switch to ``pyproject.toml`` for all configuration.
* Use ``from __future__ import annotations`` and newer annotation formats,
  e.g. ``X | Y``.
* Switched Sphinx theme to ``sphinx-book-theme``.

pycgdescent 0.4.1 (December 20, 2022)
-------------------------------------

**Fixes**

* Update CI dependencies and linting.

pycgdescent 0.4.0 (November 5, 2022)
------------------------------------

**Dependencies**

* Support Python 3.11.

**Fixes**

* Maintenance work: update CI, use REUSE license, dependabot, etc.

pycgdescent 0.3.0 (July 10, 2022)
---------------------------------

**Dependencies**

* Bump minimum version to Python 3.8.

**Features**

* Use ``black`` for formatting.

pycgdescent 0.2.3 (January 2, 2022)
-----------------------------------

**Dependencies**

* Advertise support for Python 3.10.

**Features**

* Support Numpy 1.22 type annotations.

pycgdescent 0.2.2 (September 26, 2021)
--------------------------------------

**Fixes**

* Fix an uninitialized variable in the original CG_DESCENT sources.

pycgdescent 0.2.1 (June 8, 2021)
--------------------------------

**Fixes**

* Fix version bump.

pycgdescent 0.2.0 (June 8, 2021)
--------------------------------

**Features**

* Add a very sketchy patch to limit the maximum step size in the algorithm.
* Add type checking everywhere.

pycgdescent 0.1.0 (December 24, 2020)
-------------------------------------

**Features**

* Initial release!
* A wrapper around the `CG_DESCENT <https://people.clas.ufl.edu/hager/software/>`__
  using `pybind11 <https://github.com/pybind/pybind11>`__.
* Added convenience APIs to bring it closer to
  `scipy.optimize.minimize <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`__
  (although not a drop in replacement at the moment).
