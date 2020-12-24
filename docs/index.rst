Welcome to pycgdescent's documentation!
=======================================

Example
-------

A simple example used to optimize the classic
`Rosenbrock function <https://en.wikipedia.org/wiki/Rosenbrock_function>`__
is provided below. The function is given by

.. math::

    f(x, y) = a (y - x^2)^2 + b (x - 1)^2,

which has an optimum at :math:`[1, 1]` for :math:`a > 0` and :math:`b > 0`.
We start by loading the required libraries

.. code-block:: python

    import numpy as np
    import numpy.linalg as la
    import pycgdescent as cg

Then we set up the problem with default options.

.. literalinclude:: ../tests/test_pycgdescent.py
   :start-after: START_ROSENROCK_EXAMPLE
   :end-before: END_ROSENBROCK_EXAMPLE
   :language: python
   :linenos:
   :dedent: 4

The outout of ``r.pretty()`` shows that we have found the exact solution::

                fun : 3.4584372059703856e-25
                jac : 2.3353763367596195e-11
            message : 'Convergence tolerance satisfied'
               nfev : 81
                nit : 34
               njev : 49
        nsubspaceit : 0
         nsubspaces : 0
             status : 0
            success : True
                  x : array([1., 1.])

Reference
---------

.. automodule:: pycgdescent

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
