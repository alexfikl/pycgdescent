pycgdescent documentation
=========================

Example
-------

A simple example used to optimize the classic
`Rosenbrock function <https://en.wikipedia.org/wiki/Rosenbrock_function>`__
is provided below. The function is given by

.. math::

    f(x, y) = a (y - x^2)^2 + b (x - 1)^2,

which has an optimum at :math:`[1, 1]` for :math:`a > 0` and :math:`b > 0`.
The problem is setup with default parameters in the following example.

.. literalinclude:: ../examples/high-level-api.py
   :start-after: START_ROSENROCK_EXAMPLE
   :end-before: END_ROSENBROCK_EXAMPLE
   :language: python
   :linenos:

The output of ``r.pretty()`` shows that we have found the exact solution::

                fun : 6.019745113421725e-23
                jac : 3.0587488097952063e-10
            message : 'Convergence tolerance satisfied'
               nfev : 63
                nit : 28
               njev : 35
        nsubspaceit : 0
         nsubspaces : 0
             status : 0
            success : True
                  x : array([1., 1.])

The path of the optimization can be seen in the following figure.

.. image:: _static/rosenbrock-dark.png
    :class: only-dark
    :width: 75%
    :align: center
    :alt: Rosenbrock function optimization

.. image:: _static/rosenbrock-light.png
    :class: only-light
    :width: 75%
    :align: center
    :alt: Rosenbrock function optimization

Reference
---------

.. automodule:: pycgdescent

Changelog
=========

.. include:: ../CHANGELOG.rst

Acknowledgments
===============

Work on ``pycgdescent`` was sponsored, in part, by the Office of Naval Research (ONR)
as part of the Multidisciplinary University Research Initiatives (MURI)
Program, under Grant Number N00014-16-1-2617.

The views and opinions expressed herein do not necessarily reflect those of the
funding agencies.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
