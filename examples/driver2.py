r"""
Simple example using the low level bindings to CG_DESCENT.

The function and gradient are

.. math::

        \begin{aligned}
        f(\mathbf{x}) = & \sum_{i = 0}^n e^{x_i} - \sqrt{i + 1} x_i, \\
        \nabla_i f(\mathbf{x}) = & e^{x_i} - \sqrt{i + 1}
        \end{aligned}

The operation of the code is mostly controlled by the parameters
in the ``cg_parameter`` structure.  In the following example,
the parameter ``QuadStep`` is set to *False*.  When ``QuadStep`` is *True*,
the trial step in each iteration is computed as the minimizer
of a quadratic interpolant along the search direction. In
performing the quad step, we hope to find a suitable line search
point right away, completely by-passing the secant iteration.

However, as the iterates approach a minimizer, the numerical accuracy of
the minimizer of the quadratic interpolant becomes worse. When the relative
change in the function values for two consecutive iterations reach
``QuadCutOff``, then the code completely turns off the quad step. The user
can turn off the quad step by setting ``QuadStep`` to *False*. By leaving
``QuadStep`` *True*, but increasing ``QuadCutOff`` (default ``1.0e-12``), the
code turns off the ``QuadStep`` sooner.

Below, we run the code twice, first with the ``QuadStep`` turned off,
then with the ``QuadStep`` turned on. Notice that the performance improves
with the ``QuadStep`` is on. This behavior is typical.
"""

import time
from contextlib import contextmanager

import numpy as np
import pycgdescent._private as _cg


@contextmanager
def timer():
    t_start = time.time()
    yield
    t_end = time.time()
    print("elapsed: ", t_end - t_start)


def fn(x):
    f = np.sum(np.exp(x) - t * x)
    return f


def grad(g, x):
    g[...] = np.exp(x) - t


def fngrad(g, x):
    y = np.exp(x)
    f = np.sum(y - t * x)
    g[...] = y - t
    return f


def main(n=100):
    # {{{ parameters

    global t
    x0 = np.ones(n, dtype=np.float64)
    t = np.sqrt(1 + np.arange(n))

    param = _cg.cg_parameter()
    # param.PrintParms = 1

    # }}}

    # {{{

    print("==== with QuadStep OFF ====")
    with timer():
        param.QuadStep = 0
        _, stats, _ = _cg.cg_descent(x0, 1.0e-8, param, fn, grad, fngrad, None)

    print()
    print("maximum norm for gradient: %+.16e" % stats.gnorm)
    print("function value:            %+.16e" % stats.f)
    print("cg iterations:            ", stats.iter)
    print("function evaluations:     ", stats.nfunc)
    print("gradient evaluations:     ", stats.ngrad)

    # }}}

    # {{{

    print()
    print("==== with QuadStep ON ====")
    with timer():
        param.QuadStep = 1
        _, stats, _ = _cg.cg_descent(x0, 1.0e-8, param, fn, grad, fngrad, None)

    print()
    print("maximum norm for gradient: %+.16e" % stats.gnorm)
    print("function value:            %+.16e" % stats.f)
    print("cg iterations:            ", stats.iter)
    print("function evaluations:     ", stats.nfunc)
    print("gradient evaluations:     ", stats.ngrad)

    # }}}


if __name__ == "__main__":
    t = None
    main()

# vim: fdm=marker
