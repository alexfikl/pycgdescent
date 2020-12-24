r"""
Simple example using the low level bindings to CG_DESCENT.

The function and gradient are

.. math::

        \begin{aligned}
        f(\mathbf{x}) = & \sum_{i = 0}^n e^{x_i} - \sqrt{i + 1} x_i, \\
        \nabla_i f(\mathbf{x}) = & e^{x_i} - \sqrt{i + 1}
        \end{aligned}

When the line search routine tries to bracket a local minimizer in
the search direction, it may expand the initial line search interval.
The default expansion factor is :math:`5`. You can modify this factor using
the parameter ``rho``.  In the following example, we choose a small initial
step size (initial step is ``1.0e-5``), ``QuadStep`` is *False*, and ``rho`` is
``1.5``.

The code has to do a number of expansions to reach a suitable
interval bracketing the minimizer in the initial search direction.
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
    param.QuadStep = 0
    param.step = 1.0e-5
    # param.PrintParms = 1

    # }}}

    # {{{

    print("==== with rho 1.5 ====")
    with timer():
        param.rho = 1.5
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
    print("==== with rho 5.0 ====")
    with timer():
        param.rho = 5.0
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
