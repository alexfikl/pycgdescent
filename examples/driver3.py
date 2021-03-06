r"""
Simple example using the low level bindings to CG_DESCENT.

The function and gradient are

.. math::

        \begin{aligned}
        f(\mathbf{x}) = & \sum_{i = 0}^n e^{x_i} - \sqrt{i + 1} x_i, \\
        \nabla_i f(\mathbf{x}) = & e^{x_i} - \sqrt{i + 1}
        \end{aligned}

In the line search for first iteration, there is very little information
available for choosing a suitable step size. By default, the code employs
very low order approximations to the function to estimate a suitable
step size. In some cases, this initial step size can be problematic.

For example, if the cost function contains a :math:`\log` function, the initial
step might cause the code to try to compute the :math:`\log` of a negative number.
If the cost function contains an exponential, then the initial step size
might lead to an overflow. In either case, ``NaN``\ s are potentially generated.

If the default step size is unsuitable, you can input the starting
step size using the parameter ``step``. In the following example, the initial
step size is set to 1.
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
    param.step = 1.0
    # param.PrintParms = 1

    # }}}

    # {{{ different step size

    with timer():
        _, stats, status = _cg.cg_descent(x0, 1.0e-8, param, fn, grad, fngrad,
                callback=None, work=None)

    print()
    print("status:  ", status)
    print("message: ", _cg.STATUS_TO_MESSAGE[status])

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
