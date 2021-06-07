r"""
Simple example using the low level bindings to CG_DESCENT.

The function and gradient are

.. math::

        \begin{aligned}
        f(\mathbf{x}) = & \sum_{i = 0}^n e^{x_i} - \sqrt{i + 1} x_i, \\
        \nabla_i f(\mathbf{x}) = & e^{x_i} - \sqrt{i + 1}
        \end{aligned}

Although there is a rigorous theory justifying a Wolfe line search,
the performance of the Approximate Wolfe line search is often superior.
Nonetheless, the user can turn off the Approximate Wolfe line search
by setting ``AWolfe`` to *False* and ``AWolfeFac`` to :math:`0`. Since
``AWolfe`` is *False* by default, we only need to adjust ``AWolfeFac``. When
the code detects that the Wolfe line search fails, then it will
automatically attempt the approximate Wolfe line search.

To see that the Wolfe line search failed, we also need to set the
``PrintLevel`` to at least ``1``.
"""

import time
from contextlib import contextmanager

import numpy as np
import pycgdescent._cg_descent as _cg


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
    param.AWolfe = 0
    param.AWolfeFac = 0.0
    # param.PrintLevel = 1
    # param.PrintParms = 1

    # }}}

    # {{{

    print("==== with tol 1.0e-8 ====")
    with timer():
        _, stats, _ = _cg.cg_descent(x0, 1.0e-8, param, fn, grad, fngrad,
                callback=None, work=None)

    print()
    print("maximum norm for gradient: %+.16e" % stats.gnorm)
    print("function value:            %+.16e" % stats.f)
    print("cg iterations:            ", stats.iter)
    print("function evaluations:     ", stats.nfunc)
    print("gradient evaluations:     ", stats.ngrad)

    # }}}

    # {{{

    print()
    print("==== with tol 1.0e-6 ====")
    with timer():
        _, stats, _ = _cg.cg_descent(x0, 1.0e-6, param, fn, grad, fngrad,
                callback=None, work=None)

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
