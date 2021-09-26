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

from contextlib import contextmanager
from functools import partial
from typing import Iterator

import numpy as np
import pycgdescent as cg
import pycgdescent._cg_descent as _cg

import logging
logger = logging.getLogger()


@contextmanager
def timer() -> Iterator[None]:
    import time
    t_start = time.time()
    yield
    t_end = time.time()
    logger.info("elapsed: %.3fs", t_end - t_start)


def fn(x: cg.ArrayType, t: float = 1.0) -> float:
    f: float = np.sum(np.exp(x) - t * x)
    return f


def grad(g: cg.ArrayType, x: cg.ArrayType, t: float = 1.0) -> None:
    g[...] = np.exp(x) - t


def fngrad(g: cg.ArrayType, x: cg.ArrayType, t: float = 1.0) -> float:
    y: cg.ArrayType = np.exp(x)
    f: float = np.sum(y - t * x)
    g[...] = y - t
    return f


def main(n: int = 100) -> None:
    # {{{ parameters

    x0 = np.ones(n, dtype=np.float64)
    t = np.sqrt(1 + np.arange(n))

    param = _cg.cg_parameter()
    # param.PrintParms = 1

    # }}}

    # {{{

    logger.info("==== with QuadStep OFF ====")
    with timer():
        param.QuadStep = 0
        _, stats, _ = _cg.cg_descent(x0, 1.0e-8, param,
                partial(fn, t=t), partial(grad, t=t), partial(fngrad, t=t),
                callback=None, work=None)

    logger.info("\n")
    logger.info("maximum norm for gradient: %+.16e", stats.gnorm)
    logger.info("function value:            %+.16e", stats.f)
    logger.info("cg iterations:             %d", stats.iter)
    logger.info("function evaluations:      %d", stats.nfunc)
    logger.info("gradient evaluations:      %d", stats.ngrad)

    # }}}

    # {{{

    logger.info("\n")
    logger.info("==== with QuadStep ON ====")
    with timer():
        param.QuadStep = 1
        _, stats, _ = _cg.cg_descent(x0, 1.0e-8, param,
                partial(fn, t=t), partial(grad, t=t), partial(fngrad, t=t),
                callback=None, work=None)

    logger.info("\n")
    logger.info("maximum norm for gradient: %+.16e", stats.gnorm)
    logger.info("function value:            %+.16e", stats.f)
    logger.info("cg iterations:             %d", stats.iter)
    logger.info("function evaluations:      %d", stats.nfunc)
    logger.info("gradient evaluations:      %d", stats.ngrad)

    # }}}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()

# vim: fdm=marker
