# SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

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
``logger.infoLevel`` to at least ``1``.
"""

from __future__ import annotations

from functools import partial

import numpy as np

import pycgdescent as cg

logger = cg.get_logger()


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

    x0: cg.ArrayType = np.ones(n, dtype=np.float64)
    t = np.sqrt(1 + np.arange(n))

    param = cg.cg_parameter()
    param.AWolfe = 0
    param.AWolfeFac = 0.0
    # param.logger.infoLevel = 1
    # param.logger.infoParms = 1

    # }}}

    # {{{

    logger.info("==== with tol 1.0e-8 ====")
    with cg.Timer() as time:
        _, stats, _ = cg.cg_descent(
            x0,
            1.0e-8,
            partial(fn, t=t),
            partial(grad, t=t),
            valgrad=partial(fngrad, t=t),
            param=param,
        )

    logger.info("timing: %s\n", time)

    logger.info("maximum norm for gradient: %+.16e", stats.gnorm)
    logger.info("function value:            %+.16e", stats.f)
    logger.info("cg iterations:             %d", stats.iter)
    logger.info("function evaluations:      %d", stats.nfunc)
    logger.info("gradient evaluations:      %d\n", stats.ngrad)

    # }}}

    # {{{

    x0 = np.ones(n, dtype=np.float64)

    logger.info("==== with tol 1.0e-6 ====")
    with cg.Timer() as time:
        _, stats, _ = cg.cg_descent(
            x0,
            1.0e-6,
            partial(fn, t=t),
            partial(grad, t=t),
            valgrad=partial(fngrad, t=t),
            param=param,
        )

    logger.info("timing: %s\n", time)

    logger.info("maximum norm for gradient: %+.16e", stats.gnorm)
    logger.info("function value:            %+.16e", stats.f)
    logger.info("cg iterations:             %d", stats.iter)
    logger.info("function evaluations:      %d", stats.nfunc)
    logger.info("gradient evaluations:      %d", stats.ngrad)

    # }}}


if __name__ == "__main__":
    main()

# vim: fdm=marker
