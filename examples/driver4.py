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

When the line search routine tries to bracket a local minimizer in
the search direction, it may expand the initial line search interval.
The default expansion factor is :math:`5`. You can modify this factor using
the parameter ``rho``.  In the following example, we choose a small initial
step size (initial step is ``1.0e-5``), ``QuadStep`` is *False*, and ``rho`` is
``1.5``.

The code has to do a number of expansions to reach a suitable
interval bracketing the minimizer in the initial search direction.
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
    param.QuadStep = 0
    param.step = 1.0e-5
    # param.logger.infoParms = 1

    # }}}

    # {{{

    logger.info("==== with rho 1.5 ====")
    with cg.Timer() as time:
        param.rho = 1.5
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

    logger.info("==== with rho 5.0 ====")
    with cg.Timer() as time:
        param.rho = 5.0
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
    logger.info("gradient evaluations:      %d", stats.ngrad)

    # }}}


if __name__ == "__main__":
    main()

# vim: fdm=marker
