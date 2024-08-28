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
    # param.PrintParms = 1

    # }}}

    # {{{ without fngrad

    logger.info("==== without fngrad ====")
    with cg.Timer() as time:
        x1, stats, status = cg.cg_descent(
            x0,
            1.0e-8,
            partial(fn, t=t),
            partial(grad, t=t),
            param=param,
        )

    logger.info("timing: %s\n", time)

    # }}}

    # {{{ with fngrad

    x0 = np.ones(n, dtype=np.float64)

    logger.info("==== with fngrad ====")
    with cg.Timer() as time:
        x2, stats, status = cg.cg_descent(
            x0,
            1.0e-8,
            partial(fn, t=t),
            partial(grad, t=t),
            valgrad=partial(fngrad, t=t),
            param=param,
        )

    logger.info("timing: %s\n", time)

    # }}}

    from pycgdescent import STATUS_TO_MESSAGE

    logger.info("status:  %d", status)
    logger.info("message: %s\n", STATUS_TO_MESSAGE[status])

    logger.info("maximum norm for gradient: %+.16e", stats.gnorm)
    logger.info("function value:            %+.16e", stats.f)
    logger.info("cg iterations:             %d", stats.iter)
    logger.info("function evaluations:      %d", stats.nfunc)
    logger.info("gradient evaluations:      %d", stats.ngrad)

    assert np.linalg.norm(x1 - x2) / np.linalg.norm(x2) < 1.0e-15


if __name__ == "__main__":
    main()

# vim: fdm=marker
