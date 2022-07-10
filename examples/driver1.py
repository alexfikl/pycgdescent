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

from functools import partial

import numpy as np
import pycgdescent as cg
import pycgdescent._cg_descent as _cg

import logging

logger = logging.getLogger()


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

    # param = _cg.cg_parameter()
    # param.logger.infoParms = 1

    # }}}

    # {{{ without fngrad

    logger.info("==== without fngrad ====")
    with cg.timer():
        x1, stats, status = _cg.cg_descent(
            x0,
            1.0e-8,
            None,
            partial(fn, t=t),
            partial(grad, t=t),
            None,
            callback=None,
            work=None,
        )

    # }}}

    # {{{ with fngrad

    x0 = np.ones(n, dtype=np.float64)

    logger.info("==== with fngrad ====")
    with cg.timer():
        x2, stats, status = _cg.cg_descent(
            x0,
            1.0e-8,
            None,
            partial(fn, t=t),
            partial(grad, t=t),
            partial(fngrad, t=t),
            callback=None,
            work=None,
        )

    # }}}

    from pycgdescent import STATUS_TO_MESSAGE

    logger.info("\n")
    logger.info("status:  %d", status)
    logger.info("message: %s", STATUS_TO_MESSAGE[status])

    logger.info("\n")
    logger.info("maximum norm for gradient: %+.16e", stats.gnorm)
    logger.info("function value:            %+.16e", stats.f)
    logger.info("cg iterations:             %d", stats.iter)
    logger.info("function evaluations:      %d", stats.nfunc)
    logger.info("gradient evaluations:      %d", stats.ngrad)

    assert np.linalg.norm(x1 - x2) / np.linalg.norm(x2) < 1.0e-15


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()

# vim: fdm=marker
