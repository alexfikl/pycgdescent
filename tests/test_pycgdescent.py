# SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.linalg as la
import pytest

import pycgdescent as cg

logger = cg.get_logger(__name__)


# {{{ test_optimize_options


def test_optimize_options() -> None:
    """Test options are immutable."""
    options = cg.OptimizeOptions()

    with pytest.raises(AttributeError):
        options.printLevel = 2

    options = options.replace(printLevel=2)
    assert options.printLevel == 2  # ty: ignore[unresolved-attribute]

    options2 = options.replace(step=1.0)
    logger.info("\n%s", options2.pretty())
    assert (options2.step - 1.0) < 1.0e-15
    assert options2.printLevel == 2  # ty: ignore[unresolved-attribute]

    logger.info("\n%s", options)
    logger.info("\n")
    logger.info("\n%s", options2)
    logger.info("\n")
    logger.info("\n%s", options2.pretty())


# }}}


# {{{ test_quadratic


@dataclass(init=False)
class Quadratic:
    A: cg.Matrix
    b: cg.Array
    x: cg.Array

    def __init__(self) -> None:
        self.A = np.array([[4.0, 1.0], [1.0, 3.0]])
        self.b = np.array([1.0, 2.0])
        self.x = np.array([1.0 / 11.0, 7.0 / 11.0])
        self.x0 = np.array([2.0, 1.0])

    def fun(self, x: cg.Array) -> float:
        return (x @ (self.A @ x) - x @ self.b).item()

    def jac(self, g: cg.Array, x: cg.Array) -> None:
        g[...] = self.A @ x - self.b

    def funjac(self, g: cg.Array, x: cg.Array) -> float:
        g[...] = self.A @ x - self.b
        return (x @ g).item()


@pytest.mark.parametrize("tol", [1.0e-8])
def test_quadratic(tol: float) -> None:
    """Test optimization of a quadratic function with default options."""

    # {{{ optimize

    def callback(info: cg.CallbackInfo) -> int:
        logger.info(
            "[%4d] x %.5e %.5e f %.5e g %.5e %.5e", info.it, *info.x, info.f, *info.g
        )

        return 1

    model = Quadratic()
    options = cg.OptimizeOptions(PrintLevel=3)

    r = cg.minimize(
        fun=model.fun,
        x0=model.x0,
        jac=model.jac,
        funjac=model.funjac,
        tol=tol,
        callback=callback,
        options=options,
    )
    logger.info("\n%s", r.pretty())

    # }}}

    # {{{ check

    error = la.norm(r.x - model.x) / la.norm(model.x)

    logger.info("\n%s", r.pretty())
    logger.info("\n")
    logger.info("Solution:  %s", model.x)
    logger.info("Error:     %.16e", error)

    assert r.jac < tol
    assert error < tol

    # }}}


# }}}


# {{{ test_rosenbrock


@pytest.mark.parametrize(("a", "b", "tol"), [(100.0, 1.0, 1.0e-8)])
def test_rosenbrock(a: float, b: float, tol: float) -> None:
    """Test optimization of the Rosenbrock function with default options."""

    if a < 0.0 or b < 0.0:
        raise ValueError("'a' and 'b' must be positive")

    # {{{ setup

    # https://en.wikipedia.org/wiki/Rosenbrock_function
    x0: cg.Array = np.array([-2.0, 1.0])
    x_exact: cg.Array = np.array([1.0, 1.0])

    def fun(x: cg.Array) -> float:
        return a * (x[1] - x[0] ** 2) ** 2 + b * (x[0] - 1.0) ** 2

    def jac(g: cg.Array, x: cg.Array) -> None:
        g[0] = -4.0 * a * x[0] * (x[1] - x[0] ** 2) + 2.0 * b * (x[0] - 1.0)
        g[1] = 2.0 * a * (x[1] - x[0] ** 2)

    # }}}

    # {{{ optimize

    def callback(info: cg.CallbackInfo) -> int:
        logger.info(
            "[%4d] x %.5e %.5e f %.5e g %.5e %.5e", info.it, *info.x, info.f, *info.g
        )

        return 1

    options = cg.OptimizeOptions()
    r = cg.minimize(
        fun=fun,
        x0=x0,
        jac=jac,
        tol=tol,
        callback=callback,
        options=options,
    )
    logger.info("\n%s", r.pretty())

    # }}}

    # {{{ check

    error = la.norm(r.x - x_exact) / la.norm(x_exact)

    logger.info("\n%s", r.pretty())
    logger.info("\n")
    logger.info("Solution:  %s", x_exact)
    logger.info("Error:     %.16e", error)

    assert r.jac < tol
    assert error < tol

    # }}}


# }}}


# {{{ test_exceptions


def test_exceptions() -> None:
    """Test that exceptions raised by the callbacks propagate cleanly."""

    model = Quadratic()

    # }}}

    # {{{ value

    def fun_raise(x: cg.Array) -> float:
        raise RuntimeError("value boom")

    with pytest.raises(RuntimeError, match="value boom"):
        cg.minimize(fun=fun_raise, x0=model.x0, jac=model.jac, tol=1.0e-8)

    # }}}

    # {{{ grad

    def jac_raise(g: cg.Array, x: cg.Array) -> None:
        raise RuntimeError("grad boom")

    with pytest.raises(RuntimeError, match="grad boom"):
        cg.minimize(fun=model.fun, x0=model.x0, jac=jac_raise, tol=1.0e-8)

    # }}}

    # {{{ funjac

    def funjac_raise(g: cg.Array, x: cg.Array) -> float:
        raise RuntimeError("funjac boom")

    with pytest.raises(RuntimeError, match="funjac boom"):
        cg.minimize(
            fun=model.fun, x0=model.x0, jac=model.jac, funjac=funjac_raise, tol=1.0e-8
        )

    # }}}

    # {{{ callback

    def callback_raise(info: cg.CallbackInfo) -> int:
        raise RuntimeError("callback boom")

    with pytest.raises(RuntimeError, match="callback boom"):
        cg.minimize(
            fun=model.fun,
            x0=model.x0,
            jac=model.jac,
            tol=1.0e-8,
            callback=callback_raise,
        )

    # }}}

    # {{{ mid-run

    calls: dict[str, int] = {"n": 0}

    def fun_midraise(x: cg.Array) -> float:
        calls["n"] += 1
        if calls["n"] > 2:
            raise RuntimeError("mid-run boom")

        return model.fun(x)

    with pytest.raises(RuntimeError, match="mid-run boom"):
        cg.minimize(fun=fun_midraise, x0=model.x0, jac=model.jac, tol=1.0e-8)

    # }}}


# }}}


# {{{ test_status


def test_status() -> None:
    """Test that the termination status is returned as an exact int code."""

    # {{{ setup

    A: cg.Matrix = np.array([[4.0, 1.0], [1.0, 3.0]])  # ruff:ignore[non-lowercase-variable-in-function]
    b: cg.Array = np.array([1.0, 2.0])
    x0: cg.Array = np.array([2.0, 1.0])

    def fun(x: cg.Array) -> float:
        return (x @ (A @ x) - x @ b).item()

    def jac(g: cg.Array, x: cg.Array) -> None:
        g[...] = A @ x - b

    # }}}

    # {{{ success

    r = cg.minimize(fun=fun, x0=x0, jac=jac, tol=1.0e-8)
    assert r.success
    assert r.status == 0
    assert "Convergence" in r.message

    # }}}

    # {{{ maxit

    options = cg.OptimizeOptions(maxit=1)
    r = cg.minimize(fun=fun, x0=x0, jac=jac, tol=1.0e-12, options=options)
    assert not r.success
    assert r.status == 2
    assert "Maximum number of iterations" in r.message

    # }}}

    # {{{ nan

    def fun_nan(x: cg.Array) -> float:
        return float("nan")

    r = cg.minimize(fun=fun_nan, x0=x0, jac=jac, tol=1.0e-8)
    assert not r.success
    assert r.status == 11
    assert "NaN or Inf" in r.message

    # }}}

    # {{{ callback stop

    def callback_stop(info: cg.CallbackInfo) -> int:
        return 0

    r = cg.minimize(fun=fun, x0=x0, jac=jac, tol=1.0e-8, callback=callback_stop)
    assert not r.success
    assert r.status == 13
    assert "Stopped by user callback" in r.message

    # }}}


# }}}


# {{{ test_ndarray_strided


def test_ndarray_strided() -> None:
    """Test that non-contiguous (strided) input arrays are handled correctly."""

    # {{{ setup

    model = Quadratic()
    x0: cg.Array = np.array([1.0, 42.0, 2.0, 43.0])[::2]
    assert not x0.flags.c_contiguous
    assert np.array_equal(x0, [1.0, 2.0])

    first: dict[str, cg.Array] = {}

    def fun_capture(x: cg.Array) -> float:
        if "x" not in first:
            first["x"] = x.copy()

        return model.fun(x)

    # }}}

    # {{{ optimize

    r = cg.minimize(fun=fun_capture, x0=x0, jac=model.jac, tol=1.0e-8)
    logger.info("\n%s", r.pretty())

    # }}}

    # {{{ check

    # the value function must see the logical (converted) input, not raw memory
    assert np.array_equal(first["x"], [1.0, 2.0])
    assert np.allclose(r.x, model.x)
    assert r.x.dtype == np.float64
    assert r.x.flags.c_contiguous

    # }}}


# }}}


# {{{ test_ndarray_dtype_conversion


@pytest.mark.parametrize("dtype", [np.float32, np.int64])
def test_ndarray_dtype_conversion(dtype: np.dtype) -> None:
    """Test that non-float64 input arrays are converted to float64."""

    # {{{ setup

    model = Quadratic()
    x0: cg.Array = np.array([2.0, 1.0], dtype=dtype)

    # }}}

    # {{{ optimize

    r = cg.minimize(fun=model.fun, x0=x0, jac=model.jac, tol=1.0e-8)
    logger.info("\n%s", r.pretty())

    # }}}

    # {{{ check

    assert r.x.dtype == np.float64
    assert np.allclose(r.x, model.x)

    # }}}


# }}}


# {{{ test_ndarray_list_rejected


def test_ndarray_list_rejected() -> None:
    """Test that plain Python sequences are rejected instead of converted."""

    model = Quadratic()

    with pytest.raises(TypeError):
        cg.minimize(fun=model.fun, x0=[2.0, 1.0], jac=model.jac, tol=1.0e-8)  # ty: ignore[invalid-argument-type]


# }}}


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        exec(sys.argv[1])
    else:
        raise SystemExit(pytest.main([__file__]))

# vim: fdm=marker
