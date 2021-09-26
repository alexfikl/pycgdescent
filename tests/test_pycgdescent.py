import numpy as np
import numpy.linalg as la
import pycgdescent as cg

import pytest

import logging
logger = logging.getLogger()


def test_optimize_options() -> None:
    """Test options are immutable."""
    options = cg.OptimizeOptions()

    with pytest.raises(AttributeError):
        options.printLevel = 2

    options = options.replace(printLevel=2)
    assert options.printLevel == 2      # type: ignore[attr-defined]

    options2 = options.replace(step=1.0)
    logger.info("\n%s", options2.pretty())
    assert (options2.step - 1.0) < 1.0e-15
    assert options2.printLevel == 2     # type: ignore[attr-defined]

    logger.info("\n%s", options)
    logger.info("\n")
    logger.info("\n%s", options2)
    logger.info("\n")
    logger.info("\n%s", options2.pretty())


def test_quadratic(tol: float = 1.0e-8) -> None:
    """Test optimization of a quadratic function with default options."""

    # {{{ setup

    # https://en.wikipedia.org/wiki/Conjugate_gradient_method#Numerical_example
    A = np.array([[4.0, 1.0], [1.0, 3.0]])      # noqa: N806
    b = np.array([1.0, 2.0])

    x0 = np.array([2.0, 1.0])
    x_exact = np.array([1.0 / 11.0, 7.0 / 11.0])

    def fun(x: cg.ArrayType) -> float:
        f: float = x.dot(A @ x) - x.dot(b)
        return f

    def jac(g: cg.ArrayType, x: cg.ArrayType) -> None:
        g[...] = A @ x - b

    def funjac(g: cg.ArrayType, x: cg.ArrayType) -> float:
        Ax: cg.ArrayType = A @ x - b                            # noqa: N806
        f: float = x @ g
        g[...] = Ax
        return f

    # }}}

    # {{{ optimize

    options = cg.OptimizeOptions()
    r = cg.minimize(
            fun=fun,
            x0=x0,
            jac=jac,
            funjac=funjac,
            tol=tol,
            options=options,
            )

    # }}}

    # {{{ check

    error = la.norm(r.x - x_exact) / la.norm(x_exact)           # type: ignore

    logger.info("\n%s", r.pretty())
    logger.info("\n")
    logger.info("Solution:  %s", x_exact)
    logger.info("Error:     %.16e", error)

    assert r.jac < tol
    assert error < tol

    # }}}


def test_rosenbrock(a: float = 100.0, b: float = 1.0, tol: float = 1.0e-8) -> None:
    """Test optimization of the Rosenbrock function with default options."""

    if a < 0.0 or b < 0.0:
        raise ValueError("'a' and 'b' must be positive")

    # {{{ setup

    # https://en.wikipedia.org/wiki/Rosenbrock_function
    x0 = np.array([-2.0, 1.0])
    x_exact = np.array([1.0, 1.0])

    def fun(x: cg.ArrayType) -> float:
        f: float = a * (x[1] - x[0]**2)**2 + b * (x[0] - 1.0)**2
        return f

    def jac(g: cg.ArrayType, x: cg.ArrayType) -> None:
        g[0] = -4.0 * a * x[0] * (x[1] - x[0]**2) + 2.0 * b * (x[0] - 1.0)
        g[1] = 2.0 * a * (x[1] - x[0]**2)

    # }}}

    # {{{ optimize

    options = cg.OptimizeOptions()
    r = cg.minimize(
            fun=fun,
            x0=x0,
            jac=jac,
            tol=tol,
            options=options,
            )

    # }}}

    # {{{ check

    error = la.norm(r.x - x_exact) / la.norm(x_exact)           # type: ignore

    logger.info("\n%s", r.pretty())
    logger.info("\n")
    logger.info("Solution:  %s", x_exact)
    logger.info("Error:     %.16e", error)

    assert r.jac < tol
    assert error < tol

    # }}}


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) > 1:
        exec(sys.argv[1])
    else:
        from pytest import main
        main([__file__])

# vim: fdm=marker
