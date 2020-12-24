import numpy as np
import numpy.linalg as la

import pytest
import pycgdescent as cg


def test_optimize_options():
    options = cg.OptimizeOptions()

    with pytest.raises(AttributeError):
        options.PrintLevel = 2

    options2 = options.replace(PrintLevel=2)
    assert options2.PrintLevel == 2

    print(options)
    print()
    print(options2)
    print()
    print(options2.pretty())


def test_quadratic(tol=1.0e-8):
    # {{{ setup

    # https://en.wikipedia.org/wiki/Conjugate_gradient_method#Numerical_example
    A = np.array([[4.0, 1.0], [1.0, 3.0]])      # noqa: N806
    b = np.array([1.0, 2.0])

    x0 = np.array([2.0, 1.0])
    x_exact = np.array([1.0 / 11.0, 7.0 / 11.0])

    def fun(x):
        return x.dot(A @ x) - x.dot(b)

    def jac(g, x):
        g[...] = A @ x - b

    def funjac(g, x):
        Ax = A @ x -b                           # noqa: N806
        g[...] = Ax
        return x.dot(g)

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

    error = la.norm(r.x - x_exact) / la.norm(x_exact)

    print(r.pretty())
    print()
    print("Solution: ", x_exact)
    print("Error:     %.16e" % error)

    assert r.jac < tol
    assert error < tol

    # }}}


def test_rosenbrock(tol=1.0e-8):
    # {{{ setup

    # https://en.wikipedia.org/wiki/Rosenbrock_function
    a = 100.0
    b = 1.0
    x0 = np.array([-2.0, 1.0])
    x_exact = np.array([1.0, 1.0])

    def fun(x):
        return a * (x[1] - x[0]**2)**2 + b * (x[0] - 1.0)**2;

    def jac(g, x):
        g[0] = -4.0 * a * x[0] * (x[1] - x[0]**2) + 2.0 * b * (x[0] - 1.0)
        g[1] = 2.0 * a * (x[1] - x[0]**2)

    def funjac(g, x):
        jac(g, x)
        return fun(x)

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

    error = la.norm(r.x - x_exact) / la.norm(x_exact)

    print(r.pretty())
    print()
    print("Solution: ", x_exact)
    print("Error:     %.16e" % error)

    assert r.jac < tol
    assert error < tol

    # }}}


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        exec(sys.argv[1])
    else:
        from pytest import main
        main([__file__])

# vim: fdm=marker
