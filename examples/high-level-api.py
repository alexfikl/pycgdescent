"""
Example of the high-level API with callbacks and everything!

Uses the classic Rosenbrock function as an example.
"""

# START_ROSENROCK_EXAMPLE
import numpy as np
import numpy.linalg as la
import pycgdescent as cg

from functools import partial
from dataclasses import dataclass, field
from typing import List


@dataclass(frozen=True)
class CallbackCache:
    alpha: List[float] = field(default_factory=list)
    x: List[np.ndarray] = field(default_factory=list)
    f: List[float] = field(default_factory=list)
    g: List[float] = field(default_factory=list)

    def __call__(self, info):
        self.alpha.append(info.alpha)
        self.x.append(info.x.copy())
        self.f.append(info.f)
        self.g.append(la.norm(info.g, np.inf))

        return 1


def fun(x, a=100.0, b=1.0):
    return a * (x[1] - x[0]**2)**2 + b * (x[0] - 1.0)**2


def jac(g, x, a=100.0, b=1.0):
    g[0] = -4.0 * a * x[0] * (x[1] - x[0]**2) + 2.0 * b * (x[0] - 1.0)
    g[1] = 2.0 * a * (x[1] - x[0]**2)


def main(a=100.0, b=1.0, tol=1.0e-8, visualize=False):
    callback = CallbackCache()
    x0 = np.array([-3.5, -4.0])

    options = cg.OptimizeOptions()
    r = cg.minimize(
            fun=partial(fun, a=a, b=b),
            x0=x0,
            jac=partial(jac, a=a, b=b),
            tol=tol,
            options=options,
            callback=callback)

    print(r.pretty())

    if visualize:
        plot_rosenbrock_solution(r, callback, a=a, b=b)
# END_ROSENBROCK_EXAMPLE


def savefig(fig, suffix):
    import os
    filename = os.path.join(os.path.dirname(__file__), f"rosenbrock_{suffix}")

    fig.savefig(filename)
    print("output: ", filename)

    fig.clf()


def plot_rosenbrock_solution(r, cache, a=100.0, b=1.0):
    x = np.array(cache.x).T
    alpha = np.array(cache.alpha)
    f = np.array(cache.f)
    gnorm = np.array(cache.g)

    import matplotlib.pyplot as pt
    fig = pt.figure()

    # {{{ alpha

    ax = fig.gca()
    ax.plot(alpha)
    ax.grid(True)
    savefig(fig, "alpha")

    # }}}

    # {{{ value

    ax = fig.gca()
    ax.semilogy(f)
    ax.grid(True)
    savefig(fig, "value")

    # }}}

    # {{{ grad

    ax = fig.gca()
    ax.semilogy(gnorm)
    ax.grid(True)
    savefig(fig, "gnorm")

    # }}}

    # {{{

    x1d = np.linspace(-4.0, 4.0, 128)
    xy = np.stack(np.meshgrid(x1d, x1d))
    z = fun(xy, a=a, b=b)

    ax = fig.gca()
    c = ax.contourf(xy[0], xy[1], z, levels=48, linestyles="dashed")
    ax.contour(xy[0], xy[1], z, levels=48, colors="k")
    ax.plot(x[0], x[1], "wo-")
    ax.plot(r.x[0], r.x[1], "ro")
    ax.set_aspect("equal")
    fig.colorbar(c, shrink=0.75)

    savefig(fig, "convergence")

    # }}}

    pt.close(fig)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("a", type=float, default=100.0, nargs="?",
            help="Rosenbrock function parameter")
    parser.add_argument("b", type=float, default=1.0, nargs="?",
            help="Rosenbrock function parameter")
    parser.add_argument("--tol", type=float, default=1.0e-8,
            help="stopping condition gradient tolerance")
    parser.add_argument("--visualize", action="store_true")
    args = parser.parse_args()

    main(a=args.a, b=args.b, tol=args.tol, visualize=args.visualize)
