# SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

"""
Example of the high-level API with callbacks and everything!

Uses the classic Rosenbrock function as an example.
"""

# START_ROSENROCK_EXAMPLE
from dataclasses import dataclass, field
from functools import partial
from typing import Any, List, Optional

import numpy as np
import numpy.linalg as la
import pycgdescent as cg


@dataclass(frozen=True)
class CallbackCache:
    alpha: List[cg.ScalarType] = field(default_factory=list)
    x: List[cg.ArrayType] = field(default_factory=list)
    f: List[cg.ScalarType] = field(default_factory=list)
    g: List[cg.ScalarType] = field(default_factory=list)

    def __call__(self, info: cg.CallbackInfo) -> int:
        self.alpha.append(info.alpha)
        self.x.append(info.x.copy())
        self.f.append(info.f)
        self.g.append(la.norm(info.g, np.inf))

        return 1


def fun(x: cg.ArrayType, a: float = 100.0, b: float = 1.0) -> cg.ScalarType:
    return a * (x[1] - x[0] ** 2) ** 2 + b * (x[0] - 1.0) ** 2


def jac(g: cg.ArrayType, x: cg.ArrayType, a: float = 100.0, b: float = 1.0) -> None:
    g[0] = -4.0 * a * x[0] * (x[1] - x[0] ** 2) + 2.0 * b * (x[0] - 1.0)
    g[1] = 2.0 * a * (x[1] - x[0] ** 2)


def main(
    a: float = 100.0, b: float = 1.0, tol: float = 1.0e-8, visualize: bool = False
) -> None:
    callback = CallbackCache()
    x0: cg.ArrayType = np.array([-3.5, -4.0])

    options = cg.OptimizeOptions()
    r = cg.minimize(
        fun=partial(fun, a=a, b=b),
        x0=x0,
        jac=partial(jac, a=a, b=b),
        tol=tol,
        options=options,
        callback=callback,
    )

    print(r.pretty())

    if visualize:
        plot_rosenbrock_solution(r, callback, a=a, b=b)
    # END_ROSENBROCK_EXAMPLE


def savefig(fig: Any, suffix: str, ext: Optional[str] = None) -> None:
    import pathlib

    if ext is None:
        import matplotlib.pyplot as pt

        ext = pt.rcParams["savefig.format"]

    filename = pathlib.Path(__file__).parent / f"rosenbrock_{suffix}"
    filename = filename.with_suffix(f".{ext}")

    fig.savefig(filename)
    print("output: ", filename)

    fig.clf()


def plot_rosenbrock_solution(
    r: cg.OptimizeResult, cache: CallbackCache, a: float = 100.0, b: float = 1.0
) -> None:
    x: cg.ArrayType = np.array(cache.x).T
    alpha: cg.ArrayType = np.array(cache.alpha)
    f: cg.ArrayType = np.array(cache.f)
    gnorm: cg.ArrayType = np.array(cache.g)

    import matplotlib.pyplot as pt

    pt.style.use("seaborn")
    fig = pt.figure(figsize=(10, 10), dpi=300)

    # {{{ alpha

    ax = fig.gca()
    ax.plot(alpha)
    ax.set_xlabel("$Iteration$")
    ax.set_ylabel("$Step~ Size$", fontsize="large")
    savefig(fig, "alpha")

    # }}}

    # {{{ value

    ax = fig.gca()
    ax.semilogy(f)
    ax.set_xlabel("$Iteration$")
    ax.set_ylabel("$f$")
    savefig(fig, "value")

    # }}}

    # {{{ grad

    ax = fig.gca()
    ax.semilogy(gnorm)
    ax.set_xlabel("$Iteration$")
    ax.set_ylabel("$Gradient~ Magnitude$", fontsize="large")
    savefig(fig, "gnorm")

    # }}}

    # {{{

    x1d = np.linspace(-4.0, 4.0, 128)
    xy: cg.ArrayType = np.stack(np.meshgrid(x1d, x1d))
    z = fun(xy, a=a, b=b)

    ax = fig.gca()
    c = ax.contourf(xy[0], xy[1], z, levels=48, linestyles="dashed", cmap="viridis")
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
    parser.add_argument(
        "a",
        type=float,
        default=100.0,
        nargs="?",
        help="Rosenbrock function parameter",
    )
    parser.add_argument(
        "b", type=float, default=1.0, nargs="?", help="Rosenbrock function parameter"
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1.0e-8,
        help="stopping condition gradient tolerance",
    )
    parser.add_argument("--visualize", action="store_true")
    args = parser.parse_args()

    main(a=args.a, b=args.b, tol=args.tol, visualize=args.visualize)
