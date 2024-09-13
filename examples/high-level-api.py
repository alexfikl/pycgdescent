# SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

"""
Example of the high-level API with callbacks and everything!

Uses the classic Rosenbrock function as an example.
"""

from __future__ import annotations

from functools import partial

import matplotlib.figure
import matplotlib.pyplot as mp

# START_ROSENROCK_EXAMPLE
import numpy as np
import numpy.linalg as la

import pycgdescent as cg

Array = np.ndarray[tuple[int], np.dtype[np.float64]]
Float = np.float64


class CallbackCache:
    def __init__(self) -> None:
        self.alpha: list[float] = []
        self.x: list[Array] = []
        self.f: list[float] = []
        self.g: list[Float] = []

    def __call__(self, info: cg.CallbackInfo) -> int:
        self.alpha.append(info.alpha)
        self.x.append(info.x.copy())
        self.f.append(info.f)
        self.g.append(la.norm(info.g, np.inf))

        return 1


def fun(x: Array, *, a: float, b: float) -> float:
    return a * (x[1] - x[0] ** 2) ** 2 + b * (x[0] - 1.0) ** 2  # type: ignore


def jac(g: Array, x: Array, *, a: float, b: float) -> None:
    g[0] = -4.0 * a * x[0] * (x[1] - x[0] ** 2) + 2.0 * b * (x[0] - 1.0)
    g[1] = 2.0 * a * (x[1] - x[0] ** 2)


def main(*, a: float = 100.0, b: float = 1.0, tol: float = 1.0e-8) -> None:
    callback = CallbackCache()
    x0 = np.array([-3.5, -4.0])

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
    # END_ROSENBROCK_EXAMPLE

    return r, callback  # type: ignore[return-value]


def savefig(fig: matplotlib.figure.Figure, suffix: str, ext: str | None = None) -> None:
    import pathlib

    if ext is None:
        ext = mp.rcParams["savefig.format"]

    filename = pathlib.Path(__file__).parent / f"rosenbrock_{suffix}.{ext}"

    fig.tight_layout()
    fig.savefig(filename, bbox_inches="tight")
    print("output: ", filename)

    fig.clf()


def plot_rosenbrock_solution(
    r: cg.OptimizeResult,
    cache: CallbackCache,
    *,
    a: float = 100.0,
    b: float = 1.0,
    ext: str = "png",
    dark: bool = False,
) -> None:
    x = np.array(cache.x).T
    alpha = np.array(cache.alpha)
    f = np.array(cache.f)
    gnorm = np.array(cache.g)

    facecolor = "#121212" if dark else "#FFFFFF"
    fontcolor = "#FFFFFF" if dark else "#000000"

    mp.style.use("seaborn")
    mp.rc("text", usetex=True)
    mp.rc("figure", facecolor=facecolor)
    mp.rc("axes", labelsize=32, titlesize=32, labelcolor=fontcolor)
    mp.rc("xtick", labelsize=18, color=fontcolor)
    mp.rc("ytick", labelsize=18, color=fontcolor)

    fig = mp.figure(figsize=(10, 10), dpi=300, constrained_layout=True)

    # NOTE: these are the background colors for sphinx-book-theme

    # {{{ alpha

    ax = fig.gca()
    ax.plot(alpha)
    ax.set_xlabel("$Iteration$")
    ax.set_ylabel("$Step~ Size$")
    savefig(fig, "alpha", ext=ext)

    # }}}

    # {{{ value

    ax = fig.gca()
    ax.semilogy(f)
    ax.set_xlabel("$Iteration$")
    ax.set_ylabel("$f$")
    savefig(fig, "value", ext=ext)

    # }}}

    # {{{ grad

    ax = fig.gca()
    ax.semilogy(gnorm)
    ax.set_xlabel("$Iteration$")
    ax.set_ylabel("$Gradient~ Magnitude$")
    savefig(fig, "gnorm", ext=ext)

    # }}}

    # {{{

    x1d = np.linspace(-4.0, 4.0, 128)
    xy = np.stack(np.meshgrid(x1d, x1d))
    z = fun(xy, a=a, b=b)

    ax = fig.gca()
    c = ax.contourf(xy[0], xy[1], z, levels=48, linestyles="dashed", cmap="viridis")
    ax.contour(xy[0], xy[1], z, levels=48, colors="k")
    ax.plot(x[0], x[1], "wo-")
    ax.plot(r.x[0], r.x[1], "ro")

    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    fig.colorbar(c, shrink=0.73)

    savefig(fig, "convergence", ext=ext)

    # }}}

    mp.close(fig)


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
    parser.add_argument("--dark", action="store_true")
    parser.add_argument("--visualize", action="store_true")
    args = parser.parse_args()

    result, callback = main(a=args.a, b=args.b, tol=args.tol)  # type: ignore
    if args.visualize:
        plot_rosenbrock_solution(result, callback, a=args.a, b=args.b, dark=args.dark)
