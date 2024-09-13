# SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

"""
Functions
^^^^^^^^^

.. autofunction:: minimize

.. autofunction:: min_work_size
.. autofunction:: allocate_work_for

.. autoclass:: OptimizeOptions

.. autoclass:: OptimizeResult
    :members:

.. autoclass:: CallbackInfo
    :members:

Type Aliases
^^^^^^^^^^^^

.. autodata:: FunType

.. autodata:: GradType

.. autodata:: FunGradType

.. autodata:: CallbackType

"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from importlib import metadata
from typing import TYPE_CHECKING, Any, Callable, ClassVar, TypeAlias, cast

import numpy as np

import pycgdescent._cg_descent as _cg

__version__ = metadata.version("pycgdescent")

if TYPE_CHECKING:
    ArrayType: TypeAlias = np.ndarray[Any, np.dtype[Any]]
else:
    ArrayType: TypeAlias = np.ndarray


# {{{ wrap

# NOTE: These are mostly here for the examples so they can use the low-level
# wrappers directly and not the 'minimize' fancy wrapper.

cg_parameter = _cg.cg_parameter


def cg_descent(
    x: ArrayType,
    grad_tol: float,
    value: Callable[[ArrayType], float],
    grad: Callable[[ArrayType, ArrayType], None],
    *,
    valgrad: Callable[[ArrayType, ArrayType], float] | None = None,
    callback: Callable[[_cg.cg_iter_stats], int] | None = None,
    work: ArrayType | None = None,
    param: _cg.cg_parameter | None = None,
) -> tuple[ArrayType, _cg.cg_stats, bool]:
    """A thin wrapper around the original ``cg_descent`` implementation."""
    if param is None:
        param = _cg.cg_parameter()

    return _cg.cg_descent(x, grad_tol, param, value, grad, valgrad, callback, work)


# }}}


# {{{ options


def _getmembers(obj: Any) -> list[str]:
    import inspect

    return [
        m
        for m in dir(obj)
        if not m.startswith("__") and not inspect.ismethod(getattr(obj, m))
    ]


def _stringify_dict(d: dict[str, Any]) -> str:
    width = len(max(d, key=len))
    fmt = "{:" + str(width) + "} : {}"

    items = sorted({k: repr(v) for k, v in d.items()}.items())

    from itertools import starmap

    return "\n".join(["\t" + "\n\t".join(starmap(fmt.format, items))])


class OptimizeOptions(_cg.cg_parameter):
    r"""Optimization options for the ``CG_DESCENT`` algorithm. A description
    of some of the more technical parameters can be found in the paper
    [HagerZhang2006]_ and [HagerZhang2013]_.

    .. [HagerZhang2006] W. W. Hager, H. Zhang,
        *Algorithm 851: CG_DESCENT, a Conjugate Gradient Method With Guaranteed
        Descent*,
        ACM Transactions on Mathematical Software, Vol. 32, pp. 113--137, 2006,
        `DOI <https://dx.doi.org/10.1145/1132973.1132979>`__.

    .. [HagerZhang2013] W. W. Hager, H. Zhang,
        *The Limited Memory Conjugate Gradient Method*,
        SIAM Journal on Optimization, Vol. 23, pp. 2150--2168, 2013,
        `DOI <https://dx.doi.org/10.1137/120898097>`__.

    .. automethod:: replace
    .. automethod:: pretty

    The attribute names here follow those of the original ``CG_DESCENT`` code.

    .. attribute:: PrintLevel
        :type: int

        Display level as an integer in `[0, 1, 2, 3]`.

    .. attribute:: LBFGS
        :type: bool

        Boolean flag to handle use of LBFGS. If *False*, LBFGS is only used
        when the :attr:`memory` is larger than the input size :math:`n`.

    .. attribute:: memory
        :type: int

        Number of vectors stored in memory.

    .. attribute:: SubCheck
        :type: int

    .. attribute:: SubSkip
        :type: int

        Together with :attr:`SubCheck`, it controls the frequency with which
        the subspace condition is checked. It is checked at every
        ``SubCheck * memory`` iterations and, if not satisfied, then
        it is skipped for ``SubSkip * memory`` iterations and :attr:`SubSkip`
        is doubled. Whenever the subspace condition is satisfied, :attr:`SubSkip`
        is returned to the original value.

    .. attribute:: eta0
        :type: float

        Controls relative distance from current gradient to the subspace.
        If the distance is ``<= eta0`` and the subspace dimension is
        :attr:`memory`, then the subspace is entered. This is used as
        :math:`\eta_0^2` in Equation 3.4 from [HagerZhang2013]_.

    .. attribute:: eta1
        :type: float

        Controls relative distance from current gradient to the subspace.
        If the distance is ``>= eta1``, the subspace is exited. This is used
        as :math:`\eta_1^2` in Equation 3.4 in [HagerZhang2013]_.

    .. attribute:: eta2
        :type: float

        Controls relative distance from current descent direction to the subspace.
        If the distance is ``<= eta2``, the subspace is always entered.

    .. attribute:: AWolfe
        :type: bool

    .. attribute:: AWolfeFac
        :type: float

        If :attr:`AWolfe` is *True*, then the approximate Wolfe condition is
        used when :math:`|f_{k + 1} - f_k| < \omega C_k`, for
        :math:`\omega` = :attr:`AWolfeFac`. See discussion surrounding
        Equation 25 in [HagerZhang2006]_.

    .. attribute:: Qdecay
        :type: float

        Factor in :math:`[0, 1]` used to compute the average cost magnitude
        in the Wolfe condition.

    .. attribute:: nslow
        :type: int

        Maximum number of "slow" iterations without strict improvement in
        either function values or gradient.

    .. attribute:: StopRule
        :type: bool

        If *True*, a gradient-based stopping condition is used. Otherwise,
        a function value-based stopping condition is used. They are

        .. math::

            \begin{aligned}
            \|\mathbf{g}_k\|_\infty \le &
                \max (\epsilon_g, \delta \|\mathbf{g}_0\|_\infty), \\
            \|\mathbf{g}_k\|_\infty \le & \epsilon_g (1 + |f_k|),
            \end{aligned}

        where :math:`\epsilon_g` is the tolerance in :func:`minimize` and
        :math:`\delta` = :attr:`StopFac`.

    .. attribute:: StopFac
        :type: float

        Factor used in stopping condition.

    .. attribute:: PertRule
        :type: bool

        Estimate error in function values. If *False*, use just :attr:`eps`,
        otherwise use :math:`\epsilon C_k`, where :math:`C_k` is defined in
        Equation 26 in [HagerZhang2006]_.

    .. attribute:: eps
        :type: float

        Factor used in :attr:`PertRule`.

    .. attribute:: egrow
        :type: float

        Factor by which :attr:`eps` grows when the line search fails during
        contraction.

    .. attribute:: QuadStep
        :type: bool

        If *False*, do not use a quadratic interpolation step in the line
        search. If *True*, attempt a step based on
        :math:`\epsilon` equal to :attr:`QuadCutOff` when

        .. math::

            \frac{|f_{k + 1} - f_k|}{|f_k|} > \epsilon.

    .. attribute:: QuadCutOff
        :type: float

        Factor used when :attr:`QuadStep` is *True*.

    .. attribute:: QuadSafe
        :type: float

        Maximum factor by which a quadratic step can reduce the step size.

    .. attribute:: UseCubic
        :type: bool

        Boolean flag that enables a cubic step in the line search.

    .. attribute:: CubicCutOff
        :type: float

        Factor used in the cubic step same as :attr:`QuadCutOff`.

    .. attribute:: SmallCost
        :type: float

        Tolerance for which the quadratic interpolation step can be skipped,
        checks :math:`|f_k| < \epsilon |f_0|`.

    .. attribute:: debug
        :type: bool

        Flag to control checks for decreasing function values.
        If *True*, checks that :math:`f_{k + 1} - f_k \le \delta C_k`, where
        :math:`\delta` = :attr:`debugtol`.

    .. attribute:: debugtol
        :type: float

    .. attribute:: step
        :type: float

        Initial step used in the initial line search.

    .. attribute:: max_step
        :type: float

        Maximum step size used in the descent. This is a very crude choice,
        as it mostly sidesteps the line search if large values are encountered.

    .. attribute:: maxit
        :type: int

        Maximum number of iterations.

    .. attribute:: ntries
        :type: int

        Maximum number of times the bracketing interval grows during expansion.

    .. attribute:: ExpandSafe
        :type: float

        Maximum factor by which the secand step increases in the expansion
        phase.

    .. attribute:: SecantAmp
        :type: float

        Factor by which the secant step is amplified during the expansion
        phase.

    .. attribute:: RhoGrow
        :type: float

        Factor by which :math:`rho` grows during the expansion phase.

    .. attribute:: neps
        :type: int

        Maximum number of times that :attr:`eps` can be updated.

    .. attribute:: nshrink
        :type: int

        Maximum number of times the bracketing interval shrinks.

    .. attribute:: nline
        :type: int

        Maximum number of iterations in the line search.

    .. attribute:: restart_fac
        :type: int

        Restart method after ``n * restart_fac`` iterations, where :math:`n`
        is the size of the input.

    .. attribute:: feps
        :type: float

        Tolerance for change in function values.

    .. attribute:: nan_rho
        :type: float

        Growth factor :attr:`RhoGrow` is reset to this value after
        encountering ``nan``.

    .. attribute:: nan_decay
        :type: float

        Decay factor :attr:`Qdecay` is reset to this value after
        encountering ``nan``.

    The following parameters are meant mostly for internal use. They are
    optimized or chosen specifically following results from [HagerZhang2013]_,
    so should only be modified with knowledge.

    .. attribute:: delta
        :type: float

        Parameter for the Wolfe line search in :math:`[0, 0.5]`.

    .. attribute:: sigma
        :type: float

        Parameter for the Wolfe line search in :math:`[\delta, 1]`, where
        :math:`\delta` = :attr:`delta`.

    .. attribute:: gamma
        :type: float

        Decay factor for bracket interval width in line search in :math:`(0, 1)`.

    .. attribute:: rho
        :type: float

        Growth factor in search for initial bracket interval.

    .. attribute:: psi0
        :type: float

        Factor used in starting guess for the line search.

    .. attribute:: psi_lo
        :type: float

    .. attribute:: psi_hi
        :type: float

    .. attribute:: psi2
        :type: float

        In a quadratic step, the bracket interval is given by
        :math:`[\psi_{lo}, \psi_{hi}] \times \psi_2 \times \alpha_{k - 1}`, where
        :math:`\alpha_{k - 1}` is the previous step size.

    .. attribute:: psi1
        :type: float

        If the function is approximately quadratic, this is used to estimate
        the initial step size by :math:`\psi_1 \psi_2 \alpha_{k - 1}`.

    .. attribute:: AdaptiveBeta
        :type: bool

        If *True*, :math:`\theta` is chosen adaptively.

    .. attribute:: BetaLower
        :type: float

        Lower bound for :math:`\beta`.

    .. attribute:: theta
        :type: float

        Describes the family of ``CG_DESCENT`` methods, as described in
        [HagerZhang2006]_.

    .. attribute:: qeps
        :type: float

        Parameter used in cost error estimation for the quadratic restart
        criterion.

    .. attribute:: qrestart
        :type: int

        Number of iterations the function is nearly quadratic before a restart.

    .. attribute:: qrule
        :type: float

        Tolerance used to determine if the cost can be treated as quadratic.
    """

    _changes: ClassVar[dict[str, Any]] = {}

    def __init__(self, **kwargs: Any) -> None:
        super().__init__()

        for k, v in kwargs.items():
            object.__setattr__(self, k, v)

        object.__setattr__(self, "_changes", kwargs)

    def __setattr__(self, k: str, v: Any) -> None:
        raise AttributeError(f"Cannot assign to {k!r}.")

    def replace(self, **changes: Any) -> "OptimizeOptions":
        """Creates a new instance of the same type as *self*, replacing the
        fields with values from *changes*.

        :arg changes: a :class:`dict` of new option values.
        """
        # FIXME: this should just do a deep copy object instead of keeping track
        # of the changes in all instances :(
        kwargs = self._changes.copy()
        kwargs.update(changes)

        return type(self)(**kwargs)

    def __repr__(self) -> str:
        attrs = {k: getattr(self, k) for k in _getmembers(self) if k != "_changes"}
        return f"{type(self).__name__}<{attrs}>"

    def pretty(self) -> str:
        """
        :returns: a string representation of the options in a table.
        """
        attrs = {k: getattr(self, k) for k in _getmembers(self) if k != "_changes"}
        return _stringify_dict(attrs)


# }}}


# {{{ info


@dataclass(frozen=True)
class CallbackInfo:
    it: int
    """Current iteration."""
    alpha: float
    """Step size at current iteration."""
    x: ArrayType
    """Point at which the function and gradient are evaluated."""
    f: float
    """Function value at current iteration."""
    g: ArrayType
    """Gradient (Jacobian) value at the current iteration."""
    d: ArrayType
    """Descent direction at the current iteration. This will usually not be
    the same as the gradient and can be used for debugging.
    """


def _info_from_stats(stats: _cg.cg_iter_stats) -> "CallbackInfo":
    return CallbackInfo(
        it=stats.iter,
        alpha=stats.alpha,
        x=stats.x,
        f=stats.f,
        g=stats.g,
        d=stats.d,
    )


# }}}


# {{{ result


@dataclass(frozen=True)
class OptimizeResult:
    """Based on :class:`scipy.optimize.OptimizeResult`."""

    x: ArrayType
    """Solution of the optimization to the given tolerances."""
    success: bool
    """Flag to denote a successful exit."""
    status: int
    """Termination status of the optimization."""
    message: str
    """Description of the termination status in :attr:`status`."""
    fun: float
    """Function value at the end of the optimization."""
    jac: float
    """Norm of the gradient (Jacobian) at the end of the optimization."""
    nfev: int
    """Number of function evaluations."""
    njev: int
    """Number of gradient (Jacobian) evaluations."""
    nit: int
    """Number of iterations performed by the optimizer."""

    nsubspaceit: int = field(repr=False)
    """Number of subspace iterations (valid if
    :attr:`OptimizeOptions.LBFGS` is *True*).
    """
    nsubspaces: int = field(repr=False)
    """Number of subspace (valid if :attr:`OptimizeOptions.LBFGS` is *True*)."""

    def pretty(self) -> str:
        """Aligned stringify of the results."""
        return _stringify_dict(self.__dict__)


# }}}


# {{{ status messages

# NOTE: these are similar to results with PrintFinal == True
STATUS_TO_MESSAGE = {
    0: "Convergence tolerance for gradient satisfied",
    1: "Change in function value smaller than tolerance: lower 'feps'?",
    2: "Maximum number of iterations limit exceeded",
    3: "Slope is always negative in line search: error in provided functions?",
    4: "Maximum number of line search iterations exceeded: increase 'tol'?",
    5: "Search direction is not a descent direction",
    6: (
        "Excessive updating of estimated error in function values: "
        "increase 'neps'? increase 'tol'?"
    ),
    7: "Wolfe conditions are never satisfied: increase 'eps'?",
    8: "Function values are not improving (with 'debug')",
    9: "No cost or gradient improvement: increase 'nslow'?",
    10: "Out of memory",
    11: "Function value NaN or Inf and cannot be repaired",
    12: "Invalid choice of 'memory' parameter",
    13: "Stopped by user callback",
}

# }}}


# {{{ minimize

FunType: TypeAlias = Callable[[ArrayType], float]
"""This callable takes the current guess ``x`` and returns the function value."""
GradType: TypeAlias = Callable[[ArrayType, ArrayType], None]
"""This callable takes ``(g, x)`` as arguments. The array ``g`` is updated in
place with the gradient at ``x``.
"""
FunGradType: TypeAlias = Callable[[ArrayType, ArrayType], float]
"""This callable takes ``(g, x)`` as arguments and returns the function value.
The array ``g`` needs to be updated in place with the gradient at ``x``.
"""
CallbackType: TypeAlias = Callable[[CallbackInfo], int]
"""Setting the return value to `0` will stop the iteration."""


def wrap_callback(cb: CallbackType | None) -> Callable[[_cg.cg_iter_stats], int] | None:
    from functools import wraps

    if cb is None:
        return cb

    @wraps(cb)
    def wrapper(s: _cg.cg_iter_stats) -> int:
        return cb(_info_from_stats(s))

    return wrapper


def min_work_size(options: OptimizeOptions, n: int) -> int:
    """
    Get recommended size of a *work* array.

    :param options: options used for the optimization.
    :param n: input size.
    """
    m = min(int(options.memory), n)

    if options.memory == 0:
        # original CG_DESCENT without memory
        return 4 * n

    if options.LBFGS or options.memory >= n:
        # LBFGS-based CG_DESCENT
        return 2 * m * (n + 1) + 4 * n

    # limited memory CG_DESCENT
    return (m + 6) * n + (3 * m + 9) * m + 5


def allocate_work_for(options: OptimizeOptions, n: int, dtype: Any = None) -> ArrayType:
    """
    Allocate a *work* array of a recommended size.

    :param options: options used for the optimization
    :param n: input size.
    """
    if dtype is None:
        dtype = np.dtype(np.float64)

    return cast(ArrayType, np.empty(min_work_size(options, n), dtype=dtype))


def minimize(
    fun: FunType,
    x0: ArrayType,
    *,
    jac: GradType,
    funjac: FunGradType | None = None,
    tol: float | None = None,
    options: OptimizeOptions | dict[str, Any] | None = None,
    callback: CallbackType | None = None,
    work: ArrayType | None = None,
    args: tuple[Any, ...] = (),
) -> OptimizeResult:
    """
    :param fun: a :class:`~FunType` that returns the function value at ``x``.
    :param x0: initial guess.
    :param jac: a :class:`~GradType` that computes the gradient of *fun* at ``x``.
        The gradient is stored in place.
    :param funjac: a :class:`~FunGradType` that computes both
        the function value and gradient at ``x``. This function can be used
        to improve efficiency by evaluating common parts just once, but is not
        required.
    :param tol: tolerance used to check convergence. Exact meaning depends on
        :attr:`OptimizeOptions.StopRule`.
    :param options: options used by the algorithm.
    :param callback: a :class:`~CallbackType` called at the end of each
        (successful) iteration.
    :param work: additional work array (see :func:`allocate_work_for`).
    :param args: additional arguments passed to the callables.
    """

    # {{{ setup

    if args:
        raise ValueError("Using 'args' is not supported.")

    if tol is None:
        tol = 1.0e-8

    if options is not None:
        if isinstance(options, dict):
            param = OptimizeOptions(**options)
        elif isinstance(options, OptimizeOptions):
            param = options
        else:
            raise TypeError(f"Unknown 'options' type: {type(options).__name__!r}.")
    else:
        param = OptimizeOptions()

    if work is not None:
        m = min_work_size(param, x0.size)
        if work.size >= m:
            raise ValueError(f"'work' must have size >= {m}.")

    # }}}

    # {{{ optimize

    x, stats, status = _cg.cg_descent(
        x0,
        tol,
        param,
        fun,
        jac,
        funjac,
        wrap_callback(callback),
        work,
    )

    # }}}

    return OptimizeResult(
        x=x,
        success=status == 0,
        status=int(status),
        message=STATUS_TO_MESSAGE[status],
        fun=stats.f,
        jac=stats.gnorm,
        nfev=stats.nfunc,
        njev=stats.ngrad,
        nit=stats.iter,
        nsubspaceit=stats.IterSub,
        nsubspaces=stats.NumSub,
    )


# }}}


# {{{


class Timer:
    """A simpler timer context manager.

    .. code:: python

        with Timer() as time:
            # perform some operations
            ...

        print(time)

    .. autoattribute:: t_start
    .. autoattribute:: t_end
    .. autoattribute:: t_wall
    """

    def __init__(self) -> None:
        self.t_start: float = -1.0
        """
        A time stamp (obtained from :func:`time.perf_counter`) for when the
        context manager was entered. This is not updated again until the context
        manager is entered again.
        """
        self.t_end: float = -1.0
        """
        A time stamp (obtained from :func:`time.perf_counter`) for when the
        context manager has exited. This is not valid in the ``with`` block
        itself.
        """
        self.t_wall: float = -1.0
        """The total wall time between the context manager enter and exit. This
        is only valid after exiting the context manager.
        """

    def __enter__(self) -> Timer:
        import time

        self.t_start = time.perf_counter()
        self.t_end = -1.0
        self.t_wall = -1.0

        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        import time

        self.t_end = time.perf_counter()
        self.t_wall = self.t_end - self.t_start

    def __str__(self) -> str:
        # NOTE: this matches the output of MATLAB's toc
        return f"Elapsed time is {self.t_wall:.5f} seconds."

    def __repr__(self) -> str:
        return f"{type(self).__name__}(t_start={self.t_start}, t_end={self.t_end})"


# }}}


# {{{ logger


def get_logger(
    module: str | None = None,
    *,
    level: int | str | None = None,
) -> logging.Logger:
    if isinstance(level, str):
        try:
            level = getattr(logging, level.upper())
        except AttributeError:
            level = None

    if level is None:
        level = logging.INFO

    try:
        from rich.logging import RichHandler
    except ImportError:
        try:
            # NOTE: rich is vendored by pip, so try and get it from there
            from pip._vendor.rich.logging import RichHandler  # type: ignore[assignment]
        except ImportError:
            from logging import StreamHandler as RichHandler  # type: ignore[assignment]

    logger = logging.getLogger()
    logger.setLevel(level)
    logger.addHandler(RichHandler())

    return logger


# }}}
