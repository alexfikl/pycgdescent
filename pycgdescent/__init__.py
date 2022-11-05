# SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

from contextlib import contextmanager
from dataclasses import dataclass, field
from importlib import metadata
from typing import (
    Any,
    Callable,
    Dict,
    Iterator,
    List,
    Optional,
    Tuple,
    Union,
    TYPE_CHECKING,
)

import numpy as np

import pycgdescent._cg_descent as _cg

import logging

logger = logging.getLogger()

__version__ = metadata.version("pycgdescent")

if TYPE_CHECKING:
    ArrayType = np.ndarray[Any, np.dtype[np.float64]]
else:
    ArrayType = np.ndarray

__doc__ = """
Functions
^^^^^^^^^

.. autofunction:: minimize

.. autofunction:: min_work_size
.. autofunction:: allocate_work_for

.. autoclass:: OptimizeOptions
    :no-show-inheritance:
    :exclude-members: __init__, __new__

.. autoclass:: OptimizeResult
    :no-show-inheritance:
    :exclude-members: __init__, __new__

.. autoclass:: CallbackInfo
    :no-show-inheritance:
    :exclude-members: __init__, __new__

Type Aliases
^^^^^^^^^^^^

.. class:: FunType

    ``Callable[[numpy.ndarray], float]``. This callable takes the current
    guess ``x`` and returns the function value.

.. class:: GradType

    ``Callable[[numpy.ndarray, numpy.ndarray], None]``. This callable takes
    ``(g, x)`` as arguments. The array ``g`` needs to be updated in place
    with the gradient at ``x``.

.. class:: FunGradType

    ``Callable[[numpy.ndarray, numpy.ndarray], float]``. This callable takes
    ``(g, x)`` as arguments and returns the function value. The array ``g``
    needs to be updated in place with the gradient at ``x``.

.. class:: CallbackType

    ``Callable[[CallbackInfo], int]``. Setting the return value to `0` will
    stop the iteration.
"""


# {{{ options


def _getmembers(obj: object) -> List[str]:
    import inspect

    return [
        m
        for m in dir(obj)
        if not m.startswith("__") and not inspect.ismethod(getattr(obj, m))
    ]


def _stringify_dict(d: Dict[str, Any]) -> str:
    width = len(max(d, key=len))
    fmt = f"{{:{width}}} : {{}}"

    items = sorted({k: repr(v) for k, v in d.items()}.items())

    return "\n".join(["\t" + "\n\t".join(fmt.format(k, v) for k, v in items)])


class OptimizeOptions(_cg.cg_parameter):
    r"""Optimization options for the ``CG_DESCENT`` algorithm. A description
    of some of the more technical parameters can be found in the paper
    [HagerZhang2006]_ and [HagerZhang2013]_.

    .. [HagerZhang2006] W. W. Hager, H. Zhang,
        *Algorithm 851: CG_DESCENT, a Conjugate Gradient Method With Guaranteed
        Descent*,
        ACM Transactions on Mathematical Software, Vol. 32, pp. 113--137, 2006,
        `DOI <http://dx.doi.org/10.1145/1132973.1132979>`__.

    .. [HagerZhang2013] W. W. Hager, H. Zhang,
        *The Limited Memory Conjugate Gradient Method*,
        SIAM Journal on Optimization, Vol. 23, pp. 2150--2168, 2013,
        `DOI <http://dx.doi.org/10.1137/120898097>`__.

    .. automethod:: replace
    .. automethod:: pretty

    The attribute names here follow those of the original ``CG_DESCENT`` code.

    .. attribute:: PrintLevel

        Display level as an integer in `[0, 1, 2, 3]`.

    .. attribute:: LBFGS

        Boolean flag to handle use of LBFGS. If *False*, LBFGS is only used
        when the :attr:`memory` is larger than the input size :math:`n`.

    .. attribute:: memory

        Number of vectors stored in memory.

    .. attribute:: SubCheck
    .. attribute:: SubSkip

        Together with :attr:`SubCheck`, it controls the frequency with which
        the subspace condition is checked. It is checked at every
        ``SubCheck * memory`` iterations and, if not satisfied, then
        it is skipped for ``SubSkip * memory`` iterations and :attr:`SubSkip`
        is doubled. Whenever the subspace condition is satisfied, :attr:`SubSkip`
        is returned to the original value.

    .. attribute:: eta0

        Controls relative distance from current gradient to the subspace.
        If the distance is ``<= eta0`` and the subspace dimension is
        :attr:`memory`, then the subspace is entered. This is used as
        :math:`\eta_0^2` in Equation 3.4 from [HagerZhang2013]_.

    .. attribute:: eta1

        Controls relative distance from current gradient to the subspace.
        If the distance is ``>= eta1``, the subspace is exited. This is used
        as :math:`\eta_1^2` in Equation 3.4 in [HagerZhang2013]_.

    .. attribute:: eta2

        Controls relative distance from current descent direction to the subspace.
        If the distance is ``<= eta2``, the subspace is always entered.

    .. attribute:: AWolfe
    .. attribute:: AWolfeFac

        If :attr:`AWolfe` is *True*, then the approximate Wolfe condition is
        used when :math:`|f_{k + 1} - f_k| < \omega C_k`, for
        :math:`\omega` = :attr:`AWolfeFac`. See discussion surrounding
        Equation 25 in [HagerZhang2006]_.

    .. attribute:: Qdecay

        Factor in :math:`[0, 1]` used to compute the average cost magnitude
        in the Wolfe condition.

    .. attribute:: nslow

        Maximum number of "slow" iterations without strict improvement in
        either function values or gradient.

    .. attribute:: StopRule

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

        Factor used in stopping condition.

    .. attribute:: PertRule

        Estimate error in function values. If *False*, use just :attr:`eps`,
        otherwise use :math:`\epsilon C_k`, where :math:`C_k` is defined in
        Equation 26 in [HagerZhang2006]_.

    .. attribute:: eps

        Factor used in :attr:`PertRule`.

    .. attribute:: egrow

        Factor by which :attr:`eps` grows when the line search fails during
        contraction.

    .. attribute:: QuadStep

        If *False*, do not use a quadratic interpolation step in the line
        search. If *True*, attempt a step based on
        :math:`\epsilon = ` :attr:`QuadCutOff` when

        .. math::

            \frac{|f_{k + 1} - f_k|}{|f_k|} > \epsilon.

    .. attribute:: QuadCutOff

        Factor used when :attr:`QuadStep` is *True*.

    .. attribute:: QuadSafe

        Maximum factor by which a quadratic step can reduce the step size.

    .. attribute:: UseCubic

        Boolean flag that enables a cubic step in the line search.

    .. attribute:: CubicCutOff

        Factor used in the cubic step same as :attr:`QuadCutOff`.

    .. attribute:: SmallCost

        Tolerance for which the quadratic interpolation step can be skipped,
        checks :math:`|f_k| < \epsilon |f_0|`.

    .. attribute:: debug

        Flag to control checks for decreasing function values.
        If *True*, checks that :math:`f_{k + 1} - f_k \le \delta C_k`, where
        :math:`\delta` = :attr:`debugtol`.

    .. attribute:: debugtol

    .. attribute:: step

        Initial step used in the initial line search.

    .. attribute:: max_step

        Maximum step size used in the descent. This is a very crude choice,
        as it mostly sidesteps the line search if large values are encountered.

    .. attribute:: maxit

        Maximum number of iterations.

    .. attribute:: ntries

        Maximum number of times the bracketing interval grows during expansion.

    .. attribute:: ExpandSafe

        Maximum factor by which the secand step increases in the expansion
        phase.

    .. attribute:: SecantAmp

        Factor by which the secant step is amplified during the expansion
        phase.

    .. attribute:: RhoGrow

        Factor by which :math:`rho` grows during the expansion phase.

    .. attribute:: neps

        Maximum number of times that :attr:`eps` can be updated.

    .. attribute:: nshrink

        Maximum number of times the bracketing interval shrinks.

    .. attribute:: nline

        Maximum number of iterations in the line search.

    .. attribute:: restart_fac

        Restart method after ``n * restart_fac`` iterations, where :math:`n`
        is the size of the input.

    .. attribute:: feps

        Tolerance for change in function values.

    .. attribute:: nan_rho

        Growth factor :attr:`RhoGrow` is reset to this value after
        encountering ``nan``.

    .. attribute:: nan_decay

        Decay factor :attr:`Qdecay` is reset to this value after
        encountering ``nan``.

    The following parameters are meant mostly for internal use. They are
    optimized or chosen specifically following results from [HagerZhang2013]_,
    so should only be modified with knowledge.

    .. attribute:: delta

        Parameter for the Wolfe line search in :math:`[0, 0.5]`.

    .. attribute:: sigma

        Parameter for the Wolfe line search in :math:`[\delta, 1]`, where
        :math:`\delta` = :attr:`delta`.

    .. attribute:: gamma

        Decay factor for bracket interval width in line search in :math:`(0, 1)`.

    .. attribute:: rho

        Growth factor in search for initial bracket interval.

    .. attribute:: psi0

        Factor used in starting guess for the line search.

    .. attribute:: psi_lo
    .. attribute:: psi_hi
    .. attribute:: psi2

        In a quadratic step, the bracket interval is given by
        :math:`[\psi_{lo}, \psi_{hi}] \times \psi_2 \times \alpha_{k - 1}`, where
        :math:`\alpha_{k - 1}` is the previous step size.

    .. attribute:: psi1

        If the function is approximately quadratic, this is used to estimate
        the initial step size by :math:`\psi_1 \psi_2 \alpha_{k - 1}`.

    .. attribute:: AdaptiveBeta

        If *True*, :math:`\theta` is chosen adaptively.

    .. attribute:: BetaLower

        Lower bound for :math:`\beta`.

    .. attribute:: theta

        Describes the family of ``CG_DESCENT`` methods, as described in
        [HagerZhang2006]_.

    .. attribute:: qeps

        Parameter used in cost error estimation for the quadratic restart
        criterion.

    .. attribute:: qrestart

        Number of iterations the function is nearly quadratic before a restart.

    .. attribute:: qrule

        Tolerance used to determine if the cost can be treated as quadratic.
    """

    _changes: Dict[str, Any] = {}

    def __init__(self, **kwargs: Any) -> None:
        super().__init__()

        for k, v in kwargs.items():
            object.__setattr__(self, k, v)

        object.__setattr__(self, "_changes", kwargs)

    def __setattr__(self, k: str, v: Any) -> None:
        raise AttributeError(f"cannot assign to '{k}'")

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
    """
    .. attribute:: it

        Current iteration.

    .. attribute:: alpha

        Step size at the current iteration.

    .. attribute:: x

        Point at which the function and gradient are evaluated.

    .. attribute:: f

        Function value at the current iteration.

    .. attribute:: g

        Gradient (Jacobian) value at the current iteration.

    .. attribute:: d

        Descent direction at the current iteration. This will usually not
        be the same as the gradient and can be used for debugging.
    """

    it: int
    alpha: float
    x: ArrayType
    f: float
    g: ArrayType
    d: ArrayType


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
    """Based on :class:`scipy.optimize.OptimizeResult`.

    .. attribute:: x

        Solution of the optimization.

    .. attribute:: success

        Flag to denote a successful exit.

    .. attribute:: status

        Termination status of the optimize.

    .. attribute:: message

        Description of the termination status in :attr:`status`.

    .. attribute:: fun

        Function value at the end of the optimization.

    .. attribute:: jac

        Norm of the gradient (Jacobian) at the end of the optimization.

    .. attribute:: nfev

        Number of function evaluations.

    .. attribute:: njev

        Number of gradient (Jacobian) evaluation.

    .. attribute:: nit

        Number of iterations performed by the optimizer.
    """

    x: ArrayType
    success: bool
    status: int
    message: str
    fun: float
    jac: float
    nfev: int
    njev: int
    nit: int

    nsubspaceit: int = field(repr=False)
    nsubspaces: int = field(repr=False)

    def pretty(self) -> str:
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

FunType = Callable[[ArrayType], float]
GradType = Callable[[ArrayType, ArrayType], None]
FunGradType = Callable[[ArrayType, ArrayType], float]
CallbackType = Callable[[CallbackInfo], int]


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


def allocate_work_for(
    options: OptimizeOptions, n: int, dtype: Any = np.float64
) -> ArrayType:
    """
    Allocate a *work* array of a recommended size.

    :param options: options used for the optimization
    :param n: input size.
    """
    return np.empty(min_work_size(options, n), dtype=dtype)


def minimize(
    fun: "FunType",
    x0: ArrayType,
    *,
    jac: "GradType",
    funjac: Optional["FunGradType"] = None,
    tol: Optional[float] = None,
    options: Optional[Union[OptimizeOptions, Dict[str, Any]]] = None,
    callback: Optional["CallbackType"] = None,
    work: Optional[ArrayType] = None,
    args: Tuple[Any, ...] = (),
) -> OptimizeResult:
    """
    :param fun: a :class:`~collections.abc.Callable` that returns the
        function value at ``x``.
    :param x0: initial guess.
    :param jac: a :class:`~collections.abc.Callable` that computes the
        gradient of *fun* at ``x``. The gradient is stored in place.
    :param funjac: a :class:`~collections.abc.Callable` that computes both
        the function value and gradient at ``x``. This function can be used
        to improve efficiency by evaluating common parts just once, but is not
        required.
    :param tol: tolerance used to check convergence. Exact meaning depends on
        :attr:`OptimizeOptions.StopRule`.
    :param options: options used by the algorithm.
    :param callback: a :class:`~collections.abc.Callable` called at the end
        of each (successful) iteration.
    :param work: additional work array (see :func:`allocate_work_for`).
    :param args: additional arguments passed to the callables.
    """

    # {{{ setup

    if args:
        raise ValueError("using 'args' is not supported")

    if tol is None:
        tol = 1.0e-8

    if options is not None:
        if isinstance(options, dict):
            param = OptimizeOptions(**options)
        elif isinstance(options, OptimizeOptions):
            param = options
        else:
            raise TypeError(f"unknown 'options' type: {type(options).__name__}")
    else:
        param = OptimizeOptions()

    if work is not None:
        m = min_work_size(param, x0.size)
        if work.size >= m:
            raise ValueError(f"'work' must have size >= {m}")

    wrapped_callback: Optional[Callable[[Any], int]] = None
    if callback is not None:

        def wrapped_callback(s: _cg.cg_iter_stats) -> int:
            return callback(_info_from_stats(s))  # type: ignore

    # }}}

    # {{{ optimize

    x, stats, status = _cg.cg_descent(
        x0, tol, param, fun, jac, funjac, wrapped_callback, work
    )

    # }}}

    return OptimizeResult(
        x=x,
        success=status == 0,
        status=status,
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


@contextmanager
def timer(name: str = "timer") -> Iterator[None]:
    import time

    t_start = time.time()
    yield
    t_end = time.time()
    logger.info("%s: %gs", name, t_end - t_start)


# }}}
