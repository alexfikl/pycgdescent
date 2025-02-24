from __future__ import annotations
import numpy
import typing

__all__ = ["cg_default", "cg_descent", "cg_iter_stats", "cg_parameter", "cg_stats"]

class cg_iter_stats:
    @property
    def alpha(self) -> float: ...
    @property
    def d(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]: ...
    @property
    def f(self) -> float: ...
    @property
    def g(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]: ...
    @property
    def iter(self) -> int: ...
    @property
    def x(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]: ...

class cg_parameter:
    AWolfe: int
    AWolfeFac: float
    AdaptiveBeta: int
    BetaLower: float
    CubicCutOff: float
    ExpandSafe: float
    LBFGS: int
    PertRule: int
    PrintFinal: int
    PrintLevel: int
    PrintParms: int
    Qdecay: float
    QuadCutOff: float
    QuadSafe: float
    QuadStep: int
    SecantAmp: float
    SmallCost: float
    StopFac: float
    StopRule: int
    SubCheck: int
    SubSkip: int
    UseCubic: int
    debug: int
    debugtol: float
    delta: float
    egrow: float
    eps: float
    eta0: float
    eta1: float
    eta2: float
    feps: float
    gamma: float
    max_step: float
    maxit: int
    memory: int
    nan_decay: float
    nan_rho: float
    neps: int
    nline: int
    nshrink: int
    nslow: int
    ntries: int
    psi0: float
    psi1: float
    psi2: float
    psi_hi: float
    psi_lo: float
    qeps: float
    qrestart: int
    qrule: float
    restart_fac: float
    rho: float
    sigma: float
    step: float
    theta: float
    def __init__(self) -> None: ...

class cg_stats:
    @property
    def IterSub(self) -> int: ...
    @property
    def NumSub(self) -> int: ...
    @property
    def f(self) -> float: ...
    @property
    def gnorm(self) -> float: ...
    @property
    def iter(self) -> int: ...
    @property
    def nfunc(self) -> int: ...
    @property
    def ngrad(self) -> int: ...

def cg_default(arg0: typing.Any) -> None: ...
def cg_descent(
    x: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    grad_tol: float,
    param: cg_parameter | None,
    value: typing.Callable[
        [numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]], float
    ],
    grad: typing.Callable[
        [
            numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
            numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        ],
        None,
    ],
    valgrad: typing.Callable[
        [
            numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
            numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        ],
        float,
    ]
    | None,
    callback: typing.Callable[[cg_iter_stats], int] | None,
    work: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]] | None,
) -> tuple[numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], cg_stats, bool]: ...
