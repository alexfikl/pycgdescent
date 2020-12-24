from _cg_descent import (   # noqa: F401
        cg_stats,
        cg_parameter,
        cg_default,
        cg_descent)


STATUS_TO_MESSAGE = {
    0: "Convergence tolerance satisfied",
    1: "Change in function value smaller than tolerance (see 'feps')",
    2: "Maximum number of iterations exceeded",
    3: "Slope is always negative in line search",
    4: "Maximum number of line search iterations exceeded",
    5: "Search direction is not a descent direction",
    6: "Excessive updating of estimated error in function values (see 'neps')",
    7: "Wolfe conditions are never satisfied",
    8: "Debugger is on and the function value increases",
    9: "No cost or gradient improvement (see 'nslow')",
    10: "Out of memory",
    11: "Function value NaN or Inf and cannot be repaired",
    12: "Invalid choice of 'memory' parameter",
}
