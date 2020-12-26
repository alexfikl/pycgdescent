from pycgdescent._cg_descent import (   # noqa: F401
        cg_stats,
        cg_iter_stats,
        cg_parameter,

        cg_default,
        cg_descent)


# NOTE: these are similar to results with PrintFinal == True
STATUS_TO_MESSAGE = {
    0: "Convergence tolerance for gradient satisfied",
    1: "Change in function value smaller than tolerance: lower 'feps'?",
    2: "Maximum number of iterations limit exceeded",
    3: "Slope is always negative in line search: error in provided functions?",
    4: "Maximum number of line search iterations exceeded: increase 'tol'?",
    5: "Search direction is not a descent direction",
    6: ("Excessive updating of estimated error in function values: "
        "increase 'neps'? increase 'tol'?"),
    7: "Wolfe conditions are never satisfied: increase 'eps'?",
    8: "Function values are not improving (with 'debug')",
    9: "No cost or gradient improvement: increase 'nslow'?",
    10: "Out of memory",
    11: "Function value NaN or Inf and cannot be repaired",
    12: "Invalid choice of 'memory' parameter",
    13: "Stopped by user callback",
}
