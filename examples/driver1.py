import time
from contextlib import contextmanager

import numpy as np
import pycgdescent._private as _cg

@contextmanager
def timer():
    t_start = time.time()
    yield
    t_end = time.time()
    print("elapsed: ", t_end - t_start)


def fn(x):
    f = np.sum(np.exp(x) - t * x)
    return f


def grad(g, x):
    g[...] = np.exp(x) - t


def fngrad(g, x):
    y = np.exp(x)
    f = np.sum(y - t * x)
    g[...] = y - t
    return f


n = 100
x0 = np.ones(n, dtype=np.float64)
t = np.sqrt(1 + np.arange(n))

with timer():
    x1, stats, flag = _cg.cg_descent(x0, 1.0e-8, None, fn, grad, None, None)

with timer():
    x2, stats, status = _cg.cg_descent(x0, 1.0e-8, None, fn, grad, fngrad, None)

print("status:  ", status)
print("error:   ", np.linalg.norm(x1 - x2) / np.linalg.norm(x2))

print()
print("maximum norm for gradient: %+.16e" % stats.gnorm)
print("function value:            %+.16e" % stats.f)
print("cg iterations:            ", stats.iter)
print("function evaluations:     ", stats.nfunc)
print("gradient evaluations:     ", stats.ngrad)
