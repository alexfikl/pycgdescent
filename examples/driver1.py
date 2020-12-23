import numpy as np
import pycgdescent as cg


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


n = 10
x0 = np.ones(n, dtype=np.float64)
t = np.sqrt(1 + np.arange(n))

param = cg.cg_parameter()
x, stats = cg.cg_descent(x0, param, 1.0e-8, fn, grad, fngrad)
print(stats)
print(x)
