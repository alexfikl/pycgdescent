import numpy as np
import pycgdescent as cg


def fn(x):
    t = np.arange(x.size)
    return np.sum(np.exp(x) - np.sqrt(t + 1) * x)


def grad(g, x):
    t = np.arange(x.size)
    g = np.exp(x) - np.sqrt(t + 1)


def fngrad(g, x):
    f = fn(x)
    g = grad(g, x)
    return f


n = 10
x0 = np.ones(n, dtype=np.float64)

param = cg.cg_parameter()
x, stats = cg.cg_descent(x0, param, 1.0e-8, fn, grad, fngrad)
print(stats)
print(x)
1/0
