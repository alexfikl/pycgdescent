// SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
//
// SPDX-License-Identifier: MIT

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>

#include <cstring>
#include <functional>
#include <iostream>
#include <optional>

#include "cg_user.h"

namespace nb = nanobind;

// {{{ macros

#define WRAP_RAW_POINTER(NAME, RAWNAME, SIZE)                                   \
    auto NAME = nb::ndarray<nb::numpy, double>(                                 \
        RAWNAME, {(size_t)(SIZE)}, nb::capsule(RAWNAME, [](void *) noexcept {}) \
    );

#define DEF_RO_PROPERTY(NAME) def_prop_ro(#NAME, &cl::get_##NAME)

#define DEF_RO_ARRAY_PROPERTY(NAME) def_prop_ro(#NAME, &cl::get_##NAME, nb::rv_policy::automatic)

#define DEF_PROPERTY(NAME) def_prop_rw(#NAME, &cl::get_##NAME, &cl::set_##NAME)

#define CLASS_PROPERTY(NAME, TYPE)     \
    TYPE get_##NAME() const noexcept { \
        return obj.NAME;               \
    };                                 \
    void set_##NAME(TYPE v) {          \
        obj.NAME = v;                  \
    };

#define CLASS_RO_PROPERTY(NAME, TYPE)  \
    TYPE get_##NAME() const noexcept { \
        return obj.NAME;               \
    };

#define CLASS_RO_ARRAY_PROPERTY(NAME, TYPE)                                          \
    nb::ndarray<nb::numpy, double> get_##NAME() const {                              \
        return nb::ndarray<nb::numpy, double>(                                       \
            obj.NAME, {(size_t)obj.n}, nb::capsule(obj.NAME, [](void *) noexcept {}) \
        );                                                                           \
    };

// }}}

// {{{ cg_parameter wrapper

class cg_parameter_wrapper {
   public:
    cg_parameter_wrapper(cg_parameter * p) {
        std::memcpy(&obj, p, sizeof(obj));
    }
    cg_parameter_wrapper() {
        std::memset(&obj, 0, sizeof(obj));
        cg_default(&obj);
    }
    cg_parameter_wrapper(const cg_parameter_wrapper & w) {
        memcpy(&obj, &w.obj, sizeof(obj));
    }

    CLASS_PROPERTY(PrintFinal, int)
    CLASS_PROPERTY(PrintLevel, int)
    CLASS_PROPERTY(PrintParms, int)
    CLASS_PROPERTY(LBFGS, int)
    CLASS_PROPERTY(memory, int)
    CLASS_PROPERTY(SubCheck, int)
    CLASS_PROPERTY(SubSkip, int)
    CLASS_PROPERTY(eta0, double)
    CLASS_PROPERTY(eta1, double)
    CLASS_PROPERTY(eta2, double)
    CLASS_PROPERTY(AWolfe, int)
    CLASS_PROPERTY(AWolfeFac, double)
    CLASS_PROPERTY(Qdecay, double)
    CLASS_PROPERTY(nslow, int)
    CLASS_PROPERTY(StopRule, int)
    CLASS_PROPERTY(StopFac, double)
    CLASS_PROPERTY(PertRule, int)
    CLASS_PROPERTY(eps, double)
    CLASS_PROPERTY(egrow, double)
    CLASS_PROPERTY(QuadStep, int)
    CLASS_PROPERTY(QuadCutOff, double)
    CLASS_PROPERTY(QuadSafe, double)
    CLASS_PROPERTY(UseCubic, int)
    CLASS_PROPERTY(CubicCutOff, double)
    CLASS_PROPERTY(SmallCost, double)
    CLASS_PROPERTY(debug, int)
    CLASS_PROPERTY(debugtol, double)
    CLASS_PROPERTY(step, double)
    CLASS_PROPERTY(max_step, double)
    CLASS_PROPERTY(maxit, INT)
    CLASS_PROPERTY(ntries, int)
    CLASS_PROPERTY(ExpandSafe, double)
    CLASS_PROPERTY(SecantAmp, double)
    CLASS_PROPERTY(RhoGrow, double)
    CLASS_PROPERTY(neps, int)
    CLASS_PROPERTY(nshrink, int)
    CLASS_PROPERTY(nline, int)
    CLASS_PROPERTY(restart_fac, double)
    CLASS_PROPERTY(feps, double)
    CLASS_PROPERTY(nan_rho, double)
    CLASS_PROPERTY(nan_decay, double)

    CLASS_PROPERTY(delta, double)
    CLASS_PROPERTY(sigma, double)
    CLASS_PROPERTY(gamma, double)
    CLASS_PROPERTY(rho, double)
    CLASS_PROPERTY(psi0, double)
    CLASS_PROPERTY(psi_lo, double)
    CLASS_PROPERTY(psi_hi, double)
    CLASS_PROPERTY(psi1, double)
    CLASS_PROPERTY(psi2, double)
    CLASS_PROPERTY(AdaptiveBeta, int)
    CLASS_PROPERTY(BetaLower, double)
    CLASS_PROPERTY(theta, double)
    CLASS_PROPERTY(qeps, double)
    CLASS_PROPERTY(qrule, double)
    CLASS_PROPERTY(qrestart, int)

    cg_parameter obj;
};

// }}}

// {{{ cg_stats wrapper

class cg_stats_wrapper {
   public:
    cg_stats_wrapper() {
        memset(&obj, 0, sizeof(obj));
    };

    CLASS_RO_PROPERTY(f, double)
    CLASS_RO_PROPERTY(gnorm, double)
    CLASS_RO_PROPERTY(iter, INT)
    CLASS_RO_PROPERTY(IterSub, INT)
    CLASS_RO_PROPERTY(NumSub, INT)
    CLASS_RO_PROPERTY(nfunc, INT)
    CLASS_RO_PROPERTY(ngrad, INT)

    cg_stats obj;
};

// }}}

// {{{ cg_iter_stats wrapper

class cg_iter_stats_wrapper {
   public:
    cg_iter_stats_wrapper(cg_iter_stats * stats) {
        std::memcpy(&obj, stats, sizeof(obj));
    };
    cg_iter_stats_wrapper() {
        std::memset(&obj, 0, sizeof(obj));
    };

    CLASS_RO_PROPERTY(iter, INT)
    CLASS_RO_PROPERTY(alpha, double)
    CLASS_RO_ARRAY_PROPERTY(x, double)
    CLASS_RO_PROPERTY(f, double)
    CLASS_RO_ARRAY_PROPERTY(g, double)
    CLASS_RO_ARRAY_PROPERTY(d, double)

    cg_iter_stats obj;
};

// }}}}

// {{{ cg_descent wrapper

namespace cg {

typedef nb::ndarray<nb::numpy, double> ndarray;

typedef std::function<double(cg::ndarray)> value_fn;
typedef std::function<void(cg::ndarray, cg::ndarray)> grad_fn;
typedef std::function<double(cg::ndarray, cg::ndarray)> valgrad_fn;
typedef std::function<int(cg_iter_stats_wrapper &)> callback_fn;

class FnWrapper {
   public:
    FnWrapper(value_fn * value, grad_fn * grad, valgrad_fn * valgrad, callback_fn * callback)
        : m_value(value), m_grad(grad), m_valgrad(valgrad), m_callback(callback) {};

    ~FnWrapper() {};

    value_fn * m_value;
    grad_fn * m_grad;
    valgrad_fn * m_valgrad;
    callback_fn * m_callback;
};

};  // namespace cg

double user_value(double * _x, INT n, void * User) {
    cg::FnWrapper * w = static_cast<cg::FnWrapper *>(User);
    WRAP_RAW_POINTER(x, _x, n);

    return (*w->m_value)(x);
}

void user_grad(double * _g, double * _x, INT n, void * User) {
    cg::FnWrapper * w = static_cast<cg::FnWrapper *>(User);
    WRAP_RAW_POINTER(g, _g, n);
    WRAP_RAW_POINTER(x, _x, n);

    (*w->m_grad)(g, x);
}

double user_valgrad(double * _g, double * _x, INT n, void * User) {
    cg::FnWrapper * w = static_cast<cg::FnWrapper *>(User);
    WRAP_RAW_POINTER(g, _g, n);
    WRAP_RAW_POINTER(x, _x, n);

    return (*w->m_valgrad)(g, x);
}

int user_callback(cg_iter_stats * IterStats, void * User) {
    cg::FnWrapper * w = static_cast<cg::FnWrapper *>(User);
    cg_iter_stats_wrapper wi(IterStats);

    return (*w->m_callback)(wi);
}

std::tuple<cg::ndarray, cg_stats_wrapper, bool> cg_descent_wrapper(
    cg::ndarray x,
    double grad_tol,
    std::optional<cg_parameter_wrapper *> param,
    cg::value_fn & value,
    cg::grad_fn & grad,
    std::optional<cg::valgrad_fn> valgrad,
    std::optional<cg::callback_fn> callback,
    std::optional<cg::ndarray> work
) {
    int status = 0;
    cg_stats_wrapper stats = cg_stats_wrapper();
    cg_parameter * p = param.has_value() ? &param.value()->obj : nullptr;
    double * workptr = (work.has_value() ? static_cast<double *>(work.value().data()) : nullptr);

    int n = (int)x.shape(0);
    double * ptr = new double[n];

    auto xptr = static_cast<double *>(x.data());
    std::memcpy(ptr, xptr, n * sizeof(double));

    cg::FnWrapper w(
        &value,
        &grad,
        valgrad.has_value() ? &valgrad.value() : nullptr,
        callback.has_value() ? &callback.value() : nullptr
    );
    auto * user_valgrad_p = valgrad.has_value() ? user_valgrad : nullptr;
    auto * user_callback_p = callback.has_value() ? user_callback : nullptr;

    status = cg_descent(
        ptr,
        n,
        &stats.obj,
        p,
        grad_tol,
        user_value,
        user_grad,
        user_valgrad_p,
        user_callback_p,
        workptr,
        &w
    );

    nb::capsule owner(ptr, [](void * p) noexcept { delete[] static_cast<double *>(p); });
    return std::make_tuple(cg::ndarray(ptr, {(size_t)n}, owner), std::move(stats), status);
}

// }}}

// {{{ cg_default wrapper

void cg_default_wrapper(nb::object param) {
    cg_default(&nb::cast<cg_parameter_wrapper *>(param)->obj);
}

// }}}

NB_MODULE(_cg_descent, m) {
    {
        typedef cg_parameter_wrapper cl;
        nb::class_<cl>(m, "cg_parameter")
            .def(nb::init<>())
            .DEF_PROPERTY(PrintFinal)
            .DEF_PROPERTY(PrintLevel)
            .DEF_PROPERTY(PrintParms)
            .DEF_PROPERTY(LBFGS)
            .DEF_PROPERTY(memory)
            .DEF_PROPERTY(SubCheck)
            .DEF_PROPERTY(SubSkip)
            .DEF_PROPERTY(eta0)
            .DEF_PROPERTY(eta1)
            .DEF_PROPERTY(eta2)
            .DEF_PROPERTY(AWolfe)
            .DEF_PROPERTY(AWolfeFac)
            .DEF_PROPERTY(Qdecay)
            .DEF_PROPERTY(nslow)
            .DEF_PROPERTY(StopRule)
            .DEF_PROPERTY(StopFac)
            .DEF_PROPERTY(PertRule)
            .DEF_PROPERTY(eps)
            .DEF_PROPERTY(egrow)
            .DEF_PROPERTY(QuadStep)
            .DEF_PROPERTY(QuadCutOff)
            .DEF_PROPERTY(QuadSafe)
            .DEF_PROPERTY(UseCubic)
            .DEF_PROPERTY(CubicCutOff)
            .DEF_PROPERTY(SmallCost)
            .DEF_PROPERTY(debug)
            .DEF_PROPERTY(debugtol)
            .DEF_PROPERTY(step)
            .DEF_PROPERTY(max_step)
            .DEF_PROPERTY(maxit)
            .DEF_PROPERTY(ntries)
            .DEF_PROPERTY(ExpandSafe)
            .DEF_PROPERTY(SecantAmp)
            .DEF_PROPERTY(neps)
            .DEF_PROPERTY(nshrink)
            .DEF_PROPERTY(nline)
            .DEF_PROPERTY(restart_fac)
            .DEF_PROPERTY(feps)
            .DEF_PROPERTY(nan_rho)
            .DEF_PROPERTY(nan_decay)
            // NOTE: these are not recommended to be played with, but we're making
            // them readwrite anyway for "power users"
            .DEF_PROPERTY(rho)
            .DEF_PROPERTY(delta)
            .DEF_PROPERTY(sigma)
            .DEF_PROPERTY(gamma)
            .DEF_PROPERTY(psi0)
            .DEF_PROPERTY(psi_lo)
            .DEF_PROPERTY(psi_hi)
            .DEF_PROPERTY(psi1)
            .DEF_PROPERTY(psi2)
            .DEF_PROPERTY(AdaptiveBeta)
            .DEF_PROPERTY(BetaLower)
            .DEF_PROPERTY(theta)
            .DEF_PROPERTY(qeps)
            .DEF_PROPERTY(qrule)
            .DEF_PROPERTY(qrestart);
    }

    {
        typedef cg_stats_wrapper cl;
        nb::class_<cl>(m, "cg_stats")
            .DEF_RO_PROPERTY(f)
            .DEF_RO_PROPERTY(gnorm)
            .DEF_RO_PROPERTY(iter)
            .DEF_RO_PROPERTY(IterSub)
            .DEF_RO_PROPERTY(NumSub)
            .DEF_RO_PROPERTY(nfunc)
            .DEF_RO_PROPERTY(ngrad);
    }

    {
        typedef cg_iter_stats_wrapper cl;
        nb::class_<cl>(m, "cg_iter_stats")
            .DEF_RO_PROPERTY(iter)
            .DEF_RO_PROPERTY(alpha)
            .DEF_RO_ARRAY_PROPERTY(x)
            .DEF_RO_PROPERTY(f)
            .DEF_RO_ARRAY_PROPERTY(g)
            .DEF_RO_ARRAY_PROPERTY(d);
    }

    m.def("cg_default", &cg_default_wrapper);
    m.def(
        "cg_descent",
        &cg_descent_wrapper,
        nb::arg("x").none(false),
        nb::arg("grad_tol").none(false),
        nb::arg("param").none(true),
        nb::arg("value").none(false),
        nb::arg("grad").none(false),
        nb::arg("valgrad").none(true),
        nb::arg("callback").none(true),
        nb::arg("work").none(true)
    );
}
