#include "cg_user.h"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <optional>
#include <functional>
#include <iostream>

#define DEF_RO_PROPERTY(NAME) \
    def_property_readonly(#NAME, &cl::get_##NAME)

#define DEF_PROPERTY(NAME) \
    def_property(#NAME, &cl::get_##NAME, &cl::set_##NAME)

#define CLASS_PROPERTY(NAME, TYPE) \
    TYPE get_##NAME() const noexcept { return obj->NAME; }; \
    void set_##NAME(TYPE v) { obj->NAME = v; };

#define CLASS_RO_PROPERTY(NAME, TYPE) \
    TYPE get_##NAME() const noexcept { return obj->NAME; };

#define WRAP_RAW_POINTER(NAME, RAWNAME, SIZE) \
    auto NAME = py::array(SIZE, RAWNAME, py::capsule(RAWNAME, [](void *p) {})); \
    assert(!NAME.owndata())

namespace py = pybind11;

// {{{ cg_parameter wrapper

class cg_parameter_wrapper
{
public:
    cg_parameter_wrapper()
    {
        obj = new cg_parameter;
        cg_default(obj);
    }
    ~cg_parameter_wrapper() { delete obj; };
    cg_parameter *data() { return obj; };

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

    CLASS_PROPERTY(rho, double)

    CLASS_RO_PROPERTY(delta, double)
    CLASS_RO_PROPERTY(sigma, double)
    CLASS_RO_PROPERTY(gamma, double)
    CLASS_RO_PROPERTY(psi0, double)
    CLASS_RO_PROPERTY(psi_lo, double)
    CLASS_RO_PROPERTY(psi_hi, double)
    CLASS_RO_PROPERTY(psi1, double)
    CLASS_RO_PROPERTY(psi2, double)
    CLASS_RO_PROPERTY(AdaptiveBeta, int)
    CLASS_RO_PROPERTY(BetaLower, double)
    CLASS_RO_PROPERTY(theta, double)
    CLASS_RO_PROPERTY(qeps, double)
    CLASS_RO_PROPERTY(qrule, double)
    CLASS_RO_PROPERTY(qrestart, int)

    cg_parameter *obj;
};

// }}}

// {{{ cg_stats wrapper

class cg_stats_wrapper
{
public:
    cg_stats_wrapper(): obj(new cg_stats) {};
    ~cg_stats_wrapper() { delete obj; };
    cg_stats *data() { return obj; };

    CLASS_RO_PROPERTY(f, double)
    CLASS_RO_PROPERTY(gnorm, double)
    CLASS_RO_PROPERTY(iter, INT)
    CLASS_RO_PROPERTY(IterSub, INT)
    CLASS_RO_PROPERTY(NumSub, INT)
    CLASS_RO_PROPERTY(nfunc, INT)
    CLASS_RO_PROPERTY(ngrad, INT)

    cg_stats *obj;
};

// }}}

// {{{ cg_descent wrapper

namespace cg {

typedef py::array_t<double, py::array::c_style | py::array::forcecast> array_t;

typedef std::function<double(array_t)> fn_t;
typedef std::function<void(array_t,array_t)> grad_t;
typedef std::function<double(array_t,array_t)> fngrad_t;

fn_t *fn;
grad_t *grad;
fngrad_t *fngrad;

double fn_trampoline(double *_x, INT n)
{
    WRAP_RAW_POINTER(x, _x, n);
    return (*fn)(x);
}

void grad_trampoline(double *_g, double *_x, INT n)
{
    WRAP_RAW_POINTER(g, _g, n);
    WRAP_RAW_POINTER(x, _x, n);
    (*grad)(g, x);
}

double fngrad_trampoline(double *_g, double *_x, INT n)
{
    WRAP_RAW_POINTER(g, _g, n);
    WRAP_RAW_POINTER(x, _x, n);
    return (*fngrad)(g, x);
}

};


py::tuple cg_descent_wrapper(
        cg::array_t x,
        double grad_tol,
        std::optional<cg_parameter_wrapper*> param,
        cg::fn_t &fn,
        cg::grad_t &grad,
        std::optional<cg::fngrad_t> fngrad,
        std::optional<cg::array_t> work
        )
{
    int status = 0;
    cg_stats_wrapper *stats = new cg_stats_wrapper;
    cg_parameter *p = param.has_value() ? param.value()->data() : nullptr;
    double *workptr = (work.has_value() ?
            static_cast<double*>(work.value().request().ptr) :
            nullptr);

    int n = x.shape(0);
    double *ptr = new double[n];
    auto xptr = x.unchecked();
    for (int i = 0; i < n; ++i) {
        ptr[i] = xptr(i);
    }

    // FIXME: figure out a nicer way to pass std::functions
    cg::fn = &fn;
    cg::grad = &grad;
    cg::fngrad = fngrad.has_value() ? &fngrad.value() : nullptr;
    auto *fngrad_trampoline = fngrad.has_value() ? cg::fngrad_trampoline : nullptr;

    status = cg_descent(
            ptr,
            x.shape(0),
            stats->data(),
            p,
            grad_tol,
            cg::fn_trampoline,
            cg::grad_trampoline,
            fngrad_trampoline,
            workptr);

    return py::make_tuple(
            py::array(n, ptr),
            py::cast(
                stats,
                py::return_value_policy::take_ownership),
            status
            );
}

// }}}

// {{{ cg_default wrapper

void cg_default_wrapper(py::object param)
{
    cg_default(param.cast<cg_parameter_wrapper*>()->data());
}

// }}}

PYBIND11_MODULE(_cg_descent, m)
{
    {
        typedef cg_parameter_wrapper cl;
        py::class_<cl>(m, "cg_parameter")
            .def(py::init())
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
            // NOTE: these are not recommended to play with, but are used
            // in some of the low level examples, so they're readwrite
            .DEF_PROPERTY(rho)
            // NOTE: these are not recommended to be play with, so keep them
            // read-only for now
            .DEF_RO_PROPERTY(delta)
            .DEF_RO_PROPERTY(sigma)
            .DEF_RO_PROPERTY(gamma)
            .DEF_RO_PROPERTY(psi0)
            .DEF_RO_PROPERTY(psi_lo)
            .DEF_RO_PROPERTY(psi_hi)
            .DEF_RO_PROPERTY(psi1)
            .DEF_RO_PROPERTY(psi2)
            .DEF_RO_PROPERTY(AdaptiveBeta)
            .DEF_RO_PROPERTY(BetaLower)
            .DEF_RO_PROPERTY(theta)
            .DEF_RO_PROPERTY(qeps)
            .DEF_RO_PROPERTY(qrule)
            .DEF_RO_PROPERTY(qrestart)
        ;
    }

    {
        typedef cg_stats_wrapper cl;
        py::class_<cl>(m, "cg_stats")
            .DEF_RO_PROPERTY(f)
            .DEF_RO_PROPERTY(gnorm)
            .DEF_RO_PROPERTY(iter)
            .DEF_RO_PROPERTY(IterSub)
            .DEF_RO_PROPERTY(NumSub)
            .DEF_RO_PROPERTY(nfunc)
            .DEF_RO_PROPERTY(ngrad)
        ;
    }

    m.def("cg_default", &cg_default_wrapper);
    m.def("cg_descent", &cg_descent_wrapper,
            py::arg("x").none(false),
            py::arg("grad_tol").none(false),
            py::arg("param").none(true),
            py::arg("fn").none(false),
            py::arg("grad").none(false),
            py::arg("fngrad").none(true),
            py::arg("work").none(true));
}
