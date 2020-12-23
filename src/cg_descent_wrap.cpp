#include "cg_user.h"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include <functional>
#include <iostream>

#define DEF_RO_PROPERTY(NAME) \
    def_property_readonly(#NAME, &cl::get_##NAME)

#define DEF_PROPERTY(NAME) \
    def_property(#NAME, &cl::get_##NAME, &cl::set_##NAME)

// FIXME: this probably does a copy
#define TO_NUMPY_ARRAY(ptr, n, type) \
    py::array_t<type>( \
        { n }, { sizeof(type) * n }, ptr)

namespace py = pybind11;

// {{{ cg_parameter wrapper

#define CLASS_PROPERTY(NAME, TYPE) \
    TYPE get_##NAME() const { return obj->NAME; }; \
    void set_##NAME(TYPE v) { obj->NAME = v; };

#define CLASS_RO_PROPERTY(NAME, TYPE) \
    TYPE get_##NAME() const { return obj->NAME; };

class cg_parameter_wrapper
{
public:
    cg_parameter_wrapper(cg_parameter *param): obj(std::move(param)) {};
    cg_parameter_wrapper()
    {
        obj = new cg_parameter;
        cg_default(obj);
    }
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

    CLASS_RO_PROPERTY(delta, double)
    CLASS_RO_PROPERTY(sigma, double)
    CLASS_RO_PROPERTY(gamma, double)
    CLASS_RO_PROPERTY(rho, double)
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
    cg_stats_wrapper(cg_stats *stats): obj(std::move(stats)) {}
    cg_stats_wrapper(): obj(new cg_stats) {}
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

typedef std::function<double(py::array_t<double>)> cg_fn_t;
typedef std::function<void(py::array_t<double>,py::array_t<double>)> cg_grad_t;
typedef std::function<double(py::array_t<double>,py::array_t<double>)> cg_valgrad_t;

cg_fn_t *cg_fn;
cg_grad_t *cg_grad;
cg_valgrad_t *cg_valgrad;


double fn_trampoline(double *_x, INT n)
{
    auto x = TO_NUMPY_ARRAY(_x, n, double);
    return (*cg_fn)(x);
}

void grad_trampoline(double *_g, double *_x, INT n)
{
    auto g = TO_NUMPY_ARRAY(_g, n, double);
    auto x = TO_NUMPY_ARRAY(_x, n, double);
    (*cg_grad)(g, x);
}

double valgrad_trampoline(double *_g, double *_x, INT n)
{
    auto g = TO_NUMPY_ARRAY(_g, n, double);
    auto x = TO_NUMPY_ARRAY(_x, n, double);
    return (*cg_valgrad)(g, x);
}

py::tuple cg_descent_wrapper(
        py::array_t<double> x,
        py::object UParm,
        double grad_tol,
        cg_fn_t &value,
        cg_grad_t &grad,
        cg_valgrad_t &valgrad
        )
{
    int ret = 0;
    cg_stats_wrapper *stats = new cg_stats_wrapper;
    cg_parameter *param = UParm.cast<cg_parameter_wrapper*>()->data();

    int n = x.shape(0);
    double *ptr = new double[n];
    auto xptr = x.unchecked();
    for (int i = 0; i < n; ++i) {
        ptr[i] = xptr[i];
    }

    // FIXME: figure out a nicer way to pass std::functions
    cg_fn = &value;
    cg_grad = &grad;
    cg_valgrad = &valgrad;

    ret = cg_descent(
            ptr,
            x.shape(0),
            stats->data(),
            param,
            grad_tol,
            fn_trampoline,
            grad_trampoline,
            valgrad_trampoline,
            NULL);

    return py::make_tuple(
            TO_NUMPY_ARRAY(ptr, x.shape(0), double),
            py::cast(
                stats,
                py::return_value_policy::take_ownership)
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
            // NOTE: these are not recommended to be play with, so keep them
            // read-only for now
            .DEF_RO_PROPERTY(delta)
            .DEF_RO_PROPERTY(sigma)
            .DEF_RO_PROPERTY(gamma)
            .DEF_RO_PROPERTY(rho)
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
    m.def("cg_descent", &cg_descent_wrapper);
}
