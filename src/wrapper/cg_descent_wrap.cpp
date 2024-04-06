// SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
//
// SPDX-License-Identifier: MIT

#include "cg_descent.h"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <optional>
#include <functional>
#include <iostream>

namespace py = pybind11;

// {{{ macros

#define ADD_ATTR(NAME) cls.attr(#NAME) = NAME

#define WRAP_RAW_POINTER(NAME, RAWNAME, SIZE) \
    auto NAME = py::array(SIZE, RAWNAME, py::capsule(RAWNAME, [](void *p) { (void)p; })); \
    assert(!NAME.owndata())

#define DEF_RO_PROPERTY(NAME) \
    def_property_readonly(#NAME, &cl::get_##NAME)

#define DEF_PROPERTY(NAME) \
    def_property(#NAME, &cl::get_##NAME, &cl::set_##NAME)

#define CLASS_PROPERTY(NAME, TYPE) \
    TYPE get_##NAME() const noexcept { return obj->NAME; }; \
    void set_##NAME(TYPE v) { obj->NAME = v; };

#define CLASS_RO_PROPERTY(NAME, TYPE) \
    TYPE get_##NAME() const noexcept { return obj->NAME; };

#define CLASS_RO_ARRAY_PROPERTY(NAME, TYPE) \
    py::array get_##NAME() const { \
        WRAP_RAW_POINTER(NAME, obj->NAME, obj->n); \
        return NAME; \
    };

// }}}

// {{{ cg_parameter wrapper

class cg_parameter_wrapper
{
public:
    cg_parameter_wrapper(CGparm *p): obj(p) { }
    cg_parameter_wrapper(): obj(new CGparm) { cg_default(obj); }
    cg_parameter_wrapper(const cg_parameter_wrapper &w)
    {
        obj = new CGparm;
        memcpy(obj, w.obj, sizeof(CGparm));
    }

    ~cg_parameter_wrapper() { delete obj; };
    CGparm *data() { return obj; };

    CLASS_PROPERTY(grad_tol, CGFLOAT)

    CLASS_PROPERTY(PrintStatus, int)
    CLASS_PROPERTY(PrintStat, int)
    CLASS_PROPERTY(PrintParm, int)
    CLASS_PROPERTY(PrintLevel, int)

    CLASS_PROPERTY(QuadCost, int)
    CLASS_PROPERTY(FastLA, int)

    CLASS_PROPERTY(deriv_mode, int)
    CLASS_PROPERTY(dense_Hessian, CGFLOAT)
    CLASS_PROPERTY(CG_window, CGFLOAT)
    CLASS_PROPERTY(Newton_window, CGFLOAT)
    CLASS_PROPERTY(min_error_iter, int)
    CLASS_PROPERTY(Newton_cutoff, int)
    CLASS_PROPERTY(RhoDecay, CGFLOAT)
    CLASS_PROPERTY(err_decay, CGFLOAT)
    CLASS_PROPERTY(cg_ok, int)
    CLASS_PROPERTY(Hessian_err_decay, CGFLOAT)
    CLASS_PROPERTY(Hessian_CG_solver, int)
    CLASS_PROPERTY(ss_diag_pert, CGFLOAT)
    CLASS_PROPERTY(cg_diag_pert, CGFLOAT)
    CLASS_PROPERTY(HessianSparsityFixed, int)
    CLASS_PROPERTY(TrustIterLimit, int)
    CLASS_PROPERTY(big_grad, CGFLOAT)
    CLASS_PROPERTY(QPshift, CGFLOAT)
    CLASS_PROPERTY(QPgReset_factor, CGFLOAT)
    CLASS_PROPERTY(StopFac, CGFLOAT)

    CLASS_PROPERTY(debug, int)
    CLASS_PROPERTY(debugtol, CGFLOAT)
    CLASS_PROPERTY(CheckMatrix, int)

    CLASS_PROPERTY(step, CGFLOAT)
    CLASS_PROPERTY(LBFGS, int)
    CLASS_PROPERTY(LBFGSmemory, int)
    CLASS_PROPERTY(maxit, CGINT)
    CLASS_PROPERTY(restart_fac, CGFLOAT)
    CLASS_PROPERTY(Qdecay, CGFLOAT)
    CLASS_PROPERTY(nslow, int)
    CLASS_PROPERTY(QuadStep, int)
    CLASS_PROPERTY(QuadCutOff, CGFLOAT)
    CLASS_PROPERTY(QuadSafe, CGFLOAT)
    CLASS_PROPERTY(psi_lo, CGFLOAT)
    CLASS_PROPERTY(psi_hi, CGFLOAT)
    CLASS_PROPERTY(psi1, CGFLOAT)
    CLASS_PROPERTY(qeps, CGFLOAT)
    CLASS_PROPERTY(qrule, CGFLOAT)
    CLASS_PROPERTY(qrestart, int)
    CLASS_PROPERTY(UseCubic, int)
    CLASS_PROPERTY(CubicCutOff, CGFLOAT)
    CLASS_PROPERTY(SmallCost, CGFLOAT)
    CLASS_PROPERTY(ExpandSafe, CGFLOAT)
    CLASS_PROPERTY(SecantAmp, CGFLOAT)
    CLASS_PROPERTY(approxstep, int)
    CLASS_PROPERTY(ApproxSwitchFactor, CGFLOAT)
    CLASS_PROPERTY(CostConverge, CGFLOAT)
    CLASS_PROPERTY(cgdelta, CGFLOAT)
    CLASS_PROPERTY(cgsigma, CGFLOAT)
    CLASS_PROPERTY(maxsteps, int)
    CLASS_PROPERTY(stepdecay, CGFLOAT)
    CLASS_PROPERTY(cg_infdecay, CGFLOAT)
    CLASS_PROPERTY(cg_infdecay_rate, CGFLOAT)
    CLASS_PROPERTY(cg_ninf_tries, int)
    CLASS_PROPERTY(rho, CGFLOAT)
    CLASS_PROPERTY(RhoGrow, CGFLOAT)
    CLASS_PROPERTY(BigDfactor, CGFLOAT)
    CLASS_PROPERTY(PertRule, int)
    CLASS_PROPERTY(pert_eps, CGFLOAT)
    CLASS_PROPERTY(eps_grow, CGFLOAT)
    CLASS_PROPERTY(neps, int)
    CLASS_PROPERTY(psi0, CGFLOAT)
    CLASS_PROPERTY(psi2, CGFLOAT)
    CLASS_PROPERTY(BetaLower, CGFLOAT)
    CLASS_PROPERTY(theta, CGFLOAT)
    CLASS_PROPERTY(AdaptiveTheta, int)

    CLASS_PROPERTY(SubCheck, int)
    CLASS_PROPERTY(SubSkip, int)
    CLASS_PROPERTY(eta0, CGFLOAT)
    CLASS_PROPERTY(eta1, CGFLOAT)
    CLASS_PROPERTY(eta2, CGFLOAT)

private:

    CGparm *obj;
};

// }}}

// {{{ cg_stats wrapper

class cg_stats_wrapper
{
public:
    cg_stats_wrapper(): obj(new CGstat) {};
    ~cg_stats_wrapper() { delete obj; };
    CGstat *data() { return obj; };

    CLASS_RO_PROPERTY(status, int)
    CLASS_RO_PROPERTY(f, CGFLOAT)
    CLASS_RO_PROPERTY(err, CGFLOAT)
    CLASS_RO_PROPERTY(grad_tol, CGFLOAT)
    CLASS_RO_PROPERTY(tol, CGFLOAT)
    CLASS_RO_PROPERTY(maxit, CGINT)
    CLASS_RO_PROPERTY(cg_ninf_tries, CGINT)
    CLASS_RO_PROPERTY(oldf, CGFLOAT)
    CLASS_RO_PROPERTY(newf, CGFLOAT)
    CLASS_RO_PROPERTY(maxsteps, int)
    CLASS_RO_PROPERTY(NegDiag, int)
    CLASS_RO_PROPERTY(iter, CGINT)
    CLASS_RO_PROPERTY(nfunc, CGINT)
    CLASS_RO_PROPERTY(ngrad, CGINT)
    CLASS_RO_PROPERTY(nexpand, CGINT)
    CLASS_RO_PROPERTY(nforward, CGINT)
    CLASS_RO_PROPERTY(nback, CGINT)
    CLASS_RO_PROPERTY(nCG, CGINT)
    CLASS_RO_PROPERTY(nSS, CGINT)
    CLASS_RO_PROPERTY(PRP, CGINT)
    CLASS_RO_PROPERTY(IterSub, int)
    CLASS_RO_PROPERTY(NumSub, int)

private:

    CGstat *obj;
};

// }}}

// {{{ cg_iter_stats_wrapper

class cg_iter_stats_wrapper
{
public:
    cg_iter_stats_wrapper(CGiter *stats): obj(stats) {};
    cg_iter_stats_wrapper(): obj(nullptr) {};
    ~cg_iter_stats_wrapper() { };

    CLASS_RO_PROPERTY(iter, CGINT)
    CLASS_RO_PROPERTY(n, int)
    CLASS_RO_PROPERTY(alpha, CGFLOAT)
    CLASS_RO_ARRAY_PROPERTY(x, CGFLOAT)
    CLASS_RO_PROPERTY(f, CGFLOAT)
    CLASS_RO_ARRAY_PROPERTY(g, CGFLOAT)
    CLASS_RO_ARRAY_PROPERTY(d, CGFLOAT)

    CGiter *obj;
};

// }}}

// {{{ status wrapper

class cg_status {};

// }}}

// {{{ cg_descent wrapper

namespace cg {

typedef py::array_t<double, py::array::c_style | py::array::forcecast> array;

typedef std::function<double(array)> value_fn;
typedef std::function<void(array,array)> grad_fn;
typedef std::function<double(array,array)> valgrad_fn;
typedef std::function<int(cg_iter_stats_wrapper&)> callback_fn;

}

py::tuple cg_descent_wrapper(
        cg::array x,
        double grad_tol,
        std::optional<cg_parameter_wrapper*> param,
        cg::value_fn &value,
        cg::grad_fn &grad,
        std::optional<cg::valgrad_fn> valgrad,
        std::optional<cg::callback_fn> callback,
        std::optional<cg::array> work
        )
{
    auto wrapper = [=]() -> py::tuple {
        int status = 0;

        static cg::value_fn _value;
        static cg::grad_fn _grad;
        static cg::valgrad_fn _valgrad;
        static cg::callback_fn _callback;

        int n = x.shape(0);
        CGFLOAT *ptr = new CGFLOAT[n];
        auto xptr = x.unchecked();
        for (int i = 0; i < n; ++i) {
            ptr[i] = xptr(i);
        }

        CGdata *cgdata = cg_setup();
        cgdata->n = n;
        cgdata->x = ptr;

        if (work.has_value()) {
            cgdata->Work = static_cast<CGFLOAT*>(work.value().request().ptr);
        }

        if (param.has_value()) {
            memcpy(cgdata->Parm, param.value()->data(), sizeof(CGparm));
        }
        cgdata->Parm->grad_tol = grad_tol;

        _value = std::move(value);
        cgdata->value = [](CGFLOAT *f, CGFLOAT *x, CGINT n) {
            WRAP_RAW_POINTER(_x, x, n);
            *f = _value(_x);
        };

        _grad = std::move(grad);
        cgdata->grad = [](CGFLOAT *g, CGFLOAT *x, CGINT n) {
            WRAP_RAW_POINTER(_g, g, n);
            WRAP_RAW_POINTER(_x, x, n);
            _grad(_g, _x);
        };

        if (valgrad.has_value()) {
            _valgrad = std::move(valgrad.value());
            cgdata->valgrad = [](CGFLOAT *f, CGFLOAT *g, CGFLOAT *x, CGINT n) {
                WRAP_RAW_POINTER(_g, g, n);
                WRAP_RAW_POINTER(_x, x, n);
                *f = _valgrad(_g, _x);
            };
        }

        if (callback.has_value()) {
            _callback = std::move(callback.value());
            cgdata->callback = [](CGiter *stats) -> int {
                cg_iter_stats_wrapper w(stats);
                return _callback(w);
            };
        }

        cg_stats_wrapper *cgstats = new cg_stats_wrapper;

        status = cg_descent(cgdata);
        memcpy(cgstats->data(), cgdata->Stat, sizeof(CGstat));
        cg_terminate(&cgdata);

        return py::make_tuple(
                py::array(n, ptr),
                py::cast(
                    cgstats,
                    py::return_value_policy::take_ownership),
                status
                );
    };

    return wrapper();
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
            .DEF_PROPERTY(grad_tol)
            .DEF_PROPERTY(PrintStatus)
            .DEF_PROPERTY(PrintStat)
            .DEF_PROPERTY(PrintParm)
            .DEF_PROPERTY(PrintLevel)
            .DEF_PROPERTY(QuadCost)
            .DEF_PROPERTY(FastLA)
            .DEF_PROPERTY(deriv_mode)
            .DEF_PROPERTY(dense_Hessian)
            .DEF_PROPERTY(CG_window)
            .DEF_PROPERTY(Newton_window)
            .DEF_PROPERTY(min_error_iter)
            .DEF_PROPERTY(Newton_cutoff)
            .DEF_PROPERTY(RhoDecay)
            .DEF_PROPERTY(err_decay)
            .DEF_PROPERTY(cg_ok)
            .DEF_PROPERTY(Hessian_err_decay)
            .DEF_PROPERTY(Hessian_CG_solver)
            .DEF_PROPERTY(ss_diag_pert)
            .DEF_PROPERTY(cg_diag_pert)
            .DEF_PROPERTY(HessianSparsityFixed)
            .DEF_PROPERTY(TrustIterLimit)
            .DEF_PROPERTY(big_grad)
            .DEF_PROPERTY(QPshift)
            .DEF_PROPERTY(QPgReset_factor)
            .DEF_PROPERTY(StopFac)
            .DEF_PROPERTY(debug)
            .DEF_PROPERTY(debugtol)
            .DEF_PROPERTY(CheckMatrix)
            .DEF_PROPERTY(step)
            .DEF_PROPERTY(LBFGS)
            .DEF_PROPERTY(LBFGSmemory)
            .DEF_PROPERTY(maxit)
            .DEF_PROPERTY(restart_fac)
            .DEF_PROPERTY(Qdecay)
            .DEF_PROPERTY(nslow)
            .DEF_PROPERTY(QuadStep)
            .DEF_PROPERTY(QuadCutOff)
            .DEF_PROPERTY(QuadSafe)
            .DEF_PROPERTY(psi_lo)
            .DEF_PROPERTY(psi_hi)
            .DEF_PROPERTY(psi1)
            .DEF_PROPERTY(qeps)
            .DEF_PROPERTY(qrule)
            .DEF_PROPERTY(qrestart)
            .DEF_PROPERTY(UseCubic)
            .DEF_PROPERTY(CubicCutOff)
            .DEF_PROPERTY(SmallCost)
            .DEF_PROPERTY(ExpandSafe)
            .DEF_PROPERTY(SecantAmp)
            .DEF_PROPERTY(approxstep)
            .DEF_PROPERTY(ApproxSwitchFactor)
            .DEF_PROPERTY(CostConverge)
            .DEF_PROPERTY(cgdelta)
            .DEF_PROPERTY(cgsigma)
            .DEF_PROPERTY(maxsteps)
            .DEF_PROPERTY(stepdecay)
            .DEF_PROPERTY(cg_infdecay)
            .DEF_PROPERTY(cg_infdecay_rate)
            .DEF_PROPERTY(cg_ninf_tries)
            .DEF_PROPERTY(rho)
            .DEF_PROPERTY(RhoGrow)
            .DEF_PROPERTY(BigDfactor)
            .DEF_PROPERTY(PertRule)
            .DEF_PROPERTY(pert_eps)
            .DEF_PROPERTY(eps_grow)
            .DEF_PROPERTY(neps)
            .DEF_PROPERTY(psi0)
            .DEF_PROPERTY(psi2)
            .DEF_PROPERTY(BetaLower)
            .DEF_PROPERTY(theta)
            .DEF_PROPERTY(AdaptiveTheta)
            .DEF_PROPERTY(SubCheck)
            .DEF_PROPERTY(SubSkip)
            .DEF_PROPERTY(eta0)
            .DEF_PROPERTY(eta1)
            .DEF_PROPERTY(eta2)
        ;
    }

    {
        typedef cg_stats_wrapper cl;
        py::class_<cl>(m, "cg_stats")
            .DEF_RO_PROPERTY(status)
            .DEF_RO_PROPERTY(f)
            .DEF_RO_PROPERTY(err)
            .DEF_RO_PROPERTY(grad_tol)
            .DEF_RO_PROPERTY(tol)
            .DEF_RO_PROPERTY(maxit)
            .DEF_RO_PROPERTY(cg_ninf_tries)
            .DEF_RO_PROPERTY(oldf)
            .DEF_RO_PROPERTY(newf)
            .DEF_RO_PROPERTY(maxsteps)
            .DEF_RO_PROPERTY(NegDiag)
            .DEF_RO_PROPERTY(iter)
            .DEF_RO_PROPERTY(nfunc)
            .DEF_RO_PROPERTY(ngrad)
            .DEF_RO_PROPERTY(nexpand)
            .DEF_RO_PROPERTY(nforward)
            .DEF_RO_PROPERTY(nback)
            .DEF_RO_PROPERTY(nCG)
            .DEF_RO_PROPERTY(nSS)
            .DEF_RO_PROPERTY(PRP)
            .DEF_RO_PROPERTY(IterSub)
            .DEF_RO_PROPERTY(NumSub)
        ;
    }

    {
        typedef cg_iter_stats_wrapper cl;
        py::class_<cl>(m, "cg_iter_stats")
            .DEF_RO_PROPERTY(iter)
            .DEF_RO_PROPERTY(alpha)
            .DEF_RO_PROPERTY(x)
            .DEF_RO_PROPERTY(f)
            .DEF_RO_PROPERTY(g)
            .DEF_RO_PROPERTY(d)
        ;
    }

    {
        py::class_<cg_status> cls(m, "status_code");

        ADD_ATTR(CG_ERROR_TOLERANCE_SATISFIED);
        ADD_ATTR(CG_ITERATIONS_EXCEED_MAXITS);
        ADD_ATTR(CG_SLOPE_ALWAYS_NEGATIVE);
        ADD_ATTR(CG_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS);
        ADD_ATTR(CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION);
        ADD_ATTR(CG_WOLFE_CONDITIONS_NOT_SATISFIED);
        ADD_ATTR(CG_DEBUGGER_IS_ON_AND_FUNCTION_VALUE_INCREASES);
        ADD_ATTR(CG_NO_COST_OR_GRADIENT_IMPROVEMENT);
        ADD_ATTR(CG_OUT_OF_MEMORY);
        ADD_ATTR(CG_QUADRATIC_OBJECTIVE_NO_LOWER_BOUND);
        ADD_ATTR(CG_STARTING_FUNCTION_VALUE_INFINITE_OR_NAN);
        ADD_ATTR(CG_EXCESSIVE_UPDATING_OF_PERT_EPS);
        ADD_ATTR(CG_FUNCTION_NAN_OR_INF);
        ADD_ATTR(CG_QP_LINEAR_TERM_GIVEN_BUT_HPROD_MISSING);
        ADD_ATTR(CG_N_IS_EMPTY);
        ADD_ATTR(CG_ERROR_IN_INPUT_MATRIX);
        ADD_ATTR(CG_MISSING_HESSIAN_FOR_QUADCOST);
        ADD_ATTR(CG_INVALID_DERIV_MODE_PARAMETER);
        ADD_ATTR(CG_DERIV_MODE_USES_HESSIAN_BUT_NO_HESSIAN_PROVIDED);
        ADD_ATTR(CG_SYMMETRIC_SOLVER_FAILS);
        ADD_ATTR(CG_HESSIAN_NOT_COMPUTED);
        ADD_ATTR(CG_HPROD_PLUS_HESSIAN);
        ADD_ATTR(CG_VALUE_OR_GRAD_MISSING);
        ADD_ATTR(CG_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE);
        ADD_ATTR(CG_TRIPLES_FORMAT_ERROR);
        ADD_ATTR(CG_MULTI_SOLVERS);
        ADD_ATTR(CG_USER_CALLBACK);
    }

    m.def("cg_default", &cg_default_wrapper);
    m.def("cg_descent", &cg_descent_wrapper,
            py::arg("x").none(false),
            py::arg("grad_tol").none(false),
            py::arg("param").none(true),
            py::arg("value").none(false),
            py::arg("grad").none(false),
            py::arg("valgrad").none(true),
            py::arg("callback").none(true),
            py::arg("work").none(true));
}
