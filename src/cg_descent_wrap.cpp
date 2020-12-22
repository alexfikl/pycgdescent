#include "cg_user.h"
#include "cg_descent.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

#define DEF_RO_MEMBER(NAME) \
    def_readonly(#NAME, &cl::NAME)

PYBIND11_MODULE(_cg_descent, m)
{
    {
        typedef cg_parameter cl;
        py::class_<cl>(m, "cg_parameter")
            .def(py::init([]() {
                cl* param = new cl();
                cg_default(param);
                return std::unique_ptr<cl>(param);
            }))
            .DEF_RO_MEMBER(PrintFinal)
            .DEF_RO_MEMBER(PrintLevel)
            .DEF_RO_MEMBER(PrintParms)
            .DEF_RO_MEMBER(LBFGS)
            .DEF_RO_MEMBER(memory)
            .DEF_RO_MEMBER(SubCheck)
            .DEF_RO_MEMBER(SubSkip)
            .DEF_RO_MEMBER(eta0)
            .DEF_RO_MEMBER(eta1)
            .DEF_RO_MEMBER(AWolfe)
            .DEF_RO_MEMBER(AWolfeFac)
            .DEF_RO_MEMBER(Qdecay)
            .DEF_RO_MEMBER(nslow)
            .DEF_RO_MEMBER(StopRule)
            .DEF_RO_MEMBER(StopFac)
            .DEF_RO_MEMBER(PertRule)
            .DEF_RO_MEMBER(eps)
            .DEF_RO_MEMBER(egrow)
            .DEF_RO_MEMBER(QuadStep)
            .DEF_RO_MEMBER(QuadCutOff)
            .DEF_RO_MEMBER(QuadSafe)
            .DEF_RO_MEMBER(UseCubic)
            .DEF_RO_MEMBER(CubicCutOff)
            .DEF_RO_MEMBER(SmallCost)
            .DEF_RO_MEMBER(debug)
            .DEF_RO_MEMBER(debugtol)
            .DEF_RO_MEMBER(step)
            .DEF_RO_MEMBER(maxit)
            .DEF_RO_MEMBER(ntries)
            .DEF_RO_MEMBER(ExpandSafe)
            .DEF_RO_MEMBER(SecantAmp)
            .DEF_RO_MEMBER(neps)
            .DEF_RO_MEMBER(nshrink)
            .DEF_RO_MEMBER(nline)
            .DEF_RO_MEMBER(restart_fac)
            .DEF_RO_MEMBER(feps)
            .DEF_RO_MEMBER(nan_rho)
            .DEF_RO_MEMBER(nan_decay)

            .DEF_RO_MEMBER(delta)
            .DEF_RO_MEMBER(sigma)
            .DEF_RO_MEMBER(gamma)
            .DEF_RO_MEMBER(rho)
            .DEF_RO_MEMBER(psi0)
            .DEF_RO_MEMBER(psi_lo)
            .DEF_RO_MEMBER(psi_hi)
            .DEF_RO_MEMBER(psi1)
            .DEF_RO_MEMBER(psi2)
            .DEF_RO_MEMBER(AdaptiveBeta)
            .DEF_RO_MEMBER(BetaLower)
            .DEF_RO_MEMBER(theta)
            .DEF_RO_MEMBER(qeps)
            .DEF_RO_MEMBER(qrule)
            .DEF_RO_MEMBER(qrestart)
        ;
    }

    {
        typedef cg_stats cl;
        py::class_<cl>(m, "cg_stats")
            .DEF_RO_MEMBER(f)
            .DEF_RO_MEMBER(gnorm)
            .DEF_RO_MEMBER(iter)
            .DEF_RO_MEMBER(IterSub)
            .DEF_RO_MEMBER(NumSub)
            .DEF_RO_MEMBER(nfunc)
            .DEF_RO_MEMBER(ngrad)
        ;
    }

    m.def("cg_default", cg_default);
    m.def("cg_descent", cg_descent);
}
