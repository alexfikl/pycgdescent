/* ==========================================================================
   === cg_print_status ======================================================
   ==========================================================================
    Print the status at termination of the run
   ========================================================================== */
#include "cg_descent.h"
#ifdef MATLAB_MEX_FILE
#define LPAREN "("
#define RPAREN ")"
#define OFFSET 1
#else
#define LPAREN "["
#define RPAREN "]"
#define OFFSET 0
#endif

void cg_print_status
(
    CGdata *Data /* pointer to cgdata structure */
)
{
    CGstat *Stat = Data->Stat ;
    int status = Stat->status ;
    printf ("\nCG_DESCENT (Version %d.%d.%d, %s):\n\n",
        CG_MAIN_VERSION, CG_SUB_VERSION, CG_SUBSUB_VERSION, CG_DATE) ;

    /* CG error message strings */
     const char mess1 [] = "Possible causes of this error message:" ;
     const char mess2 [] = "   - your tolerance may be too strict: "
                           "grad_tol = " ;
     const char mess3 [] = "   - your gradient routine has an error" ;

    if ( status == CG_ERROR_TOLERANCE_SATISFIED )
    {
        printf ("Success: Error %e satisfies error "
                "tolerance %e.\n", Stat->err, Stat->tol) ;
    }
    else if ( status == CG_ITERATIONS_EXCEED_MAXITS )
    {
        printf ("The number of iterations exceed "
                 "specified limit of %ld.\n\n%s\n%s %e\n",
                (LONG) Stat->maxit, mess1, mess2, Stat->tol) ;
    }
    else if ( status == CG_SLOPE_ALWAYS_NEGATIVE )
    {
        printf ("The slope is always negative in line search.\n%s\n"
                "   - your cost function has an error\n%s\n", mess1, mess3) ;
    }
    else if ( status == CG_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS )
    {
        printf ("Unable to find an acceptable step in the\n"
                "line search before hitting the maxsteps limit %i\n"
                "in the parameter structure.\n%s\n%s %e\n",
                Stat->maxsteps, mess1, mess2, Stat->tol) ;
    }
    else if ( status == CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION )
    {
        printf ("The search direction was not a descent direction.\n") ;
    }
    else if ( status == CG_EXCESSIVE_UPDATING_OF_PERT_EPS )
    {
        printf ("The line search fails due to excessive "
                "updating of the parameter pert_eps.\n%s\n%s %e\n%s\n",
                mess1, mess2, Stat->tol, mess3) ;
    }
    else if ( status == CG_WOLFE_CONDITIONS_NOT_SATISFIED )
    {
        /* line search fails */
        printf ("The line search fails.\n%s\n%s %e\n%s\n"
                "   - Parm->pert_eps may be too small.\n",
                mess1, mess2, Stat->tol, mess3) ;
    }
    else if ( status == CG_DEBUGGER_IS_ON_AND_FUNCTION_VALUE_INCREASES )
    {
        printf ("The debugger in CG_DESCENT was turned on and the function "
                "value did not improve.\n"
                "new value: %25.16e old value: %25.16e\n",
                 Stat->newf, Stat->oldf) ;
    }
    else if ( status == CG_NO_COST_OR_GRADIENT_IMPROVEMENT )
    {
        printf ("Parm->nslow iterations performed in CG_DESCENT without "
                "strict improvement in cost or gradient.\n") ;
    }
    else if ( status == CG_OUT_OF_MEMORY )
    {
        printf ("CG_DESCENT ran out of memory.\n") ;
    }
    else if ( status == CG_QUADRATIC_OBJECTIVE_NO_LOWER_BOUND )
    {
        printf ("The quadratic objective has no lower bound "
                "over the feasible region. If the\n"
                "Hessian is approximately positive semidefinite, "
                "you could try to regularize\n"
                "the Hessian by adding a small positive number "
                "to the diagonal. This adjustment\n"
                "is achieved by setting the parameter QPshift to a "
                "small positive number.\n"
                "If PASA was called, then adjust PASA's parameter "
                "QPshift.\n"
                "If CG_DESCENT was called, then adjust CG_DESCENT's "
                "parameter QPshift.\n") ;
    }
    else if ( status == CG_STARTING_FUNCTION_VALUE_INFINITE_OR_NAN )
    {
        printf ("The function value is nan at the starting point.\n") ;
    }
    else if ( status == CG_FUNCTION_NAN_OR_INF )
    {
        printf ("The line search could not locate a finite objective value\n"
                "after Parm->cg_ninf_tries attempts.\n"
                "-- currently Parm->cg_ninf_tries = %i\n", Stat->cg_ninf_tries);
    }
    else if ( status == CG_QP_LINEAR_TERM_GIVEN_BUT_HPROD_MISSING )
    {
        printf ("The linear cost term for a quadratic was given to CG_DESCENT\n"
                "but not the rule hprod for multiplying a vector by the\n"
                "Hessian of the objective.\n") ;
    }
    else if ( status == CG_N_IS_EMPTY )
    {
        printf ("In CG_DESCENT, the problem dimension was not provided in\n"
                "cgdata->n\n") ;
    }
    else if ( status == CG_ERROR_IN_INPUT_MATRIX )
    {
        printf ("In CG_DESCENT, the Hessian of a quadratic was provided using\n"
                "one of the three possible input formats: sparse matrix,\n"
                "triples, or dense matrix.  An error was detected in the\n"
                "input matrix which was described above in more detail.\n") ;
    }
    else if ( status == CG_MISSING_HESSIAN_FOR_QUADCOST )
    {
        printf ("The CG objective was specified as quadratic,\n"
                "however, the user did not provide the Hessian, either\n"
                "as a sparse matrix, as a routine cgdata->hessian\n"
                "for evaluating the objective Hessian, or as a routine\n"
                "cgdata->hprod for evaluating the product between the\n"
                "Hessian and a vector.\n") ;
    }
    else if ( status == CG_INVALID_DERIV_MODE_PARAMETER )
    {
        printf ("In CG_DESCENT, the parameter deriv_mode = %i is invalid.\n",
                Data->Parm->deriv_mode) ;
    }
    else if (status == CG_DERIV_MODE_USES_HESSIAN_BUT_NO_HESSIAN_PROVIDED)
    {
        printf ("In CG_DESCENT, the parameter deriv_mode = %i implies\n"
                "that the Hessian of the objective should be exploited,\n"
                "however, neither the function to evaluate the Hessian nor\n"
                "the Hessian of a quadratic were given in the input cgdata\n"
                "structure.\n",
                Data->Parm->deriv_mode) ;
    }
    else if (status == CG_SYMMETRIC_SOLVER_FAILS)
    {
        printf ("In CG_DESCENT, deriv_mode set to 2, which implies that a\n"
                "Newton Hessian-based implementation is to be used, however,\n"
                "the symmetric solver failed. Note that the symmetric solver\n"
                "is normally specified through OPTFLAGS in the file\n"
                "SuiteOPTconfig/Userconfig.mk\n") ;
    }
    else if (status == CG_HESSIAN_NOT_COMPUTED)
    {
        printf ("In CG_DESCENT, the user specified the Hessian-based\n"
                "implementation of CG_DESCENT, however, the user's Hessian\n"
                "function did not return the Hessian matrix when it was\n"
                "invoked.\n") ;
    }
    else if ( status == CG_HPROD_PLUS_HESSIAN )
    {
        printf ("In CG, user provides the routine cgdata->hprod to\n"
                    "evaluate the product between the Hessian of a quadratic\n"
                    "and a vector, and the routine cgdata->hessian for\n"
                    "evaluating the Hessian of a matrix. Only one routine\n"
                    "should be provided.\n") ;
    }
    else if ( status == CG_VALUE_OR_GRAD_MISSING )
    {
        printf ("In CG, the objective is not quadratic and the user failed\n"
                "to provide the routine to evaluate either the objective\n"
                "or its gradient.\n") ;
    }
    else if ( status == CG_TRIPLES_FORMAT_ERROR )
    {
        printf ("If the CG user returns the Hessian at a given point\n"
                "using triples format, then only the elements on the main\n"
                "diagaonal and on one side of the diagonal should be\n" 
                "given, and the indexing should employ Fortran format\n"
                "where the first row and column are one. Since the\n"
                "user set the parameter sym = FALSE, it appears that\n"
                "this manditory requirement is being violated.\n") ;
    }
    else if ( status == CG_MULTI_SOLVERS )
    {
        printf ("The compiled version of PASA contains both the solver MUMPS\n"
                "and Harwell routine MA57. PASA can employ at most one of\n"
                "these solvers. In the file SuiteOPTconfig/Userconfig.mk,\n"
                "remove -DUSE_MUMPS to use MA57, remove -DUSE_HARWELL to\n"
                "use MUMPS, and recompile.\n") ;
    }
    else
    {
        if ( (status > 0) && (Stat->NegDiag == TRUE) )
        {
            printf ("NOTE: A negative diagonal element was encountered in a QR "
                    "factorization.  The parameter eta2 may be too small.\n") ;
        }
    }
}

void cg_print_stat
(
    CGdata *Data  /* pointer to cgdata structure */
)
{
    CGstat *Stat = Data->Stat ;
    if ( Stat == NULL ) return ;

    printf ("\nCG_DESCENT (Version %d.%d.%d, %s) run statistics:\n\n",
        CG_MAIN_VERSION, CG_SUB_VERSION, CG_SUBSUB_VERSION, CG_DATE) ;

    int status = Stat->status ;
    if ( (status != CG_OUT_OF_MEMORY) &&
         (status != CG_STARTING_FUNCTION_VALUE_INFINITE_OR_NAN) &&
         (status != CG_QP_LINEAR_TERM_GIVEN_BUT_HPROD_MISSING) &&
         (status != CG_N_IS_EMPTY) &&
         (status != CG_MULTI_SOLVERS) )
    {
        if ( Stat->f < CGINF )
        {
            printf("Final f in CG       : %22.15e\n", Stat->f);
        }
        if ( Stat->err < CGINF )
        {
            printf("sup-norm of gradient: %22.15e\n", Stat->err);
        }

        printf("Number of  iterations: %-10ld\n", (LONG) Stat->iter);
        printf("Function  evaluations: %-10ld\n", (LONG) Stat->nfunc);
        printf("Gradient  evaluations: %-10ld\n", (LONG) Stat->ngrad) ;
        printf("Times cgexpand called: %-10ld\n", (LONG) Stat->nexpand) ;
        printf("Number forward expand: %-10ld\n", (LONG) Stat->nforward) ;
        printf("Number backward moves: %-10ld\n", (LONG) Stat->nback) ;
        printf("Newton steps (PRP_CG): %-10ld\n", (LONG) Stat->nCG) ;
        printf("Newton steps     (SS): %-10ld\n", (LONG) Stat->nSS) ;
        printf("PRP  iterations in CG: %-10ld\n", (LONG) Stat->PRP) ;
        if ( Stat->IterSub > 0 )
        {
            printf ("Subspace iterations : %-10d\n", Stat->IterSub) ;
            printf ("Number of subspaces : %-10d\n", Stat->NumSub) ;
        }
        printf ("\n") ;
    }
    else
    {
        printf ("There were no CG statistics due to the following error:\n") ;
        cg_print_status (Data) ;
    }
}

/* ==========================================================================
   === cg_print_parm ========================================================
   ==========================================================================
    Print data in the CGparm structure
   ========================================================================== */
void cg_print_parm
(
    CGdata *Data /* pointer to cgdata structure */
)
{
    CGparm *Parm = Data->Parm ;
    printf ("\nCG_DESCENT parameter settings (Version %d.%d.%d, %s):\n",
        CG_MAIN_VERSION, CG_SUB_VERSION, CG_SUBSUB_VERSION, CG_DATE) ;

    printf ("(see cg_default for definitions)\n\n") ;

    printf ("grad_tol ..............................: %e\n",
             Parm->grad_tol) ;
    printf ("PrintStatus ...........................: ") ;
             cg_print_TF (Parm->PrintStatus) ;
    printf ("PrintStat .............................: ") ;
             cg_print_TF (Parm->PrintStat) ;
    printf ("PrintParm .............................: ") ;
             cg_print_TF (Parm->PrintParm) ;
    printf ("Print level (0 = none, 3 = maximum) ...: %i\n",
             Parm->PrintLevel) ;
    printf ("QuadCost ..............................: ") ;
             cg_print_TF (Parm->QuadCost) ;
    printf ("FastLA ................................: %i\n",
             Parm->FastLA) ;
    printf ("deriv_mode ............................: %i\n",
             Parm->deriv_mode) ;
    printf ("dense_Hessian .........................: %e\n",
             Parm->dense_Hessian) ;
    printf ("CG_window .............................: %i\n",
             Parm->CG_window) ;
    printf ("Newton_window .........................: %i\n",
             Parm->Newton_window) ;
    printf ("min_error_iter ........................: %i\n",
             Parm->min_error_iter) ;
    printf ("Newton_cutoff .........................: %i\n",
             Parm->Newton_cutoff) ;
    printf ("RhoDecay ..............................: %e\n",
             Parm->RhoDecay) ;
    printf ("err_decay .............................: %e\n",
             Parm->err_decay) ;
    printf ("cg_ok .................................: %i\n",
             Parm->cg_ok) ;
    printf ("Hessian_err_decay .....................: %e\n",
             Parm->Hessian_err_decay) ;
    printf ("Hessian_CG_solver .....................: %i\n",
             Parm->Hessian_CG_solver) ;
    printf ("ss_diag_pert ..........................: %e\n",
             Parm->ss_diag_pert) ;
    printf ("cg_diag_pert ..........................: %e\n",
             Parm->cg_diag_pert) ;
    printf ("HessianSparsityFixed ..................: ") ;
             cg_print_TF (Parm->HessianSparsityFixed) ;
    printf ("TrustIterLimt .........................: %i\n",
             Parm->TrustIterLimit) ;
    printf ("big_grad ..............................: %e\n",
             Parm->big_grad) ;
    printf ("QPshift ...............................: %e\n",
             Parm->QPshift) ;
    printf ("QPgReset_factor .......................: %e\n",
             Parm->QPgReset_factor) ;
    printf ("StopFac ...............................: %e\n",
             Parm->StopFac) ;
    printf ("debug .................................: ") ;
             cg_print_TF (Parm->debug) ;
    printf ("debugtol ..............................: %e\n",
             Parm->debugtol) ;
    printf ("CheckMatrix ...........................: ") ;
             cg_print_TF (Parm->CheckMatrix) ;
    printf ("step ..................................: %e\n",
             Parm->step) ;
    printf ("LBFGS (0 = cg, 1 = lbfgs, 2 = either) .: %i\n",
             Parm->LBFGS) ;
    printf ("LBFGSmemory ...........................: %i\n",
             Parm->LBFGSmemory) ;
    printf ("maxit .................................: %ld\n",
             (LONG) Parm->maxit) ;
    printf ("restart_fac (restart_fac*n iterations).: %e\n",
             Parm->restart_fac) ;
    printf ("Qdecay (factor for averaging cost) ....: %e\n",
             Parm->Qdecay) ;
    printf ("nslow .................................: %i\n",
             Parm->nslow) ;
    printf ("QuadStep ..............................: ") ;
             cg_print_TF (Parm->QuadStep) ;
    printf ("QuadCutOff.............................: %e\n",
             Parm->QuadCutOff) ;
    printf ("QuadSafe ..............................: %e\n",
             Parm->QuadSafe) ;
    printf ("psi_lo ................................: %e\n",
             Parm->psi_lo) ;
    printf ("psi_hi ................................: %e\n",
             Parm->psi_hi) ;
    printf ("psi1 ..................................: %e\n",
             Parm->psi1) ;
    printf ("qeps ..................................: %e\n",
             Parm->qeps) ;
    printf ("qrule .................................: %e\n",
             Parm->qrule) ;
    printf ("qrestart ..............................: %i\n",
             Parm->qrestart) ;
    printf ("UseCubic ..............................: ") ;
             cg_print_TF (Parm->UseCubic) ;
    printf ("CubicCutOff ...........................: %e\n",
             Parm->CubicCutOff) ;
    printf ("SmallCost .............................: %e\n",
             Parm->SmallCost) ;
    printf ("ExpandSafe ............................: %e\n",
             Parm->ExpandSafe) ;
    printf ("SecantAmp .............................: %e\n",
             Parm->SecantAmp) ;
    printf ("RhoGrow ...............................: %e\n",
             Parm->RhoGrow) ;
    printf ("BigDfactor ............................: %e\n",
             Parm->BigDfactor) ;
    printf ("approxstep ............................: ") ;
             cg_print_TF (Parm->approxstep) ;
    printf ("ApproxSwitchFactor ....................: %e\n",
             Parm->ApproxSwitchFactor) ;
    printf ("CostConverge ..........................: %e\n",
             Parm->CostConverge) ;
    printf ("cgdelta ...............................: %e\n",
             Parm->cgdelta) ;
    printf ("cgsigma ...............................: %e\n",
             Parm->cgsigma) ;
    printf ("maxsteps ..............................: %i\n",
             Parm->maxsteps) ;
    printf ("stepdecay .............................: %e\n",
             Parm->stepdecay) ;
    printf ("cg_infdecay ...........................: %e\n",
             Parm->cg_infdecay) ;
    printf ("cg_infdecay_rate ......................: %e\n",
             Parm->cg_infdecay_rate) ;
    printf ("cg_ninf_tries .........................: %i\n",
             Parm->cg_ninf_tries) ;
    printf ("rho ...................................: %e\n",
             Parm->rho) ;
    printf ("RhoGrow ...............................: %e\n",
             Parm->RhoGrow) ;
    printf ("BigDfactor ............................: %e\n",
             Parm->BigDfactor) ;
    printf ("PertRule ..............................: %i\n",
             Parm->PertRule) ;
    printf ("pert_eps ..............................: %e\n",
             Parm->pert_eps) ;
    printf ("ncontract .............................: %i\n",
             Parm->ncontract) ;
    printf ("eps_grow ..............................: %e\n",
             Parm->eps_grow) ;
    printf ("neps (max # of times eps is updated) ..: %i\n",
             Parm->neps) ;
    printf ("psi0 ..................................: %e\n",
             Parm->psi0) ;
    printf ("psi2 ..................................: %e\n",
             Parm->psi2) ;
    printf ("BetaLower .............................: %e\n",
             Parm->BetaLower) ;
    printf ("theta .................................: %e\n",
             Parm->theta) ;
    printf ("AdaptiveTheta .........................: ") ;
             cg_print_TF (Parm->AdaptiveTheta) ;
    /* limited memory CG parameters */
    printf ("SubCheck ..............................: %i\n",
             Parm->SubCheck) ;
    printf ("SubSkip ...............................: %i\n",
             Parm->SubSkip) ;
    printf ("eta0 ..................................: %e\n",
             Parm->eta0) ;
    printf ("eta1 ..................................: %e\n",
             Parm->eta1) ;
    printf ("eta2 ..................................: %e\n",
             Parm->eta2) ;
}
