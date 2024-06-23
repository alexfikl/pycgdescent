/* ==========================================================================
   === SSM_print_status ====================================================
   ==========================================================================
    Print the status at termination of the run
   ========================================================================== */

#include "SSM.h"

void SSM_print_status
(
    int status  /* status from SSM */
)
{
    printf ("\nSSM run status (Version %d.%d.%d, %s):\n\n",
               SSM_MAIN_VERSION, SSM_SUB_VERSION, SSM_SUBSUB_VERSION,
               SSM_DATE) ;

    if ( status == SSM_TOO_MANY_ITERATIONS_IN_MINRES )
    {
        printf ("SSM used too many iterations in refinement (SSMminres)\n") ;
    }
    else if ( status == SSM_VIOLATES_SPHERE_CONSTRAINT )
    {
        printf ("SSM violates sphere constraint\n") ;
    }
    else if ( status == SSM_RESTARTS_EXCEEDS_LIMIT )
    {
        printf ("SSM exceeds limit on number of restarts\n") ;
    }
    else if ( status == SSM_DIMENSION_NONPOSITIVE )
    {
        printf ("SSM input problem dimension was not positive\n") ;
    }
    else if ( status == SSM_QR_DIAGONALIZATION_FAILS )
    {
        printf ("In SSM, unable to diagonalize the matrix\n") ;
    }
    else if ( status == SSM_NOT_ENOUGH_SPACE_FOR_QR_DIAGONALIZATION )
    {
        printf ("In SSM, unable to diagonalize the matrix\n") ;
    }
    else if ( status == SSM_OUT_OF_MEMORY )
    {
        printf ("Out of memory in SSM\n") ;
    }
}

/* ==========================================================================
   === SSMprint_parms =======================================================
   ==========================================================================
   print values in SSMparm structure
   ========================================================================== */
void SSMprint_parms
(
    SSMParm *Parm /* SSMparm structure to be printed */
)
{
    /* numerical parameters */
    printf ("print level (0 = none, 3 = maximum) .......... PrintLevel: %i\n",
             Parm->PrintLevel) ;
    printf ("machine epsilon ..................................... eps: %e\n",
             Parm->eps) ;
    printf ("subdiagonal tolerance in Lanczos iteration .. Lanczos_tol: %e\n",
             Parm->Lanczos_tol) ;
    printf ("change from Householder to Lanczos at dim ... House_limit: %i\n",
             (int) Parm->House_limit) ;
    printf ("relative eigenvalue convergence tol .......... qr_tol_rel: %e\n",
             Parm->qr_tol_rel) ;
    printf ("absolute eigenvalue convergence tol .......... qr_tol_abs: %e\n",
             Parm->qr_tol_abs) ;
    printf ("lower bound for QR scaling factor .............. qr_lower: %e\n",
             Parm->qr_lower) ;
    printf ("max number QR algorithm iterations qr_its * dim .. qr_its: %i\n",
             (int) Parm->qr_its) ;
    printf ("multiplier accuracy for diagonal matrix .. diag_optim_tol: %e\n",
             Parm->diag_optim_tol) ;
    printf ("accuracy of diag in tridiagonal matrix .... diag_diag_eps: %e\n",
             Parm->diag_eps) ;
    printf ("tolerance used when checking cost decay ....... check_tol: %e\n",
             Parm->check_tol) ;
    printf ("tolerance used when checking KKT error ........ check_kkt: %e\n",
             Parm->check_kkt) ;
    printf ("radius_flex times radius ok in interior ..... radius_flex: %e\n",
             Parm->radius_flex) ;
    printf ("error decay factor in SQP iteration ........... sqp_decay: %e\n",
             Parm->sqp_decay) ;
    printf ("error decay factor in inverse power iteration . eig_decay: %e\n",
             Parm->eig_decay) ;
    printf ("all eigenvalues of matrix are >= .............. eig_lower: %e\n",
             Parm->eig_lower) ;
    printf ("error decay factor in SSM iteration ........... ssm_decay: %e\n",
             Parm->ssm_decay) ;
    printf ("mu multiplication factor if mu = 0 possible ...... shrink: %e\n",
             Parm->shrink) ;
    printf ("max number MINRES iterations is fac*n .... minres_its_fac: %e\n",
             Parm->minres_its_fac) ;
    printf ("max number SSM iterations is fac*n .......... ssm_its_fac: %e\n",
             Parm->ssm_its_fac) ;
    printf ("upper bound on number of Lanczos iterations . Lanczos_bnd: %i\n",
             (int) Parm->Lanczos_bnd) ;
    printf ("restart Lanczos, iterations grow by factor . grow_Lanczos: %e\n",
             Parm->grow_Lanczos) ;
    printf ("number of Lanczos restarts attempted ........... nrestart: %i\n",
             Parm->nrestart) ;

    printf ("\n") ;
    printf ("The number of Lanczos iterations L is given by the formula:\n"
            "    if      ( n <= dim1 ) L = n-1\n"
            "    else if ( n <= dim2 ) L = L1\n"
            "    else                  L = L1 + L2*(log10 (n/dim2)\n") ;
    printf ("...................................................  dim1: %i\n",
             (int) Parm->Lanczos_dim1) ;
    printf ("...................................................  dim2: %i\n",
             (int) Parm->Lanczos_dim2) ;
    printf ("...................................................    L1: %i\n",
             (int) Parm->Lanczos_L1) ;
    printf ("...................................................    L2: %i\n",
             (int) Parm->Lanczos_L2) ;

    /* logical parameters */
    printf ("\nLogical parameters:\n") ;
    if ( Parm->BndLessThan == SSMTRUE )
        printf ("    Constraint is ||x|| <= r\n") ;
    else
        printf ("    Constraint is ||x|| = r\n") ;
    if ( Parm->PrintParms == SSMTRUE )
        printf ("    Print the parameter structure\n") ;
    else
        printf ("    Do not print parameter structure\n") ;
    if ( Parm->PrintFinal == SSMTRUE )
        printf ("    Print final status and statistics\n") ;
    else
        printf ("    Do not print final status and statistics\n") ;
    if ( Parm->IPM == SSMTRUE )
        printf ("    Use inverse power method to estimate smallest "
                     "eigenvalue\n") ;
    else
        printf ("    Use SQP method to estimate smallest eigenvalue\n") ;
}
