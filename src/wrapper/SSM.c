/* =========================================================================
   ================================ SSM ====================================
   =========================================================================
       ________________________________________________________________
      |    Solve a sphere constrained quadratic program of the form    |
      |                                                                |
      |           1                                                    |
      |       min - x'Ax + b'x  subject to ||x|| <= r or ||x|| = r     |
      |           2                                                    |
      |                                                                |
      |    using an implementation of the sequential subspace method   |
      |    which includes:                                             |
      |                                                                |
      |       1. SQP acceleration                                      |
      |       2. Solution of SQP system using MINRES and SSOR          |
      |          preconditioning                                       |
      |                                                                |
      |                     Version 2.0 (May 26, 2023)                 |
      |                  Version 1.1 (September 25, 2009)              |
      |                     Version 1.0 (May 5, 2009)                  |
      |                                                                |
      |                         William W. Hager                       |
      |                        hager@math.ufl.edu                      |
      |                     Department of Mathematics                  |
      |                       University of Florida                    |
      |                     Gainesville, Florida 32611                 |
      |                         352-392-0281 x 244                     |
      |                   http://www.math.ufl.edu/~hager               |
      |                                                                |
      |                   Copyright by William W. Hager                |
      |________________________________________________________________|

       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|*/

/* Contents:

    1.  SSM         - solve sphere constrained optimization problem
    2.  SSMdefault  - sets default parameter values in the SSMParm structure
    3.  SSMallocate - allocates the arrays needed in the SSM algorithm
    4.  SSMdestroy  - free the allocated memory
    5.  SSMinitProb - store problem, order rows in col, locate 1st super diag
    6.  SSMprint_parms - print the parameter array
    7.  SSMballdense- solve a dense sphere constrained QP
    8.  SSMsubspace - solve the SSM subspace problem
    9.  SSMorth     - generate orthonormal vectors using Householder
    10. SSMmineig   - eigenvector associated with smallest eigenvalue
    11. SSMdiag     - diagonalize a tridiagonal matrix
    12. SSMcheckKKT - check KKT conditions (only used in debug mode) */

#include "SSM.h"

int SSM /* return 0 (error tolerance satisfied)
                        -1 (min residual convergence failure in SQP)
                        -2 (SSM failed to converge in specified iteration)
                        -3 (number of SSM restarts exceeded limit)
                        -4 (dimension <= 0)
                        -5 (failure of QR diagonalization)
                        -6 (insufficient space allocated for QR diagonalization)
                        -7 (starting Lanczos vector vanishes) */
(
/* output: */
    SSMFLOAT     *x, /* solution (output) */
    ssm_stat  *Stat, /* NULL means do not return statistics */

/* input: */
    SSMINT        n, /* size of x */
    SSMFLOAT    *Ax, /* numerical entries in matrix, packed array */
    SSMINT      *Ai, /* row indices for each column */
    SSMINT      *Ap, /* points into Ax and Ai, Ap [j] = start of column j */
    SSMFLOAT     *b, /* linear term in objective function */
    SSMFLOAT      r, /* radius of ball */
    SSMFLOAT    tol, /* 2-norm KKT error tolerance for solution */
    SSMFLOAT *guess, /* starting, NULL means no guess given */
    SSMParm  *Uparm/* user parameters, NULL means use default parameters */
)
{
    int it, Allocate, PrintLevel ;
    SSMINT j ;
    SSMFLOAT normx, t, *xnew, *v, *vnew ;
    SSMDiagOpt *DO ;
    SSMDiag *DT ;
    SSMLanczos *LH ;
    SSMParm *Parm, ParmStruc ;
    SSMProblem PB ;
    SSMcom Com ;
    int status = 0 ;

    /*print_matrix (Ai, Ap, Ax, "A", n) ;*/
    /* initialize counters for statistics */
    Com.ssm_its = Com.minres_its = 0 ;
    Com.mults = SSMZERO ;
    Com.cost_old = SSMINF ;
    Com.emin_old = SSMINF ;
    Com.error = SSMINF ;
    it = 0 ;

    /* no memory has been allocated */
    Allocate = SSMFALSE ;
    if ( n <= 0 )
    {
        status = -4 ;
        goto Exit ;
    }

    /* ---------------------------------------------------------------------- */
    /* set parameter values */
    /* ---------------------------------------------------------------------- */
    Parm = Uparm ;
    if ( Uparm == NULL )
    {
        Parm = &ParmStruc ;
        SSMdefault (Parm) ;
    }
    PrintLevel = Parm->PrintLevel ;

    /* ---------------------------------------------------------------------- */
    /* if Parm->PrintParms is true, then print parameter values */
    /* ---------------------------------------------------------------------- */
    if ( Parm->PrintParms ) SSMprint_parms (Parm) ;

    /* trivial case n = 1 */
    if ( n == 1 )
    {
        if ( Parm->BndLessThan ) /* ||x|| <= r */
        {
            if ( Ap [0] )
            {
                t = -b [0]/Ax [0] ;
                if ( fabs (t) > r ) t = r*t/fabs (t) ;
            }
            else if ( b [0] ) t = -r*b [0]/fabs (b [0]) ;
            else              t = SSMZERO ;
        }
        else
        {
            if ( Ap [0] )
            {
                t = -b [0]/Ax [0] ;
                if ( t != SSMZERO ) t = r*t/fabs (t) ;
                else                t = r ;
            }
            else if ( b [0] ) t = -r*b [0]/fabs (b [0]) ;
            else              t = r ;
        }
        x [0] = t ;
        Com.error = SSMZERO ;
        goto Exit ;
    }

    /* check if x = 0 is a solution */
    if ( Parm->BndLessThan )
    {
        t = SSMZERO ;
        for (j = 0; j < n; j++) t += b [j]*b [j] ;
        t = sqrt (t) ;
        if ( t <= tol )
        {
            for (j = 0; j < n; j++) x [j] = SSMZERO ;
            Com.error = t ;
            goto Exit ;
        }
    }

    Com.tol = tol ; /* KKT error tolerance for solution */

    /* set srand, if requested (default is 1) */
    if ( Parm->seed != 0 )
    {
        srand (Parm->seed) ;
    }

    /* allocate memory */
    SSMallocate (&Com, &PB, Parm, n) ;
    /* memory has been allocated */
    Allocate = SSMTRUE ;
    v = Com.v ;
    vnew = Com.MINRES+(5*n) ;
    xnew = Com.MINRES+(6*n) ;
    /* NOTE: Lanczos starting guess after a restart is computed by
             SSMreallocate and stored in Com.MINRES+(8*n) */

    /* allocate memory for PB and setup problem structure
       (sort cols, extract diagonal D and last nonzero in each column) */
    SSMinitProb (Ax, Ai, Ap, n, b, r, &PB, Parm, &Com) ;
    if ( PB.Amax == SSMZERO ) /* the matrix is zero */
    {
        t = SSMZERO ;
        for (j = 0; j < n; j++) t += b [j]*b [j] ;
        if ( t != SSMZERO )
        {
            t = r/sqrt (t) ;
            for (j = 0; j < n; j++) x [j] = -t*b [j] ;
        }
        else /* anything is optimal */
        {
            for (j = 0; j < n; j++) x [j] = SSMZERO ;
            x [0] = r ;
        }
        goto Exit ;
    }

    DO = PB.DO ;
    DT = PB.DT ;
    LH = PB.LH ;

    /* if SSM does not converge, then restart the tridiagonalization
       process in order to obtain a better starting guess for SSM.
       Below the variable "it" is the number of times the
       tridiagonalization has been restarted */

    for (it = 0; it < Parm->nrestart; it++)
    {
        /* tridiagonalize the matrix A, populate LH structure */
        if ( PrintLevel < 0 )
        {
            printf ("before tridiag\n") ;
            fflush (stdout) ;
        }
        if ( it == 0 ) status = SSMtridiag (guess, it, LH, &PB, &Com) ;
        else           status = SSMtridiag (    x, it, LH, &PB, &Com) ;
        if ( PrintLevel < 0 )
        {
            printf ("after tridiag\n") ;
            fflush (stdout) ;
        }
        if ( status ) goto Exit ; /* Lanczos starting vector vanishes */
 
        /* diagonalize the tridiagonal matrix, populate the DT structure */
        if ( PrintLevel < 0 )
        {
            printf ("before diag\n") ;
            fflush (stdout) ;
        }
        status = SSMdiag (LH->d, LH->u, LH->ncols, DT,
                          Com.wi1, Com.wx1, Com.wx2, Com.wx3) ;
        if ( PrintLevel < 0 )
        {
            printf ("after diag\n") ;
            fflush (stdout) ;
        }
        if ( status ) goto Exit ; /* failure of QR diagonalization */

        /* copy eigenvalues and dimension into the SSMDiagOpt structure*/
        DO->n = LH->ncols ;
        for (j = 0; j < LH->ncols; j++) DO->d [j] = DT->e [j] ;

        if ( PrintLevel < 0 )
        {
            printf ("before mult\n") ;
            fflush (stdout) ;
        }
        /* Compute (b)'V and store it in the SSMDiagOpt structure */
        SSMtDenseMult (DO->f, b, LH->V, n, LH->ncols) ;
        if ( PrintLevel < 0 )
        {
            printf ("after mult\n") ;
            fflush (stdout) ;
        }

        /* Multiply DO->f by the Givens rotations x' G_1 G_2 ... */
        if ( PrintLevel < 0 )
        {
            printf ("before G multmult\n") ;
            fflush (stdout) ;
        }
        SSMGivensMult (DO->f, DT) ;
        if ( PrintLevel < 0 )
        {
            printf ("after G mult\n") ;
            fflush (stdout) ;
        }

        /* Solve reduced problem */
        if ( PrintLevel < 0 )
        {
            printf ("before diagopt\n") ;
            fflush (stdout) ;
        }
        SSMdiagopt (Com.wx3, r, Parm->BndLessThan, DO, &Com) ;
        if ( PrintLevel < 0 )
        {
            printf ("after diagopt\n") ;
            fflush (stdout) ;
        }
        Com.mu = DO->mu ;
        if ( DO->mu > SSMZERO ) Com.Active = SSMTRUE ;
        else                    Com.Active = SSMFALSE ;

        /* Multiply solution of (WP) by rotations G_1 G_2 ... x */
        if ( PrintLevel < 0 )
        {
            printf ("before tG multmult\n") ;
            fflush (stdout) ;
        }
        SSMtGivensMult (Com.wx3, DT) ;
        if ( PrintLevel < 0 )
        {
            printf ("after tG mult\n") ;
            fflush (stdout) ;
        }

        /* Multiply by V and solve the subspace problem */
        if ( it == 0 ) /* no prior eigenvector estimate exists */
        {
            /* evaluate starting estimate of solution to sphere problem */
            SSMDenseMult (x, Com.wx3, LH->V, n, LH->ncols) ;

            /* estimate smallest eigenvalue and associated eigenvector */
            SSMmineig (&Com.emin, v, Com.wx1, &PB) ;

            /* evaluate kkt error and eigenvector residual */
            status = SSMsubspace (x, t, Com.Ax, v, NULL, NULL, NULL, 2,
                                  SSMTRUE, &Com);
            if ( status ) goto Exit ; /* failure of QR diagonalization */
        }
        else           /* preserve prior eigenvector estimate */
        {
            if ( PrintLevel < 0 )
            {
                printf ("before dense mult\n") ;
                fflush (stdout) ;
            }
            SSMDenseMult (xnew, Com.wx3, LH->V, n, LH->ncols) ;
            if ( PrintLevel < 0 )
            {
                printf ("after dense mult\n") ;
                fflush (stdout) ;
            }
            /* new eigenvector estimate */
            SSMmineig (&Com.emin, vnew, Com.wx1, &PB) ;

            /* compute ||x|| */
            normx = SSMZERO ;
            for (j = 0; j < n; j++)
            {
                t = x [j] ;
                normx += t*t ;
            }
            normx = sqrt (normx) ;

            /* solve subspace problem, basis: x, xnew, v, vnew */
            if ( PrintLevel < 0 )
            {
                printf ("before subspace\n") ;
                fflush (stdout) ;
            }
            status = SSMsubspace (x, normx, Com.Ax, xnew, v, vnew,
                                  NULL, 4, SSMFALSE, &Com) ;
            if ( PrintLevel < 0 )
            {
                printf ("after subspace\n") ;
                fflush (stdout) ;
            }
            if ( status ) goto Exit ; /* QR method convergence failure */
        }
        if ( PrintLevel >= 1 )
        {
            printf ("\nbegin restart #: %i ball_err: %e ball_tol: %e\n",
                 it, Com.error, Com.tol) ;
            printf ("emin: %e mu: %e eig_err: %e\n",
                  Com.emin, Com.mu, Com.eig_error) ;
        }

        /* for a small problem, we stop immediately */
        status = 0 ;
        if ( n <= 5 ) goto Exit ;
        if ( Com.error <= tol ) goto Exit ;

        /* refine solution */
        if ( Parm->BndLessThan && !Com.Active ) status = -8 ;
        else                                    status = -9 ;
        while ( status <= -8 )
        {
            if ( (status == -9) || (PB.Dmin < SSMZERO) )
            {
                if ( PrintLevel >= 1 ) printf ("\nBOUNDARY SOLUTION\n") ;
                status = SSMboundary (x, &Com) ;
            }
            else
            {
                if ( PrintLevel >= 1 ) printf ("\nINTERIOR SOLUTION\n") ;
                status = SSMinterior (x, &Com) ;
            }
            if ( PrintLevel >= 1 )
            {
                if ( status > -8 ) 
                {
                    printf ("\ndone with restart #: %i error: %e status: %i\n",
                             it, Com.error, status) ;
                }
                else
                {
                    printf ("\ncontinue restart #: %i error: %e status: %i\n",
                             it, Com.error, status) ;
                }
            }
            /* check if SSM failed to converge or QR method failed */
            if ( (status == -2) || (status == -5) || (status == -6) ) goto Exit;
        }

#ifndef NDEBUG
        /* check KKT conditions for x */
        if ( status == 0 ) SSMcheckKKT (x, r, &PB, &Com) ;
#endif
        /* convergence tolerance satisfied */
        if ( status == 0 ) goto Exit ;

        /* else SSM failed to converge, generate a larger tridiagonal matrix */
        SSMreallocate (x, &Com) ;
    }

    /* convergence not achieved within specified number of restarts */
    status = -3 ;

    Exit:
    if ( Stat != NULL )
    {
        Stat->mults = Com.mults ;
        Stat->ssm_its = Com.ssm_its ;
        Stat->minres_its = Com.minres_its ;
        Stat->restart_its = it ;
        Stat->error = Com.error ;
    }
    if ( Parm->PrintFinal || (PrintLevel >= 1) )
    {
        printf ("\n") ;
        if ( status == 0 )
        {
            printf ("Convergence tolerance statisfied\n" ) ;
        }
        else if ( status == -1 )
        {
            printf ("Minimum residual algorithm failed to converge.\n") ;
            printf ("The iteration limit %i was reached\n",
                     (int) Com.minres_limit) ;
        }
        else if ( status == -2 )
        {
            printf ("SSM failed to converge in %i iterations\n",
                     (int) Com.ssm_limit);
        }
        else if ( status == -3 )
        {
            printf ("SSM failed to converge within %i restart\n",
                     Parm->nrestart) ;
        }
        else if ( status == -4 )
        {
            printf ("The specified problem dimension %i < 0\n", (int) n) ;
        }
        else if ( status == -5 )
        {
            printf ("The QR method was unable to diagonalize the "
                    "tridiagonal matrix\n") ;
        }
        else if ( status == -6 )
        {
            printf ("In the QR method, there was not enough space to store\n") ;
            printf ("the Givens rotations. The parameter Parm->qr_its\n") ;
            printf ("should be increased so that more space will be "
                    "allocated.\n") ;
        }
        else if ( status == -7 )
        {
            printf ("Starting vector in the Lanczos process vanishes\n") ;
        }
        printf ("\n") ;
        printf ("Number multiplications by matrix:%10i\n", (int) Com.mults) ;
        printf ("SSM iterations .................:%10i\n", (int) Com.ssm_its);
        printf ("Minimum residual iterations ....:%10i\n",
                (int) Com.minres_its);
        printf ("Number of Lanczos restarts .....:%10i\n", it);
        printf ("2-norm of KKT error ............:%10.3e\n", Com.error);
    }
    /* destroy allocated memory */
    if ( Allocate ) SSMdestroy (&Com) ;
    return (status) ;
}

/* ==========================================================================
   === SSMdefault =========================================================== 
   ========================================================================== 
    Set default parameter values.
   ========================================================================== */

void SSMdefault
(
    SSMParm  *Parm /* pointer to parameter structure */
)
{
    SSMFLOAT t, eps ;

    /* random number generator seed for srand */
    Parm->seed = 1 ;

    /* compute machine epsilon for SSMFLOAT */
    eps = SSMONE ;
    t = SSMONE ;
    while ( t > 0 )
    {
        eps /= 2. ;
        t = SSMONE + eps ;
        t -= SSMONE ;
    }
    eps *= 2 ;
    Parm->eps = eps ;

    /* the number of Lanczos iterations L is given by the formula:
       if      ( n <= dim1 ) L = n-1 ;
       else if ( n <= dim2 ) L = L1 ;
       else                  L = L1 + L2*(log10 (n/dim2)) ; */
    Parm->Lanczos_dim1 = 30 ;
    Parm->Lanczos_dim2 = 100 ;
    Parm->Lanczos_L1 = 30 ;
    Parm->Lanczos_L2 = 80 ;
    Parm->grow_Lanczos = 1.3e0;    /*factor to grow size of Lanczos startup*/
    Parm->Lanczos_bnd = 10000 ;    /* upper bound on number Lanczos iterations*/
    Parm->Lanczos_tol = eps ;      /* stop Lanczos when |u_j| <= tol */
    Parm->House_limit = 30 ;       /* use Householder when n <= limit */
    Parm->qr_its = 4 ;             /* mult qr_its by dim for max # qr steps*/
    Parm->qr_tol_rel = eps ;       /* relative eigenvalue convergence tol*/
    Parm->qr_tol_abs = eps*1.e-3  ;/* absolute eigenvalue convergence tol*/

    /* lower bound for qr scaling factor; if this is made large, then
       allocation for DT->gi and DT->gs should be increased so that
       -log_2(qr_lower) < the factor in these allocation */
    Parm->qr_lower = 1.e-9 ;

    Parm->diag_optim_tol = 1.e5*eps ;/* accuracy of mult. for diag. matrix */

    /* if |d [i]-d [j]| <= diag_eps * n * absolute maximum diagonal element
       in diagopt, then d [i] and d [j] considered equal */
    Parm->diag_eps = 100*eps ;

    Parm->check_tol = 1.e-8 ;  /* tolerance used when checking SSM structure */
    Parm->check_kkt = 1.e-6 ;  /* KKT error tolerance is check_kkt*n */
    Parm->eig_refine = 5.e-1 ; /* factor determines when to refine eig in SSM */

    /* immediately switch from interior routine to boundary routine when
       interior iterate has norm >= radius_flex * sphere radius */
    Parm->radius_flex = 1.1 ;

    /* Factor which determines when to stop SQP iteration. Error in SQP
       system <= sqp_decay times KKT error */
    Parm->sqp_decay = 1.e-2 ;

    /* Factor which determines when to stop inverse power iteration.
       Norm of the residual <= eig_decay */
    Parm->eig_decay = 1.e-2 ;

    /* eig_lower is an actual lower bound for the smallest eigenvalue,
       not an estimate. Such a lower bound is useful in the safeguarding
       process */
    Parm->eig_lower = -SSMINF ;

    /* Factor which determines when to restart iteration. Error in SSM
       should decay by ssm_decay in one iteration, otherwise perform
       more Lanczos iterations */
    Parm->ssm_decay = 8.e-1 ;

    /* when optimal multiplier mu in diagopt is near zero, shrink to zero
       by factor shrink in each iteration */
    Parm->shrink = 1.e-1 ;     /* mu multiplied by shrink in diagopt */
    Parm->BndLessThan = SSMTRUE ;/* constraint is ||x|| <= r */

    /* TRUE  (use inverse power method to estimate smallest eigenvalue)
       FALSE (use SQP method to estimate smallest eigenvalue) */
    Parm->IPM = SSMFALSE ;

    Parm->minres_its_fac = 3 ;  /* max number MINRES iterations is fac*n */
    Parm->ssm_its_fac = 1 ;     /* max number SSM iterations is fac*n */
    Parm->nrestart = 40 ;       /* number of Lanczos restarts attempted */

    /* PrintLevel = 0, 1, 2, 3
       PrintLevel = 3 gives maximum printing of iteration data while
       PrintLevel = 0 gives no printing of iteration data */
    Parm->PrintLevel = 0 ;

    /* PrintParms = TRUE  means to print parameter values */
    Parm->PrintParms = SSMFALSE ;

    /* PrintFinal = TRUE means to print final status and statistics */
    Parm->PrintFinal = SSMFALSE ;
}

/* ==========================================================================
   === SSMallocate ==========================================================
   ==========================================================================
   Allocate memory for the SSM structures and copy parameter values into SSMcom 
   ==========================================================================*/

int SSMallocate
(
    SSMcom    *Com, /* pointer to SSMcom structure */
    SSMProblem *PB, /* problem specification */
    SSMParm  *Parm, /* parameters, needed to determine allocation */
    SSMINT       n  /* problem dimension */
)
{
    SSMINT L, N, qr_its ;
    SSMDiag *DT ;
    SSMLanczos *LH ;
    SSMDiagOpt *DO ;
    SSMProblem *P ;
    int status = 0 ;

    /* malloc work arrays */
    Com->wi1 = sopt_malloc (&status, n, sizeof (SSMINT)) ;
    Com->wi2 = sopt_malloc (&status, n, sizeof (SSMINT)) ;
    Com->wx1 = sopt_malloc (&status, 3*n, sizeof (SSMFLOAT)) ;
    /* estimate eigenvector */
    Com->v = sopt_malloc (&status, n, sizeof (SSMFLOAT)) ;
    Com->Av = sopt_malloc (&status, n, sizeof (SSMFLOAT)) ; /* A times v */
    Com->Ax = sopt_malloc (&status, n, sizeof (SSMFLOAT)) ; /* A times x */
    Com->r0 = sopt_malloc (&status, n, sizeof (SSMFLOAT)) ; /* r0 = Ax + b */

    /* create Lanczos/Householder structure */
    LH = sopt_malloc (&status, (SSMINT) 1, sizeof (SSMLanczos)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;
    Com->wx2 = (Com->wx1)+n ;
    Com->wx3 = (Com->wx2)+n ;
    PB->LH = LH ;
    /* max number of Lanczos iterations, needed for mallocs */
    if ( n <= Parm->House_limit ) L = n ; /* need space for full matrix */
    else if ( n <= Parm->Lanczos_dim1 ) L = n-1 ;
    else if ( n <= Parm->Lanczos_dim2 ) L = Parm->Lanczos_L1 ;
    else L = Parm->Lanczos_L1 + Parm->Lanczos_L2*
             (log10 ((double) n/ (double) Parm->Lanczos_dim2)) ;
    /* allocate space for Lanczos vectors/Householder vectors */
    LH->V = sopt_malloc (&status, L*n, sizeof (SSMFLOAT)) ;
    /* diagonal of tridiagonal matrix */
    LH->d = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
    /* superdiagonal of tridiag. matrix */
    LH->u = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;

    /* create the diagonalization structure for tridiagonal matrix */
    DT = sopt_malloc (&status, (SSMINT) 1, sizeof (SSMDiag)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;
    PB->DT = DT ;
    DT->e = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;  /* eigenvalues */

    /* max number of qr iteration */
    qr_its = (Parm->qr_its)*L ;          /* max QR algorithm iterations */
    DT->nalloc = L ;
    DT->qr_lower = Parm->qr_lower ;
    DT->gx = sopt_malloc (&status, qr_its*L, sizeof (SSMFLOAT)) ;
    DT->gy = sopt_malloc (&status, qr_its*L, sizeof (SSMFLOAT)) ;
    DT->gz = sopt_malloc (&status, qr_its*L, sizeof (short)) ;
    /* space for ending -1 in each iteration, the L final rescaling at end,
       and the in-iteration recaling which could occur at most every
       30 steps when qr_lower = 1.e-9 */
    DT->gs = sopt_malloc (&status, (qr_its + L + 5 + L*(qr_its/30)),
                                                            sizeof (SSMFLOAT)) ;
    DT->gi = sopt_malloc (&status, (qr_its + L + 5 + L*(qr_its/30)),
                                                              sizeof (SSMINT)) ;
    DT->gj = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;
    DT->gk = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;
    DT->gf = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;

    /* create the diagonal matrix optimization structure */
    DO = sopt_malloc (&status, (SSMINT) 1, sizeof (SSMDiagOpt)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;
    PB->DO = DO ;
    DO->d = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
    DO->f = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
    DO->f2= sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
    DO->dshift = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;

    /* storage for d, s, p, q arrays used in SSOR iteration*/
    Com->SSOR = sopt_malloc (&status, 4*n, sizeof (SSMFLOAT)) ;

    /*for AV (5n) or MINRES y, Ay, r, aj, vj, vj1, wj, pr, rj1, rj2, zj1, zj2*/
    Com->MINRES = sopt_malloc (&status, 12*n, sizeof (SSMFLOAT)) ;

    Com->W = sopt_malloc (&status, 5*n, sizeof (SSMFLOAT)) ;
    Com->V = sopt_malloc (&status, 5*n, sizeof (SSMFLOAT)) ;
    Com->VAV = sopt_malloc (&status, 25, sizeof (SSMFLOAT)) ;

    /* set up a dense problem structure for matrices of size up to 5 by 5 */
    P = sopt_malloc (&status, (SSMINT) 1, sizeof (SSMProblem)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;
    Com->PBdense = P ;
    N = 5 ;
    P->b = sopt_malloc (&status, N, sizeof (SSMFLOAT)) ;

    /* create Lanczos/Householder structure */
    LH = sopt_malloc (&status, (SSMINT) 1, sizeof (SSMLanczos)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;
    P->LH = LH ;
    /* default number of Lanczos iterations, needed for mallocs */
    /* orthonormal vectors */
    LH->V = sopt_malloc (&status, N*N, sizeof (SSMFLOAT)) ;
    /* diagonal of tridiagonal matrix*/
    LH->d = sopt_malloc (&status, N, sizeof (SSMFLOAT)) ;
    /* superdiag. of tridiag. matrix */
    LH->u = sopt_malloc (&status, N, sizeof (SSMFLOAT)) ;

    /* create the diagonalization structure for tridiagonal matrix */
    DT = sopt_malloc (&status, (SSMINT) 1, sizeof (SSMDiag)) ;
    P->DT = DT ;
    DT->e = sopt_malloc (&status, N, sizeof (SSMFLOAT)) ;/* eigenvalues */
    qr_its = 30 ;
    DT->nalloc = N ;
    DT->qr_lower = Parm->qr_lower ;
    DT->gx = sopt_malloc (&status, qr_its*N, sizeof (SSMFLOAT)) ;
    DT->gy = sopt_malloc (&status, qr_its*N, sizeof (SSMFLOAT)) ;
    DT->gz = sopt_malloc (&status, qr_its*N, sizeof (short)) ;
    DT->gs = sopt_malloc (&status, (SSMINT) (2*qr_its + N + 5 + N*(qr_its/30)),
                                                            sizeof (SSMFLOAT)) ;
    DT->gi = sopt_malloc (&status, (SSMINT) (2*qr_its + N + 5 + N*(qr_its/30)),
                                                            sizeof (SSMINT)) ;
    DT->gj = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;
    DT->gk = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;
    DT->gf = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;

    /* create the diagonal matrix optimization structure */
    DO = sopt_malloc (&status, (SSMINT) 1, sizeof (SSMDiagOpt)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;
    P->DO = DO ;
    DO->d = sopt_malloc (&status, N, sizeof (SSMFLOAT)) ;
    DO->f = sopt_malloc (&status, N, sizeof (SSMFLOAT)) ;
    DO->f2= sopt_malloc (&status, N, sizeof (SSMFLOAT)) ;
    DO->dshift = sopt_malloc (&status, N, sizeof (SSMFLOAT)) ;
    return (sopt_convert_error (LSSM, status)) ;
}

/* ==========================================================================
   === SSMreallocate ========================================================
   ==========================================================================
   SSM failed to converge. Tridiagonalize a larger matrix in order to
   obtain a better starting guess. Need to reallocate storage for the
   Lanczos process.
   ==========================================================================*/

int SSMreallocate
(
    SSMFLOAT    *x, /* current estimate of solution */
    SSMcom    *Com  /* SSMcom structure */
)
{
    SSMINT j, L, n, qr_its ;
    SSMFLOAT mu, *start_vector ;
    SSMDiag *DT ;
    SSMLanczos *LH ;
    SSMDiagOpt *DO ;
    SSMProblem *PB ;
    SSMParm *Parm ;
    int status = 0 ;

    PB = Com->PB ;
    Parm = Com->Parm ;
    n = PB->n ;
    LH = PB->LH ;
    DT = PB->DT ;
    DO = PB->DO ;
    L = DT->nalloc ;

    /* max number of Lanczos iterations, needed for mallocs */
    L *= Parm->grow_Lanczos ;
    if ( L >= n ) L = n ;
    else /* compute starting guess for Lanczos process */
    {
        L = SSMMIN (L, Parm->Lanczos_bnd) ;
        start_vector = Com->MINRES+(8*n) ;
        mu = Com->mu ;
        for (j = 0; j < n; j++) start_vector [j] = Com->r0 [j] + mu*x [j] ;
        /* remove projections on prior Lanczos vectors */
/*      ncols = LH->ncols ;
        Vj = LH->V ;

        for (j = 0; j < ncols; j++)
        {
            t = SSMZERO ;
            for (i = 0; i < n; i++) t += start_vector [i]*Vj [i] ;
            for (i = 0; i < n; i++) start_vector [i] -= t*Vj [i] ;
        } */
    }
    if ( DT->nalloc < L )
    {
        DT->nalloc = L ;
        /* free the prior Lanczo/Householder structure */
        ssm_free (LH->d) ;
        ssm_free (LH->u) ;
        ssm_free (LH->V) ;

        /* create new Lanczos structure */
        /* orthonormal vectors */
        LH->V = sopt_malloc (&status, L*n, sizeof (SSMFLOAT)) ;
        /*diagonal of tridiagonal matrix*/
        LH->d = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
        /* superdiagonal of tri. matrix */
        LH->u = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
        if ( status ) return (sopt_convert_error (LSSM, status)) ;
        LH->max_its = L ;

        /* free the prior tridiagonalization structure */
        ssm_free (DT->e) ;
        ssm_free (DT->gx) ;
        ssm_free (DT->gy) ;
        ssm_free (DT->gz) ;
        ssm_free (DT->gs) ;
        ssm_free (DT->gi) ;
        ssm_free (DT->gj) ;
        ssm_free (DT->gk) ;
        ssm_free (DT->gf) ;
        /* create the diagonalization structure for tridiagonal matrix */
        DT->e = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;  /* eigenvalues */
        qr_its = (Parm->qr_its)*L ;          /* max QR algorithm iterations */
        DT->max_its = SSMMAX (30, qr_its*L) ;
        DT->gx = sopt_malloc (&status, qr_its*L, sizeof (SSMFLOAT)) ;
        DT->gy = sopt_malloc (&status, qr_its*L, sizeof (SSMFLOAT)) ;
        DT->gz = sopt_malloc (&status, qr_its*L, sizeof (short)) ;
        DT->gs = sopt_malloc (&status, (SSMINT) (qr_its + L + 5 +L*(qr_its/30)),
                                                            sizeof (SSMFLOAT)) ;
        DT->gi = sopt_malloc (&status, (SSMINT) (qr_its + L + 5 +L*(qr_its/30)),
                                                              sizeof (SSMINT)) ;
        DT->gj = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;
        DT->gk = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;
        DT->gf = sopt_malloc (&status, qr_its, sizeof (SSMINT)) ;
        if ( status ) return (sopt_convert_error (LSSM, status)) ;

        /* free the prior diagonal optimization structure */
        ssm_free (DO->d) ;
        ssm_free (DO->f) ;
        ssm_free (DO->f2) ;
        ssm_free (DO->dshift) ;
        /* create the diagonal matrix optimization structure */
        DO->d = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
        DO->f = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
        DO->f2= sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
        DO->dshift = sopt_malloc (&status, L, sizeof (SSMFLOAT)) ;
        return (sopt_convert_error (LSSM, status)) ;
    }
}

/* ==========================================================================
   === SSMdestroy ===========================================================
   ==========================================================================
    Free the storage memory
   ========================================================================== */
void SSMdestroy
(
    SSMcom *Com  /* SSMcom structure to free */
)
{
    SSMProblem *PB, *P ;

    ssm_free (Com->wi1) ;
    ssm_free (Com->wi2) ;
    ssm_free (Com->wx1) ;
    ssm_free (Com->v) ;
    ssm_free (Com->Av) ;
    ssm_free (Com->Ax) ;
    ssm_free (Com->r0) ;

    PB = Com->PB ;
    ssm_free (PB->D) ;
    ssm_free (PB->i) ;
    ssm_free (PB->p) ;
    ssm_free (PB->u) ;
    ssm_free (PB->x) ;
    ssm_free (PB->LH->V) ;
    ssm_free (PB->LH->d) ;
    ssm_free (PB->LH->u) ;
    ssm_free (PB->LH) ;

    ssm_free (PB->DT->e) ;
    ssm_free (PB->DT->gx) ;
    ssm_free (PB->DT->gy) ;
    ssm_free (PB->DT->gz) ;
    ssm_free (PB->DT->gs) ;
    ssm_free (PB->DT->gi) ;
    ssm_free (PB->DT->gj) ;
    ssm_free (PB->DT->gk) ;
    ssm_free (PB->DT->gf) ;
    ssm_free (PB->DT) ;

    ssm_free (PB->DO->d) ;
    ssm_free (PB->DO->f) ;
    ssm_free (PB->DO->f2) ;
    ssm_free (PB->DO->dshift) ;
    ssm_free (PB->DO) ;

    P = Com->PBdense ;
    ssm_free (P->b) ;
    /* destroy Lanczos/Householder structure */
    ssm_free (P->LH->V) ;
    ssm_free (P->LH->d) ;
    ssm_free (P->LH->u) ;
    ssm_free (P->LH) ;

    /* destroy the diagonalization structure for tridiagonal matrix */
    ssm_free (P->DT->e) ;
    ssm_free (P->DT->gx) ;
    ssm_free (P->DT->gy) ;
    ssm_free (P->DT->gz) ;
    ssm_free (P->DT->gs) ;
    ssm_free (P->DT->gi) ;
    ssm_free (P->DT->gj) ;
    ssm_free (P->DT->gk) ;
    ssm_free (P->DT->gf) ;
    ssm_free (P->DT) ;

    /* destroy the diagonal matrix optimization structure */
    ssm_free (P->DO->d) ;
    ssm_free (P->DO->f) ;
    ssm_free (P->DO->f2) ;
    ssm_free (P->DO->dshift) ;
    ssm_free (P->DO) ;
    ssm_free (P) ;

    ssm_free (Com->MINRES) ;
    ssm_free (Com->W) ;
    ssm_free (Com->V) ;
    ssm_free (Com->VAV) ;
    if ( Com->SSOR != NULL ) ssm_free (Com->SSOR) ;
}

/* ==========================================================================
   === SSMinitProb ==========================================================
   ==========================================================================
    Copy matrix, order column, extract diagonal, locate first element
    beneath diagonal
   ========================================================================== */
int SSMinitProb
(
    SSMFLOAT   *Ax, /* numerical entries in matrix, packed array */
    SSMINT     *Ai, /* row indices for each column */
    SSMINT     *Ap, /* points into Ax and Ai, Ap [0], ... Ap [n], packed */
    SSMINT       n, /* problem dimension */
    SSMFLOAT    *b, /* linear term in objective function */
    SSMFLOAT   rad, /* radius of sphere */
    SSMProblem *PB, /* problem structure */
    SSMParm  *Parm, /* parameter structure */
    SSMcom    *Com  /* SSMcom structure */
)
{
    SSMINT k, l, nnz, p, *Ap1, *Bp, *Bp1, *Bu ;
    SSMINT ai, j, *Bi, *wi ;
    SSMFLOAT ax, Amax, Dmin, dj, *Bx, *D ;
    SSMLanczos *LH ;
    SSMDiag *DT ;
    SSMDiagOpt *DO ;
    int status = 0 ;

    PB->n = n ;
    PB->rad = rad ;
    PB->b = b ;
    Com->minres_limit = (int) (((SSMFLOAT) n)*Parm->minres_its_fac) ;
    Com->ssm_limit = (int) (((SSMFLOAT) n)*Parm->ssm_its_fac) ;

    /* flag nonzero diagonal elements in matrix, diagonal stored in
       separate array PB->D */
    wi = Com->wi1 ;
    Ap1 = Ap+1 ;
    k = 0 ;
    for (j = 0; j < n; j++)
    {
        l = Ap1 [j] ;
        wi [j] = 0 ;
        for (; k < l; k++)
        {
            if ( Ai [k] == j )
            {
                wi [j] = 1 ;
                k = l ;
                break ;
            }
        }
    }

    Bp = PB->p = sopt_malloc (&status, (n+1), sizeof (SSMINT)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;
    Bp1 = Bp+1 ;
    p = 0 ;
    nnz = 0 ;
    for (j = 0; j < n; j++)
    {
        /* number of nonzero diagonal elements up to end of column j */
        p += wi [j] ;
        Bp [j] = nnz ;
        /* number nonzeros up to end of column j excluding diagonal elements*/
        nnz = Ap1 [j] - p ;
    }
    Bp [n] = nnz ;

    /* allocate memory for matrix storage */
    D  = PB->D = sopt_malloc (&status, n, sizeof (SSMFLOAT)) ;
    Bu = PB->u = sopt_malloc (&status, n, sizeof (SSMINT)) ;
    Bi = PB->i = sopt_malloc (&status, nnz, sizeof (SSMINT)) ;
    Bx = PB->x = sopt_malloc (&status, nnz, sizeof (SSMFLOAT)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;

    /* ======================================================================
       === B = A'============================================================
       ======================================================================
          - this transpose operation ensures that columns of matrix are sorted
          - store diagonal in PB->D
          - compute absolute largest element of A
       ====================================================================== */

    Amax = SSMZERO ;
    Dmin = SSMINF ;
    k = 0 ;
    for (j = 0; j < n; j++)
    {
        l = Ap1 [j] ;
        dj = SSMZERO ;
        for (; k < l; k++)
        {
            ai = Ai [k] ;
            ax = Ax [k] ;
            if ( fabs (ax) > Amax ) Amax = fabs (ax) ;
            if ( ai == j )  /* store the diagonal in D, not  Bx */
            {
                dj = ax ;
            }
            else
            {
                p = Bp [ai]++ ;
                Bi [p] = j ;
                Bx [p] = ax ;
            }
        }
        D [j] = dj ;
        if ( dj < Dmin ) Dmin = dj ;
    }
    PB->Amax = Amax ;
    PB->Dmin = Dmin ;
    for (j = n; j > 0; j--) Bp [j] = Bp [j-1] ;
    Bp [0] = 0 ;

    /* Bu points just past last nonzero above diagonal in each column */
    Bu = PB->u ;
    k = 0 ;
    for (j = 0; j < n; j++)
    {
        l = Bp1 [j] ;
        Bu [j] = l ;
        for (; k < l; k++)
        {
            ai = Bi [k] ;
            if ( Bi [k] > j )
            {
                Bu [j] = k ;
                k = l ;
                break ;
            }
        }
    }

    /* store Parm and PB in Com */
    Com->Parm = Parm ;
    Com->PB = PB ;

    /* store Lanczos parameters in LH structure */
    LH = PB->LH ;
    LH->tol = Amax*Parm->Lanczos_tol ;   /*u_j termination tol in Lanczos iter*/
    LH->House_limit = Parm->House_limit ;
    LH->nrows = n ;

    /* store diagonalization and tridiagonalization parameters in DT structure*/
    DT = PB->DT ;
    LH->max_its = DT->nalloc ;
    DT->tol1 = Parm->qr_tol_rel ;        /* relative tolerance for QR method */
    DT->tol2 = Parm->qr_tol_abs ;        /* absolute tolerance for QR method */
    DT->max_its = SSMMAX (30, Parm->qr_its*DT->nalloc) ;

    /* store diagonal optimization parameter in DO structure*/
    DO = PB->DO ;
    DO->tol = Parm->diag_optim_tol ;

    /* now store parameters for the 5 by 5 matrices arising in refinement */
    DT = Com->PBdense->DT ;
    DT->tol1 = Parm->qr_tol_rel ;      /* relative tolerance for QR method*/
    DT->tol2 = Parm->qr_tol_abs ;      /* absolute tolerance for QR method*/
    DT->max_its = SSMMAX (30, Parm->qr_its*5) ;

    DO = Com->PBdense->DO ;
    DO->tol = Parm->diag_optim_tol ;
    return (sopt_convert_error (LSSM, status)) ;
}

/* ==========================================================================
   === SSMballdense ==========================================================
   ==========================================================================

   Solve a dense sphere constrained quadratic program of the form

       min  x'Ax + 2b'x  subject to ||x|| = r
   ========================================================================== */

int SSMballdense /*return 0 (error tolerance satisfied)
                                 -5 (failure of QR diagonalization)
                                 -6 (insufficient space in QR diagonalization)*/
(
    SSMFLOAT     *x, /* n-by-1 solution vector (output) */
    SSMProblem  *PB, /* problem structure associated with A */
    SSMcom     *Com /* SSMcom structure */
)
{
    int status ;
    SSMINT j, n ;
    SSMFLOAT *b, *A ;
    SSMDiag *DT ;
    SSMDiagOpt *DO ;
    SSMLanczos *LH ;

    DT = PB->DT ;
    DO = PB->DO ;
    LH = PB->LH ;
    b = PB->b ;
    A = LH->V ;
    n = PB->n ;
    SSMtriHouse (A, LH->d, LH->u, Com->wx1, n) ;
    LH->ncols = n ;
    status = SSMdiag (LH->d, LH->u, n, DT,
                      Com->wi1, Com->wx1, Com->wx2, Com->wx3) ;

    /* copy the eigenvalues and the dimension into the SSMDiagOpt structure*/
    DO->n = n ;
    for (j = 0; j < n; j++) DO->d [j] = DT->e [j] ;
    SSMtDenseMult (DO->f, b, LH->V, n, n) ;
    SSMGivensMult (DO->f, DT) ;
    SSMdiagopt (Com->wx3, Com->PB->rad, Com->Parm->BndLessThan, DO, Com) ;
    Com->Active = SSMTRUE ;
    if ( Com->Parm->BndLessThan && (DO->mu <= SSMZERO) ) Com->Active = SSMFALSE;
    SSMtGivensMult (Com->wx3, DT) ;
    SSMDenseMult (x, Com->wx3, LH->V, n, n) ;
    return (status) ;
}

/* ==========================================================================
   === SSMsubspace ==========================================================
   ==========================================================================
   Minimize the quadratic over a subspace spanned by m vectors, m <= 5.
   Since ||v1|| and the product A*v1 have generally been computed already,
   they are provided as input parameters to avoid their recomputation.
   NOTE: if only_residual is TRUE, it is assumed that v1 contains x
   ========================================================================== */
int SSMsubspace /* return 0 (error tolerance satisfied)
                                -5 (failure of QR diagonalization)
                                -6 (insufficient space in QR diagonalization)*/
(
    SSMFLOAT      *v1, /* subspace vector 1, solution is returned in v1 */
    SSMFLOAT   normv1, /* norm of vector 1 */
    SSMFLOAT     *Av1, /* A*v1 */
    SSMFLOAT      *v2, /* subspace vector 2 */
    SSMFLOAT      *v3, /* subspace vector 3 */
    SSMFLOAT      *v4, /* subspace vector 4 */
    SSMFLOAT      *v5, /* subspace vector 5 */
    int             m, /* number of vectors */
    int only_residual, /* only compute kkt error and eigenvector residual */
    SSMcom       *Com
)
{
    int BndLessThan, status ;
    SSMINT i, j, n ;
    SSMFLOAT rad, r, s, t, u, X [5], evec [5], mu, ball_err, normx,
          *V, *W, *AV, *VAV, *VAVsave, *b, *v, *Av, *Ax, *r0, *emin, *x ;
    SSMProblem *PB, *P ;
#ifndef NDEBUG
    SSMFLOAT cost, VAVcompare [25], kkt_res [5] ;
    if ( (m > 5) || (m < 1) )
    {
        DPRINT (("in SSMsubpace, the dimension %i outside range [1, 5]\n", m)) ;
        SOPTERROR ("stop") ;
    }
#endif

    PB = Com->PB ;
    b = PB->b ;
    n = PB->n ;
    rad = PB->rad ;
    P = Com->PBdense ;
    P->n = m ;
    V = Com->V ;        /* storage for orthogonal vectors */
    W = Com->W ;        /* storage for Householder vectors */
    AV= Com->MINRES ;   /* storage for A*V */
    VAV= P->LH->V ;     /* storage for V'AV */
    VAVsave = Com->VAV ;/* save copy of V'AV in case eigenvector refined */
    Ax = Com->Ax ;
    Av = Com->Av ;
    v  = Com->v ;
    r0 = Com->r0 ;
    emin = &Com->emin ;
    BndLessThan = Com->Parm->BndLessThan ;
    x = v1 ;            /* solution is return in array x = v1 */
    if ( only_residual == SSMTRUE )
    {
        /* if constraint active or ||x|| > rad, then normalize x
           NOTE: x = v1 */
        normx = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = x [j] ;
            normx += t*t ;
        }
        normx = sqrt (normx) ;
        if ( (PB->DO->mu != SSMZERO) || !BndLessThan || (normx > rad) )
        {
            if ( normx != SSMZERO )
            {
                t = rad/normx ;
                for (j = 0; j < n; j++) x [j] *= t ;
            }
        }
        SSMmult (Ax, x, PB->x, PB->D, PB->i, PB->p, n, Com) ;
        SSMmult (Av, v, PB->x, PB->D, PB->i, PB->p, n, Com) ;
        P->DO->mu = PB->DO->mu ;
        goto Residual_computation ;
    }
    if ( m < 5 )
    {

        /* compute orthonormal basis for given basis vectors */
        t = SSMONE/sqrt (normv1*(normv1 + fabs (v1 [0]))) ;
        if ( v1 [0] > SSMZERO ) W [0] = t*(v1 [0] + normv1) ;
        else                    W [0] = t*(v1 [0] - normv1) ;
        for (j = 1; j < n; j++) W [j] = v1 [j]*t ; /* H = I - w*w' */
        if ( m >= 2 ) SSMorth (V+(1*n), W, v2, 1, n) ;
        if ( m >= 3 ) SSMorth (V+(2*n), W, v3, 2, n) ;
        if ( m >= 4 ) SSMorth (V+(3*n), W, v4, 3, n) ;

        /* normalize v1 and A*v1, they are first column of V and AV */
        t = SSMONE/normv1 ;
        if ( t != SSMONE )
        {
            for (j = 0; j < n; j++) AV [j] = t*Av1 [j] ;
            for (j = 0; j < n; j++)  V [j] = t* v1 [j] ;
        }
        else
        {
            for (j = 0; j < n; j++) AV [j] =   Av1 [j] ;
            for (j = 0; j < n; j++)  V [j] =    v1 [j] ;
        }

        /* multiply remaining columns of V by A */
        for (i = 1; i < m; i++)
            SSMmult (AV+(i*n), V+(i*n), PB->x, PB->D, PB->i, PB->p, n, Com) ;
    
        /* compute elements of V'AV on or below diagonal */
        for (i = 0; i < m; i++)
        {
            SSMtDenseMult (VAV+(i*m + i), V+(i*n), AV+(i*n), n, m-i) ;
        }
    
        /* elements above diagonal known from symmetry */
        for (j = 0; j < m-1; j++)
        {
            for (i = j+1; i < m; i++)
            {
                /* (i,j) element inserted in location (j,i) above diagonal */
                VAV [m*i+j] = VAV [m*j+i] ;
            }
        }
#ifndef NDEBUG
        /* compare with directly computed product */
        for (i = 0; i < m; i++)
        {
            SSMtDenseMult (VAVcompare+(i*m), V+(i*n), AV, n, m) ;
        }
        t = SSMZERO ;
        for (i = 0; i < m*m; i++) t += fabs (VAVcompare [i] - VAV [i]) ;
        if ( t/PB->Amax > Com->Parm->check_kkt )
        {
            DPRINT (("PV product in refine/restart, error: %e\n", t)) ;
            SOPTERROR ("stop") ;
        }
#endif
    
        /* save VAV in case we also compute a refined eigenvector */
        for (i = 0; i < m*m; i++) VAVsave [i] = VAV [i] ;
    
        /* compute linear term */
        SSMtDenseMult (P->b, b, V, n, m) ; /* linear term P->b = V'b */
    
        /* subproblem dimension is m */
        P->n = m ;
    
        /* X = solution of dense QP */
        status = SSMballdense (X, P, Com) ;
        if ( status ) return (status) ;
#ifndef NDEBUG
        /* check solution of dense subproblem */
        t = SSMZERO ;
        for (i = 0; i < m; i++) t += X [i]*X [i] ;
        t = sqrt (t) ;
        /* check radius and complementary slackness */
        if ( BndLessThan )
        {
            if ( (t - rad)/rad > Com->Parm->check_kkt )
            {
                DPRINT (("||x|| = %e > rad: %e in SSMsubpace\n", t, rad)) ;
                SOPTERROR ("stop") ;
            }
            mu = P->DO->mu ;
            if ( mu < SSMZERO )
            {
                DPRINT (("multiplier (%e) has wrong sign in SSMsubspace\n",
                          mu)) ;
                SOPTERROR ("stop") ;
            }
            s = SSMMIN (rad - t, mu) ;
            if ( s/rad > Com->Parm->check_kkt )
            {
                DPRINT (("Complementary slackness violated in SSMsubspace\n"));
                DPRINT (("||x||: %e radius: %e mu: %e\n", t, rad, mu)) ;
                SOPTERROR ("stop") ;
            }
        }
        else
        {
            if ( fabs (t-rad)/rad > Com->Parm->check_kkt )
            {
                DPRINT (("||x|| = %e != rad: %e in SSMsubpace\n", t, rad)) ;
                SOPTERROR ("stop") ;
            }
        }

        /* check kkt error */
        SSMtDenseMult (kkt_res, X, VAVcompare, m, m) ;
        s = SSMZERO ;
        for (i = 0; i < m; i++)
        {
            s = SSMMAX (s, fabs (kkt_res [i] + P->b [i] + mu*X [i])) ;
        }
        if ( s > PB->Amax*Com->Parm->check_kkt )
        {
            DPRINT (("kkt error %e exceeds check tolerance in SSMsubspace\n",
                      s)) ;
            SOPTERROR ("stop") ;
        }
#endif
        /*x = VX, solution x is returned in array v1 (x = v1)*/
        SSMDenseMult (x, X, V, n, m) ;
        /* compute Ax */
        SSMmult (Ax, x, PB->x, PB->D, PB->i, PB->p, n, Com) ;
    
        /* smallest eigenpair of V'AV */
        SSMmineig (emin, evec, Com->wx1, P) ;
        SSMDenseMult (v, evec, V, n, m) ;/* v = V*evec */
        SSMmult (Av, v, PB->x, PB->D, PB->i, PB->p, n, Com) ; /* Av */
    }
    else /* special case where the eigenvector was refined */
    {
        /* vector v5 used to generate 5th orthonormal column of V */
        SSMorth (V+(4*n), W, v5, 4, n) ;
#ifndef NDEBUG
        /* check orthogonality of columns of V */
        for (j = 0; j < 5; j++)
        {
             SSMINT k ;
             for (i = 0; i < 5; i++)
             {
                 t = SSMZERO ;
                 for (k = 0; k < n; k++) t += V [k+j*n]*V [k+i*n] ;
                 if ( i != j )
                 {
                     if ( fabs (t) > Com->Parm->check_tol*PB->Amax )
                     {
                         DPRINT (("in refine,columns of V not orth\n"));
                         DPRINT (("i: %i j: %i err: %e\n", i, j, t)) ;
                     }
                 }
                 else
                 {
                     if ( fabs (t-SSMONE) > Com->Parm->check_tol*PB->Amax )
                     {
                         DPRINT (("in refine, cols of V not normal\n"));
                         DPRINT (("i: %i norm: %e\n", i, t)) ;
                     }
                 }
             }
        }
#endif
        /* multiply 5th column of V by A to obtain 5th column of AV */
        SSMmult (AV+(4*n), V+(4*n), PB->x, PB->D, PB->i, PB->p, n, Com) ;
        /* scatter 4 by 4 matrix in VAVsave into 5 by 5 matrix VAV */
        for (j = 3; j >= 0; j--)
        {
            for (i = 3; i >= 0; i--) VAV [5*j+i] = VAVsave [4*j+i] ;
        }
        /* compute 5th column of 5 by 5 matrix of VAV */
        SSMtDenseMult (VAV+20, AV+(4*n), V, n, 5) ;
    
        /* use symmetry to generate 5th row of VAV */
        VAV [4]  = VAV [20] ; VAV [9]  = VAV [21] ;
        VAV [14] = VAV [22] ; VAV [19] = VAV [23] ;

#ifndef NDEBUG
        /* compare with directly computed product. Since the first
           4 columns of AV were wiped out, reconstruct them */
        for (i = 0; i < 4; i++)
            SSMmult (AV+(i*n), V+(i*n), PB->x, PB->D, PB->i, PB->p, n, Com) ;
        for (i = 0; i < m; i++)
        {
            SSMtDenseMult (VAVcompare+(i*m), V+(i*n), AV, n, m) ;
        }
        t = SSMZERO ;
        for (i = 0; i < m*m; i++) t += fabs (VAVcompare [i] - VAV [i]) ;
        if ( t/PB->Amax > Com->Parm->check_kkt )
        {
            DPRINT (("PV product in refine/restart, error: %e\n", t)) ;
            SOPTERROR ("stop") ;
        }
#endif
        /* due to the 5th column of V, there is a 5th component of P->b */
        SSMtDenseMult ((P->b)+4, b, V+(4*n), n, 1) ;
        P->n = 5 ;
        /* X = solution of dense QP */
        status = SSMballdense (X, P, Com) ;
        if ( status ) return (status) ;
#ifndef NDEBUG
        /* check solution of dense subproblem */
        t = SSMZERO ;
        for (i = 0; i < m; i++) t += X [i]*X [i] ;
        t = sqrt (t) ;
        /* check radius and complementary slackness */
        if ( BndLessThan )
        {
            if ( (t - rad)/rad > Com->Parm->check_kkt )
            {
                DPRINT (("||x|| = %e > rad: %e in SSMsubpace\n", t, rad)) ;
                SOPTERROR ("stop") ;
            }
            mu = P->DO->mu ;
            if ( mu < SSMZERO )
            {
                DPRINT (("multiplier (%e) has wrong sign in SSMsubspace\n",
                          mu)) ;
                SOPTERROR ("stop") ;
            }
            s = SSMMIN (rad - t, mu) ;
            if ( s/rad > Com->Parm->check_kkt )
            {
                DPRINT (("Complementary slackness violated in SSMsubspace\n"));
                DPRINT (("||x||: %e radius: %e mu: %e\n", t, rad, mu)) ;
                SOPTERROR ("stop") ;
            }
        }
        else
        {
            if ( fabs (t-rad)/rad > Com->Parm->check_kkt )
            {
                DPRINT (("||x|| = %e != rad: %e in SSMsubpace\n", t, rad)) ;
                SOPTERROR ("stop") ;
            }
        }

        /* check kkt error */
        SSMtDenseMult (kkt_res, X, VAVcompare, m, m) ;
        s = SSMZERO ;
        for (i = 0; i < m; i++)
        {
            s = SSMMAX (s, fabs (kkt_res [i] + P->b [i] + mu*X [i])) ;
        }
        if ( s > PB->Amax*Com->Parm->check_kkt )
        {
            DPRINT (("kkt error %e except check tolerance in SSMsubspace\n",
                      s)) ;
            SOPTERROR ("stop") ;
        }
#endif

        /* x = VX, return solution in v1 (x = v1) */
        SSMDenseMult (x, X, V, n, 5) ;
        /* compute Ax */
        SSMmult (Ax, x, PB->x, PB->D, PB->i, PB->p, n, Com) ;
    
        /* smallest eigenpair of V'AV */
        SSMmineig (emin, evec, Com->wx1, P) ;
        /* v = V*evec */
        SSMDenseMult (v, evec, V, n, m) ;
        SSMmult (Com->Av, v, PB->x, PB->D, PB->i, PB->p, n, Com) ; /* Av */
    }

    Residual_computation:
#ifndef NDEBUG
    /* check for cost function decay */
    cost = SSMZERO ;
    for (j = 0; j < n; j++) cost += b [j]*x [j] ;
    for (j = 0; j < n; j++) cost += 0.5*Ax [j]*x [j] ;
    DPRINT (("cost: %30.15e cost_old: %30.15e\n", cost, Com->cost_old)) ;
    if ( cost > Com->cost_old + Com->Parm->check_tol*PB->Amax )
    {
        DPRINT (("cost does not decay in SSM\n")) ;
        DPRINT (("old: %30.15e new: %30.15e\n", Com->cost_old, cost)) ;
        SOPTERROR ("stop") ;
    }
    Com->cost_old = cost ;
    if ( *emin > Com->emin_old + Com->Parm->check_tol*PB->Amax )
    {
        DPRINT (("eigen estimate does not decay\n")) ;
        DPRINT (("old: %30.15e new: %30.15e\n", Com->emin_old, *emin)) ;
        SOPTERROR ("stop") ;
    }
    Com->emin_old = *emin ;
#endif

    /* compute residual r0 and estimate multiplier for the
       original problem if subproblem mu != 0 */
    if ( (P->DO->mu != SSMZERO) || !BndLessThan )
    {
        mu = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = b [j] + Ax [j] ;
            mu -= t*x [j] ;
            r0 [j] = t ;
        }
        mu /= rad*rad ;  /* mu = r0'*x/rad^2 */
        Com->mu = mu ;

        /* estimate error */
        ball_err = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = r0 [j] + x [j]*mu ;
            ball_err += t*t ;
        }
        ball_err = sqrt (ball_err) ;
        if ( (mu < SSMZERO) && BndLessThan )
        {
            ball_err = SSMMAX (ball_err, fabs (mu)*rad) ;
        }
        Com->error = ball_err ;
    }
    else
    {
        /* compute residual */
        Com->mu = SSMZERO ;
        ball_err = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = b [j] + Ax [j] ;
            r0 [j] = t ;
            ball_err += t*t ;
        }
        ball_err = sqrt (ball_err) ;
        Com->error = ball_err ;
    }

    /* compute eigenvector residual */
    t = SSMZERO ;
    u = SSMZERO ;
    for (j = 0; j < n; j++)
    {
        r = v [j] ;
        s = Av [j] - *emin*r ;    /* residual */
        t += s*s ;               /* norm of residual */
        u += r*r ;               /* norm of v */
    }
    Com->eig_error = sqrt (t/u) ;/* norm (residual)/norm (eigenvector) */
    return (0) ;
}

/* ==========================================================================
   === SSMorth ==============================================================
   ==========================================================================
    One step in the Householder orthogonalization of a collection of vectors
    or equivalently, one step in the QR factorization of the matrix whose
    columns are the vectors to be orthogonalized.  W stores the Householder
    vectors h0, h1, ... , hk-1 in its columns. The columns are packed
    to squeeze out the zeros.  Given a new vector x,
    the goal is to generate the next orthonormal vector. As in the QR
    factorization of a matrix, we first compute
    y = (I - hk-1 hk-1') ... (I - h1 h1')(I - h0 h0')x
    Then we evaluate the Householder vector hk which annihilates components
    k+1 through n-1 in y. Finally, we multiply the product
    (I - h0 h0') ... (I - hk hk') by the k-th column of the identity to
    obtain the desired orthonormal vector w
   ========================================================================== */

void SSMorth
(
    SSMFLOAT  *w, /* k-th orthonormal vector */
    SSMFLOAT  *W, /* packed matrix storing Householder vectors */
    SSMFLOAT  *x, /* k-th new vector */
    SSMINT     k,
    SSMINT     n  /* dimension of x */
)
{
    SSMINT i, j, jp1, K ;
    SSMFLOAT xk, s, t, u, *Wp ;
    if ( k > n-1 )
    {
        DPRINT (("in SSMorth, k = %i > n-1 = %i "
                "(the orthogonal vector is zero)\n", k, n-1)) ;
        SOPTERROR ("stop") ;
    }
    K = k ;
    if ( k == n-1 )
    {
        Wp = W+((n*(n-1) - 2)/2) ;
        t = Wp [k] ;
        w [k] = SSMONE - t*t ;
        K-- ;
        w [K] = -Wp [K]*t ;
        goto compute_w ;
    }

    /* compute y = (I - hk-1 hk-1') ... (I - h1 h1')(I - h0 h0')x, store in w */
    for (j = 0; j < n; j++) w [j] = x [j] ;
    Wp = W ;
    for (j = 0; j < k;)
    {
        s = SSMZERO ;
        for (i = j; i < n; i++) s += Wp [i]*w [i] ;
        for (i = j; i < n; i++) w [i] -= Wp [i]*s ;
        j++ ;
        Wp += n - j ;  /* advance to next column, columns are packed */
    }

    u = SSMZERO ;
    for (j = k; j < n; j++) u = SSMMAX (u, fabs (w [j])) ;
    if ( u == SSMZERO )
    {
        for (j = k; j < n; j++) Wp [j] = SSMZERO ;
        goto init_w ;
    }
    s = SSMZERO ;
    for (j = k; j < n; j++)
    {
        t = w [j]/u ;
        Wp [j] = t ;
        s += t*t ;
    }
    xk = Wp [k] ;
    s = sqrt (s) ;
    if ( xk >= SSMZERO ) Wp [k] = xk + s ;
    else              Wp [k] = xk - s ;
    t = SSMONE/sqrt (s*(s + fabs (xk))) ;
    for (j = k; j < n; j++) Wp [j] *= t ;

    /* compute column k of the product (I - h0 h0') ... (I - hk hk') */
    init_w:
    t = Wp [k] ;
    for (i = k; i < n; i++) w [i] = -t*Wp [i] ;
    w [k] += SSMONE ;

    compute_w:
    for (j = K-1; j >= 0; j--)
    {
        jp1 = j + 1 ;
        Wp -= n - jp1 ;
        s = SSMZERO ;
        for (i = jp1; i < n; i++) s += w [i]*Wp [i] ;
        w [j] = -s*Wp [j] ;
        for (i = jp1; i < n; i++) w [i] -= s*Wp [i] ;
    }
}

/* ==========================================================================
   === SSMmineig ============================================================
   ==========================================================================
    Estimate a smallest eigenvalue and associated eigenvector for the
    matrix associated with PB. It is assumed that the diagonal optimization
    structure, the diagonlization, and the Lanczos/Householder
    tridiagonalization of the matrix have already been performed and
    stored in PB->DO, PB->DT, PB->LH
   ==========================================================================*/

void SSMmineig
(
    SSMFLOAT *emin, /* estimated smallest eigenvalue */
    SSMFLOAT    *v, /* associated eigenvector */
    SSMFLOAT    *w, /* work array of size n */
    SSMProblem *PB  /* problem specification */
)
{
    SSMINT imin, j, n ;
    SSMDiagOpt *DO ;
    SSMDiag *DT ;
    SSMLanczos *LH ;

    DO = PB->DO ;
    DT = PB->DT ;
    LH = PB->LH ;
    n = PB->n ;
    *emin = DO->dmin ;  /* smallest eigenvalue */
    imin  = DO->imin ;  /* index of smallest eigenvalue */
    for (j = 0; j < n; j++) w [j] = SSMZERO ;
    w [imin] = SSMONE ;

    /* Multiply by rotations ... G_2' G_1' diagonalizing tridiagonal matrix*/
    SSMtGivensMult (w, DT) ;

    /* Multiply by V */
    SSMDenseMult (v, w, LH->V, n, LH->ncols) ; /* v = estimated eigenvector */
    return ;
}

/* ========================================================================== */
/* === SSMdiag ============================================================== */
/* ========================================================================== */
/*
        Input: An n-dimensional tridiagonal matrix stored in arrays d and u

       Output: The eigenvalues along with the Givens rotations needed to
               diagonalize the matrix

    Algorithm: Inverse shifted QR algorithm with fast Givens rotations.
               The number of Givens rotations needed for the diagonalization
               is related to the number of QR iterations. This algorithm
               is an explicit implicit version of the QR method. It is
               explicit since the rotations are computed as in a
               QR factorization. It is implicit since each rotation
               is applied to both the right and left side of the matrix
               simultaneously, as is done with an implicit implementation
               of the QR method. An advantage of this explicit implicit
               algorithm is that the loss of shift phenomenon associated
               with the implicit QR algorithm does not occur since the
               rotations are computed from the explicit QR factorization.
               The speed associated with the implicit QR method is
               retained since each iteration entails a single pass over
               the matrix from upper left to lower right corner. In
               the fast Givens implementation, the matrix is stored as
               a product DTD where D is diagonal and T is tridiagonal.
               In the QR algorithm, the Givens rotations are applied to
               DTD - sigma I = D(T - sigma/D^2)D. In the explicit QR
               algorithm, the partially factored matrix has the form:

                         x  x  f  0  0  0  0
                         0  x  x  f  0  0  0
                         0  0  p  x  0  0  0
                         0  0  q  d  u  0  0
                         0  0  0  u  d  u  0
                         0  0  0  0  u  d  u

               Here f denotes fill, d and u are the original diagonal,
               subdiagonal, and superdiagonal elements, and p and q
               are the elements on the diagonal and subdiagonal which
               are used to generate the next Given rotation G. The
               rotation is designed to annihilate q while replacing
               p by sqrt(p*p + q*q). In performing the left and right
               multiplications by G, one needs to recall the structure
               of the matrix that is operated on:

                         d  u  0  0  0
                         u  d  u  a  0
                         0  u  d  u  0
                         0  a  u  d  u
                         0  0  0  u  d

               The matrix is tridiagonal with a bulge corresponding to
               the element a. In the code below, the iteration starts
               below the comment "perform one pass". The variables p, q,
               and a inside the subsequent code are identical to the
               p, q, and a above.

               An advantage of the fast Givens rotation is that several
               multiplications associated with the usual Givens rotation
               are eliminated as well as the square root. For details
               concerning fast Givens, see W. W. Hager, Applied Numerical
               Linear Algebra, pages 221-223. In the code that follows, the
               squares of the diagonal of the diagonal scaling matrix D are
               stored in the array sc (scaling array). As the algorithm
               progresses, the elements of sc tend to 0 since they are being
               multiplied by factors between 1/2 and 1.  When sc[i] drops
               below qr_lower, we multiply row i of T by sc[i] and reset
               sc[i] = 1. The code below rectifies an unflow problem
               which could occur in the analogous NAPACK code tdg.f.
               In tdg.f, we check sc at the start of each qr iteration.
               However, for a matrix with 5000 or more rows, an underflow
               could potentially occur within a QR iteration. Hence, in
               the new code, we check for small scales within the iteration.
               The code also checks to see if the matrix can be split due
               to a small off-diagonal element.
   ========================================================================== */
int SSMdiag /* return: 0 (convergence tolerance satisfied)
                              -5 (algorithm did not converge)
                              -6 (insufficient space allocated for Givens) */
(
    SSMFLOAT *din, /* diagonal of tridiagonal matrix */
    SSMFLOAT *uin, /* superdiagonal of tridiagonal matrix, use uin (0:n-2) */
    SSMINT      n, /* dimensional of tridiagonal matrix */
    SSMDiag   *DT, /* structure for storing diagonalization */
    SSMINT   *wi1, /* work array of size n */
    SSMFLOAT *wx1, /* work array of size n */
    SSMFLOAT *wx2, /* work array of size n */
    SSMFLOAT *wx3  /* work array of size n */
)
{
    SSMFLOAT a, b, c, delta, lb, tol1, tol2, break_tol,
         *gx, *gy, *gs, *e, *d, *u, *sc ;
    SSMFLOAT p, q, r, s, scale, sigma, t, w, x, y, z0, z1 ;
    SSMINT *gj ;
    SSMINT i, ip1, im1, k, l, l1, nm1, istart, nstart,
       *gi, *gk, *gf, *start_index ;
    short *gz ;
    SSMINT ngivens,     /* number of fast Givens */
        nscale ;        /* number of scalings performed */
    int npass,          /* number of passes */
        max_npass ;

    /* output, not defined on input: */
    DT->n = n ;
    e  = DT->e,        /* eigenvalues */
    gx = DT->gx ;      /* fast Givens x factor, SSMFLOAT */
    gy = DT->gy ;      /* fast Givens y factor, SSMFLOAT */
    gz = DT->gz ;      /* rotation type 1 or 2, short */
    gs = DT->gs ;      /* scaling factors, SSMFLOAT */
    gi = DT->gi ;      /* row indices of scale factors, SSMINT */
    gj = DT->gj ;      /* indexed by QR pass, points into Givens factors, INT */
    gk = DT->gk ;      /* index of last diagonal element in QR iteration, INT */
    gf = DT->gf ;      /* index of first diagonal element in QR iteration */

    /* input, not modified: */
    nm1 = n - 1 ;
    tol1 = DT->tol1 ;  /* stop when u_{k-1}^2 <= tol1*|d_k| */
    tol2 = DT->tol2 ;  /* stop when |u_{k-1}| <= tol2 */
    max_npass = DT->max_its ;/* maximum number of QR algorithm passes */

    /* workspace vectors of size n */
    d = wx1 ;
    u = wx2 ;
    sc = wx3 ;
    start_index = wi1 ;
    nstart = 0 ;
    start_index [0] = 0 ;

    /* ---------------------------------------------------------------------- */

    DT->ngivens = 0 ;      /* number of fast Givens */
    DT->its = 0 ;          /* number of QR iterations */
    DT->nscale = 0 ;       /* number of scalings */

    /* ---------------------------------------------------------------------- */
    /* check for 1-by-1 case */
    /* ---------------------------------------------------------------------- */

    if (n == 1)
    {
        e [0] = din [0] ;
        return (0) ;
    }
    /* ---------------------------------------------------------------------- */
    /* determine the problem scaling */
    /* ---------------------------------------------------------------------- */
    scale = fabs (din [nm1]) ;
    for (i = 0 ; i < nm1 ; i++)
    {
        scale = SSMMAX (scale, fabs (din [i])) ;
        t = uin [i] ;
        scale = SSMMAX (scale, fabs (t)) ;
        if ( t == SSMZERO )
        {
            nstart++ ;
            start_index [nstart] = i + 1 ;
        }
    }
    if (scale == 0 )
    {
        for (i = 0 ; i < n ; i++) e[i] = 0 ;
        return (0) ;
    }
    scale = 1/scale ;
    /* ---------------------------------------------------------------------- */
    /* make a scaled copy of the input problem */
    /* ---------------------------------------------------------------------- */
    for (i = 0; i < nm1; i++)
    {
        d [i] = scale*din [i] ; /* scaled diagonal is needed in next QR phase */
        u [i] = scale*uin [i] ;
        sc[i] = 1. ;
    }
    t = scale*din [nm1] ;
    d [nm1] = t ; /* scaled diagonal is needed in next QR phase */
    sc[nm1] = 1. ;
    u [nm1] = 0 ;
    /* ---------------------------------------------------------------------- */
    /* inverse shifted QR method */
    /* ---------------------------------------------------------------------- */
    lb = DT->qr_lower ; /* threshold in fast Givens for rescaling of rows */
    break_tol = lb*tol2 ;
    ngivens = 0 ;
    npass = 0 ;
    nscale = 1 ;
    gi [0] = 0 ;
    k = n ;
    l = n-1 ;

    while (k > 1)
    {
        l1 = 0 ;
        k = l ;
        l-- ;
        while ( k == start_index [nstart] )
        {
            k = l ;
            l-- ;
            nstart-- ;
            if ( k == 0 ) break ;
        }
        if ( k == 0 ) break ;

        while ( 1 )
        {
            if ( l1 >= 30 ) return (-5) ;
            t = sqrt(sc [l]*sc [k])*fabs(u [l]) ;
            if ( t <= tol2 ) break ;
            if ( t <= tol1*fabs(sc [k]*d [k]) ) break ;
            /* -------------------------------------------------------------- */
            /* log the start of a pass */
            /* -------------------------------------------------------------- */
            l1++ ;
            /* check if enough space was allocated for Givens rotations,
               if not, need to increase Parm->qr_its */
            if (npass >= max_npass) return (-6) ;
            gk [npass] = k ;
            istart = gf [npass] = start_index [nstart] ;
            gj [npass] = ngivens ;

            npass++ ;
            /* -------------------------------------------------------------- */
            /* compute shift sigma */
            /* -------------------------------------------------------------- */
            /* note that a, b, c, t, and delta */
            /* are not used outside this scope */
            a = d [k]*sc [k] ;
            b = d [l]*sc [l] ;
            c = sqrt(sc [l]*sc [k])*fabs(u [l]) ;
            delta = .5*(b-a) ;
            if (delta >= 0)
            {
                if (delta <= c)
                {
                    t = delta/c ;
                    sigma = a - c/(t+sqrt(t*t+1)) ;
                }
                else
                {
                    t = c/delta ;
                    sigma = a - c*t/(1+sqrt(t*t+1)) ;
                }
            }
            else
            {
                if (-delta <= c)
                {
                    t = -delta/c ;
                    sigma = a + c/(t+sqrt(t*t+1)) ;
                }
                else
                {
                    t = -c/delta ;
                    sigma = a + c*t/(1+sqrt(t*t+1)) ;
                }
            }

            /* -------------------------------------------------------------- */
            /* compute the initial p and q */
            /* -------------------------------------------------------------- */

            p = d [istart] - sigma/sc [istart] ;
            q = u [istart] ;
            /* -------------------------------------------------------------- */
            /* perform one pass */
            /* -------------------------------------------------------------- */
            for (i = istart; i < k; i++)
            { 
                /* test for break up of matrix */
                if ( i > istart )
                {
                    if ( fabs (u [i-1])*z0*z1 < break_tol )
                    {
                        nstart++ ;
                        start_index [nstart] = i ;
                    }
                }
	        /* ---------------------------------------------------------- */
	        /* scale the problem, if necessary */
	        /* ---------------------------------------------------------- */
                im1 = i - 1 ;
                z0 = sc [i] ;
                if ( z0 < lb )
                {
                    gi [nscale] = (i+1) ;
                    d [i] *= z0 ;
                    z0 = sqrt (z0) ;
                    gs [nscale] = z0 ;
                    nscale++ ;
                    p *= z0 ;
                    u [i] *= z0 ;
                    if ( i > 0 ) u [im1] *= z0 ;
                    z0 = SSMONE ;
                    sc [i] = SSMONE ;
                }
                ip1 = i + 1 ;
                z1 = sc [ip1] ;
                if ( z1 < lb )
                {
                    gi [nscale] = -(i+1) ;
                    d [ip1] *= z1 ;
                    z1 = sqrt (z1) ;
                    gs [nscale] = z1 ;
                    nscale++ ;
                    q *= z1 ;
                    u [i] *= z1 ;
                    u [ip1] *= z1 ;
                    if ( i > 0 ) a *= z1 ;
                    z1 = SSMONE ;
                    sc [ip1] = SSMONE ;
                }

                if (fabs (p) > fabs (q))
                {
                    r = z1/z0 ;
                    s = q/p ;
                    t = r*s*s ;

                    if (t <= 1)
                    {
                        t = 1/(1+t) ;
                        sc [i] = z0*t ;
                        sc [ip1]=z1*t ;
                        y = s ;
                        x = r*s ;
                        gx [ngivens] = x ;
                        gy [ngivens] = y ;
                        gz [ngivens] = 1 ;
                        ngivens++ ;

                        if (i > 0) u [im1] = u [im1] + x*a ;
                        t = d [i] ;
                        p = u [i] ;
                        q = t + x*p ;
                        r = p - y*t ;
                        w = d [ip1] - y*p ;
                        s = p + x*d [ip1] ;
                        u [i] = s - y*q ;
                        d [i] = q + x*s ;
                        d [ip1] = w - r*y ;
                        q = u [ip1] ;
                        a = x*q ;
                        p = w - sigma/z1 ;
                    }
                    else
                    {

                        t = t/(1+t) ;
                        r = z0/z1 ;
                        s = p/q ;
                        sc [i] = z1*t ;
                        sc [ip1] = z0*t ;
                        y = s ;
                        x = r*s ;
                        gx [ngivens] = x ;
                        gy [ngivens] = y ;
                        gz [ngivens] = 2 ;
                        ngivens++ ;

                        if (i > 0) u [im1] = x*u [im1] + a ;
                        t = d [i] ;
                        p = u [i] ;
                        q = x*t+p ;
                        r = y*p-t ;
                        w = y*d [ip1]-p ;
                        s = x*p+d [ip1] ;
                        u [i] = y*s-q ;
                        d [i] = x*q+s ;
                        d [ip1] = y*w-r ;
                        a = u [ip1] ;
                        q = a ;
                        u [ip1] = y*a ;
                        p = w - y*sigma/z1 ;
                    }

                }
                else
                {
                    if (q == 0)
                    {
                        gx [ngivens] = 0 ;
                        gy [ngivens] = 0 ;
                        gz [ngivens] = 1 ;
                        ngivens++ ;

                        p = d [ip1] - sigma/z1 ;
                        q = u [ip1] ;
                    }
                    else
                    {
                        r = z0/z1 ;
                        s = p/q ;
                        t = r*s*s ;

                        if (t < 1)
                        {
                            t = 1/(1+t) ;
                            sc [i] = z1*t ;
                            sc [ip1] = z0*t ;
                            y = s ;
                            x = r*s ;
                            gx [ngivens] = x ;
                            gy [ngivens] = y ;
                            gz [ngivens] = 2 ;
                            ngivens++ ;

                            if (i > 0) u [im1] = x*u [im1] + a ;
                            t = d [i] ;
                            p = u [i] ;
                            q = x*t+p ;
                            r = y*p-t ;
                            w = y*d [ip1] - p ;
                            s = x*p + d [ip1] ;
                            u [i] = y*s - q ;
                            d [i] = x*q + s ;
                            d [ip1] = y*w - r ;
                            a = u [ip1] ;
                            q = a ;
                            u [ip1] = y*a ;
                            p = w - y*sigma/z1 ;
                        }
                        else
                        {
                            t = t/(1+t) ;
                            r = z1/z0 ;
                            s = q/p ;
                            sc [i] = z0*t ;
                            sc [ip1] = z1*t ;
                            y = s ;
                            x = r*s ;
                            gx [ngivens] = x ;
                            gy [ngivens] = y ;
                            gz [ngivens] = 1 ;
                            ngivens++ ;

                            if (i > 0) u [im1] = u [im1] + x*a ;
                            t = d [i] ;
                            p = u [i] ;
                            q = t + x*p ;
                            r = p - y*t ;
                            w = d [ip1] - y*p ;
                            s = p + x*d [ip1] ;
                            u [i] = s - y*q ;
                            d [i] = q + x*s ;
                            d [ip1] = w - r*y ;
                            q = u [ip1] ;
                            a = x*q ;
                            p = w - sigma/z1 ;
                        }
                    }
                }
            }
            gi [nscale] = 0 ;
            nscale++ ;
        }
    }

    /* ---------------------------------------------------------------------- */

    scale = 1/scale ;
    for (i = 0 ; i < n ; i++ )
    {
        e [i] = sc [i]*d [i]*scale ;
    }
    DT->nscale = nscale ;

    if ( ngivens > 0 )
    {
        for (i = 0 ; i < n ; i++)
        {
            gs [nscale] = sqrt (sc [i]) ;

            nscale++ ;
        }
    }

    DT->ngivens = ngivens ;
    DT->its = npass ;
    return (0) ;
}

#ifndef NDEBUG
/* ==========================================================================
   === SSMcheckKKT ==========================================================
   ==========================================================================
    Check that x satisfies KKT conditions
   ========================================================================== */

int SSMcheckKKT
(
    SSMFLOAT    *x,  /* n-by-1 solution vector (output) */
    SSMFLOAT     r,  /* radius of sphere */
    SSMProblem *PB,  /* problem specification */
    SSMcom    *Com   /* pointer to SSMcom structure */
)
{
    SSMINT i, n ;
    SSMFLOAT err, s, t, check_kkt, check_rad, mu, *p ;
    int status = 0 ;

    n = PB->n ;
    mu = Com->mu ;
    check_kkt = r*PB->Amax*Com->Parm->check_kkt*(SSMFLOAT) n ;
    check_rad = r*Com->Parm->check_kkt*(SSMFLOAT) n ;

    /* check ||x|| <= r */
    s = 0 ;
    for (i = 0; i < n; i++)
    {
        t = x [i] ;
        s += t*t ;
    }
    if ( Com->Parm->BndLessThan ) err = SSMMAX (sqrt (s) - r, SSMZERO)/r ;
    else                          err = fabs (sqrt (s) - r)/r ;
    DPRINT (("constraint error: %e  check_rad: %e\n", err, check_rad)) ;
    if ( err > check_rad )
    {
         DPRINT (("||x||-r: %e > %e\n", err, check_rad)) ;
         SOPTERROR ("stop") ;
    }

    /* check complementary slackness */
    if ( Com->Parm->BndLessThan )
    {
        err = SSMZERO ;
        if ( sqrt (s) < r ) err = SSMMIN (r-sqrt (s), fabs (mu)) ;
        else                err = SSMMAX (SSMZERO, -mu) ;
        err /= r ;
        DPRINT (("complementary slackness error:  %e check_rad: %e\n",
            err, check_rad)) ;
        if ( err > check_rad )
        {
             DPRINT (("complementary slackness err: %e  > %e\n",
                 err, check_rad)) ;
             SOPTERROR ("stop") ;
        }
    }

    /* check 1st-order optimality condition
       Ax + b + mu*x = 0 */
    
    /* allocate workspace [ */
    p = sopt_malloc (&status, n, sizeof (SSMFLOAT)) ;
    if ( status ) return (sopt_convert_error (LSSM, status)) ;

    /* compute (A+D)x */
    SSMmult (p, x, PB->x, PB->D, PB->i, PB->p, n, Com) ;

    /* compute Ax + b + mu*x */
    for (i = 0; i < n; i++) p [i] += PB->b [i] + mu*x [i] ;

    /* sup norm error */
    err = 0. ;
    for (i = 0; i < n; i++) err = SSMMAX (err, fabs (p [i])) ;
    DPRINT (("KKT error (%e) check_kkt: %e\n", err, check_kkt)) ;
    if ( err > check_kkt )
    {
         DPRINT (("KKT error (%e) > %e\n", err, check_kkt)) ;
         SOPTERROR ("stop") ;
    }

    /* free workspace ] */
    ssm_free (p) ;
    return (sopt_convert_error (LSSM, status)) ;
}
#endif
