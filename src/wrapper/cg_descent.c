      /* Code parts:
         CGSTART start of main cg_descent loop
         AUTO    Autodetect before start of Newton
         NEWTON  use_hessian = TRUE
         DCG     dCG = TRUE (use CG direction)
         #NCG    while loop for Newton CG
         CGTRUST trust region in Newton CG
         #NSS    start of SS Newton
         */
/* =========================================================================
   ======================== CG_DESCENT =====================================
   =========================================================================
       ________________________________________________________________
      |      A conjugate gradient method with guaranteed descent       |
      |                Based on cg_descent version 6.8                 |
      |                                                                |
      |           William W. Hager    and   Hongchao Zhang             |
      |             hager@ufl.edu         hozhang@math.lsu.edu         |
      |                   Department of Mathematics                    |
      |                     University of Florida                      |
      |                 Gainesville, Florida 32611 USA                 |
      |                      352-392-0281 x 244                        |
      |                                                                |
      |        Copyright by William W. Hager and Honchao Zhang         |
      |                       November 1, 2019                         |
      |                                                                |
      |          http://www.math.ufl.edu/~hager/papers/CG              |
      |                                                                |
      |  Disclaimer: The views expressed are those of the authors and  |
      |              do not reflect the official policy or position of |
      |              the Department of Defense or the U.S. Government. |
      |                                                                |
      |      Approved for Public Release, Distribution Unlimited       |
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
      |________________________________________________________________|
      Alternative licenses are also available.  Please contact William Hager
      for details.

      References:
      1. W. W. Hager and H. Zhang, A new conjugate gradient method
         with guaranteed descent and an efficient line search,
         SIAM Journal on Optimization, 16 (2005), 170-192.
      2. W. W. Hager and H. Zhang, Algorithm 851: CG_DESCENT,
         A conjugate gradient method with guaranteed descent,
         ACM Transactions on Mathematical Software, 32 (2006), 113-137.
      3. W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
         methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58.
      4. W. W. Hager and H. Zhang, Limited memory conjugate gradients,
         SIAM Journal on Optimization, 23 (2013), 2150-2168. */

#ifdef PASA
#include "pasa.h"
#else
#include "cg_descent.h"
#endif

int XXCG(descent)
(
    CGdata    *cgdata /* CG data structure */
#ifdef PASA
   ,PASAcom  *pasacom  /* common variables from pasa */
#endif
)
{
    int     AutoChoice, CG_window, Newton_window, dCG, dCG_old, find_max_iter,
            max_d_iter_bump, deriv_mode, IterQuad, QuadF, mem, memk, mlast, mp,
            mpp, spp, LBFGS, qrestart, SSsetupDone, smallq, trustfails,
            use_hessian, skipHessian, HessianSparsityFixed, sF, location ;
    int     status = 0 ;
    int    *CGNfail, CGNfail_int, *SSfail, SSfail_int ;
    CGINT   Cdim, d_iter, FirstIter, hnnz, i, j, k,
            IterRestart, maxit, ncol, nnzt, nrestart, nslow, slowlimit,
            start_it ;
    CGINT  *max_d_iter ;
    CGFLOAT *speed ;
    char    *speed_info ;
#ifdef PASA
    CGINT l, p ;
#endif
    CGFLOAT alpha, alphaold, alphabound, Amax, ApproxSwitchFactor, beta, Ck, Qk,
            delta2, denom, dnorm2, dphi, dphi0, dq, dQd, err, err_mark,
            f, fbest, fnew, gbest, Cmax, Hmax, maxstep, rtr, rtv, scale,
            ykyk, ykPyk, gkPyk, dkyk, QuadTrust, s, t, trial,
            *d, *dense, *g, *gH, *gnew, *gproj, *gnewproj,
            *Hwork, *Hd, *r, *rhs, *v, *Sk, *SkYk,
            *tau, *work, *Yk, *xnew, *Xp ;
    CGcom   *cgcom, cgcom_struc ;
    LONG    li, le ;

    /* method control variables */
    CGFLOAT maxerr, minerr, htic, method_tic, method_toc ;
    int     end_it, OldMethod, NewMethod ;
    int     use_hessian_old, RestartAfterHessian, RestartSetupDone, SetTimer ;

    /* limited memory variables */
    int     cg_ok, l1, l2,  memsq, memk_begin, mlast_sub, mp_begin, nCG, nSS,
            NewtonStep, nsub, spp1, SkFstart, SkFlast, Subspace,
            use_memory, Restart, InvariantSpace,
            IterSubRestart, FirstFull, SubSkip, SubCheck,
            StartSkip, StartCheck, DenseCol1, memk_is_mem, d0isg ;
    CGINT   iter_since_err_update, iter_check, nrestartsub ;
    CGFLOAT err_decay, Hessian_err_decay, gnorm2, gsubnorm2,  ratio,
            stgkeep, zeta, yty, ytg, t1, t2, t3, t4,
           *D, *Rk, *Re, *SkF, *stemp,
           *dsub, *gsub, *gsubtemp, *gkeep, *vsub, *wsub ;

    SOPT_Cmatrix *C ;
    void (*hessian) (SOPT_matrix *, CGFLOAT *, CGINT) ;

    cgcom = &cgcom_struc ;
    cgcom->cgdata = cgdata ;
    cgcom->Stat = cgdata->Stat ;
    cgcom->Parm = cgdata->Parm ;
    int SSinitDone = FALSE ;

#ifdef USE_HARWELL
    int sI ;
    CGINT  nrow = 0 ;
    sI = sizeof (CGINT) ;
    LBL_Data *lbl ;
    lbl = NULL ;
    cgcom->lbl = lbl ;
    CGINT *Cparent, *Cmember ;
    FILE *lblfile ;
    cholmod_sparse *Achol ;
    cholmod_common *cmm ;
    cgcom->cmm = NULL ;
#endif

#ifdef USE_MUMPS
    DMUMPS_STRUC_C *id ;
    CGINT nrow = 0 ;
#ifdef USE_HARWELL
    /* print an error message if both MUMPS and HARWELL are active */
    status = cg_wrapup (CG_MULTI_SOLVERS, TRUE, cgcom) ;
    return (status) ;
#endif
#endif

#ifdef TIME
    CGFLOAT tic ;
    cgcom->loop_hessian = CGZERO ;
    cgcom->loop_init = CGZERO ;
    cgcom->loop_prp = CGZERO ;
    cgcom->loop_SSmumps = CGZERO ;
    cgcom->loop_factor = CGZERO ;
    cgcom->loop_sstrust = CGZERO ;
#endif
    
    sF = sizeof (CGFLOAT) ;
    /* extract variables from cgdata structure */
    CGFLOAT           *x = cgdata->x ;
    CGINT const        n = cgdata->n ;
    CGparm         *Parm = cgdata->Parm ;
    CGstat         *Stat = cgdata->Stat ;
    CGFLOAT        *Work = cgdata->Work ;
    int const     FastLA = Parm->FastLA ;
    int const PrintLevel = Parm->PrintLevel ;

    /* if diag_pert used for Hessian in Newton step, initialize Cmax */
    if ( Parm->cg_diag_pert ) /* initialization if diag_pert used for Hessian */
    {
        Cmax = -CGONE ;
    }
    else /* no diagonal perturbation when Newton direction computed */
    {
        Cmax = CGZERO ;
    }

    /* if diag_pert used for Hessian in symmetric solves, initialize Hmax */
    if ( Parm->ss_diag_pert )
    {
        Hmax = -CGONE ;
    }
    else /* no diagonal perturbation when Newton direction computed */
    {
        Hmax = CGZERO ;
    }

    /* By default, when the optimization problem is unconstrained and if the
       Hessian is provided, always try first an iterative PRP+ solver
       applied to the Hessian-based model. */
    htic = CGZERO ;
    dCG = TRUE ;
    minerr = maxerr = CGZERO ;
    RestartAfterHessian = FALSE ;
#ifdef PASA
    int UnitNewtonStep = FALSE ;
    CGINT ni, *order ;
    CGFLOAT fp, maxeqerr, penalty, *AAd, *gpen, *lambda ;
    location = LPASA ;         /* version of CG compiled with PASA */
    SSinitDone = pasacom->SSinitDone ;
    hnnz = pasacom->hnnz ;
    hessian = pasacom->hessian ;
    speed = pasacom->speed ;
    int const QuadCost = pasacom->QP ;
    SOPT_matrix *H = pasacom->pasadata->H ;
    int const use_napheap = pasacom->use_napheap ;
    /* Even though use_napheap is true, the constraint may be inactive.
       nap_present is TRUE if the constraint exists and is active. */
    int const nap_present = ( use_napheap && pasacom->nap_constraint ) ? TRUE :
                                                                        FALSE ;
    PASAparm *pasaparm = pasacom->pasaparm ;
    int const TrustIterLimit = pasaparm->TrustIterLimit ;
    int const use_pproj = pasacom->use_pproj ;
    int const use_restoration = pasacom->use_restoration ;
    int const Aexists = pasacom->Aexists ;
    /* cg does optimization over the free variables */
    pasacom->maxbndstep = PASAINF ;
    pasacom->maxconstep = PASAINF ;

    /* pasa has precedence over CG in choice of parameters */
    deriv_mode = pasacom->deriv_mode ;
    /* if the problem is small and deriv_mode is 2, no need to time or to
       compare speed with CG Newton */
    if ( (deriv_mode == 2) && (n <= Parm->Newton_cutoff) )
    {
        NewMethod = NewtonSS ;
        SetTimer = FALSE ;
    }

    /* err_mark is updated to err_decay times the current value of the local
       error, either pasacom->e (PASA) or gnorm (CG), whenever the current
       error <= err_mark.  iter_since_err_update counts the number
       of CG descent iterations since the last update of err_mark.
       If the number of iterations exceeds cg_ok, then we consider starting
       the Hessian-based iterations. See cg_default for more details. */
    err_decay = pasaparm->err_decay ;
    cg_ok     = pasaparm->cg_ok ;

    /* solution accuracy of the CG Hessian optimizer is equal to
       Hessian_err_decay times the starting error */
    Hessian_err_decay = pasaparm->Hessian_err_decay ;

    order = NULL ;
    if ( QuadCost )
    {
        order = pasacom->order ;
    }
    else
    {
        if ( pasacom->order_with_ifree )
        {
            order = pasacom->ifree ;
        }
        else
        {
            if ( pasacom->order_with_colperm )
            {
                order = pasacom->colperm ;
            }
        }
    }

    CGNfail = &pasacom->CGNfail ; /* *CGNfail = T => CG Newton solver failed */
    SSfail = &pasacom->SSfail ;   /* *SSfail  = T => SS failed or not exist */
    /* if both CG and SS have failed, or CG has failed and SS is not available,
       then need to use ordinary cg_descent */
    if ( *CGNfail && *SSfail )
    {
        pasacom->deriv_mode = deriv_mode = 1 ;
    }
    /* if one or both solvers work and using deriv-based implementation, then
       choose the solver */
    else if ( deriv_mode > 1 )
    {
        if ( !*SSfail ) /* the solvers are available and have not failed */
        {
           /* if the solvers were previously initialized, grab some info */
            if ( SSinitDone )
            {
#ifdef USE_MUMPS
                id = pasacom->id ;
#endif
#ifdef USE_HARWELL
                Cparent = pasacom->Cparent ;
                Cmember = pasacom->Cmember ;
                lbl = pasacom->lbl ;
                lblfile = pasacom->lblfile ;
                Achol = pasacom->Achol ;
                cmm = pasacom->cmm ;
#endif
            }
        }
        /* if either the symmetric solver or the CG Newton methods have not
           failed, then perform initialization */
        if ( !(*SSfail) || !(*CGNfail) )
        {
            C = pasacom->C ;
            C->Aexists = Aexists ;
            HessianSparsityFixed = C->SparsityFixed ;
            /* In pasa, ncol is the problem dimension before removing fixed or
               bound variables. This is the dimension of the Hessian returned
               by the user. */
            ncol = pasacom->ucol ;
            max_d_iter = &pasacom->max_d_iter ;
            C->order = order ;
            if ( order != NULL )
            {
                if ( C->col == NULL )
                {
                    C->col = sopt_malloc (&status, ncol, sizeof (SOPTINT)) ;
                    if ( status ) return (sopt_convert_error (LCG, status)) ;
                }
                sopt_initi (C->col, EMPTY, ncol) ;
                for (j = 0; j < n; j++)
                {
                    C->col [order [j]] = j ;
                }
            }
            else C->col = NULL ;
        }

        /* by default dCG = TRUE is set at start of code.
           If both the symmetric solver and Newton CG have worked,
           then determine whether to start with the CG Newton or the symmetric
           solver direction. */
        if ( !(*SSfail) && !(*CGNfail) )
        {
            /* If cgparm->Hessian_CG_solver = TRUE, then will initially
               use the PR+ version of CG to solve the Hessian-based model.
               Once it is determined that CG is slow, we switch to the
               factorization based solver and use it subsequently by toggeling
               the corresponding pasacom parameter to FALSE */
            if ( pasacom->Hessian_CG_solver ) dCG = TRUE ;  /* CG */
            else                              dCG = FALSE ; /* SS */

            /* modify this rule depending on the measured speed of methods,
               note that the speed is -infinity if not yet tried */
            if      ( speed [NewtonCG] > speed [NewtonSS] ) dCG = TRUE ;
            /* use NewtonSS if it is strictly faster than NewtonCG */
            else if ( speed [NewtonSS] > speed [NewtonCG] ) dCG = FALSE ;
        }
        /* otherwise, use the method that did not fail */
        else if ( !*CGNfail ) dCG = TRUE ; /* only Newton CG worked */
        else if ( !*SSfail  ) dCG = FALSE; /* only SS worked */
    }
#else
    /* standalone CG */
    CGFLOAT cgtol, gnorm, gReset ;
    int const Aexists = FALSE ;
    location = LCG ;

    hessian = cgdata->hessian ;
    /* if n is unassigned, then terminate with error */
    if ( n <= 0 )
    {
        status = cg_wrapup (CG_N_IS_EMPTY, TRUE, cgcom) ;
        return (status) ;
    }
    /* in cg, ncol is the problem dimension */
    ncol = n ;
    /* in cg, row = 0 since there are no constraints */

    if ( (cgdata->hprod != NULL) && (hessian != NULL) )
    {
        status = cg_wrapup (CG_HPROD_PLUS_HESSIAN, TRUE, cgcom) ;
        return (status) ;
    }
    /* is the problem a QP? */
    SOPT_matrix *H = cgdata->H ;
    cgcom->H = H ;

    /* if the Hessian for a quadratic is input as either a dense matrix or
       as a triple, then it is converted to sparse matrix format */
    status = (H->p != NULL) || (H->rows != NULL) || (H->by_rows != NULL)
                                                 || (H->by_cols != NULL) ;
    if ( status || (Parm->QuadCost == 1) || cgdata->hprod )
    {
        /* Objective is quadratic */
        if ( (cgdata->hprod == NULL) && (status == FALSE) )
        {
            /* objective quadrataic since Parm->QuadCost = TRUE. Try to
               generate the sparse matrix by calling the Hessian routine.
               If this is not possible, then report an error. */
            if ( hessian != NULL )
            {
                /* evaluate Hessian at starting point */
                hessian (H, x, ncol) ;
                status = TRUE ;
            }
            else
            {
                status = CG_MISSING_HESSIAN_FOR_QUADCOST ;
                cg_wrapup (status, TRUE, cgcom) ;
                return (status) ;
            }
        }

            /* if the Hessian for a quadratic is input as either a dense matrix
              or as a triple, then convert to sparse matrix format */
        if ( status )
        {
            status = sopt_convert_to_sparse (H, ncol, ncol, TRUE, TRUE, LCG) ;
            if ( status > 0 )
            {
                cg_wrapup (status, TRUE, cgcom) ;
                return (status) ;
            }
            if ( status == SOPT_SPARSE_MATRIX_CREATED )
            {
                cgdata->H_created = TRUE ;
            }

            if ( Parm->CheckMatrix )
            {
                status = sopt_check_matrix (H->p, H->i, H->x, ncol) ;
                if ( status )
                {
                    status = sopt_convert_error (LCG, status) ;
                    cg_wrapup (status, TRUE, cgcom) ;
                    return (status) ;
                }
            }
        }
    }

    /* check to see if the objective is quadratic */
    int const QuadCost = ((H->p == NULL) && (cgdata->hprod == NULL)) ?
                                                         FALSE : TRUE ;
    hnnz = cgdata->hnnz ;
    deriv_mode = Parm->deriv_mode ;
    if ( deriv_mode == -1 )
    {
        if ( (hessian != NULL) || (H->p != NULL) )
        {
            deriv_mode = 12 ;
        }
        else /* no Hessian or cgdata->hprod != NULL */
        {
            deriv_mode = 1 ; /* gradient-based algorithm */
        }
    }
    else if ( deriv_mode == -2 )
    {
        /* If no Hessian is provided, then deriv_mode = 1 */
        if ( (!QuadCost && (cgdata->hessian == NULL)) || (cgdata->hprod!=NULL) )
        {
            deriv_mode = 1 ;
        }
        else
        {
            /* if possible, determine nnz's in Hessian */
            if ( (hnnz == EMPTY) && (H->p != NULL) )
            {
                hnnz = H->p [ncol] ;
            }
            /* use deriv_mode = 1 if Hessian is mostly dense */
            t = CGZERO ;
            if ( hnnz > CGZERO )
            {
                /* t = density = fraction of matrix elements that are nonzero */
                t = hnnz/((CGFLOAT) ncol * (CGFLOAT) ncol) ;
            }
            if ( t >= Parm->dense_Hessian )
            {
                deriv_mode = 1 ;
            }
            else if ( QuadCost )
            {
                deriv_mode = 2 ;
            }
            else
            {
                deriv_mode = 12 ;
            }
        }
    }
    else /* deriv_mode is positive */
    {
        if ( (deriv_mode != 1) && (deriv_mode != 2) && (deriv_mode != 12) )
        {
            status = cg_wrapup (CG_INVALID_DERIV_MODE_PARAMETER, TRUE, cgcom) ;
            return (status) ;
        }
        /* check for Hessian if deriv_mode >= 2 */
        if ( deriv_mode >= 2 )
        {
            if ( (hessian == NULL) && (H->p == NULL) )
            {
                status = CG_DERIV_MODE_USES_HESSIAN_BUT_NO_HESSIAN_PROVIDED ;
                cg_wrapup (status, TRUE, cgcom) ;
                return (status) ;
            }
        }
    }
    if ( deriv_mode > 1 )
    {
        H->nrow = H->ncol = ncol ;
    }

    if ( PrintLevel )
    {
        if ( deriv_mode == 1 )
        {
            printf ("Based on the input data, cg_descent employs a gradient-\n"
                    "based (not Hessian-based) implementation\n") ;
        }
        else
        {
            printf ("Based on the input data, cg_descent employs a Hessian-\n"
                    "based implementation\n") ;
        }
    }

    cgcom->deriv_mode = deriv_mode ;
    int const TrustIterLimit = Parm->TrustIterLimit ;

    if ( deriv_mode > 1 )
    {
        SSfail = &SSfail_int ;
        CGNfail = &CGNfail_int ;
        *SSfail = FALSE ;
#if ( !defined (USE_MUMPS) && !defined (USE_HARWELL) )
        /* If no symmetric solver is available, then *SSfail = TRUE */
        *SSfail = TRUE ;
#endif
        /* Since CG has not yet been tried, it cannot have failed */
        *CGNfail = FALSE ;
        
        max_d_iter = &cgcom->max_d_iter ;
        *max_d_iter = EMPTY ;
        speed = cg_malloc (&status, (CGINT) 3, sizeof (CGFLOAT)) ;
        cgcom->speed = speed ;
        speed [0] = speed [1] = speed [2] = -CGINF ;
#ifdef USE_MUMPS
        cgcom->id = NULL ;
#endif
#ifdef USE_HARWELL
        cgcom->cmm = NULL ;
        cgcom->Achol = NULL ;
        cgcom->lbl = NULL ;
        cgcom->lblfile = NULL ;
        cgcom->Cparent = NULL ;
        cgcom->Cmember = NULL ;
#endif
        cgcom->basisWork = NULL ;
        cgcom->basisiWork = NULL ;
        C = cgcom->C = cg_malloc (&status, (CGINT) 1, sizeof(SOPT_Cmatrix)) ;
        if ( status )
        {
            status = cg_wrapup (sopt_convert_error (LCG, status), TRUE, cgcom) ;
            return (status) ;
        }
        sopt_Cmatrix_default (C) ;
        if ( QuadCost ) HessianSparsityFixed = TRUE ;
        else            HessianSparsityFixed  = Parm->HessianSparsityFixed ;
        C->SparsityFixed = HessianSparsityFixed ;
        C->Aexists = FALSE ;

        /* If the Hessian-based quadratic model is used to obtain the search
       direction in cg_descent, then the optimum of the quadratic model is
       initially obtained by the PRP+ CG method when Hessian_CG_solver is TRUE.
       Otherwise, the optimum is obtained by the symmetric solver. */
        if ( Parm->Hessian_CG_solver ) dCG = TRUE ;  /* use CG */
        else                           dCG = FALSE ; /* use SS */
    }
    else
    {
        cgcom->C = NULL ;
    }

    /* err_mark is updated to err_decay times the current value of the local
       error [either pasacom->e (PASA) or gnorm (CG)] whenever the current
       error <= err_mark.  iter_since_err_update counts the number
       of CG descent iterations since the last update of err_mark.
       If the number of iterations exceeds cg_ok, then we consider starting
       the Hessian-based iterations. See cg_default for more details. */
    err_decay = Parm->err_decay ;
    cg_ok     = Parm->cg_ok ;

   /* solution accuracy of the CG Hessian optimizer is equal to
      Hessian_err_decay times the starting error */
    Hessian_err_decay = Parm->Hessian_err_decay ;
    if ( Parm->PrintParm )
    {
        cg_print_parm (cgdata) ;
    }

    if ( QuadCost )
    {
        if ( H->p != NULL )
        {
            cgcom->hprod_format = 1 ; /* builtin product */
        }
        else /* ( cgdata->hprod != NULL )*/ cgcom->hprod_format = 0 ;
        /* if ( cgdata->c == NULL ), then linear term treated as zero */
        cgcom->QPshift = Parm->QPshift ; /* regularization 0.5*QPshift*x'x */
    }

    /* initialize the statistics */
    cgcom->FirstIter = Stat->iter  = 0 ;/* total number of iterations */
    Stat->nfunc    = 0 ; /* number of function evaluations */
    Stat->ngrad    = 0 ; /* number of gradient evaluations */
    Stat->nexpand  = 0 ; /* number of times cgexpand was called */
    Stat->nforward = 0 ; /* number of forward (expansion) steps */
    Stat->nback    = 0 ; /* number of backward steps */
    Stat->IterSub  = 0 ; /* number of iterations in subspace */
    Stat->NumSub   = 0 ; /* total number of subspaces */
    Stat->nCG      = 0 ; /* number Newton steps computed by CG */
    Stat->nSS      = 0 ; /* number Newton steps computed by SS */
    Stat->PRP      = 0 ; /* number of PRP iterations in Newton CG */
    cgcom->f = CGINF ;
    Stat->err = CGINF ;
    /* end of stand alone CG */
#endif

    cgcom->QuadCost = QuadCost ;
    use_hessian = FALSE ;
    NewtonStep = FALSE ; /* change to TRUE in Newton step */
    if ( deriv_mode > 1 ) /* set parameters for Hessian-based implementation */
    {
        Newton_window = Parm->Newton_window ;
        CG_window = Parm->CG_window ;
        dCG_old = dCG ;
        if ( (deriv_mode == 2) && !(*SSfail) )
        {
            use_hessian = TRUE ;
            /* if n <= Newton_cutoff, then set dCG = FALSE */
            if ( n <= Parm->Newton_cutoff )
            {
                dCG = FALSE ; /* start with SS */
            }
        }

        C->n = n ;
        /* if hnnz is empty, then extract from H->p if available */
        if (  H->p != NULL )  hnnz = H->p [ncol] ;
        if (  hnnz == EMPTY ) iter_check = -1 ;
        else                  iter_check = hnnz / ncol ;
        if ( PrintLevel )
        {
            printf ("hnnz: %li ncol: %li iter_check: %li\n",
                    (LONG) hnnz, (LONG) ncol, (LONG) iter_check);
        }
        iter_since_err_update = 0 ;
        nCG = 0 ; /* # search directions computed by CG applied to Hessian */
        nSS = 0 ; /* # search directions computed by direct symmetric solver */
        skipHessian = FALSE ; /* changes to TRUE when Hessian computed */
        SSsetupDone = FALSE ;
        /* if SS has not failed, then find max_d_iter, the max # CG iterations
           if SS has failed, then max_d_iter is set to n */
        if ( *SSfail )
        {
            find_max_iter = FALSE ;
            if ( *max_d_iter == EMPTY ) *max_d_iter = n ;
        }
        else /* symmetric solver did not fail */
        {
            /* find max # CG iterations if SS did not fail, max_d_iter is
               EMPTY, deriv_mode = 2 or 12, and CG Newton has not failed */
            if ( (*max_d_iter == EMPTY) && (deriv_mode >= 2) && !(*CGNfail) )
            {
                find_max_iter = TRUE ;
            }
            else
            {
                /* use the provided value for max_d_iter */
                find_max_iter = FALSE ;
            }
        }
        trustfails = 0 ;
    }

    /* if x is NULL, then allocate memory and initialize x to zero */
    if ( x == NULL )
    {
        x = (CGFLOAT *) cg_malloc (&status, n, sizeof (CGFLOAT)) ;
        if ( status )
        {
            status = cg_wrapup (sopt_convert_error (LCG, status), TRUE, cgcom) ;
            return (status) ;
        }
        cg_initx (x, CGZERO, n) ;
        cgdata->x = cgdata->x_created = x ;
    }

    if ( PrintLevel >= 1 )
    {
        printf ("\nSTART CG\n") ;
    }
    CGFLOAT const grad_tol = Parm->grad_tol ;

    /* initialize statistics structure */
    cgcom->work_created = NULL ;

    cgcom->neps = 0 ;
    maxit = Parm->maxit ; /* abort when # iterations reaches maxit */

    Stat->maxit = maxit ;
    Stat->NegDiag = FALSE ; /* no negative diagonal element in QR factor */
    Stat->maxsteps = Parm->maxsteps ;
    Stat->cg_ninf_tries = Parm->cg_ninf_tries ;

    /* cgcom structure initialization */
    cgcom->x = x ;
    cgcom->n = n ; /* problem dimension */
    cgcom->approxstep = Parm->approxstep ;
    ApproxSwitchFactor = Parm->ApproxSwitchFactor ;

    fbest = CGINF ; /* stores best function value */
    gbest = CGINF ; /* stores best norm of gradient */
    nslow = 0 ;     /* number of slow iterations */
    slowlimit = 2*n + Parm->nslow ;
    memk = 0 ;      /* number of vectors in subspace current memory */
    cgcom->AvoidFeval = FALSE ;

    Ck = CGZERO ; /* average cost magnitude */
    Qk = CGZERO ; /* used in estimate of average cost magnitude */

    /* limited memory initializations */
    use_memory = FALSE ;/* => ignore CG memory */
    Subspace = FALSE ; /* => work in full space */
    FirstFull = FALSE ;/* => not first full iteration after leaving subspace */
    Restart = FALSE ;  /* => no restart in limited memory algorithm */
    IterRestart = 0 ; /* counts number of iterations since last restart */
    IterQuad = 0 ;    /* counts number of iterations that function change
                         is close to that of a quadratic */

    /* copy more values from Parm */
    cgcom->pert_eps = Parm->pert_eps ;
    cgcom->PertRule = Parm->PertRule ;
    cgcom->Wolfe    = FALSE ; /* initially a Wolfe line search not performed */
    delta2          = 2*Parm->cgdelta - CGONE ;

    /*  LBFGS = 0 => use cg_descent
                1 => use L-BFGS
                2 => use L-BFGS when LBFGSmemory >= n,
                     use cg_descent when LBFGSmemory < n
                3 => use L-BFGS when LBFGSmemory >= n, use limited
                     memory cg_descent when LBFGSmemory < n */
    LBFGS        = Parm->LBFGS ;

    qrestart     = CGMIN (n, Parm->qrestart) ;
    /* Initially the function is assumed to be not approximately quadratic.
       However, if QuadCost is TRUE, then the value of QuadF is ignored. */
    QuadF = FALSE ;

    /* the conjugate gradient algorithm is restarted every nrestart iteration */
    nrestart = (CGINT) (((CGFLOAT) n)*Parm->restart_fac) ;

    /* adjust LBFGS depending on value of mem */
    mem = Parm->LBFGSmemory ;
    mem = CGMIN (mem, n) ;
    if ( (mem == n) && (LBFGS >= 2) )
    {
        LBFGS = 1 ; /* use LBFGS */
        /* switch the deriv_mode to 1 unless deriv_mode = 2 and dimension <=
           Newton_cutoff (a small problem that will be solved by Newton) */
        if ( (deriv_mode != 2) || (n > Parm->Newton_cutoff) )
        {
            deriv_mode = 1 ;
            use_hessian = FALSE ;
        }
    }
    /* if mem <= 2 and limited memory CG is requested, then switch to
       ordinary CG since the code breaks in Yk update. Note that
       L-BFGS still works with memory = 2, so LBFGS = 1 could be specified. */
    if ( (mem <= 2) && (LBFGS >= 2) )
    {
        mem = 0 ;
        LBFGS = 0 ; /* use cg_descent */
    }

    /* memory allocation */
    if ( Work == NULL )
    {
        /* memory for g, d, xnew, gnew */
        i = 4*n ;
        if ( QuadCost == TRUE )
        {
            i += n ;
        }
        if ( LBFGS == 1 ) /* use L-BFGS */
        {
            /* Sk (mem*n), Yk (mem*n), SkYk (mem), tau (mem) */
            i += 2*mem*(n+1) ;
        }
        else if ( LBFGS == 3 )
        {
            /* SkF (mem*n), stemp (n), gkeep (n), Sk (mem*mem), Rk (mem*mem),
               Re (mem+1) Yk (mem*(mem+1) +2), SkYk (mem), tau (mem),
               dsub (mem), gsub (mem+1), wsub (mem+1), vsub (mem),
               gsubtemp (mem) */
            i += (mem+2)*n + (3*mem+9)*mem + 5 ;
        }
        if ( deriv_mode >= 2 )
        {
            /* Hd, r, h_dir */
            i += 3*n ;
        }
        work = (CGFLOAT *) cg_malloc (&status, i, sizeof (CGFLOAT)) ;
        if ( status )
        {
            status = cg_wrapup (sopt_convert_error (LCG, status), TRUE, cgcom) ;
            return (status) ;
        }
        cgcom->work_created = work ;
    }
    else
    {
        work = Work ;
    }

#ifndef PASA
    /* running CG by itself, not through PASA, no A */
    cgcom->g    = g    = work ; work += n ;
    cgcom->d    = d    = work ; work += n ;
    cgcom->xnew = xnew = work ; work += n ;
    cgcom->gnew = gnew = work ; work += n ;

    /* with the follows definitions, some formulas work with either
       cg_descent or pasa */
    gproj = g ;
    gnewproj = gnew ;
    D = d ;
    if ( QuadCost )
    {
        cgcom->Qd = work ; work += n ;
        cgcom->hprod = cgdata->hprod ;
        cgcom->c = cgdata->c ;
    }
    else /* not a quadratic */
    {
        if ( (cgdata->value == NULL) || (cgdata->grad == NULL) )
        {
            status = cg_wrapup (CG_VALUE_OR_GRAD_MISSING, TRUE, cgcom) ;
            return (status) ;
        }
    }
    cgcom->value   = cgdata->value ;
    cgcom->grad    = cgdata->grad ;
    cgcom->valgrad = cgdata->valgrad ;
    t = CGZERO ;
    XXCG(evaluate) (t, &t, "fg", PASA_CG_COM) ;
    f = cgcom->f ;
    if ( (f != f) || (f >= CGINF) || (f <= -CGINF) )
    {
        status = CG_STARTING_FUNCTION_VALUE_INFINITE_OR_NAN ;
        cg_wrapup (status, TRUE, cgcom) ;
        return (status) ;
    }

    /* set d = -g, compute gnorm  = infinity norm of g and
                           gnorm2 = square of 2-norm of g */
    gnorm = cg_update_inf2 (g, g, d, &gnorm2, n) ;
    Stat->tol = cgtol = CGMAX (gnorm*Parm->StopFac, grad_tol) ;
    Stat->err = gnorm ;
    if ( gnorm <= cgtol )
    {
        status = cg_wrapup (CG_ERROR_TOLERANCE_SATISFIED, FALSE, cgcom) ;
        return (status) ;
    }

    /* For a QP, the gradient calculation above was done from scratch.
       Subsequently, the gradient is updated during each iteration.
       gReset is the gradient norm that triggers a fresh computation of
       the gradient from scratch. */
    if ( QuadCost )
    {
        gReset = Parm->QPgReset_factor * gnorm ;
    }
    dnorm2 = gnorm2 ;
    dphi0 = -gnorm2 ;
    err_mark = err_decay*gnorm ;

    /* determine starting stepsize */
    alpha = Parm->step ;
    if ( alpha == CGZERO )
    {
        CGFLOAT xnorm ;
        xnorm = cg_sup_normx (x, n) ;
        if ( xnorm == CGZERO )
        {
            if ( f != CGZERO ) alpha = 2.*fabs (f)/gnorm2 ;
            else               alpha = CGONE ;
        }
        else alpha = Parm->psi0*xnorm/gnorm ;
    }
    cgcom->maxstep = maxstep = CGINF ; /* there are no bounds or constraints */
#else
    /* CG is run from PASA */
#ifndef USE_MUMPS
#ifndef USE_HARWELL
    CGINT nrow ;
#endif
#endif
    nrow = pasacom->nrow ;
    Stat->tol = grad_tol ;
    cgcom->FirstIter = Stat->iter ; /* store first iteration */
    pasacom->location = PASA_CG_DESCENT ;
    int use_penalty = pasacom->use_penalty ;

    /* cg_bb_est is initialized to be -1. Its value is reset to an estimate
       for the BB parameter at various places in the code where it is
       easy to generate an estimate */
    pasacom->cg_bb_est = -1.0 ;

    /* If use_penalty is TRUE, then A must exist.
       Check to see if there are active constraints */
    if ( use_penalty )
    {
        /* check for an active equality */
        if ( use_pproj )
        {
            PPINT *RLinkUp ;
            PPwork *W ;
            W = pasacom->ppcom->Work ;
            RLinkUp = W->RLinkUp ;
            status = ( nrow == RLinkUp [nrow] ) ; /* T => no active equality */
        }
        else /* use_napheap = TRUE */
        {
            status = (pasacom->nap_constraint == 0) ||
                     (pasacom->nap_a2 == CGZERO) ; /* T => no active equality */
        }

        if ( status == TRUE )
        {
            use_penalty = FALSE ;
            pasacom->use_penalty = FALSE ;
        }
        else
        {
            penalty = pasaparm->penalty ;
            pasacom->penalty = penalty ;
        }
    }

    /* array allocation */
    g    = cgcom->g    = pasacom->g ;    /* gradient at x */
    d    = cgcom->d    = pasacom->d ;    /* search direction */
    xnew = cgcom->xnew = pasacom->xnew ; /* new x at x + alpha*d */
    gnew = cgcom->gnew = pasacom->gnew ; /* new g at x + alpha*d */
    AAd = NULL ;
    gpen = NULL ;
    if ( Aexists == TRUE )
    {
        if ( use_penalty == TRUE )
        {
             if ( PrintLevel >= 1 )
             {
                 printf ("include penalty term in objective\n") ;
             }
            lambda = pasacom->lambda_pen ;
            gpen = pasacom->gpen ;
            AAd  = work ;  work += n ;
        }
        gproj    = work ;  work += n ;
        gnewproj = work ;  work += n ;
        D        = work ;  work += n ;

        /* print bound rows */
        if ( PrintLevel >= 2 )
        {
            printf ("bound rows:\n") ;
            for (i = 0; i < pasacom->nr; i++)
            {
                j = pasacom->bound_rows [i] ;
                if ( j < 0 )
                {
                    printf ("%ld lower bound\n", (LONG) -(j+1)) ;
                }
                else
                {
                    printf ("%ld upper bound\n", (LONG) (j-1)) ;
                }
            }
        }
    }
    else /* only bound constraints are present */
    {
        gproj = g ;
        gnewproj = gnew ;
        D = d ;
    }

    cgcom->gpen = gpen ;
    cgcom->AAd = AAd ;
    cgcom->approxstep = pasacom->approxstep ;/* PASA overwrites CG start value*/

    /* penalty cost and derivative terms initialized to zero */
    pasacom->fp = CGZERO ; /* value of the penalty function at current iterate*/
    pasacom->dp = CGZERO ; /* derivative of penalty in search direction */
    pasacom->Ad2 = CGZERO ;/* penalty*||Ad||^2 */
    ni = 0 ;               /* no inequality constraints by default */

    /* If use_penalty = TRUE, then compute:
       (1) multiplier estimate lambda = (BB')^{-1}Bg, B = active rows and cols,
       (2) the gradient gpen of the penalty term lambda'(b-Bx) + 0.5p||b-Bx||^2
           gpen = B'(p(Bx-b) - lambda)
       (3) the value of the penalty term fp = lambda'(b-Bx) + 0.5p||b-Bx||^2 */
    fp = CGZERO ;
    maxeqerr = CGZERO ;
#ifndef NOPPROJ
/*  Begin code to compute orthonormal basis for active rows of A.
    Replace with Davis' sparse QR ASAP. */
    if ( pasaparm->use_QR )
    {
        PASAINT const *RLinkUp = pasacom->ppcom->Work->RLinkUp ;
        PPFLOAT const     *ATx = pasacom->ppcom->Work->ATx ;
        PPINT   const     *ATi = pasacom->ppcom->Work->ATi ;
        PPINT   const     *ATp = pasacom->ppcom->Work->ATp ;

        PASAFLOAT *A = work ;  work = work+(n*nrow) ;
        pasacom->Z   = work ;  work = work+(n*PASAMIN (n, nrow)+nrow) ;
        PASAFLOAT *AP = A ;
        j = 0 ;
        for (i = RLinkUp [nrow]; i < nrow; i = RLinkUp [i])
        {
            t = PASAZERO ;
            PPINT const q = ATp [i+1] ;
            for (p = ATp [i]; p < q; p++)
            {
                t += ATx [p]*ATx [p] ;
            }
            t = sqrt (t) ;
            if ( t > PASAZERO ) /* add the row to a column of A */
            {
                j++ ;
                pasa_initx (AP, PASAZERO, n) ;
                t = PASAONE/t ;
                for (p = ATp [i]; p < q; p++)
                {
                    AP [ATi [p]] = t*ATx [p] ;
                }
                AP = AP+n ;
            }
        }
#if 0
        PASAFLOAT *Acopy = (PASAFLOAT *) malloc (n*j*sizeof (PASAFLOAT)) ;
        pasa_copyx (Acopy, A, n*j) ;
#endif
        pasacom->Zncol = pasa_null (pasacom, A, j) ;
        printf ("projection matrix: %ld by %ld nonzero rows of original "
                "A: %ld\n", (LONG) n, (LONG) pasacom->Zncol, (LONG) j) ;
#if 0
        /* check Z */
        PASAINT Acols = j ;
        PASAFLOAT *Z = pasacom->Z ;
        for (j = 0; j < pasacom->Zncol; j++)
        {
            t = pasa_dot (Z+j*n, Z+j*n, n) ;
            if ( fabs (1 - t) > 1.e-12 )
            {
                printf ("error in norm of column %i, dot: %e\n", j, t) ;
                pasa_error (-1, __FILE__, __LINE__, "stop") ;
            }
        }
        for (i = 0; i < pasacom->Zncol; i++)
        {
            for (j = i+1; j < pasacom->Zncol; j++)
            {
                t = pasa_dot (Z+i*n, Z+j*n, n) ;
            }
            if ( fabs (t) > 1.e-12 )
            {
                printf ("error in column %i dot column %i: %e\n", i, j, t) ;
                pasa_error (-1, __FILE__, __LINE__, "stop") ;
            }
        }
        for (j = 0; j < Acols; j++)
        {
            for (i = 0; i < pasacom->Zncol; i++)
            {
                t = pasa_dot (Acopy+j*n, Z+i*n, n) ;
                if ( fabs (t) > 1.e-12 )
                {
                    printf ("error in row %i dot col %i: %e\n", j, i, t) ;
                    pasa_error (-1, __FILE__, __LINE__, "stop") ;
                }
            }
        }
        free (Acopy) ;
        pasa_error (-1, __FILE__, __LINE__, "stop") ;
#endif
    }
/*  end code to compute orthonormal basis for active rows of A */
    if ( (use_penalty == TRUE) && (use_pproj == TRUE) )
    {
        PPINT *Ap, *Ai, *Anz, *ATp, *ATi, *ir, *RLinkDn, *RLinkUp ;
        PPFLOAT c, *Ax, *ATx, *Axk, *b ;
        PPprob *Prob ;
        PPwork *W ;
        W = pasacom->ppcom->Work ;
        RLinkDn = W->RLinkDn ;
        RLinkUp = W->RLinkUp ;
        Prob = pasacom->ppcom->Prob ;
        b = pasacom->b ;
        Ap = Prob->Ap ;
        Anz = Prob->Anz ;
        Ai = Prob->Ai ;
        Ax = Prob->Ax ;
        pasa_initx (lambda, PASAZERO, nrow) ;
        /* this initx could be replaced by the following more precise
           code, however, valgrind errors could be generated in the back solve
           below since the factorization could have zeros that were not
           squeezed out
        k = nrow ;
        while ( (k = RLinkUp [k]) < nrow )
        {
            lambda [k] = PASAZERO ;
        } */
        for (j = 0; j < n; j++)
        {
            t = g [j] ;
            if ( t != PASAZERO )
            {
                k = Ap [j] ;
                l = k + Anz [j] ;
                for (; k < l; k++)
                {
                    lambda [Ai [k]] += t*Ax [k] ;
                }
            }
        }
        pproj_lsol (W->L, lambda, RLinkUp [nrow], nrow, RLinkUp) ;
        k = RLinkUp [nrow] ;
        /* momentarily set the initial RLinkDn to -1, this simplifies
           indexing in dltsolve */
        RLinkDn [k] = -1 ;
        l = RLinkDn [nrow] ;
        pproj_dltsol (W->L, lambda, lambda, l, k, RLinkDn) ;

        RLinkDn [k] = nrow ; /* restore RLinkDn */

        pasa_initx (gpen, PASAZERO, n) ;
        ATp = W->ATp ;
        ATi = W->ATi ;
        ATx = W->ATx ;
        Axk = pasacom->Axk ;
        ir = W->ir ;
        c = 0.5*penalty ;
        p = 0 ;
        l = 0 ;
        ni = Prob->ni ;
        for (i = 0; i < nrow; i++)
        {
            s = CGZERO ;
            PPINT const q = ATp [i+1] ;
            for (; p < q; p++)
            {
                s += ATx [p]*x [ATi [p]] ;
            }
            if ( ir [i] == 0 )
            {
                t = b [i] - s ;
                fp += t*(c*t + lambda [i]) ;
                s = -(penalty*t + lambda [i]) ;
                t = fabs (t) ;
                maxeqerr = CGMAX (t, maxeqerr) ;
                p = ATp [i] ;
                for (; p < q; p++)
                {
                    gpen [ATi [p]] += ATx [p]*s ;
                }
            }
            else
            {
                ASSERT (ir [i] > ni) ;
                l++ ;
                Axk [l] = s ;
            }
        }
        ASSERT (l == ni) ;
        cgcom->maxeqerr = maxeqerr ;
        if ( QuadCost == FALSE )
        {
            /* Save function value before adding the penalty term to it.  */
            pasacom->f_orig = pasacom->f ;
            pasacom->f += fp ;
        }
        /* else if QuadCost = TRUE, then pasacom->f is unchanged */
    }
    else if ( (use_penalty == TRUE) && (use_napheap == TRUE) )
    {
        PASAFLOAT c, *a ;

        a = pasacom->nap_a ;
        /* multiplier estimate = a'g/a'a */
        const CGFLOAT Lambda = pasa_dot (a, g, n)/pasacom->nap_a2 ;
        *lambda = Lambda ;
        c = 0.5*penalty ;
        s = pasa_dot (a, x, n) ;
        t = pasacom->nap_bl - s ;
        cgcom->maxeqerr = fabs (t) ;
/*printf ("i: %i lambda: %e Ax-b: %e b: %e\n", i, lambda [i], t, b [i]) ;*/
/*printf ("i: %i t: %25.15e lambda: %25.15e\n", i, t, lambda [i]) ;*/

        fp = t*(c*t + Lambda) ;
        s = -(penalty*t + Lambda) ;
        pasa_scale (gpen, a, s, n) ;
        if ( QuadCost == FALSE )
        {
            /* Save function value before adding the penalty term to it.  */
            pasacom->f_orig = pasacom->f ;
            pasacom->f += fp ;
        }
        /* else if QuadCost = TRUE, then pasacom->f is unchanged.
           function values for quadratic cost problems exclude penalty term */
    }
    else if ( Aexists && use_pproj )
    {
        PPprob *Prob ;
        PPwork *W ;
        PPINT row, *ineq_row, *ATp, *ATi ;
        PPFLOAT *ATx, *Axk ;
        /* Compute A*x for the inactive inequalities */
        Prob = pasacom->ppcom->Prob ;
        ni = Prob->ni ;
        if ( ni > 0 )
        {
            ineq_row = Prob->ineq_row ;
            W = pasacom->ppcom->Work ;
            ATp = W->ATp ;
            ATi = W->ATi ;
            ATx = W->ATx ;
            Axk = pasacom->Axk ;
            for (i = 1; i <= ni; i++)
            {
                row = ineq_row [i] ;
                s = CGZERO ;
                CGINT const q = ATp [row+1] ;
                for (p = ATp [row]; p < q; p++)
                {
                    s += ATx [p]*x [ATi [p]] ;
                }
                Axk [i] = s ;
            }
        }
    }
    else if ( Aexists && use_napheap )
    {
        if ( pasacom->nap_constraint == 0 )
        {
            pasacom->Axk [1] = pasa_dot (pasacom->nap_a, x, n) ;
        }
    }
#endif

    pasacom->fp = fp ;
    cgcom->f = f = pasacom->f ;

    /* Compute projected gradient and store in gproj.
       gproj is the projection of the total gradient g + gpen */
    pasa_null_project (gproj, g, gpen, TRUE, pasacom) ;

    /* see if the convergence tolerance is satisfied */
    if ( pasacom->e <= pasacom->testtol )
    {
        status = pasa_check_error (pasacom) ;
        if ( (status == PASA_ERROR_TOLERANCE_SATISFIED) ||
             (status == PASA_GRAD_PROJ) )
        {
            status = XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
            return (status) ;
        }
        Stat->tol = pasacom->testtol ;
    }
    err_mark = err_decay*pasacom->e ;

    /* D is the search direction before the final projection */
    cg_scale (D, gproj, -CGONE, n) ;

    /* d is the final search direction */
    pasa_null_project (d, D, NULL, FALSE, pasacom) ;
    dnorm2 = cg_dot (d, d, n) ;

    /* derivative in search direction without penalty term */
    dphi0 = cg_dot (g, d, n) ;
    if ( use_penalty == TRUE )
    {
        pasacom->dp = cg_dot (gpen, d, n) ;
    }

    alpha = pasacom->bbk ; /* starting step in cg is previous bb step */
    scale = pasacom->bbk ; /* scale is the approximation to inverse
                              Hessian in LBFGS */

    /* if there is no estimate for the bb parameter, use the same
       estimates for the initial stepsize that is used in the pure cg code */
    /* set d = -g, compute gnorm  = infinity norm of g and
                           gnorm2 = square of 2-norm of g */

    if ( alpha < PASAZERO )
    {
        alpha = Parm->step ;
        if ( alpha == CGZERO )
        {
            CGFLOAT xnorm ;
            xnorm = cg_sup_normx (x, n) ;
            if ( xnorm == CGZERO )
            {
                if ( f != CGZERO ) alpha = 2.*fabs (f)/dnorm2 ;
                else               alpha = CGONE ;
            }
            else alpha = Parm->psi0*xnorm/cg_sup_normx (d, n) ;
        }
        /* there are no bounds or constraints */
        cgcom->maxstep = maxstep = CGINF ;
    }
/* end of pasa initialization */
#endif

    if ( LBFGS == 1 ) /* allocate storage connected with LBFGS */
    {
        if ( PrintLevel >= 1 )
        {
            printf ("use LBFGS, dimension %ld\n", (LONG) n) ;
        }
        LBFGS = TRUE ;      /* use L-BFGS */
        mlast = -1 ;
        Sk = work ;    work += mem*n ;
        Yk = work ;    work += mem*n ;
        SkYk = work ;  work += mem ;
        tau = work ;   work += mem ;
    }

    else if ( LBFGS == 3) /* allocate storage connected with limited memory CG*/
    {
        if ( PrintLevel >= 1 )
        {
            printf ("use Limited Memory CG, dimension %ld\n", (LONG) n) ;
        }
        use_memory = TRUE ;           /* previous search direction is saved */
        SubSkip = 0 ;                 /* # iterations to skip checking memory*/
        SubCheck = mem*Parm->SubCheck;/* number of iterations to check */
        StartCheck = Stat->iter ;     /* start checking memory at initial iter*/
        InvariantSpace = FALSE ;      /* iterations not in invariant space */
        FirstFull = TRUE ;            /* first iteration in full space */
        nsub = 0 ;                    /* initial subspace dimension */
        memsq = mem*mem ;
        SkF = work ;   work += mem*n ;/* directions in memory (x_k+1 - x_k) */
        stemp = work ; work += n ;    /* stores x_k+1 - x_k */
        gkeep = work ; work += n ;    /* store grad when 1st direction != -g */
        Sk = work ;    work += memsq ;/* Sk = Rk at start LBFGS in subspace */
        Rk = work ;    work += memsq ;/* upper triangular factor in SkF=Zk*Rk*/
        Re = work ;    work += mem+1 ;/* end column of Rk, for new direction*/
        Yk = work ;    work += memsq+mem+2 ;
        SkYk = work ;  work += mem ;  /* dot products sk'yk in the subspace */
        tau = work ;   work += mem ;  /* stores alpha in Nocedal and Wright */
        dsub = work ;  work += mem ;  /* direction projection in subspace */
        gsub = work ;  work += mem+1 ;/* gradient projection in subspace */
        wsub = work ;  work += mem+1 ;/* work array for triangular solve */
        vsub = work ;  work += mem ;  /* work array for triangular solve */
        gsubtemp = work ; work += mem;/* new gsub before update */

        cg_initx (Rk, CGZERO, memsq) ; /* initialize Rk to 0 */
    }

    cgcom->f0 = f + f ;
    cgcom->SmallCost = fabs (f)*Parm->SmallCost ;
    cgcom->df0 = -2.0*fabs(f)/alpha ;
    FirstIter = cgcom->FirstIter + 1 ;

    /* When deriv_mode = 2 or 12, we compare the error decay rate to the
       speed to decide whether to switch to a different method:
       PureCG (0), NewtonCG (1), or NewtonSS (2).  */
    AutoChoice = ( (deriv_mode >= 2) &&
                   (LBFGS != 1 || (mem != n) || (n > Parm->Newton_cutoff)) ) ;
    /* Initialize timer */
    if ( AutoChoice && (FirstIter == 1) ) sopt_timer () ;

    /* Start the conjugate gradient iteration.
       alpha starts as old step, ends as final step for current iteration
       f is function value for alpha = 0
       QuadOK = TRUE means that a quadratic step was taken */

    /* CGSTART */
    while ( Stat->iter < maxit )
    {
#ifdef PASA
#ifndef NDEBUG
        pasa_check_feas (pasacom) ;
#endif
        maxstep = XXCG(maxstep) (pasacom, cgcom) ;
        /* iterate hits boundary immediately if maxstep = 0 */
        if ( maxstep == CGZERO )
        {
            status = CG_HITS_BOUNDARY ;
            /* When cg hits boundary of the feasible region,
               return to active gradient projection scheme. But first,
               perform an expansion step if the problem is only bound
               constrained or a Newton step was computed. */
            if ( QuadCost && (pasacom->OnlyBounds || NewtonStep) )
            {
                int flag = pasa_cgexpandQP (x, d, n, maxstep, PASAZERO, pasacom,
                                            NewtonStep, UnitNewtonStep) ;
                if ( flag != PASA_OK ) status = flag ;
                /* NOTE: x and gtot were updated as well as pasacom->f */

                if ( order != NULL )
                {
                    pasa_convert_to_pproj (g, pasacom->gtot, order, n) ;
                }
                else
                {
                    pasa_copyx (g, pasacom->gtot, n) ;
                }
#ifndef NDEBUG
                pasa_check_feas (pasacom) ;
#endif
                status = XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                return (status) ;
            }
            else if ( pasacom->OnlyBounds || NewtonStep ) /* problem not a QP */
            {
                int flag = pasa_cgexpand (x, d, n, maxstep, cgcom->f, pasacom,
                                          cgcom, NewtonStep, UnitNewtonStep) ;
                if ( flag != PASA_OK ) status = flag ;

                /* NOTE: x = xnew and g = gnew */
                status = XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                /* either restart cg_descent or branch to active gradproj */
                return (status) ;
            }
            /* NOTE: In the trust region part of the Hessian-based routine,
               the objective is evaluated and the current cgcom->f is destroyed.
               The variable f contains the function value at step = 0.
               Hence, we need to restore its value in the structures. */
            if ( use_hessian )
            {
                cgcom->f = f ;
#ifdef PASA
                pasacom->f = f ;
#endif
            }
            pasacom->update_factor = TRUE ;
            status = XXCG(wrapup) (CG_HITS_BOUNDARY, FALSE, PASA_CG_COM) ;
            return (status) ;
        }
        err = pasacom->e ;
#else
        err = gnorm ;
#endif
        Stat->iter++ ;
        if ( PrintLevel >= 3 )
        {
            printf ("CG iter: %ld f: %25.15e err: %e memk: %i dphi0: %e "
                    "dCG: %i use_hessian: %i\n",
                    (LONG) Stat->iter, f, err, memk, dphi0, dCG, use_hessian) ;
        }
        if ( QuadCost )
        {
            /* function value for quadratic cost problems exclude penalty term*/
            cgcom->f0 = f ; /* save prior function value */
#ifdef PASA
            t = -CGONE ;
            /* the following evaluation computes the product between
               the Hessian and the search direction */
            status = pasa_evaluate (CGZERO, &t, pasacom, "fg") ;
            if ( status == PASA_FUNCTION_NAN_OR_INF )
            {
                XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                return (status) ;
            }
            if ( order != NULL )
            {
                dQd = PASAZERO ;
                /* dQd = pasa_dot (d, cgcom->>Qd, n) but Qd is in user's coor */
                for (j = 0; j < n; j++)
                {
                    dQd += d [j]*pasacom->Qd [order [j]] ;
                }
            }
            else /* no bounds, no linear constraints, no fixed variables */
            {
                dQd = pasa_dot (d, pasacom->Qd, n) ;
            }
            /* when dQd <= 0, we take the max step along the search direction
               include the penalty term in this computation */
            if ( dQd + pasacom->Ad2 <= CGZERO )
            {
                if ( PrintLevel >= 2 )
                {
                    printf ("dQd = %e <= 0, take maxstep\n", dQd) ;
                }
                if ( maxstep >= CGINF ) /* unbounded objective */
                {
                    status = CG_QUADRATIC_OBJECTIVE_NO_LOWER_BOUND ;
                    XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                    return (status) ;
                }
                alpha = maxstep ;
            }
            else /* compute the exact minimizer in search direction */
            {
                /* printf ("dphi0: %e dp: %e dQd: %e Ad2: %e\n",
                            dphi0, pasacom->dp, dQd, pasacom->Ad2) ;*/
                /* penalty term is included */
                alpha = -(dphi0+pasacom->dp)/(dQd + pasacom->Ad2) ;
            }
            if ( PrintLevel )
            {
                printf ("cg QP stepsize: %e maxstep: %e\n", alpha, maxstep) ;
            }
            if ( alpha >= maxstep ) /* hits boundary */
            {
                if ( PrintLevel )
                {
                    printf ("cg hits boundary in QuadCost\n") ;
                }
                /* cg hits boundary of the feasible region;
                   if boundary is hit in the first iteration, then
                   method, method_tic, ...  have not yet been set */
                if ( AutoChoice && (Stat->iter > FirstIter) )
                {
                    cg_update_speed (speed, NewMethod, minerr, maxerr,
                                     method_tic, sopt_timer ());
                }

                PASAFLOAT Df = maxstep*(0.5*dQd*maxstep + dphi0) ;
                int flag = pasa_cgexpandQP (x, d, n, maxstep, Df, pasacom,
                                            NewtonStep, UnitNewtonStep) ;
                if ( flag != PASA_OK ) status = flag ;
                else                   status = CG_HITS_BOUNDARY ;
                /* NOTE: x and gtot were updated as well as pasacom->f */

                if ( order != NULL )
                {
                    pasa_convert_to_pproj (g, pasacom->gtot, order, n) ;
                }
                else
                {
                    pasa_copyx (g, pasacom->gtot, n) ;
                }
#ifndef NDEBUG
                pasa_check_feas (pasacom) ;
#endif
                status = XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                return (status) ;
            }

            /* take the step alpha */
            cg_step (xnew, x, d, alpha, n) ;

            if ( PrintLevel )
            {
                printf ("sup-norm of x: %e\n", pasa_sup_normx (xnew, n)) ;
            }

            /* update cost without penalty */
            pasacom->f += alpha*(0.5*dQd*alpha + dphi0) ;
            cgcom->f = pasacom->f ;

            /* update gradient without penalty */
            pasa_step (pasacom->gtot, pasacom->gtot, pasacom->Qd,
                       alpha, pasacom->ncol) ;

            /* copy relevant part of gtot into gnew */
            if ( order != NULL )
            {
                pasa_convert_to_pproj (gnew, pasacom->gtot, order, n) ;
            }
            else
            {
                pasa_copyx (gnew, pasacom->gtot, n) ;
            }

            /* when the penalty term is included, the derivative in the
               search direction is zero in theory */
            cgcom->df = CGZERO ;
            /* in the subspace routines, dphi0 includes the penalty term */
            dphi0 += pasacom->dp ;

            cgcom->QuadOK = TRUE ;
            cgcom->alpha = alpha ;
            status = CG_WOLFE_OK ;
#else
            /* CG is not being used with pasa, there is no penalty.
               This is the same as the block of code above, but with all
               the penalty and constraint related statements removed.
               In the block of code below, the objective is an unconstrained
               quadratic, and the stepsize is the Cauchy step.  */
            t = -CGONE ;
            XXCG(evaluate) (t, &t, "fg", PASA_CG_COM) ;
            dQd = cg_dot (d, cgcom->Qd, n) ;
            if ( dQd <= CGZERO )
            {
                status = CG_QUADRATIC_OBJECTIVE_NO_LOWER_BOUND ;
                XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                return (status) ;
            }

            /* compute the exact minimizer in search direction */
            alpha = -dphi0/dQd ;

            /* take the step alpha */
            cg_step (xnew, x, d, alpha, n) ;

            /* update cost */
            cgcom->f -= 0.5*dphi0*dphi0/dQd ;

            /* update gradient */
            cg_step (gnew, g, cgcom->Qd, alpha, n) ;

            cgcom->df = CGZERO ; /* Cauchy step, deriv in search direction 0 */
            cgcom->QuadOK = TRUE ;
            cgcom->alpha = alpha ;
            status = CG_WOLFE_OK ;
#endif
        }
        else /* general nonlinear objective */
        {
#ifdef PASA
            dphi0 += pasacom->dp ;
#endif

            if ( use_hessian && NewtonStep )
            {
                alpha = CGONE ;
                cgcom->f_at = CGONE ;
                cgcom->df_at = -CGONE ; /* derivative not yet computed */
                cgcom->QuadOK = TRUE ;
                alphabound = CGINF ;
            }
            else
            {
                /* will store the point where f or df is evaluated */
                cgcom->f_at = -CGONE ;
                cgcom->df_at = -CGONE ;
                cgcom->QuadOK = FALSE ;

                alpha *= Parm->psi2 ;
                if ( f != CGZERO )
                {
                    t = fabs ((f-cgcom->f0)/f) ;
                }
                else
                {
                    t = CGONE ;
                }
                cgcom->UseCubic = TRUE ;
                if ( (t < Parm->CubicCutOff) || !Parm->UseCubic )
                {
                    cgcom->UseCubic = FALSE ;
                }
                if ( Parm->QuadStep )
                {
                    /* test if quadratic interpolation step should be tried */
                    if ( ((t > Parm->QuadCutOff)&&(fabs(f) >= cgcom->SmallCost))
                          || QuadF )
                    {
                        /* check if the function is approximately quadratic */
                        if ( QuadF )
                        {
                            alpha *= Parm->psi1 ;
#ifdef PASA
                            /* truncate the nominal step if it exceeds maxstep*/
                            if ( alpha >= maxstep )
                            {
                                alpha = maxstep ;
                            }
#endif
                            t = alpha ;
                            status = XXCG(evaluate)
                                     (CGZERO, &alpha, "g", PASA_CG_COM) ;
                            /* alpha may decrease if df is infinite or nan */
                            if ( status == CG_FUNCTION_NAN_OR_INF )
                            {
                                XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                                return (status) ;
                            }
                            cgcom->df_at = alpha ;

                            /* if stepsize reduced due to infinite function
                               values, store the valid stepsize */
                            if ( alpha < t )
                            {
                                alphabound = alpha ;
                            }
                            else
                            {
                                alphabound = CGINF ;
                            }

                            /* secant approximation to step */
                            if ( cgcom->df > dphi0 )
                            {
                                alpha = -dphi0/((cgcom->df-dphi0)/alpha) ;
                                cgcom->QuadOK = TRUE ;
                            }
                            else if ( LBFGS == 1 )
                            {
                                if ( memk >= n )
                                {
                                    alpha = CGONE ;
                                    cgcom->QuadOK = TRUE ;
                                }
                                else
                                {
                                    alpha = 2. ;
                                }
                            }
                            else if ( Subspace )
                            {
                                if ( memk >= nsub )
                                {
                                    alpha = CGONE ;
                                    cgcom->QuadOK = TRUE ;
                                }
                                else  alpha = 2. ;
                            }
                        }
                        else /* not approximately quadratic */
                        {
                            /* Attempt a quadratic quadratic interpolation step
                               based on the starting function value and
                               derivative, and the function value at the point
                               "trial" below.
                               If the step is successful in the sense that the
                               quadratic is convex, then keep a safeguarded
                               version of the step. Otherwise, retain the
                               original alpha. */

                            t = CGMAX (Parm->psi_lo,
                                                cgcom->df0/(dphi0*Parm->psi2)) ;
                            trial = alpha*CGMIN (t, Parm->psi_hi) ;
#ifdef PASA
                            /* truncate the nominal step if it exceeds maxstep*/
                            if ( trial >= maxstep )
                            {
                                trial = maxstep ;
                            }
#endif
                            t = trial ;
                            status = XXCG(evaluate)
                                            (CGZERO, &trial, "f", PASA_CG_COM) ;
                            /* trial may decrease if f is infinite or nan */
                            if ( status == CG_FUNCTION_NAN_OR_INF )
                            {
                                XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                                return (status) ;
                            }
                            cgcom->f_at = trial ;

                            /* if stepsize reduced due to infinite function
                               values, store the valid stepsize */
                            if ( trial < t )
                            {
                                alphabound = trial ;
                            }
                            else
                            {
                                alphabound = CGINF ;
                            }
                            fnew = cgcom->f ;
                            denom = 2.*(((fnew-f)/trial)-dphi0) ;
                            if ( denom > CGZERO ) /* convex quadratic */
                            {
                                /* t = quadratic interpolation iterate */
                                t = -dphi0*trial/denom ;
                                /* safeguard */
                                if ( fnew >= f )
                                {
                                    alpha = CGMAX (t, trial*Parm->QuadSafe) ;
                                }
                                else
                                {
                                    alpha = t ;
                                }
                                cgcom->QuadOK = TRUE ;
                            }
                            else /* quadratic is concave */
                            {
                                if ( (trial > alpha) &&
                                   (fabs(fnew-f) > ApproxSwitchFactor*Ck) )
                                {
                                    alpha = trial ;
                                }
                            }
                        }
                        if ( PrintLevel >= 2 )
                        {
                            if ( denom <= CGZERO )
                            {
                                printf("Quad step fails (denom = %14.6e)\n",
                                       denom) ;
                            }
                            else if ( cgcom->QuadOK )
                            {
                                printf ("Quad step %14.6e OK\n",
                                         CGMIN (alpha, maxstep)) ;
                            }
                            else
                            {
                                printf ("Quad step %14.6e done, but not OK\n",
                                         CGMIN (alpha, maxstep)) ;
                            }
                        }
                    }
                    else if ( PrintLevel >= 2 )
                    {
                        printf ("No quad step (chg: %14.6e, cut: %10.2e)\n",
                                 t, Parm->QuadCutOff) ;
                    }
                }
            }
            if ( alpha > maxstep )
            {
                alpha = maxstep ;
                cgcom->QuadOK = FALSE ;
            }
            if ( alpha > alphabound )
            {
                alpha = alphabound ;
                cgcom->QuadOK = FALSE ;
            }

            cgcom->f0 = f ; /* save prior value (corresponding to step = 0) */
            cgcom->df0 = dphi0 ;

            if ( cgcom->PertRule == TRUE )
            {
                cgcom->fpert = f + cgcom->pert_eps*fabs (f) ;
            }
            else
            {
                cgcom->fpert = f + cgcom->pert_eps ;
            }

            cgcom->wolfe_hi = Parm->cgdelta*dphi0 ;
            cgcom->wolfe_lo = Parm->cgsigma*dphi0 ;
            cgcom->awolfe_hi = delta2*dphi0 ;
            cgcom->alpha = alpha ;

            /* perform line search */
            status = XXCG(line) (FALSE, PASA_CG_COM) ;
            if ( status == CG_FUNCTION_NAN_OR_INF )
            {
                XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                return (status) ;
            }

#ifdef PASA
            /* cg hits boundary of the feasible region, either restart
               cg_descent or return to activeGP if Parm->use_activeGP is TRUE */
            if ( status == CG_HITS_BOUNDARY )
            {
                if ( PrintLevel >= 1 )
                {
                     printf ("cg hits boundary in line search\n") ;
                }
                /* if boundary is hit in the first iteration, then
                   method, method_tic, ...  have not yet been set */
                if ( AutoChoice && (Stat->iter > FirstIter) )
                {
                    cg_update_speed (speed, NewMethod, minerr, maxerr,
                                     method_tic, sopt_timer ());
                }

                int flag = pasa_cgexpand (x, d, n, maxstep, cgcom->f, pasacom,
                                          cgcom, NewtonStep, UnitNewtonStep) ;
                if ( flag != PASA_OK ) status = flag ;

                /* NOTE: x = xnew and g = gnew */
                status = XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                /* either restart cg_descent or branch to active gradproj */
                return (status) ;
            }
#endif
            /*try approximate Wolfe line search if ordinary Wolfe fails */
            if ( (status != CG_WOLFE_OK) && !cgcom->approxstep )
            {
                if ( PrintLevel >= 1 )
                {
                     printf ("\nWOLFE LINE SEARCH FAILS\n") ;
                }
                if ( status != CG_SLOPE_ALWAYS_NEGATIVE )
                {
                    cgcom->approxstep = TRUE ;
                    cgcom->alpha = alpha ;
                    status = XXCG(line) (TRUE, PASA_CG_COM) ;
                    if ( status == CG_FUNCTION_NAN_OR_INF )
                    {
                        XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                        return (status) ;
                    }
                }
            }
        }

#ifdef PASA
#ifndef NOPPROJ
        /* no restoration when performing either purely unconstrained or
           purely bound-constrained optimization */
        if ( use_restoration && (use_pproj || use_napheap) )
        {
            int tstatus = status ;
            status = pasa_restoration (pasacom, cgcom) ;
            if ( status == CG_HITS_BOUNDARY )
            {
                /* if boundary is hit in the first iteration, then
                   method, method_tic, ...  have not yet been set */
                if ( AutoChoice && (Stat->iter > FirstIter) )
                {
                    cg_update_speed (speed, NewMethod, minerr, maxerr,
                                     method_tic, sopt_timer ());
                }
                pasacom->update_factor = TRUE ;
                status = XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                /* either restart cg_descent or branch to active gradproj */
                return (status) ;
            }
            status = tstatus ;
        }
#endif
#endif

        /* parameters in Wolfe, approximate Wolfe conditions, and update */
        Qk = Parm->Qdecay*Qk + CGONE ;
        cgcom->Ck = Ck = Ck + (fabs (f) - Ck)/Qk ; /* average cost magnitude */

        alpha = cgcom->alpha ;
        f = cgcom->f ; /* final value of objective after line search */
        dphi = cgcom->df ;

        if ( status != CG_WOLFE_OK ) /* line search fails */
        {
            XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
            return (status) ;
        }

#ifdef PASA
        /* update the penalty gradient if penalty term present */
#ifndef NOPPROJ
        if ( use_penalty )
        {
            cg_daxpy (gpen, AAd, alpha*penalty, n) ;
            pasacom->fp += alpha*(pasacom->dp + 0.5*alpha*pasacom->Ad2);
        }

        /* update A*xk based on the current stepsize:
           new Axk = old Axk + alpha * Adk */
        if ( ni > 0 )
        {
            cg_daxpy (pasacom->Axk+1, pasacom->Adk+1, alpha, ni) ;
        }
#endif
#endif

        if ( !QuadCost && !use_hessian )
        {
            /* test how close the cost function changes are to that of a
               quadratic QuadTrust = 0 means the function change matches
               that of a quadratic*/
            t = alpha*(dphi+dphi0) ;
            if ( fabs (t) <= Parm->qeps*CGMIN (Ck, CGONE) )
            {
                QuadTrust = CGZERO ;
            }
            else
            {
                QuadTrust = fabs((2.0*(f-cgcom->f0)/t)-CGONE) ;
            }
            if ( QuadTrust <= Parm->qrule)
            {
                IterQuad++ ;
            }
            else
            {
                IterQuad = 0 ;
            }

            if ( IterQuad == qrestart )
            {
                QuadF = TRUE ;
            }

            /* as function values converge, they contain less information
               and are avoided when possible in the line search */
            cgcom->AvoidFeval = FALSE ;
            if ( fabs (f-cgcom->f0) <= Parm->CostConverge*fabs (f) )
            {
                cgcom->AvoidFeval = TRUE ;
            }

            /* check whether to change to approximate Wolfe line search */
            if ( !cgcom->approxstep )
            {
                if ( fabs (f-cgcom->f0) < ApproxSwitchFactor*Ck )
                {
                    if ( PrintLevel >= 1 )
                    {
                         printf ("change to approximate Wolfe line search\n") ;
                    }
                    cgcom->approxstep = TRUE ;
                    if ( cgcom->Wolfe == TRUE )
                    {
                        Restart = TRUE ;
                    }
                }
            }
        }

        IterRestart++ ;

        /* some bookkeeping for limited memory CG */
        if ( (LBFGS == 3) && !use_hessian )
        {
            if ( use_memory == TRUE )
            {
                if ( (Stat->iter - StartCheck > SubCheck) && !Subspace )
                {
                    StartSkip = Stat->iter ;
                    use_memory = FALSE ;
                    if ( SubSkip == 0 ) SubSkip = mem*Parm->SubSkip ;
                    else                SubSkip *= 2 ;
                    if ( PrintLevel >= 1 )
                    {
                        printf ("skip subspace %i iterations\n", SubSkip) ;
                    }
                }
            }
            else
            {
                if ( Stat->iter - StartSkip > SubSkip )
                {
                    StartCheck = Stat->iter ;
                    use_memory = TRUE ;
                    memk = 0 ;
                }
            }
        }
        if ( !use_hessian && (deriv_mode > 1) )
        {
            if ( n > mem ) /* deriv_mode = 12, if n <= mem, use L-BFGS */
            {
                /* when problem small enough, use NewtonSS before NewtonCG*/
                if ( n <= Parm->Newton_cutoff )
                {
                    dCG = FALSE ;
                }
                /* above err is set to pasacom->e (pasa) or to gnorm (cg)*/
                if ( err < err_mark )
                {
                    err_mark = err_decay*err ;
                    iter_since_err_update = 0 ;
                }
                else
                {
                    iter_since_err_update++ ;
                }
                /* possibly set use_hessian TRUE based on sparsity */
                if ( (iter_since_err_update > cg_ok) &&
                     (iter_since_err_update >= iter_check) )
                {
                    if ( iter_check < 0 ) /* iter_check not yet computed */
                    {
                        Xp = xnew ;
#ifdef PASA
                        if ( order != NULL )
                        {
                            Xp = pasacom->userx ;
                            for (j = 0; j < n; j++) Xp [order [j]] = xnew[j] ;
                        }
#endif
                        htic = sopt_timer () ;
                        hessian (H, Xp, ncol) ;
                        /* store the time to evaluate the Hessian */
                        htic = sopt_timer () - htic ;
                        /* Since the hessian was just evaluated, set
                           skipHessian to TRUE to avoid another Hessian
                           evaluation below. It is set to FALSE after
                           the search direction computation. */
                        skipHessian = TRUE ;
                        if      ( H->p    != NULL ) hnnz = H->p [ncol] ;
                        else if ( H->rows != NULL )
                        {
                            if ( H->sym ) hnnz = H->nnz ;
                            else          hnnz = 2*H->nnz ;
                        }
                        else
                        {
                            if ( (H->by_rows == NULL)&&(H->by_cols == NULL))
                            {
                                status = SOPT_HESSIAN_NOT_COMPUTED ;
                                status =sopt_convert_error(location,status) ;
                                XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                                return (status) ;
                            }
                            if ( H->by_rows != NULL ) dense = H->by_rows ;
                            else                      dense = H->by_cols ;
                            hnnz = 0 ;
                            le =  (LONG) n * (LONG) n ;
                            for (li = 0; li < le; i++)
                            {
                                if ( dense [i] ) hnnz++ ;
                            }
                        }
#ifdef PASA
                        pasacom->hnnz = hnnz ;
#endif
                        iter_check = hnnz / ncol ;
                        if ( PrintLevel )
                        {
                            printf ("hnnz: %li ncol: %li iter_check: %li\n",
                                  (LONG) hnnz, (LONG) ncol, (LONG) iter_check) ;
                        }
                        if ( iter_since_err_update < iter_check )
                        {
                            /* we are not going to switch to Hessian-based
                               implementation, so will need to compute
                               the Hessian at a new x value when doing a
                               Newton step */
                            skipHessian = FALSE ;
                        }
                    }
                    if ( iter_since_err_update >= iter_check )
                    {
                        use_hessian = TRUE ;
                    }
                }
            }
        }

        if ( PrintLevel && (deriv_mode > 1) )
        {
            printf ("use_hessian: %i iter_since_err_update: %li iter_check: "
                    "%li mem: %i dCG: %i\n", use_hessian,
                    (LONG) iter_since_err_update, (LONG) iter_check, mem, dCG) ;
            printf ("    err: %e err_mark: %e\n", err, err_mark) ;
        }

        /* if we do not plan to do a restart, put the projection of gnew
           into gnewgproj, otherwise put it into gproj */
        if ( use_hessian || ((LBFGS == 1) && (IterRestart == nrestart)) ||
             (((LBFGS == 0) || (LBFGS == 2) ||
               ((LBFGS == 3) && (use_memory == FALSE)))
                           && ((IterRestart >= nrestart) ||
                     ((IterQuad == qrestart) && (IterQuad != IterRestart)))) )
        {
            RestartSetupDone = TRUE ;
            /* project gnew into the null space of the active
               constraint gradients and store projection in gproj */
#ifdef PASA
            if ( Aexists )
            {
                pasa_null_project (gproj, gnew, gpen, TRUE, pasacom) ;
            }
            else
            {
                pasa_copyx (gproj, gnew, n) ;
                pasacom->e = pasa_sup_normx (gnew, n) ;
            }
#else
            /* gnorm2 is needed during a restart */
            gnorm2 = cg_dot (gnew, gnew, n) ;
#endif
            if ( (LBFGS == 3) && !use_hessian )
            {
                Restart = TRUE ;
            }
        }
#ifdef PASA
        else /* no restart */
        {
            RestartSetupDone = FALSE ;
            /* project gnew into the null space of the active
               constraint gradients and store projection in gnewproj */
            pasa_null_project (gnewproj, gnew, gpen, TRUE, pasacom) ;
        }

        if ( PrintLevel >= 2 )
        {
            printf ("Local error: %e testtol: %e Global error Egp: %e E1: %e\n",
                    pasacom->e, pasacom->testtol, pasacom->Egp, pasacom->E1) ;
        }
        /* test nominal stopping condition */
        err = pasacom->e ;
        if ( err <= pasacom->testtol )
        {
            status = pasa_check_error (pasacom) ;
            if ( (status == PASA_ERROR_TOLERANCE_SATISFIED) ||
                 (status == PASA_GRAD_PROJ) )
            {
                /* set x = xnew, g = gnew */
                cg_copy (x, xnew, n) ;
                cg_copy (g, gnew, n) ;
                if ( AutoChoice && (status != PASA_ERROR_TOLERANCE_SATISFIED)
                                && (Stat->iter > FirstIter) )
                {
                    if      ( err < minerr ) minerr = err ;
                    else if ( err > maxerr ) maxerr = err ;
                    cg_update_speed (speed, OldMethod, minerr, maxerr,
                                     method_tic, sopt_timer ());
                }
                status = XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                return (status) ;
            }
            Stat->tol = pasacom->testtol ;
        }
#else
        /* stand alone cg, no pasa */
        t = gnorm ;
        err = Stat->err = gnorm = cg_sup_normx (gnew, n) ;
        if ( gnorm <= cgtol )
        {
            cg_copy (x, xnew, n) ;
            status = CG_ERROR_TOLERANCE_SATISFIED ;
            status = XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
            return (status) ;
        }
        /* for a quadratic objective, both the function value and gradient
           are evaluated by an update process. If the gradient shows a big
           drop, then re-evaluate objective value and gradient to
           maintain accuracy */

        if ( QuadCost && (gnorm <= gReset) )
        {
            if ( PrintLevel >= 1 )
            {
                printf ("re-evaluate quadratic objective and gradient\n") ;
            }
            gReset = Parm->QPgReset_factor * gnorm ;
            cgcom->x = xnew ;
            cgcom->g = gnew ;
            t = CGZERO ;
            XXCG(evaluate) (t, &t, "fg", PASA_CG_COM) ;
            f = cgcom->f ;
            cgcom->x = x ;
            cgcom->g = g ;
        }
#endif

        /* When choosing between PureCG, NewtonCG, and NewtonSS, we need an
           estimate for the speed of the different methods */
        /* AUTO */
        if ( AutoChoice )
        {
            /* In first iteration, a special startup process is used */
            if ( Stat->iter == FirstIter )
            {
                start_it = Stat->iter ;
                if ( use_hessian )
                {
                    /* Newton-based algorithm, Newton timer is started after
                       Newton initialization process (SS analysis and
                       flop count comparison between NewtonSS and NewtonCG) */
                    SetTimer = TRUE ;
                    if ( dCG ) NewMethod = NewtonCG ;
                    else       NewMethod = NewtonSS ;
                    end_it = start_it + Newton_window ;
                }
                else /* use PureCG */
                {
                    NewMethod = PureCG ;
                    /* PureCG timer can be started right away since there is
                       no initialization process with PureCG */
                    method_tic = sopt_timer() ;
                    end_it = start_it + CG_window ;
                }
                OldMethod = NewMethod ;
                maxerr = minerr = err ;
                use_hessian_old = use_hessian ;
                if ( PrintLevel )
                {
                    printf ("AutoChoice: %i deriv_mode: %i LBFGS: %i mem: %i "
                            "n: %li\n",
                             AutoChoice, deriv_mode, LBFGS, mem, (LONG) n) ;
                }
            }
            /* compute speed of current method in any of the following cases:
               1. use_hessian != use_hessian_old (we were either using a
                  Newton-based iteration and switch to PureCG or we were
                  using PureCG but now switch to a Newton-based method.
               2. use_hessian = use_hessian_old = TRUE, but dCG != dCG_old
                  (we switch the Newton methods).
               3. iter >= end_it (we have performed a block of iterations
                  for a method, will estimate its speed, and will switch to
                  the fastest method) */
            if ( PrintLevel >= 3 )
            {
                speed_info = "" ;
                printf ("Old Method: %i maxerr: %e minerr: %e start_it: %li\n"
                        "    end_it: %li", OldMethod, maxerr,
                         minerr, (LONG) start_it, (LONG) end_it) ;
                printf (" speeds: %10.2e (CG) %10.2e (NCG) %10.2e (NSS)\n",
                             speed [0], speed [1], speed [2]) ;
            }
            if      ( err < minerr ) minerr = err ;
            else if ( err > maxerr ) maxerr = err ;
            if ( (use_hessian != use_hessian_old) || (dCG != dCG_old) ||
                 (Stat->iter >= end_it) )
            {
                /* new start_it if method potentially changes */
                if ( PrintLevel >= 3 ) start_it = Stat->iter ;
                /* estimate speed of the current method (OldMethod) */
                method_toc = sopt_timer () ;
                cg_update_speed (speed, OldMethod, minerr, maxerr,
                                 method_tic, method_toc) ;
                minerr = maxerr = err ;
                method_tic = method_toc ;
                if ( use_hessian != use_hessian_old )
                {
                    if ( use_hessian )
                    /* previously used PureCG, now use either NewtonCG
                       or NewtonSS */
                    {
                        if ( dCG ) NewMethod = NewtonCG ;
                        else       NewMethod = NewtonSS ;

                        /* if PureCG is faster than the new Newton method,
                           then modify iter_check to make it less likely
                           to employ the Newton method */
                        if ( (speed [NewMethod] > -CGINF) &&
                             (speed [PureCG] > speed [NewMethod]) )
                        {
                            err_mark = err_decay*err ;
                            iter_since_err_update = 0 ;
                            if ( iter_check > 1 ) iter_check *= 2 ;
                            else                  iter_check  = 2 ;
                        }
                        /* if Newton faster than PureCG, expand Newton window */
                        else if ( speed [NewMethod] > speed [PureCG] )
                        {
                            Newton_window *= 1.3 ;
                        }
                        end_it = Stat->iter + Newton_window ;
                    }
                    else /* Hessian previously used, but now use PureCG */
                    {
                        /* NewtonSS failed or the NewtonCG search direction
                           was poor */
                        NewMethod = PureCG ;
                        end_it = Stat->iter + CG_window ;
                        err_mark = err_decay*err ;
                        iter_since_err_update = 0 ;
                    }
                }
                /* Newton is being used, but type switches */
                else if ( use_hessian && (dCG != dCG_old) )
                {
                    end_it = Stat->iter + Newton_window ;
                    if ( dCG ) NewMethod = NewtonCG ;
                    else       NewMethod = NewtonSS ;
                }
                else if ( Stat->iter >= end_it ) /* use the fastest method */
                {
                    if ( PrintLevel >= 3 )
                    {
                        speed_info = "(change method for speed)" ;
                    }

                    CGFLOAT best = -CGINF ;
                    int mbest = 1 ;
                    for (i = 0; i < 3; i++)
                    {
                        if ( speed [i] > best )
                        {
                            best = speed [i] ;
                            mbest = i ;
                        }
                    }
                    /* if the fastest method is not the current (NewMethod),
                       then switch NewMethod to the fastest method */
                    if ( mbest != NewMethod )
                    {
                        NewMethod = mbest ;
                        /* if the fastest method is  Newton-based, then bias
                           future choice so as to favor Newton-based methods */
                        if ( NewMethod >= NewtonCG )
                        {
                            use_hessian = TRUE ;
                            Newton_window *= 1.3 ;
                            end_it += Newton_window ;

                            if ( NewMethod == NewtonCG )
                            {
                                /* if NewtonCG beat NewtonSS in speed, but
                                   NewtonSS was the OldMethod and it was close
                                   to the speed of NewtonCG, then continue
                                   using NewtonSS (since NewtonCG could
                                   require more iterations to converge, making
                                   it less attractive) */
                                if ( (speed [NewtonSS] > -CGINF) &&
                                     (OldMethod == NewtonSS) &&
                                     (1.3*speed [NewtonCG] < speed [NewtonSS]) )
                                {
                                    NewMethod = NewtonSS ;
                                    /* dCG = FALSE continues */
                                }
                                else
                                {
                                    dCG = TRUE ;
                                }
                            }
                            else /* NewtonSS is winner, continue using it */
                            {
                                dCG = FALSE ;
                            }
                        }
                        else /* new fastest method is PureCG */
                        {
                            /* mark the point where we start to evaluate the
                               speed of PureCG in order to decide when it would
                               be good to switch to a Newton method */
                            err_mark = err_decay*err ;
                            iter_since_err_update = 0 ;
                            use_hessian = FALSE ;
                            CG_window *= 1.3 ;
                            end_it += CG_window ;
                            /* if we return to Newton, then use the previous
                               version which was fastest */
                            if ( speed [NewtonCG] > speed [NewtonSS] )
                            {
                                dCG = TRUE ;
                            }
                            else
                            {
                                dCG = FALSE ;
                            }
                        }
                    }
                    /* mbest = current method, retain it and increase window */
                    else
                    {
                        if ( NewMethod >= NewtonCG )
                        {
                            Newton_window *= 1.3 ;
                            end_it += Newton_window ;
                        }
                        else
                        {
                            CG_window *= 1.3 ;
                            end_it += CG_window ;
                        }
                    }
                }

                /* a CG restart is needed when we change from a Newton-based
                   method to PureCG */
                RestartAfterHessian = FALSE ;
                if ( (OldMethod >= NewtonCG) && (NewMethod == PureCG) )
                {
                    if (LBFGS == 0 || LBFGS == 1 || LBFGS == 2)
                    {
                        IterRestart = nrestart ;
                    }
                    else if (LBFGS == 3 )
                    {
                        RestartAfterHessian = TRUE ;
                        Restart = TRUE ;
                    }
                    /* in the block of code above that does the setup for a
                       restart, we need to make some adjustments if the
                       setup was not done */
                    if ( !RestartSetupDone )
                    {
#ifdef PASA
                        /* the new projected gradient was stored in gnewproj,
                           now move it to gproj where the restart code expects
                           to find it */
                        sopt_copyx (gproj, gnewproj, n) ;
#else
                        /* gnorm2 is needed during a restart */
                        gnorm2 = cg_dot (gnew, gnew, n) ;
#endif
                    }
                }

                /* The time stored in method_toc above can be used as
                   the method_tic for PureCG. Otherwise, for a Newton method,
                   the timer needs to be started after any initialization
                   is complete. */
                /* For a Newton method, we delay setting the timer. */
                if ( NewMethod >= NewtonCG )
                {
                    SetTimer   = TRUE ;
                }
                /* For PureCG, the timer start is set to method_toc */
                else
                {
                    method_tic = method_toc ;
                }
            
                OldMethod = NewMethod ;
            }

            if ( PrintLevel >= 3 )
            {
                printf ("New Method: %i maxerr: %e minerr: %e start_it: "
                        "%li\n    end_it: %li ", NewMethod, maxerr, minerr,
                         (LONG) start_it, (LONG) end_it) ;
                printf ("speeds: %10.2e (CG) %10.2e (NCG) %10.2e (NSS) %s\n",
                             speed [0], speed [1], speed [2], speed_info) ;
            }
        }

        /* compute search direction, starting with a Newton-based direction
           when use_hessian is TRUE */
        /* NEWTON */
        use_hessian_old = use_hessian ;
        if ( use_hessian )
        {
            /* PASA is expecting that the new projected gradient is stored in
               gproj. If not, then we need to copy it to this location.
               The storage into gproj is done above in the section for
               Restart Setup. */
#ifdef PASA
            if ( !RestartSetupDone )
            {
                pasa_copyx (gproj, gnewproj, n) ;
            }
#endif
            /* save dCG to see if it is changed in the Newton routine */
            dCG_old = dCG ;
            /* Let H denote current Hessian. Compute a solution to the problem

               min 1/2 d' (H + sigma I) d + g'd subject to Bd = 0,
               

               where B = active constraints.
               First try to compute a solution using PRP+ version of CG.
               Let P denote the projection into the null space of B. We
               make the change of variables d = PD to obtain the following
               equivalent problem which we feed to CG:

               min 1/2 (PD)' (H + sigma I) (PD) + gproj' D, gproj = Pg

               If convergence of CG is slow, then permanently switch to the
               factorization-based solution approach.  */

            /* set x = xnew, g = gnew */
            cg_copy (x, xnew, n) ;
            cg_copy (g, gnew, n) ;

#ifdef TIME
tic = sopt_timer() ;
#endif
            /* Evaluate the Hessian at x if problem not a QP */
            if ( !QuadCost && !skipHessian )
            {
                Xp = x ;
#ifdef PASA
                if ( order != NULL )
                {
                    Xp = pasacom->userx ;
                    for (j = 0; j < n; j++) Xp [order [j]] = xnew [j] ;
                }
#endif
                /* evaluate the Hessian, and time the evaluation if htic == 0 */
                if ( htic )
                {
                    hessian (H, Xp, ncol) ;
                }
                else
                {
                    htic = sopt_timer () ;
                    hessian (H, Xp, ncol) ;
                    htic = sopt_timer () - htic ;
                }
            }
            skipHessian = FALSE ;
#ifdef TIME
cgcom->loop_hessian += sopt_timer () - tic ;
#endif

#if 0
/*sopt_printAMATLAB (ncol, H->p, H->i, NULL, H->x, "H") ;*/
sopt_printxMATLAB (cgdata->c, n, "c") ;
SOPT_matrix HH, *pHH ;
pHH = &HH ;
sopt_matrix_default (pHH) ;
i = H->nnz ;
printf ("i: %i\n", i) ;
pHH->rows = sopt_malloc (&status, i, sizeof (SOPTINT)) ;
pHH->cols = sopt_malloc (&status, i, sizeof (SOPTINT)) ;
pHH->vals = sopt_malloc (&status, i, sF) ;
printf ("XXstatus: %i\n", status) ;
sopt_copyi (pHH->rows, H->rows, i) ;
sopt_copyi (pHH->cols, H->cols, i) ;
sopt_copyx (pHH->vals, H->vals, i) ;
pHH->nnz = i ;
pHH->sym = TRUE ;
status = sopt_convert_to_sparse (pHH, ncol, ncol, TRUE, TRUE, LCG) ;
printf ("status: %i\n", status) ;
fflush(stdout) ;
sopt_printAMATLAB (ncol, pHH->p, pHH->i, NULL, pHH->x, "H") ;
sopt_free (pHH->p) ;
sopt_free (pHH->i) ;
sopt_free (pHH->x) ;
/*sopt_printxMATLAB (x, n, "x") ;*/
/*sopt_printxMATLAB (gproj, n, "g") ;*/
#endif

            /* --- DECLARATIONS --- */

            /* for SS and MUMPS: */
#ifdef USE_MUMPS
            char **argv ;
            char *name = "c_example" ;
            argv = &name ;
            int argc, myid ;
            argc = 1 ;
#endif
            CGINT Anrow, Annz, m, nt, *basisiWork, *Mi, *Mp, *perm ;
            /* end of SS */

            /* --- for CG --- */
            int NegCurv ;
            CGINT *col, *row ;
            CGFLOAT beta_k, qcost, dHd, dtol, dHdsave, rho,
                    rnorm, rnorm_start, ss_flops, step, *basisWork,
                    *Mb, *Mx, *pd, *pv, *T, *V, *X, *Xb, *Z, *Z1, *Z2 ;
#ifdef PASA
            CGFLOAT *PHd ;

            /* NewtonStep is TRUE if a Newton step was performed. If the
               trust region step terminates at the starting point (the
               approximate minimizer of the quadratic) before any
               backtracking, then UnitNewtonStep is TRUE. If in the
               subsequent line search, it is bound that a new constraint is
               activated, and UnitNewtonStep is TRUE, then we expand the
               stepsize in the line search, and project onto the feasible
               set after each expansion. The expansion process continues
               until the cost no longer decays. When only bound constraints
               are present, we always perform the expansion phase since the
               projection process is trivial. */
            UnitNewtonStep = FALSE ;
#endif
            /* NewtonStep is TRUE if an approximation to a Newton step
               was generated, and the line search terminated successfully
               at a point with a smaller objective value, and the Newton
               search direction is a descent direction */
            NewtonStep = TRUE ;

            /* --- end of CG --- */

            /* --- END OF DECLARATIONS --- */

#ifdef TIME
tic = sopt_timer() ;
#endif

            /* the first time entering the Hessian-based algorithm, set
               up the workspace for the trust region step */
            if ( nCG + nSS == 0 )
            {
#ifdef PASA
                basisWork  = pasacom->basisWork ;
                basisiWork = pasacom->basisiWork ;
#else
                basisWork  = cgcom->basisWork ;
                basisiWork = cgcom->basisiWork ;
#endif
                if ( basisWork == NULL )
                {
                    /* 3*ncol (Z), 6 (T),  3*ncol (H), 3 (Anorm), 3 (colnorm)
                       9     (Mx), 3 (Mb), 2 (Xb) */
                    basisWork = sopt_malloc (&status, 9*ncol+29, sF) ;
                    if ( status )
                    {
                        status = sopt_convert_error (location,status) ;
                        XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                        return (status) ;
                    }
                }
                if ( basisiWork == NULL )
                {
                    /* perm [3], Mi [9], Mp [4] */
                    basisiWork = sopt_malloc(&status,16,sizeof(CGINT));
                    if ( status )
                    {
                        status = sopt_convert_error (location, status) ;
                        XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                        return (status) ;
                    }
                }
#ifdef PASA
                pasacom->basisWork = basisWork ;
                pasacom->basisiWork = basisiWork ;
#else
                cgcom->basisWork = basisWork ;
                cgcom->basisiWork = basisiWork ;
#endif
                Z = basisWork ;  basisWork += 3*n ;
                V = basisWork ;  basisWork += 3*n ;
                T = basisWork ;  basisWork += 6 ;
                Mx= basisWork ;  basisWork += 9 ;
                Mb= basisWork ;  basisWork += 3 ;
                X = basisWork ;  basisWork += 3 ;
                Xb= basisWork ;  basisWork += 2 ;

                perm = basisiWork ; basisiWork += 3 ;
                Mi   = basisiWork ; basisiWork += 9 ;
                Mp   = basisiWork ; basisiWork += 4 ;

                /* Initialize the symmetric solvers if they are available
                   and they have not failed.  This initialization is done
                   only one time for either solver. */
                if ( !*SSfail && !SSinitDone )
                {
                    SSinitDone = TRUE ;
#ifdef USE_MUMPS
                    /* MPI_Init and MPI_Comm_rank return int ierr */
                    MPI_Init (&argc, &argv) ;
                    MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
                    id = sopt_malloc (&status, sizeof (DMUMPS_STRUC_C), 1) ;
                    id->comm_fortran = USE_COMM_WORLD ;
                    id->par = 1 ;
                    id->sym = 2 ;
                    id->job = JOB_INIT ;

                    dmumps_c (id) ;

                    /* the following parameter is recommended by the MUMPS
                       team to help allow more node amalgamation in
                       the elimination tree */
                    id->keep [196] = 1 ;
                    /* icntl [3] = 0 no printing, = 1 error + warning,
                                 = 2 = 1 + stats */
                    id->icntl [3] = 1 ;
                    id->icntl [13] = 100 ;
                    /* NOTE: one less than shown in manual, could use:
                                 #define ICNTL(I) icntl[(I)-1]*/
#ifdef PASA
                    pasacom->id = id ;
                    pasacom->SSinitDone = TRUE ;
#else
                    cgcom->id = id ;
#endif
                    /* end of MUMPS initialization */
#endif
#ifdef USE_HARWELL
                    lblfile = fopen ("lbl.log", "w") ;
                    cmm = sopt_malloc
                               (&status, (SOPTINT) 1, sizeof (cholmod_common)) ;
                    CHOLMOD (start) (cmm) ;
                    Cmember = Cparent = NULL ;
                    Achol = sopt_malloc
                               (&status, (SOPTINT) 1, sizeof (cholmod_sparse)) ;
                    Achol->stype = 1 ;
                    Achol->packed = TRUE ;
                    Achol->sorted = FALSE ;
                    Achol->itype = CHOLMOD_INT ;
                    Achol->dtype = CHOLMOD_DOUBLE ;
                    /* no numerical value since only using CHOLMOD for order */
                    Achol->xtype = CHOLMOD_PATTERN ;
                    Achol->p = NULL ;  /* matrix is not malloc'd */
#ifdef PASA
                    pasacom->lblfile = lblfile ;
                    pasacom->cmm = cmm ;
                    pasacom->Cmember = pasacom->Cparent = NULL ;
                    pasacom->Achol = Achol ;
#else
                    cgcom->lblfile = lblfile ;
                    cgcom->cmm = cmm ;
                    cgcom->Cmember = cgcom->Cparent = NULL ;
                    cgcom->Achol = Achol ;
#endif
                    /* end of HARWELL initialization */
#endif
                }
            }
            /* find_max_iter = TRUE => find max number of iterations in CG
                                       before switch to SS */
            /* dCG = TRUE => initially compute the Newton search
                             direction by minimizing the quadratic
                             model using PRP+ */
            /* The CG version of the Hessian-based algorithm employs a
               compressed matrix C which eliminates bound variables and
               active variables. This compression setup only needs to be
               done in the first iteration. Thereafter, at new iterates,
               we simply map the Hessian into C. */
            if ( dCG && nCG == 0 ) /* dCG = TRUE => CG search direction */
            {
                sopt_convert_H_to_C_in_CG (C, H, !nCG, location) ;

                /* if diagonal perturbation used in Newton CG, compute Cmax */
                if ( Cmax < CGZERO )
                {
                    /* if Hmax has no value, then compute Cmax from scratch */
                    if ( Hmax <= CGZERO )
                    {
                        Cmax = sopt_sup_normx (C->x, C->p [n]) ;
                        if ( Hmax < CGZERO ) /* Hmax is used */
                        {
                            Hmax = Cmax*Parm->ss_diag_pert ;
                        }
                        Cmax *= Parm->cg_diag_pert ;
                    }
                    else /* Hmax already set, express Cmax in terms of Hmax */
                    {
                        Cmax = Parm->cg_diag_pert*(Hmax/Parm->ss_diag_pert) ;
                    }
                }
                /* otherwise, either Cmax was already evaluated or diagonal
                   perturbation not used (Cmax = 0) */
#ifndef NDEBUG
                XXCG(debug_C_in_CG) (C, H) ;
#endif
            }
            /* find_max_iter = TRUE means that we will compute the number
               of CG iterations that are equivalent to a single SS iterate.
               This yields the switching point from CG to SS. We only need
               to compute this if we are going employ the Newton CG iteration
               (that is, dCG = TRUE). If SS previously failed, then max_iter
               is set to n. */
            if ( dCG && find_max_iter && !(*SSfail) )
            {
                /* estimate flops for an iteration of CG */
                SOPTFLOAT cg_flops = C->p [n] + 6*n ;
#ifdef PASA
                if ( Aexists )
                {
                    if ( use_napheap ) cg_flops += 6*n ; /* 2 projections */
                    else
                    {
                        /* two projection and 2*Lnnz ops for each solve */
                        cg_flops += 4*pasacom->ppcom->Work->Lnnz ;
                    }
                }
#endif

                /* setup AT, first determine number of nonzeros */
                Annz = 0 ;  /* nnz's in A */
                Cdim = ncol ; /* total dimension of C */
#ifdef PASA
                if ( Aexists )
                {
                    if ( use_napheap )     /* napsack constraint exists */
                    {
                        if ( nap_present ) /* napsack constraint active */
                        {
                            Cdim++ ;
                            Annz = pasacom->nap_nnz ;
                        }
                    }
                    else /* linear constraint present */
                    {
                        Cdim += nrow ;
                        Annz = pasacom->Copy->ATp [nrow] ;
                    }
                }
#endif
                /* The triples needed by the symmetric solver are extracted
                   from the H-matrix and stored in the C-matrix, which
                   becomes the input for the symmetric solver. */
                sopt_convert_H_to_C_in_SS (C, H, !nSS, Cdim, Annz, location) ;

                /* if diagonal perturbation of Hessian is used, compute Hmax
                   if not yet computed */
                if ( Hmax < CGZERO ) /* Hmax not yet computed */
                {
                    if ( Cmax <= CGZERO ) /* Cmax not yet computed */
                    {
                        Hmax = sopt_sup_normx (C->vals, C->nnz) ;
                        if ( Cmax < CGZERO )
                        {
                            Cmax = Hmax*Parm->cg_diag_pert ;
                        }
                        Hmax *= Parm->ss_diag_pert ;
                    }
                    else /* express Hmax in terms of Cmax */
                    {
                        Hmax = Parm->ss_diag_pert*(Cmax/Parm->cg_diag_pert) ;
                    }
                }
                /* otherwise, either Hmax was already evaluated or diagonal
                   perturbation not used (Hmax = 0) */
#ifndef NDEBUG
                XXCG(debug_C_in_SS) (C, H) ;
#endif

                /* Insert the A matrix in C if it exists. Space was left for
                   A when the triples arrays were malloc'd in the
                   convert routine. This only needs to be done in the
                   first iteration if the Hessian sparsity is fixed. It needs
                   to be done in every iteration when the Hessian sparsity
                   is not fixed. */
                Amax = CGZERO ;
                li = C->nnz ;
                j  = n ;
#ifdef PASA
                if ( use_napheap ) /* napsack constraint exists */
                {
                    if ( nap_present ) /* napsack constraint active */
                    {
                        j = n + 1 ; /* fortran indexing */
                        for (i = 0; i < n; i++)
                        {
                            CGFLOAT const ax = pasacom->nap_a [i] ;
                            if ( ax )
                            {
                                if ( fabs (ax) > Amax ) Amax = fabs (ax) ;
                                C->vals [li] = ax ;
                                C->rows [li] = i + 1 ;
                                C->cols [li] = j ;
                                li++ ;
                            }
                        }
                    }
                }
                else if ( Aexists )
                {
                    PPwork *W = pasacom->ppcom->Work ;
                    PASAINT   *ir  = W->ir ;
                    PASAINT   *ATp = W->ATp ;
                    PASAINT   *ATi = W->ATi ;
                    PASAFLOAT *ATx = W->ATx ;
                    for (i = 0; i < nrow; i++)
                    {
                        if ( ir [i] == 0 )
                        {
                            /* Fortran indexing so column is n + 1 */
                            j++ ;
                            PASAINT const q = ATp [i+1] ;
                            for (p = ATp [i]; p < q; p++)
                            {
                                CGFLOAT const ax = ATx [p] ;
                                C->vals [li] = ax ;
                                if ( fabs (ax) > Amax ) Amax = fabs (ax) ;
                                C->rows [li] = ATi [p] + 1 ;
                                C->cols [li] = j ;
                                li++ ;
                            }
                        }
                    }
                }
                Amax *= Parm->ss_diag_pert ;
#endif
                Anrow = j - n ; /* number of active rows */

                /* total dim, Hessian plus active rows in A */
                nt = j ;

                /* total nnz, Hessian nnz + active rows + diag regularization */
                nnzt = li + nt ;

                /* Insert the diagonal regularization terms into C */
                cg_initx (C->vals+li, Hmax, n) ;
                cg_initx (C->vals+(li+n), -Amax, Anrow) ;
                row = C->rows+(li-1) ;
                col = C->cols+(li-1) ;
                for (i = 1; i <= nt; i++)
                {
                    row [i] = i ;
                    col [i] = i ;
                }

#ifdef USE_MUMPS
                id->n   = nt ;     /* total n: active rows + Hessian */
                id->nnz = nnzt ;   /* Hessian nnz + active rows + diag reg. */
                id->irn = C->rows ;
                id->jcn = C->cols ;
                id->a   = C->vals ;

                /* factorization analysis */
                id->job = 1 ;
                dmumps_c (id) ;
                if ( id->info [0] < 0 )
                {
                    *SSfail = TRUE ;
#ifdef PASA
                    speed [NewtonSS] = -CGINF ;
#endif
                }
                else
                {
                    /* operations in factorization */
                    ss_flops = id->rinfog [0] ;
                }
                /* number of entries in factors: id->infog [19]) */
                /* end of MUMPS initialization */
#endif

#ifdef USE_HARWELL
                if ( lbl != NULL )
                {
                    lbl->irn = lbl->ma57->irn = NULL ;
                    lbl->jcn = lbl->ma57->jcn = NULL ;
                    LBL_Finalize (lbl) ;
                }
                lbl = LBL_Initialize (nnzt, nt, lblfile, MA57_SOLVER) ;
                /* no scaling */
                /*LBL_set_int_parm(lbl, LBL_I_SCALING, 0) ;*/
                /* provide ordering = 1,
                   AMD pivoting     = 2,
                   AMD/ND           = 5 */
                LBL_set_int_parm(lbl, LBL_I_PIV_SELECTION, 1) ;
                /* do pivoting */
                LBL_set_int_parm(lbl, LBL_I_PIV_NUMERICAL, 1) ;
                /* pivot threshold */
                /*LBL_set_real_parm(lbl, LBL_D_PIV_THRESH, 0.5) ;*/

                /* MA57 malloc's space for rows and cols, free this space
                   since we already have space in C */
                cg_free (lbl->ma57->irn) ;
                cg_free (lbl->ma57->jcn) ;
                lbl->irn = lbl->ma57->irn = C->rows ;
                lbl->jcn = lbl->ma57->jcn = C->cols ;

#ifdef PASA
                pasacom->lbl = lbl ;
#else
                cgcom->lbl = lbl ;
#endif
                /* CHOLMOD is used to find the ordering for the rows and
                   columns of the matrix. The ordering is stored in
                   lbl->ma57->keep. Begin by building the matrix as a CHOLMOD
                   sparse matrix. */
                if ( Achol->p == NULL ) /* malloc cholmod_sparse matrix */
                {
                    /* number of column pointers = max #cols + max #rows + 1 */
                    i = ncol + nrow + 1 ;
                    Achol->p = (CGINT *)   cg_malloc (&status, i, sI) ;
                    i = Annz + C->nnzsym ; /* most nnz in symmetric matrix */
                    Achol->i = (CGINT *)   cg_malloc (&status, i, sI) ;
                    Achol->x = NULL ;
                    if ( status )
                    {
                        status = sopt_convert_error (location, status) ;
                        XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                        return (status) ;
                    }
                }
                Achol->nrow = nt ;
                Achol->ncol = nt ;
                /* compute column counts */
                CGINT *cc = Achol->p ;
                cg_initi (cc, (CGINT) 0, nt+1) ;
                row = C->rows ;
                col = C->cols ;
                for (i = 0; i < li; i++)
                {
                    /* Fortran indexing for row and col, diagonal not needed
                       when finding the ordering */
                    if ( row [i] != col [i] ) cc [col [i]]++ ;
                }

                /* evaluate the column pointers */
                for (i = 0; i < nt; i++) cc [i+1] += cc [i] ;
                /* scatter the matrix elements in Achol */
                CGINT *Ai = Achol->i ;
                for (i = 0; i < li; i++)
                {
                    /* Fortran indexing for row and col, skip diagonal */
                    if ( row [i] != col [i] )
                    {
                        CGINT const ccj = cc [col [i]-1]++ ;
                        Ai [ccj] = row [i] - 1 ;
                    }
                }
                /* shift cc up by one and cc [0] = 0 */
                for (j = nt; j > 0; j--) cc [j] = cc [j-1] ;
                cc [0] = 0 ;
/*sopt_printAMATLAB (nt, cc, Ai, NULL, NULL, "Achol") ;*/

                CGINT *cholperm = lbl->ma57->keep ; /* column permutation */
                /* decide between AMD and Metis based on graph density */
                t = (double) nt ;
                t = 0.4*t*(t+1) ;
                if ( cc [nt] > t )
                {
                    if ( PrintLevel )
                    {
                        printf ("use AMD in Hessian-based cg_descent\n") ;
                    }
                    /* use AMD */
                    CHOLMOD (camd) (Achol, NULL, 0, NULL, cholperm, cmm);
                }
                else /* use Metis */
                {
                    if ( PrintLevel )
                    {
                        printf ("use Metis in Hessian-based cg_descent\n") ;
                    }
                    if ( Cparent == NULL )
                    {
#ifdef PASA
                        Cmember = cg_malloc (&status, ncol + nrow, sI) ;
                        Cparent = cg_malloc (&status, ncol + nrow, sI) ;
                        pasacom->Cparent = Cparent ;
                        pasacom->Cmember = Cmember ;
#else
                        Cmember = cg_malloc (&status, ncol, sI) ;
                        Cparent = cg_malloc (&status, ncol, sI) ;
                        cgcom->Cparent = Cparent ;
                        cgcom->Cmember = Cmember ;
#endif
                    }
                    cmm->current = 0 ;
                    cmm->method [0].nd_components = 1 ;
                    CHOLMOD (nested_dissection) (Achol, NULL, 0,
                                       cholperm, Cparent, Cmember, cmm) ;
                    /* CHOLMOD (metis) (Achol, NULL, 0, TRUE, cholperm, cmm) ;*/
                }
                /* convert permutation to Fortran format */
                for (i = 0; i < nt; i++) cholperm [i]++ ;

/* FILE *FP ;
long II ;
FP = fopen ("/tmp/MUMPS_perm.txt", "r") ;
for (i = 0; i < nt; i++)
{
    if ( fscanf (FP, "%ld\n", &II) == 1 ) cholperm [i] = II+1 ;
printf ("i: %i perm: %li\n", i, II) ;
}
fclose (FP) ; */

                /* iflag = 0 => automatic pivot choice in ma27 */
                status = LBL_Analyze(lbl, 0) ;
                if ( PrintLevel )
                {
                    printf ("Forecast number of reals in factors:     %i\n",
                            lbl->ma57->info [4]) ;
                    printf ("Forecast number of integers in factors:  %i\n",
                            lbl->ma57->info [5]) ;
                    printf ("Forecast max frontal size:               %i\n",
                            lbl->ma57->info [6]) ;
                    printf ("Number of nodes in assembly tree:        %i\n",
                            lbl->ma57->info [7]) ;
                    printf ("Min length of fact (no compression):     %i\n",
                            lbl->ma57->info [8]) ;
                    printf ("Min length of ifact (no compression):    %i\n",
                            lbl->ma57->info [9]) ;
                    printf ("Min length of fact (with compression):   %i\n",
                            lbl->ma57->info [10]) ;
                    printf ("Min length of ifact (with compression):  %i\n",
                            lbl->ma57->info [11]) ;
                    printf ("Forcast number of flops for elimin:      %e\n",
                            lbl->ma57->rinfo [1]) ;
                }

#if 0
                if ( PrintLevel )
                {
                    printf ("Number of entries in factors:           %i\n",
                            lbl->ma57->info [13]) ;
                    printf ("Storage for real data in factors:       %i\n",
                            lbl->ma57->info [14]) ;
                    printf ("Min length of fact:                     %i\n",
                            lbl->ma57->info [16]) ;
                    printf ("Min length of ifact:                    %i\n",
                            lbl->ma57->info [17]) ;
                    printf ("Min length of fact without compression :%i\n",
                            lbl->ma57->info [18]) ;
                    printf ("Min length of ifact without compression:%i\n",
                            lbl->ma57->info [19]) ;
                    printf ("Order of largest frontal matrix:        %i\n",
                            lbl->ma57->info [20]) ;
                    printf ("Number of 2x2 numerical pivots:         %i\n",
                            lbl->ma57->info [21]) ;
                    printf ("Number of fully-summed variables:       %i\n",
                            lbl->ma57->info [22]) ;
                    printf ("Number of negative eigenvalues:         %i\n",
                            lbl->ma57->info [23]) ;
                    printf ("Rank of the factorization:              %i\n",
                            lbl->ma57->info [24]) ;
                }
#endif
                if( status < 0 )
                {
                    printf("Error %i returned by Harwell MA57 in analysis.\n",
                            status);
                    printf ("Continue solution process using CG to solve\n"
                            "the quadratic model.\n") ;
                    *SSfail = TRUE ;
                    ss_flops = CGINF ;
#ifdef PASA
                    speed [NewtonSS] = -CGINF ;
#endif
                }
                else if ( status >= 0 )
                {
                    if ( status > 0 )
                    {
                        printf("Warning %i returned by Harwell MA57 in "
                               "analysis.\n", status);
                        printf ("Continue solution process (see MA57\n"
                                "documentation for explanation of error)\n") ;
                    }
                    ss_flops = lbl->ma57->rinfo [1]/2 ;
                }
/* end of calls to MA57 */
#endif
                if ( *SSfail )
                {
                    *max_d_iter = n ;
                    dCG = TRUE ;
                    if ( PrintLevel )
                    {
                        printf ("In Hessian-based cg_descent, factorization\n"
                                "of Hessian fails; will use an iterative\n"
                                "solver to compute the Newton search "
                                "direction.\n") ;
                    }
                }
                else
                {
                    *max_d_iter = 1.5*ss_flops/cg_flops ;
                    if ( *max_d_iter > n )
                    {
                        *max_d_iter = n ;
                    }
                    else /* max_d_iter  <= n */
                    {
                        if ( n <= 100 ) *max_d_iter = n ;
                    }
                    if ( *max_d_iter <= 3 ) *max_d_iter = 3 ;
                }
                if ( PrintLevel )
                {
                    printf ("max_d_iter: %li ss_flops: %e cg_flops: %e\n",
                            (LONG) *max_d_iter, ss_flops, cg_flops) ;
                }
                find_max_iter = FALSE ;
                SSsetupDone = TRUE ;
            }
            else if ( *SSfail )
            {
                *max_d_iter = n ;
                find_max_iter = FALSE ;
            }

            /* DCG */
            if ( dCG )
            {
                ASSERT (NewMethod == NewtonCG) ;
                /* we start the NewtonCG timer here, but this misses the
                   time htic associated with the Hessian evaluation; we
                   save that time in htic and backup to include this time */
                if ( SetTimer )
                {
                    method_tic = sopt_timer () - htic ;
                    SetTimer = FALSE ;
                }
/*  --------------  MALLOC MEMORY FOR CG if nCG = 0 ------------------------  */
                /* map Hessian from H into C if nCG > 0 (nCG = 0 done above) */
                if ( nCG )
                {
                    /* Map the current Hessian into C. The setup of C was done
                       above. For a quadratic, this step can be skipped since
                       the Hessian does not depend on x. */
                    if ( !QuadCost )
                    {
                        sopt_convert_H_to_C_in_CG (C, H, FALSE, location) ;
                    }
#ifndef NDEBUG
                    XXCG(debug_C_in_CG) (C, H) ;
#endif
                }
                else /* nCG = 0, initial iterate, allocate work memory */
                {
                    /* allocate memory for CG iteration (Hwork) */
#ifdef PASA
                    if ( C->Hwork == NULL )
                    {
                        /* malloc space for Hd, PHd, r, and v (all ncol) */
                        C->Hwork =  pasa_malloc (&status, 5*ncol, sF) ;
                        if ( status )
                        {
                            status =sopt_convert_error(location,status) ;
                            XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                            return (status) ;
                        }
                    }
                    Hwork = C->Hwork ;
                    PHd = Hwork ;  Hwork += n ;
#else
                    if ( C->Hwork == NULL )
                    {
                        /* malloc space for Hd, PHd, r, and v (all ncol) */
                        C->Hwork = cg_malloc (&status, 4*ncol, sF) ;
                        if ( status )
                        {
                            status = sopt_convert_error (location, status) ;
                            XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                            return (status) ;
                        }
                    }
                    Hwork = C->Hwork ;
#endif
                    Hd = Hwork ;  Hwork += n ;
                    r  = Hwork ;  Hwork += n ; /* gradient */
                    v  = Hwork ;  Hwork += n ; /* search direction */
                    /* PH*gproj is stored in gH */
                    gH = Hwork ;  Hwork += n ;
/*  ----------------------  END OF MALLOC   --------------------------------  */
                }

/*  ---------------   CG INITIALIZATION   ----------------------------------  */
                /* rnorm is used in stopping condition for CG */
#ifdef PASA
                rnorm = pasacom->e ;
#else
                rnorm = gnorm ;
#endif
                rnorm_start = rnorm ;
                rtr = CGZERO ;
                max_d_iter_bump = TRUE ;
                for (j = 0; j < n; j++)
                {
                    t = gproj [j] ;
                    r [j] = t ;
                    rtr += t*t ;
                    v [j] = -t ;
                }
                rtv = rtr ;
                CGFLOAT start_rtr = rtr ;
                cg_initx (D, CGZERO, n) ;
                dtol = Hessian_err_decay*rnorm ;

                d_iter = 0 ;
                NegCurv = FALSE ;
                CGFLOAT dcost = CGZERO ;
/*  ---------------   END CG INITIALIZATION   ------------------------------  */
#ifdef TIME
cgcom->loop_init += sopt_timer () - tic ;
#endif

/*  ---------------------   START CG PRP+ ITERATION   ----------------------  */
               /* while loop breaks when NegCurv, d_iter = max_d_iter, or
                  objective gradient is sufficiently small */
#ifdef TIME
tic = sopt_timer () ;
#endif
                /* #NCG */
                while ( 1 )
                {
                    d_iter++ ;
                    /* need to compute PHP * search direction v where P is
                       either the null space projection if Aexists or I.
                       pv = Pv or v */
                    pv = v ;
#ifdef PASA
                    /* start computation of PHP*v, v = search direction */
                    if ( Aexists )
                    {
                        /* d is used temporarily for storing the projection */
                        pasa_null_project (d, v, NULL, FALSE, pasacom) ;
                        pv = d ;
                    }
#endif
                    /* multiply pv by H and store in Hd */
                    cg_builtin_hprod (Hd, pv, n, C->p, C->i, C->x) ;
                    /* compute d*PHP*d, pv = P*v if Aexists */
                    dHd = cg_dot (pv, Hd, n) ;
                    if ( dHd > CGZERO ) /* objective locally convex */
                    {
                        /* in first iteration, store pd in gH for line search */
                        if ( d_iter == 1 )
                        {
#ifdef PASA
                            if ( Aexists )
                            {
                                /* compute P*H*Pv = P*Hd and store in gH for
                                   trust region computation */
                                pasa_null_project (gH, Hd, NULL, FALSE,pasacom);
                            }
                            else /* copy Hd to gH */
#endif
                                sopt_copyx (gH, Hd, n) ;
                                
                            /* gH = PHP*d = (PHP*(-g), which is needed for trust
                               region */
                            dHdsave = dHd ;
                        }

                        /* apply regularization */
                        if ( Cmax )
                        {
                            CGFLOAT const sx = Cmax ;
                            t = CGZERO ;
                            for (i = 0; i < n; i++)
                            {
                                CGFLOAT const pvi = pv [i] ;
                                Hd [i] += sx*pvi ;
                                t += pvi*pvi ;
                            }
                            dHd += sx*t ;
                        }

                        step = rtv/dHd ;
                        /* QP cost improvement due to step */
                        if ( PrintLevel )
                        {
                            dcost -= 0.5*rtr*step ; /* = -0.5||g||^4/dHd */
                        }
                        /* D = Pd is variable in transformed space */
                        cg_daxpy (D, v, step, n) ;
#if 0
SOPTFLOAT *XX = malloc (n*sizeof(SOPTFLOAT)) ;
cg_builtin_hprod (XX, D, n, C->p, C->i, C->x) ;
t = sopt_dot (D, gproj, n) + 0.5*sopt_dot (XX, D, n) ;
printf ("cost at D: %e\n", t) ;
cg_builtin_hprod (XX, gproj, n, C->p, C->i, C->x) ;
t = sopt_dot (gproj, gproj, n) ;
t = -0.5*t*t/sopt_dot (XX, gproj, n) ;
printf ("Cauchy: %e\n", t) ;
#endif

                    }
                    else /* negative curvature */
                    {
                        if ( PrintLevel )
                        {
                            printf ("negative curvature %e in CG-Hessian: "
                                    "d'g = %e\n", dHd, sopt_dot (Hd, gproj, n));
                        }
                        NegCurv = TRUE ;
                        break ;
                    }

                    /* D has been updated, break if d_iter >= max_d_iter */
                    if ( (d_iter >= *max_d_iter) )
                    {
                        /* if NewtonSS failed, then increase max_d_iter and
                           continue the iteration, but never beyond n */
                        if ( *SSfail )
                        {
                            /* if both NewtonCG and NewtonSS fail, then
                               use gradient descent */
                            if ( d_iter >= n )
                            {
                                use_hessian = FALSE ;
                                iter_check = INT_MAX ;
                                iter_since_err_update = 0 ;
                                deriv_mode = 1 ;
                                AutoChoice = FALSE ;
                                *CGNfail = TRUE ;
#ifdef PASA
                                pasacom->deriv_mode = 1 ;
#endif
                                if ( PrintLevel )
                                {
                                    printf ("\n\n## Both the symmetric solver "
                                            "NewtonSS and NewtonCG routines "
                                            "have failed.\nSwitch to ordinary "
                                            "cg_descent. ##\n") ;
                                }
                                break ;
                            }
                            else /* increase max_d_iter and continue the iter */
                            {
                               if ( *max_d_iter > 4 )
                               {
                                   *max_d_iter *= 1.3 ;
                                   if ( *max_d_iter > n ) *max_d_iter = n ;
                               }
                               else *max_d_iter = 4 ;
                           }
                        }
                        else /* SS did not fail */
                        {
                            /* allow one increase in max_d_iter if NewtonCG
                               was faster, otherwise switch to NewtonSS */
                            if ( (speed [NewtonSS] > -CGINF) &&
                                 (speed [NewtonCG] > speed [NewtonSS]) &&
                                 (*max_d_iter < n) && (max_d_iter_bump) )
                            {
                               if ( *max_d_iter > 4 )
                               {
                                   *max_d_iter *= 1.3 ;
                                   if ( *max_d_iter > n ) *max_d_iter = n ;
                               }
                               else *max_d_iter = 4 ;
                               max_d_iter_bump = FALSE ;
                            }
                            else /* otherwise switch to NewtonSS */
                            {
                                *CGNfail = TRUE ;
                                dCG = FALSE ;
                                if ( PrintLevel )
                                {
                                    printf ("\n## Iterations in CG-Hessian "
                                            "reach limit %li switch to "
                                            "symmetric solver ##\n\n",
                                            (LONG) *max_d_iter) ;
                                }
                                break ;
                            }
                        }
                    }

                    /* update r using a vector pd which is either Hd or PHd */
                    pd = Hd ;
#ifdef PASA
                    if ( Aexists )
                    {
                        /* compute P*H*Pv = P*Hd */
                        pasa_null_project (PHd, Hd, NULL, FALSE, pasacom) ;
                        pd = PHd ;
                    }
#endif
                    /* update gradient by step*matrix*d:
                       g^new = g^old + step*A*d
                       t = gnew'*(gnew-gold)
                       rtr = g^new'*g^new
                       rnorm = max (abs (rnew)) */
                    t = CGZERO ;
                    CGFLOAT const rtr_old = rtr ;
                    rtr = CGZERO ;
                    rnorm = CGZERO ;
                    CGFLOAT const Step = step ;
                    for (i = 0; i < n; i++)
                    {
                        CGFLOAT const pds = pd [i]*Step ;
                        CGFLOAT const rnew = r [i] + pds ;
                        t += rnew*pds ;
                        r [i] = rnew ;
                        rtr += rnew*rnew ;
                        if ( fabs (rnew) > rnorm ) rnorm = fabs (rnew) ;
                    }
                    if ( (rnorm <= dtol) && (d_iter > 1) )
                    {
                        if ( PrintLevel )
                        {
                            printf ("NewtonCG tolerance reached: rnorm %e "
                                    "dtol %e\n", rnorm, dtol) ;
                        }
                        break ;
                    }
                    if ( rtr >= start_rtr*Parm->big_grad )
                    {
                        if ( PrintLevel )
                        {
                            printf ("rtr: %e start_rtr: %e break and treat "
                                    "the current iterate as minimizer of QP "
                                    "in CG-Hessian\n", rtr, start_rtr);
                        }
                        /* NewtonStep = TRUE by default */
                        break ;
                    }
                    /* PRP+ CG:  beta_k = max (gnew'(gnew - g)/g'g, 0) */
                    t = CGMAX (t, CGZERO) ;
                    beta_k = t/rtr_old ;
                    /* if jamming, restart CG */
                    if ( beta_k == CGZERO )
                    {
                        rtv = rtr ;
                        cg_scale (v, r, -CGONE, n) ;
                        /* dCG = FALSE ;*/ /* use symmetric solver next time */
                    }
                    else
                    {
                        rtv = CGZERO ;
                        for (j = 0; j < n; j++)
                        {
                            v [j] = -r [j] + beta_k*v [j] ;
                            rtv -=  v [j]*r [j] ;
                        }
                    }
                }
#ifdef TIME
cgcom->loop_prp += sopt_timer () - tic ;
#endif
/*  -----------------------   END CG PRP_ ITERATION   ----------------------  */
#ifdef TIME
tic = sopt_timer () ;
#endif
                Stat->PRP += d_iter ; /* update total number of PRP+ iter */
                if ( PrintLevel )
                {
                    printf ("number PRP+ iter: %li total: %li quad dcost: %e"
                            "\n\n", (LONG) d_iter, (LONG) Stat->PRP, dcost) ;
                }
                /* Perform a search over solutions to a subspace-restricted
                   problems of the following form:

                   (H) min q(z) = (1/2) z'PHPz + (Pg)'z
                       s. t. ||z|| <= rho, z in S,

                   where S is the subspace. S and the initial
                   rho are chosen so that D, the approximate Newton step,
                   is feasible. We compare the function value change
                   f(x) - f(x+z) to the optimal value
                   H* for the quadratic (H).  If the ratio

                   (tt)    [f(x) - f(x+z)]/|H*| >= tau,

                   where 0 < tau < 1 is the acceptance threshold, then
                   we accept the step z. Otherwise, we continue to
                   decrease rho until (tt) is satisfied. As rho tends to
                   zero, the ratio (tt) should approach 1 because the
                   quadratic term in (H) becomes negligible when compared
                   compared to the linear term as as rho tends to zero.
                   Consequently, z approaches the direction gproj = Pg.
                   If the solution z* of (H) is nearly a multiple of gproj,
                   then exit Hessian-based pasa and return to ordinary
                   cg_descent.

                   The subspace S includes gproj and D. Moreover, if the
                   final search direction v is a direction of negative
                   curvature, then it is also inserted in S. Let V be
                   the n by m matrix of vectors spanning the subspace, where
                   m is currently at most 3. Column 1 of V contains gproj,
                   column 2 contains v (when NegCurv exists) or D
                   (if NegCurv does not exist), and column 3 contains D
                   or it is empty. If Z is a matrix whose columns are
                   an orthonormal basis for the columns of V, then we
                   replace (H) by the problem:

                   (H') min  (1/2)y'My +  (Pg)'Zy, ||y|| <= rho

                   where M = Z'*PHP*Z and y lies in R^m, m <= 3. To avoid
                   unnecessary multiplications between the columns of Z
                   and H, the products PHP*gproj and PHP*v (when v is a
                   direction of negative curvature) were saved.
                   This enables the formation of PHP*Z(:, 0) and
                   PHP*Z(:, 1) using the saved products.  */

                if ( PrintLevel )
                {
                    printf ("d_iter: %li NegCurv: %i\n",
                                                      (LONG) d_iter, NegCurv) ;
                }
                if ( QuadCost && !NegCurv )
                {
                    /* for a QP and no negative curvature, do a Newton step
                       along D, note that projection is still needed */
#ifdef PASA
                    if ( Aexists )
                    {
                        pasa_null_project(d, D, NULL, FALSE, pasacom) ;
                    }
#endif
                }
                else if ( d_iter == 1 && NegCurv )
                {
                    /* only option is a search along the negative gradient */
                    NewtonStep = FALSE ;
                    cg_scale (D, gproj, -CGONE, n) ;
#ifdef PASA
                    if ( Aexists )
                    {
                        pasa_null_project(d, D, NULL, FALSE, pasacom) ;
                    }
#endif
                }
                else while ( 1 )
                {
                    /* CGTRUST */
                    /* compute search direction using trust region approach */
                    cg_copy (V, gproj, n) ;
                    /* m stores number of columns of V used in trust region */
                    if ( d_iter == 1 )
                    {
                        cg_copy (V+n, D, n) ; /*    D exists */
                        m = 2 ;
                    }
                    else /* d_iter >= 2 */
                    {
                        /* if NegCurv present, D = previous direction */
                        if ( NegCurv )
                        {
                            cg_copy (V+n, v, n) ;
                            cg_copy (V+(n+n), D, n) ;
                            m = 3 ;
                        }
                        else
                        {
                            cg_copy (V+n, D, n) ;
                            m = 2 ;
                        }
                    }
                    k = cg_basis (Z, T, V, perm, PrintLevel, n, m,
                                  1.e-8, basisWork) ;
                    if ( PrintLevel )
                    {
                        printf ("trust region basis vectors: %li\n", (LONG) k) ;
                    }
                    /* if there is only one independent vector, only option
                       is negative gradient search */
                    if ( k == 1 )
                    {
                        NewtonStep = FALSE ;
                        /* negative gradient = search direction */
                        cg_scale (D, gproj, -CGONE, n) ;
#ifdef PASA
                        if ( Aexists )
                        {
                            pasa_null_project(d, D, NULL, FALSE, pasacom) ;
                        }
#endif
                        break ;
                    }
                    /* otherwise there are at least 2 independent vectors, one
                       is gproj
                       NOTE: T denotes the elements of the triangular matrix
                             generated in the basis computation.
                       p0 = PHP*Z(:, 0) = T(0)*PHP*V(:, 0)
                       if perm [1] = 1, then
                       p1 = PHP*Z(:, 1) = T(1)*PHP*V(:, 0) + T(2)*PHP*V(:, 1)
                           NOTE: PHPV(:, 0) = -gH
                       otherwise
                       p1 = PHP*Z(:, 1) (compute from scratch)
                       p2 = PHP*Z(:, 2) (if it exists, compute from scratch) */

                    /* form the matrix M of problem (H') starting with products
                       between columns of Z and PHP*Z(:, 0) */

                    Mx [0] = dHdsave*T [0]*T [0] ;
                    Z1 = Z+n ;
                    Z2 = Z1+n ;
                    Mx [1] = Mx [k] = -T [0]*cg_dot (Z1, gH, n) ;
                    if ( k > 2 ) /* k = 3 */
                    {
                        Mx [2] = Mx [6] = -T [0]*cg_dot (Z2, gH, n) ;
                        /* the first row and column of Mx have been computed
                           next, evaluate pd = PHP*Z(:, 1)
                           use stored pd & gH, use p1 formula */
                        if ( perm [1] == 1 )
                        {
                            t = -T [1]/T [2] ;
                            /* If A exists, then project Hd, store result in
                               pd. We do this since Hd = H*P*v when we need
                               P*H*P*v where v stores NegCurv direction */
                            pd = Hd ;
#ifdef PASA
                            if ( Aexists )
                            {
                                /* compute P*H*Pv = P*Hd */
                                pasa_null_project
                                           (PHd, Hd, NULL, FALSE, pasacom) ;
                                pd = PHd ;
                            }
#endif
                            cg_daxpy (pd, gH, t, n) ; /* pd += t*gH */
                            cg_scale (pd, pd, T [2], n) ;
                        }
                        else /* directly compute pd = PHP*Z1 */
                        {
                            pd = XXCG(PHPz) (Z1, r, v, Aexists, FALSE, C_PASA) ;
                        }
                        Mx [4] = cg_dot (pd, Z1, n) ;

                        /* pd = PHP*Z1 */
                        Mx [5] = Mx [7] = cg_dot (Z2, pd, n) ;
                        /* directly compute pd = PHP*Z2 */
                        pd = XXCG(PHPz) (Z2, r, v, Aexists, FALSE, C_PASA) ;
                        Mx [8] = cg_dot (pd, Z2, n) ;
                    }
                    else /* k = 2, directly compute pd = PHP*Z1 */
                    {
                        pd = XXCG(PHPz) (Z1, r, v, Aexists, FALSE, C_PASA) ;
                        Mx [2] = Mx [1] ;
                        Mx [3] = cg_dot (pd, Z1, n) ;
                    }

#if 0
pd = XXCG(PHPz) (D, Newr, Newv, Aexists, FALSE, C_PASA) ;
printf ("dcost should be: %e\n", 0.5*cg_dot(D, pd, n) + cg_dot(gproj, D, n)) ;

/* check quadratic, dcost versus value of quadartic associated with D */
CGFLOAT xD [3] ;
for (j = 0; j < k; j++)
{
    xD [j] = cg_dot (D, Z+(j*n), n) ;
}
CGFLOAT ncost ;
ncost = 0 ;
for (j = 0; j < k; j++)
{
    ncost += xD [j]*cg_dot (xD, Mx+(j*k), k) ;
}
ncost = 0.5*ncost + cg_dot (Mb, xD, k) ;
printf ("XX k: %i ncost: %e\n", k, ncost) ;
#endif

                    /* create Mi and Mp since input should be sparse */
                    CGINT *MiPointer = Mi ;
                    m = 0 ;
                    for (j = 0; j < k; j++)
                    {
                        Mp [j] = m ;
                        m += k ;
                        for (i = 0; i < k; i++)
                        {
                            MiPointer [i] = i ;
                        }
                        MiPointer += k ;
                    }
                    Mp [k] = m ;
                    /* build the linear term in (H') */
                    for (i = 0; i < k; i++)
                    {
                        Mb [i] = cg_dot (gproj, Z+(i*n), n) ;
                    }

                    /* Use SSM to minimize the quadratic over the sphere of
                       radius rho, We want D to be included in the sphere so
                       we choose the radius to be the length of D */
                    rho = sqrt (cg_dot (D, D, n)) ;
                    /* If negative curvature occurred, then rho is increased. */
                    if ( NegCurv ) rho *= 10 ;

                    /* the current function value at x is stored in f */
                    CGFLOAT tbest = CGINF ;
                    int t_iter ;
                    /* In the loop below, we compute the minimizer of the
                       quadratic model over a sphere of radius rho. We
                       then evaluate the objective function at xk plus
                       the minimizer of the quadratic. If this point
                       generates a reduction in the objective value, then
                       we break and use this point in a line search.
                       If this point increases the objective value, then
                       we continue to reduce rho until finding a point
                       where either the objective increases or the decrease in
                       the objective is at least a small fraction of the
                       decreased anticipated from the quadratic model.
                       The point used in the CG line search is always the
                       point with the smallest objective value. Since we
                       do not save the points where the objective is evaluated,
                       the variable restore_x let us know whether the
                       current evaluation point has the smallest objective
                       value. restore_x is TRUE if we need to backtrack to
                       a prior evaluation point, while restore_x is FALSE
                       if the current evaluation point has the smallest
                       objective value. */
                    int restore_x = FALSE ;
                    for (t_iter = 0; t_iter < TrustIterLimit; t_iter++)
                    {
                        /* evaluate solution to (H') */
                        SSM (X, NULL, k, Mx, Mi, Mp, Mb, rho, 1.e-10,NULL,NULL);

                        /* compare objective at D = Z*X to quadratic at X */
                        cg_scale (D, Z, X [0], n) ;
                        if ( k > 1 ) cg_daxpy (D, Z1, X [1], n) ;
                        if ( k > 2 ) cg_daxpy (D, Z2, X [2], n) ;

                        /* if Aexists, project D on null space of constraints */
#ifdef PASA
                        if ( Aexists )
                        {
                            pasa_null_project (d, D, NULL, FALSE, pasacom) ;
                        }
#endif

                        /* evaluate objective at x + d (different evaluate
                           process for a quadratic) */
                        if ( QuadCost ) trial = -2 ;
                        else            trial = CGONE ;
                        status = XXCG(evaluate)
                                        (CGZERO, &trial, "f", PASA_CG_COM) ;
                        if ( QuadCost ) dq = cgcom->dq ;
                        else            dq = cgcom->f - f ;

                        /* compute value of the Hessian-based model at X */
                        qcost = X [0]*cg_dot (X, Mx, k) ;
                        if ( k > 1 ) qcost += X [1]*cg_dot (X, Mx+k, k) ;
                        if ( k > 2 ) qcost += X [2]*cg_dot (X, Mx+(k+k), k) ;
                        qcost = 0.5*qcost + cg_dot (Mb, X, k) ;

                        if ( PrintLevel > 1 )
                        {
                            printf ("trust rho: %e qcost: %e obj. cost: %e "
                                    "NegCurv: %i\n",
                                     rho, qcost, dq, NegCurv) ;
                        }

                        if ( dq < tbest )
                        {
                            restore_x = FALSE ;
                            tbest = dq ;
                            Xb [0] = X [0] ;
                            if ( k == 2 )
                            {
                                Xb [1] = X [1] ;
                            }
                            else if ( k > 2 )
                            {
                                Xb [1] = X [1] ;
                                Xb [2] = X [2] ;
                            }
                        }
                        /* otherwise the values of Xb will be used to
                           restore the prior best evaluation point */
                        else restore_x = TRUE ;

                        if ( t_iter == 0 && tbest <= CGZERO )
                        {
#ifdef PASA
                            UnitNewtonStep = TRUE ;
#endif
                            break ;
                        }
                        else if ( (dq > tbest) && (tbest < CGZERO) )
                        {
                            /* the previous best value was smaller than
                               f at the starting point and the cost grew
                               when rho decreased */
                            break ;
                        }

                        /* Trial may decrease if f is infinite or nan
                           We continue of decrease rho as long as there is
                           some possibility that the decrease will lead to
                           a better objective value. If the recent absolute
                           objective values are much smaller than the
                           cost improvement (-qcost) based on the quadratic
                           model, stop trying to decrease rho */
                        if ( (status == CG_FUNCTION_NAN_OR_INF) ||
                             (trial < CGONE) || ((dq >= 0.1*qcost)
                             && (-qcost > 1.e-12*Ck)) )
                        {
                            /* continue but for smaller rho */
                            rho *= Parm->RhoDecay ;
                            continue ;
                        }
                        break ;
                    }

                    /* store current function value if not a quadratic and
                       t_iter > 0 */
                    if ( t_iter && !QuadCost )
                    {
                        cgcom->f = f + tbest ;
#ifdef PASA
                        pasacom->f = cgcom->f ;
#endif
                    }
                    /* If t_iter = TrustIterLimit, then the search direction
                       is essentially the negative gradient. Return to using
                       cg_descent. */
                    if ( t_iter == TrustIterLimit )
                    {
                        use_hessian = FALSE ;
                        iter_since_err_update = 0 ;
                        err_mark = err_decay*err ;
                        /* the negative gradient is the search direction */
                        cg_scale (D, gproj, -CGONE, n) ;
                        NewtonStep = FALSE ;
#ifdef PASA
                        if ( Aexists )
                        {
                            pasa_null_project (d, D, NULL, FALSE, pasacom) ;
                        }
#endif
                    }
                    else if ( restore_x )
                    {
                        cg_scale (D, Z, Xb [0], n) ;
                        if ( k > 1 ) cg_daxpy (D, Z1, Xb [1], n) ;
                        if ( k > 2 ) cg_daxpy (D, Z2, Xb [2], n) ;
                        /* if Aexists, project D on null space of constraints */
#ifdef PASA
                        if ( Aexists )
                        {
                            pasa_null_project (d, D,NULL,FALSE,pasacom);
                        }
#endif
                    }
                    /* else the move is stored in d */
                    break ;
                }
#ifdef TIME
cgcom->loop_cgtrust += sopt_timer () - tic ;
#endif
                /* The search direction in original x space is stored in d.
                   If the search direction is an ascent direction, then discard
                   it and set the search direction to be the negative gradient
                   (or the negative projected gradient) */
                dphi0  = cg_dot (g, d, n) ;
                if ( dphi0 > CGZERO ||
                     (rnorm > rnorm_start && NegCurv == FALSE) )
                {
                    /* When restarting cg_descent, the gradient must be
                       projected twice. gproj has been projected once, so
                       now perform the second projection */
                    cg_scale(D, gproj, -CGONE, n) ;
#ifdef PASA
                    if ( Aexists )
                    {
                        pasa_null_project(d, D, NULL, FALSE, pasacom) ;
                    }
#endif
                    dphi0 = cg_dot (g, d, n) ;
                    NewtonStep = FALSE ;
                    if ( dphi0 > CGZERO )
                    {
                        if ( PrintLevel )
                        {
                            printf("search direction not descent direction "
                                   "after Newton CG\n") ;
                        }
                        status = CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION ;
                        XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                        return (status) ;
                    }
                }

                nCG++ ;
                Stat->nCG++ ;
#ifdef PASA
                /* if dCG flipped from TRUE to FALSE, then also set
                   Hessian_CG_solver to FALSE so that the SS will get first
                   preference at the start */
                if ( !dCG ) pasacom->Hessian_CG_solver = FALSE ;
#endif

                if ( PrintLevel )
                {
                    printf ("dphi0: %e nCG: %li total: %li\n",
                             dphi0, (LONG) nCG, (LONG) Stat->nCG) ;
                }
                /* If pasacom->Hessian_CG_solver = TRUE, then
                   initially try to solve the Hessian-based model using
                   PRP+ CG. Otherwise, go directly to the symmetric solver. */
            }
            /* #NSS */
            else /* use symmetric solver */
            {
                ASSERT (NewMethod == NewtonSS) ;
#ifdef TIME
double tic = sopt_timer () ;
#endif
                /* The triples needed by the symmetric solver are extracted
                   from the H-matrix and stored in the C-matrix, which
                   becomes the input for the symmetric solver. */
                if ( SSsetupDone && HessianSparsityFixed )
                {
                    sopt_convert_H_to_C_in_SS (C, H, FALSE, Cdim,Annz,location);
#ifndef NDEBUG
                    status = XXCG(debug_C_in_SS) (C, H) ;
                    if ( status ) return (sopt_convert_error (LCG, status)) ;
#endif
                }
                else if ( (nSS == 0) || !HessianSparsityFixed )
                {
                   /* If the setup was not done previously, then we not only
                      strip the values needed for the SS, but also determine
                      the mapping from H to the triples in C. The processing
                      of AT, only done in first iteration nSS = 0 */
                    Annz = 0 ;  /* nnz's in A */
                    Cdim = ncol ;
#ifdef PASA
                    if ( Aexists )
                    {
                        if ( use_napheap ) /* napsack constraint exists */
                        {
                            if ( nap_present ) /* constraint is active */
                            {
                                Cdim++ ;
                                Annz = pasacom->nap_nnz ;
                            }
                        }
                        else /* linear constraints */
                        {
                            Cdim += nrow ;
                            Annz = pasacom->Copy->ATp [nrow] ;
                        }
                    }
#endif
                    sopt_convert_H_to_C_in_SS (C, H, TRUE, Cdim,Annz,location) ;
                    /* if diagonal perturbation of Hessian is used, compute Hmax
                       if not yet computed */
                    if ( Hmax < CGZERO ) /* Hmax not yet computed */
                    {
                        if ( Cmax <= CGZERO ) /* Cmax not yet computed */
                        {
                            Hmax = sopt_sup_normx (C->vals, C->nnz) ;
                            if ( Cmax < CGZERO )
                            {
                                Cmax = Hmax*Parm->cg_diag_pert ;
                            }
                            Hmax *= Parm->ss_diag_pert ;
                        }
                        else /* express Hmax in terms of Cmax */
                        {
                            Hmax = Parm->ss_diag_pert*(Cmax/Parm->cg_diag_pert);
                        }
                    }
                    /* otherwise, either Hmax was already evaluated or diagonal
                       perturbation not used (Hmax = 0) */
#ifndef NDEBUG
                    status = XXCG(debug_C_in_SS) (C, H) ;
                    if ( status ) return (sopt_convert_error (LCG, status)) ;
#endif
                }

                /* only need to install the constraints in the triples if
                   sparsity not fixed or sparsity is fixed, nSS = 0, and
                   SS setup was not done */
                if ( !HessianSparsityFixed || (!nSS && !SSsetupDone) )
                {
                    /* Insert the A matrix in C if it exists. Space was left for
                       A when the triples arrays were malloc'd in the
                       convert routine. This only needs to be done in the first
                       iteration if the Hessian sparsity is fixed. It needs
                       to be done in every iteration when the Hessian sparsity
                       is not fixed. */
                    Amax = CGZERO ;
                    li = C->nnz ;
                    j = n ;
#ifdef PASA
                    if ( use_napheap ) /* napsack constraint exists */
                    {
                        if ( nap_present ) /* napsack constraint active */
                        {
                            j = n + 1 ; /* fortran indexing */
                            for (i = 0; i < n; i++)
                            {
                                CGFLOAT const ax = pasacom->nap_a [i] ;
                                if ( ax )
                                {
                                    if ( fabs (ax) > Amax ) Amax = fabs (ax) ;
                                    C->vals [li] = ax ;
                                    C->rows [li] = i + 1 ;
                                    C->cols [li] = j ;
                                    li++ ;
                                }
                            }
                        }
                    }
                    else if ( Aexists ) /* linear constraints present */
                    {
                        PPwork *W = pasacom->ppcom->Work ;
                        PASAINT   *ir  = W->ir ;
                        PASAINT   *ATp = W->ATp ;
                        PASAINT   *ATi = W->ATi ;
                        PASAFLOAT *ATx = W->ATx ;
                        for (i = 0; i < nrow; i++)
                        {
                            if ( ir [i] == 0 )
                            {
                                /* Fortran indexing so column is n + 1 */
                                j++ ;
                                PASAINT const q = ATp [i+1] ;
                                for (p = ATp [i]; p < q; p++)
                                {
                                    CGFLOAT const ax = ATx [p] ;
                                    C->vals [li] = ax ;
                                    if ( fabs (ax) > Amax ) Amax = fabs (ax) ;
                                    C->rows [li] = ATi [p] + 1 ;
                                    C->cols [li] = j ;
                                    li++ ;
                                }
                            }
                        }
                    }
                    Amax *= Parm->ss_diag_pert ;
#endif
                    Anrow = j - n ; /* number of active rows */

                    /* total dim, Hessian plus active rows in A */
                    nt = j ;

                    /* total nnz, Hessian nnz + active rows + diag regular.. */
                    nnzt = li + nt ;

                    /* Insert the diagonal regularization terms into C */
                    cg_initx (C->vals+li, Hmax, n) ;
                    cg_initx (C->vals+(li+n), -Amax, Anrow) ;
                    row = C->rows+(li-1) ;
                    col = C->cols+(li-1) ;
                    for (i = 1; i <= nt; i++)
                    {
                        row [i] = i ;
                        col [i] = i ;
                    }

#ifdef USE_MUMPS
                    /* MUMPS initialization */
                    if ( !SSinitDone )
                    {
                        /* MPI_Init and MPI_Comm_rank return int ierr */
                        MPI_Init (&argc, &argv) ;
                        MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
                        id = sopt_malloc (&status, sizeof (DMUMPS_STRUC_C), 1) ;
#ifdef PASA
                        pasacom->id = id ;
                        pasacom->SSinitDone = TRUE ;
#else
                        cgcom->id = id ;
#endif
                        SSinitDone = TRUE ;
                        id->comm_fortran = USE_COMM_WORLD ;
                        id->par = 1 ;
                        id->sym = 2 ;
                        id->job = JOB_INIT ;
                        dmumps_c (id) ;

                        /* 0 = no printing, 1 = error + warning, 2 = 1 + stats*/
                        id->icntl [3] = 1 ;
                        /* NOTE: one less than shown in manual, could use:
                                     #define ICNTL(I) icntl[(I)-1]*/
                        /* allow extra space for factors due to pivoting */
                        /* id->icntl [13] = 100 ; (100% increase in space) */
                    }
                    else /* solver was previously initialized */
                    {
#ifdef PASA
                        id = pasacom->id ;
#else
                        id = cgcom->id ;
#endif
                    }

                    /* input to MUMPS */
                    id->n  = nt ;      /* total n, active rows + Hessian  */
                    id->nnz = nnzt ;   /* Hessian nnz + active rows + diag */
                    if ( PrintLevel )
                    {
                        printf ("CG symmetric solve: %li (Hessian) %li active "
                                "rows %li nnz\n",
                                       (LONG) n, (LONG) Anrow, (LONG) li + nt) ;
                    }
                    id->irn = C->rows ;
                    id->jcn = C->cols ;
                    id->a   = C->vals ;
/* end of MUMPS set up */
#endif

#ifdef USE_HARWELL
                    if ( lbl != NULL )
                    {
                        lbl->irn = lbl->ma57->irn = NULL ;
                        lbl->jcn = lbl->ma57->jcn = NULL ;
                        LBL_Finalize(lbl) ;
                    }
                    lbl = LBL_Initialize (nnzt, nt, lblfile, MA57_SOLVER) ;
                    /* no scaling */
                    LBL_set_int_parm(lbl, LBL_I_SCALING, 0) ;
                    /* provide ordering = 1,
                       AMD pivoting     = 2,
                       AMD/ND           = 5 */
                    LBL_set_int_parm(lbl, LBL_I_PIV_SELECTION, 1) ;
                    /* do pivoting */
                    LBL_set_int_parm(lbl, LBL_I_PIV_NUMERICAL, 1) ;
                    /* pivot threshold */
                    /*LBL_set_real_parm(lbl, LBL_D_PIV_THRESH, 0.5) ;*/

                    /* MA57 malloc's space for rows and cols, free this space
                       since we already have space in C */
                    cg_free (lbl->ma57->irn) ;
                    cg_free (lbl->ma57->jcn) ;
                    lbl->irn = lbl->ma57->irn = C->rows ;
                    lbl->jcn = lbl->ma57->jcn = C->cols ;
#ifdef PASA
                    pasacom->lbl = lbl ;
#else
                    cgcom->lbl = lbl ;
#endif
                    /* CHOLMOD is used to find the ordering for the rows and
                       columns of the matrix. The ordering is stored in
                       lbl->ma57->keep. Begin by building the matrix as a
                       CHOLMOD sparse matrix. */
                    if ( Achol->p == NULL ) /* malloc cholmod_sparse matrix */
                    {
                        i = ncol + nrow + 1 ;
                        Achol->p = (CGINT *)   cg_malloc (&status, i, sI) ;
                        i = Annz + C->nnzsym ; /* most nnz in matrix*/
                        Achol->i = (CGINT *)   cg_malloc (&status, i, sI) ;
                        Achol->x = (CGFLOAT *) cg_malloc (&status, i, sF) ;
                    }
                    /* compute column counts */
                    CGINT *cc = Achol->p ;
                    cg_initi (cc, (CGINT) 0, nt+1) ;
                    Achol->nrow = nt ;
                    Achol->ncol = nt ;
                    row = C->rows ;
                    col = C->cols ;
                    for (i = 0; i < li; i++)
                    {
                        if ( row [i] != col [i] ) cc [col [i]]++ ;
                    }
                    /* evaluate the column pointers */
                    for (i = 0; i < nt; i++) cc [i+1] += cc [i] ;
                    /* scatter the matrix elements in Achol */
                    CGINT   *Ai = Achol->i ;
                    for (i = 0; i < li; i++)
                    {
                        /* Fortran indexing for row and col */
                        if ( row [i] != col [i] )
                        {
                            CGINT const ccj = cc [col [i]-1]++ ;
                            Ai [ccj] = row [i] - 1 ;
                        }
                    }
                    /* shift cc up by one and cc [0] = 0 */
                    for (j = nt; j > 0; j--) cc [j] = cc [j-1] ;
                    cc [0] = 0 ;
                    CGINT *cholperm = lbl->ma57->keep ;
                    /* decide between AMD and Metis based on graph density */
                    t = (double) nt ;
                    if ( cc [nt] > 0.4*t*(t+1) )
                    {
                        if ( PrintLevel == 2 )
                        {
                            printf ("Order Hessian with AMD in cg_descent\n") ;
                        }
                        /* use AMD */
                        CHOLMOD (camd) (Achol, NULL, 0, NULL, cholperm, cmm) ;
                    }
                    else
                    {
                        if ( PrintLevel == 2 )
                        {
                            printf ("Order Hessian with Metis in cg_descent\n");
                        }
                        /* use Metis */
                        if ( Cparent == NULL )
                        {
#ifdef PASA
                            Cmember = cg_malloc (&status, ncol + nrow, sI) ;
                            Cparent = cg_malloc (&status, ncol + nrow, sI) ;
                            pasacom->Cparent = Cparent ;
                            pasacom->Cmember = Cmember ;
#else
                            Cmember = cg_malloc (&status, ncol, sI) ;
                            Cparent = cg_malloc (&status, ncol, sI) ;
                            cgcom->Cparent = Cparent ;
                            cgcom->Cmember = Cmember ;
#endif
                        }
                        cmm->current = 0 ;
                        cmm->method [0].nd_components = 1 ;
                        CHOLMOD (nested_dissection) (Achol, NULL, 0,
                                       cholperm, Cparent, Cmember, cmm) ;
                    }
                    /* convert permutation to Fortran format */
                    for (i = 0; i < nt; i++) cholperm [i]++ ;
/* end of Harwell setup */
#endif
                }
                if ( !nSS )
                {
                    /* malloc space for right side in symmetric solve
                       and also space for a Hessian-vector product in the
                       SSM section */
                    if ( C->rhs == NULL )
                    {
                        /* Only need the first n components of solution
                           after we solve the linear system with n + Anrow
                           equations. Also 2n extra space for the Hessian-vector
                           product after solving the linear system, but
                           Anrow is no longer needed. Thus the space is
                           ncol + MAX (ncol, nrow) */
#ifdef PASA
                        i = ncol + CGMAX (nrow, ncol) ;
#else
                        i = ncol + ncol ;
#endif
                        C->rhs = rhs = sopt_malloc (&status, i, sF) ;
                        if ( status )
                        {
                            status = sopt_convert_error (location, status) ;
                            XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                            return (status) ;
                        }
                    }
                    else rhs = C->rhs ;
                }

                /* install negative gradient top part of right side */
                sopt_scale (rhs, g, -CGONE, n) ;

                /* If Aexists, install the trailing entries of rhs
                   corresponding to the violation of the linear constraints */
#ifdef PASA
                if ( use_napheap ) /* napsack constraint exists */
                {
                    if ( nap_present ) /* napsack constraint active */
                    {
                        PASAFLOAT rowsum = PASAZERO ;
                        for (i = C->nnz; i < li; i++)
                        {
                            /* convert fortran column to C column (subtract 1)*/
                            rowsum += x [C->rows [i] - 1]*C->vals [i] ;
                        }
                        if ( pasacom->nap_constraint < 0 )
                        {
                            rhs [n] = pasacom->nap_bl - rowsum ;
                        }
                        else
                        {
                            rhs [n] = pasacom->nap_bu - rowsum ;
                        }
                    }
                }
                else if ( Aexists ) /* linear constraint present */
                {
                    /* this code puts zero on the right side of the equations */
                    pasa_initx (rhs+n, PASAZERO, Anrow) ;
#if 0
                    /* this block of code incorporates the constraint violation
                       in the right side for the linear constraint */
                    PPwork *W = pasacom->ppcom->Work ;
                    PASAINT   *ir  = W->ir ;
                    PASAINT   *ATp = W->ATp ;
                    PASAINT   *ATi = W->ATi ;
                    PASAFLOAT *ATx = W->ATx ;
                    nt = n ;
                    for (i = 0; i < nrow; i++)
                    {
                        if ( ir [i] == 0 )
                        {
                            PASAINT const q = ATp [i+1] ;
                            PASAFLOAT rowsum = PASAZERO ;
                            for (p = ATp [i]; p < q; p++)
                            {
                                rowsum += x [ATi [p]]*ATx [p] ;
                            }
                            rhs [nt] = pasacom->b [i] - rowsum ;
                            nt++ ;
                        }
                    }
#endif
                }
#endif
#ifndef NDEBUG
                /* save the right side to check solution accuracy later */
                CGFLOAT *rhsDebug = sopt_malloc (&status, nt, sF) ;
                sopt_copyx (rhsDebug, rhs, nt) ;
#endif

#ifdef TIME
t = sopt_timer () ;
#endif

#ifdef USE_MUMPS
                /* analysis of matrix only done in first iteration when
                   Hessian sparsity is fixed, otherwise done every iteration */
                id->rhs = rhs ;
                if ( (!nSS && !SSsetupDone) || !HessianSparsityFixed )
                {
                    /* analyze matrix (permute) */
                    id->job = 1 ;
                    dmumps_c (id) ;
                    if ( id->info [0] < 0 ) *SSfail = TRUE ;
                }
                if ( !HessianSparsityFixed )
                {
                    /* factor matrix and solve linear system */
                    if ( SetTimer )
                    {
                        method_tic = sopt_timer () - htic ;
                        SetTimer = FALSE ;
                    }
                    id->job = 5 ;
                    dmumps_c (id) ;
                    if ( id->info [0] < 0 ) *SSfail = TRUE ;
                }
                else /* HessianSparsityFixed */
                {
                    if ( !nSS ) /* initial iteration */
                    {
                        /* start the timer if not a quadratic */
                        if ( !QuadCost && SetTimer)
                        {
                            method_tic = sopt_timer () - htic ;
                        }
                        id->job = 2 ;
                        dmumps_c (id) ;
                        if ( id->info [0] < 0 )
                        {
                            *SSfail = TRUE ;
                        }
                        else if ( QuadCost && SetTimer )
                        {
                            method_tic = sopt_timer () - htic ;
                        }
                        SetTimer = FALSE ;
                    }
                    else if ( !QuadCost ) /* nSS > 0 and not quadratic */
                    {
                        if ( SetTimer )
                        {
                            method_tic = sopt_timer () - htic ;
                            SetTimer = FALSE ;
                        }
                        id->job = 2 ;
                        dmumps_c (id) ;
                        if ( id->info [0] < 0 ) *SSfail = TRUE ;
                    }
                    /* solve linear system if the factorization did not fail */
                    if ( !(*SSfail) )
                    {
                        id->job = 3 ;
                        dmumps_c (id) ;
                        if ( id->info [0] < 0 ) *SSfail = TRUE ;
                    }
                }
                if ( PrintLevel && *SSfail )
                {
                    printf ("MUMPS fails in Hessian-based solver\n"
                            "info [1]: %li info [2]: %li\n",
                            (LONG) id->info [0], (LONG) id->info [1]) ;
                }
/* end of MUMPS solve */
#endif

#ifdef USE_HARWELL
                if ( (!nSS && !SSsetupDone) || !HessianSparsityFixed )
                {
                    /* iflag = 0 => automatic pivot choice in ma27 */
                    status = LBL_Analyze(lbl, 0) ;
                    if( status < 0 )
                    {
                        printf ("In Hessian-based cg_descent, analysis fails.\n"
                                "Error %i returned by Harwell MA57\n", status) ;
                        *SSfail = TRUE ;
                    }
                }
                if ( !(*SSfail) ) /* Harwell analysis phase did not fail */
                {
                    status = 0 ;
                    if ( !HessianSparsityFixed )
                    {
                        /* factor matrix */
                        if ( SetTimer )
                        {
                            method_tic = sopt_timer () - htic ;
                            SetTimer = FALSE ;
                        }
                        status = LBL_Factorize (lbl, C->vals) ;
                        if ( PrintLevel && status == 0 )
                        {
                            XXCG(harwell_fact_data) (lbl) ;
                        }
                    }
                    else /* Hessian sparsity is fixed */
                    {
                        if ( !nSS ) /* initial iteration */
                        {
                            /* start the timer if not a quadratic */
                            if ( !QuadCost && SetTimer )
                            {
                               method_tic = sopt_timer () - htic ;
                            }

                            /* factor matrix */
                            status = LBL_Factorize (lbl, C->vals) ;
                            if ( QuadCost && SetTimer )
                            {
                                method_tic = sopt_timer () - htic ;
                            }
                            SetTimer = FALSE ;
                            if ( PrintLevel && status == 0 )
                            {
                                XXCG(harwell_fact_data) (lbl) ;
                            }
                        }
                        else if ( !QuadCost ) /* nSS > 0 and not quadratic */
                        {
                            if ( SetTimer )
                            {
                                method_tic = sopt_timer () - htic ;
                                SetTimer = FALSE ;
                            }

                            /* factor matrix */
                            status = LBL_Factorize (lbl, C->vals) ;
                            if ( PrintLevel && status == 0 )
                            {
                                XXCG(harwell_fact_data) (lbl) ;
                            }
                        }
                    }
                    if( status < 0 )
                    {
                        printf("Error %i returned by Harwell MA57 in "
                               "factorization\n", status) ;
                        *SSfail = TRUE ;
                    }
                    /* solve linear system if the factorization did not fail */
                    else
                    {
                        status = LBL_Solve(lbl, rhs) ;
                        if( status < 0 )
                        {
                            *SSfail = TRUE ;
                            if ( PrintLevel )
                            {
                                printf("Error %i returned by Harwell MA57 in "
                                       "solve\n", status) ;
                            }
                        }
                    }
                }
#endif
#ifdef TIME
cgcom->loop_factor += sopt_timer () - t ;
#endif
                SSsetupDone = TRUE ;
                nSS++ ;
                Stat->nSS++ ;
#ifdef TIME
t = sopt_timer () - tic ;
printf ("factor + solve time: %e\n", t) ;
cgcom->loop_SSmumps += t ;
#endif

#ifndef NDEBUG
                /* evaluate solution accuracy */
#ifdef USE_MUMPS
                status = XXCG(checkSS) (id, C->vals, rhs, rhsDebug, nnzt) ;
#elif USE_HARWELL
                status = XXCG(checkSS) (lbl, C->vals, rhs, rhsDebug, nnzt) ;
#endif
#endif
                /* If the symmetric solver did not fail, then form an
                   orthonormal basis for the space spanned by the projected
                   gradient and by the solution of the KKT system. */
#ifdef TIME
tic = sopt_timer () ;
#endif
                if ( !(*SSfail) )
                {
                    cg_copy (V, gproj, n) ;
                    cg_copy (V+n, rhs, n) ;
                    k = cg_basis (Z, T, V, perm, PrintLevel, n, 2,
                                      1.e-8, basisWork) ;
                    if ( PrintLevel )
                    {
                        printf ("trust region basis vectors: %li\n", (LONG) k) ;
                    }
                }

                /* if there is only one independent vector, or the symmetric
                   solver fails, then the only option is negative gradient
                   search */
                if ( *SSfail || k == 1 )
                {
#if 0
                    if ( *SSfail && (deriv_mode == 2) )
                    {
                        status = CG_SYMMETRIC_SOLVER_FAILS ;
                        XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                        return (status) ;
                    }
#endif
                    /* we take a gradient step */
                    NewtonStep = FALSE ;

                    /* If SS fails, then we should never return to NewtonSS.
                       Thus if use_hessian is TRUE, as it currently is, then
                       should have dCG = TRUE */
                    dCG = TRUE ;

                    /* SS failed if only one independent vector */
                    *SSfail = TRUE ;
                    speed [NewtonSS] = -CGINF ;

                    /* if both solvers fail, then switch to gradient-based
                       deriv_mode = 1 */
                    if ( *CGNfail )
                    {
                        use_hessian = FALSE ;
                        iter_check = INT_MAX ;
                        iter_since_err_update = 0 ;
                        AutoChoice = FALSE ;
                        deriv_mode = 1 ;
#ifdef PASA
                        pasacom->deriv_mode = 1 ;
#endif
                        if ( PrintLevel )
                        {
                            printf ("\n\n## Both the symmetric solver "
                                    "NewtonSS and NewtonCG routines "
                                    "have failed.\nSwitch to ordinary "
                                    "cg_descent. ##\n") ;
                        }
                    }
                    else /* CG Newton has not failed */
                    {
                        /* Continue Hessian-based routine if NewtonCG is at
                           least close to the speed of PureCG; otherwise,
                           return to PureCG */
                        if ( speed [PureCG] > 0.8*speed [NewtonCG] )
                        {
                            /* switch to PureCG */
                            use_hessian = FALSE ;
                            iter_since_err_update = 0 ;
                            err_mark = err_decay*err ;
                        }
                        /* else since dCG changed from FALSE to TRUE,
                           AutoChoice will see that it should change to
                           the CG Newton routine and set a new timer */
                    }
           
                    /* negative gradient = search direction */
                    cg_scale (D, gproj, -CGONE, n) ;
#ifdef PASA
                    if ( Aexists )
                    {
                        pasa_null_project(d, D, NULL, FALSE, pasacom) ;
                    }
#endif
                }
                /* There are at least 2 independent vectors, one
                   is the gradient p0 = PHP*Z(:, 0) p1 = PHP*Z(:, 1)
                   We perform both a trust region step and a Newton step,
                   and use the step that gives the best reduction in the
                   objective value. To begin, let us form the matrix M
                   of problem (H') starting with products between
                   columns of Z and PHP*Z(:, 0) */

                else
                {
                    /* build sparse matrix, D is only used when Aexists */
                    Z1 = Z+n ;
                    pd = XXCG(PHPz) (Z, rhs+n, D, Aexists, TRUE, C_PASA) ;
                    Mx [0] = cg_dot (Z, pd, n) ;
                    Mx [1] = Mx [2] = cg_dot (Z+n, pd, n) ;

                    pd = XXCG(PHPz) (Z1, rhs+n, D, Aexists, TRUE, C_PASA) ;
                    Mx [3] = cg_dot (Z1, pd, n) ;

                    Mi [0] = Mi [2] = 0 ; Mi [1] = Mi [3] = 1 ;
                    Mp [0] = 0 ; Mp [1] = 2 ; Mp [2] = 4 ;
                    /* build the linear term in (H') */
                    Mb [0] = cg_dot (gproj, Z, n) ;
                    Mb [1] = cg_dot (gproj, Z1, n) ;
#if 0
sopt_printx (Mx, 4, "Mx") ;
t = sopt_dot (gproj, D, n) /sqrt(sopt_dot (gproj, gproj, n)*sopt_dot (D, D, n)) ;
printf ("angle: %e\n", t) ;
printf ("DHD/DD: %e\n", cg_dot(D, g, n)/cg_dot(D, D, n)) ;

for (i = 0; i < 2; i++)
{
    X [i] = cg_dot (rhs, Z+i*n, n) ;
}
t = 0.5*(X [0]*(X [0]*Mx [0] + X [1]*Mx [2])  +
         X [1]*(X [0]*Mx [1] + X [1]*Mx [3])) +
         X [0]*Mb [0] + X [1]*Mb [1] ;
printf ("cost for hessian-based solution: %e\n", t) ;

pd = XXCG(PHPz) (rhs, rhs+n, D, Aexists, TRUE, C_PASA) ;
t = 0.5*cg_dot (pd, rhs, n) + cg_dot (rhs, gproj, n) ;
printf ("cost based on direct computation: %e\n", t) ;
#endif


                    /* Use SSM to minimize the quadratic over the sphere of
                       radius rho. We want D to be included in the sphere so
                       we choose the radius to be the length of solution
                       to the linear system */
                    rho = sqrt (cg_dot (rhs, rhs, n)) ;
                    dphi0  = cg_dot (g, rhs, n) ;
                    /* dphi0 > 0 => negative curvature => expand rho */
                    if ( dphi0 > CGZERO )
                    {
                        NegCurv = TRUE ;
                        rho *= 10 ;
                    }
                    else NegCurv = FALSE ;

                    /* the current function value at x is stored in f */
                    CGFLOAT tbest = CGINF ;
                    int t_iter ;
                    /* In the loop below, we compute the minimizer of the
                       quadratic model over a sphere of radius rho and the
                       space spanned by the projected gradient and the computed
                       search direction D. We then evaluate the objective
                       at xk plus the minimizer of the quadratic. If this point
                       generates a reduction in the objective value, then
                       we break and use this point in a line search.
                       If this point increases the objective value, then
                       we continue to reduce rho until finding a point
                       where either the objective increases or the decrease in
                       the objective is at least a small fraction of the
                       decreased anticipated from the quadratic model.
                       The point used in the CG line search is always the
                       point with the smallest objective value. If this point
                       is not near the solution to the first-order optimality
                       conditions computed by the symmetric solver, then
                       we abandon the symmetric solver.  */
                    for (t_iter = 0; t_iter < TrustIterLimit; t_iter++)
                    {
                        /* compute solution of (H') */
                        SSM (X, NULL, 2, Mx, Mi, Mp, Mb, rho, 1.e-10,NULL,NULL);

                        /* compare objective at D = Z*X to quadratic at X */
                        cg_scale (D, Z, X [0], n) ;
                        cg_daxpy (D, Z1,X [1], n) ;
#if 0
                        /* If the solution to the trust region problem is not
                           close to the solution to the solution of the
                           linearly constrained QP, then abandon the
                           symmetric solver and return to CG, but increase
                           the bound on the number of CG iterations. */
                        if ( t_iter == 0 )
                        {
                            tmax = smax = CGZERO ;
                            for (i = 0; i < n; i++)
                            {
                                if ( t = fabs (rhs [i] - D [i]) > tmax ) tmax=t;
                                if ( s = fabs (D [i]) > smax ) smax = s ;
                                s = fabs (D [i]) ;
                            }
                            if ( (smax == CGZERO) || (tmax/smax > 0.5) )
                            {
                                /* symmetric solve does not seem to yield a good
                                   search direction; return to CG */
                                dCG = TRUE ;
                                if ( find_max_iter )
                                {
                                    /* estimate flops for an iteration of CG */
                                    SOPTFLOAT cg_flops = 2*C->nnz + 6*n ;
#ifdef PASA
                                    if ( Aexists )
                                    {
                                        if ( use_napheap )
                                        {
                                            cg_flops += 6*n ; /* 2 projections*/
                                        }
                                        else
                                        {
                                            /* two projection and 2*Lnnz ops
                                               for each solve */
                                            cg_flops +=
                                                  4*pasacom->ppcom->Work->Lnnz ;
                                        }
                                    }
#endif
                                    /* operations in factorization */
                                    CGFLOAT ss_flops = id->rinfog [0] ;
                                    /* # entries in factors: id->infog [19]) */

                                    *max_d_iter = ss_flops/cg_flops ;
                                    if ( n < *max_d_iter ) *max_d_iter = n ;
                                    if ( PrintLevel )
                                    {
                                        printf ("max_d_iter: %li ss_flops: %e cg_flops: %e\n",
                                                (LONG) *max_d_iter, ss_flops, cg_flops) ;
                                    }
                                    find_max_iter = FALSE ;
                                }
                                *max_d_iter *= 2 ;
                            }
                        }
#endif
#ifdef PASA
                        if ( Aexists )
                        {
                            pasa_null_project(d, D, NULL, FALSE, pasacom) ;
                        }
#endif
/* compare qp solution to ss solution */
#if 0
CGFLOAT *diff = malloc (n*sF) ;
sopt_step (diff, rhs, D, -1, n) ;
t = sopt_sup_normx (diff, n) ;
printf ("sup norm D: %e\n", sopt_sup_normx (D, n)) ;
printf ("sup norm rhs: %e\n", sopt_sup_normx (rhs, n)) ;
printf ("relative diff: %e\n", t/sopt_sup_normx (D, n)) ;
t = sopt_dot (rhs, D, n) /sqrt(sopt_dot (rhs, rhs, n)*sopt_dot (D, D, n)) ;
printf ("angle: %e\n", t) ;
#endif
                        /* compute value of the quadratic at X */
                        qcost = X [0]*cg_dot (X, Mx, 2) ;
                        qcost += X [1]*cg_dot (X, Mx+2, 2) ;
                        qcost = 0.5*qcost + cg_dot (Mb, X, 2) ;

                        /* evaluate objective at x + d (different evaluate
                           process for a quadratic) */
                        if ( QuadCost ) trial = -2 ;
                        else            trial = CGONE ;
                        status = XXCG(evaluate)
                                        (CGZERO, &trial, "f", PASA_CG_COM) ;
                        if ( QuadCost ) dq = cgcom->dq ;
                        else            dq = cgcom->f - f ;

                        if ( PrintLevel > 1 )
                        {
                            printf ("trust rho: %e ssqcost: %e obj. cost: %e "
                                    "NegCurv: %i\n", rho, qcost, dq, NegCurv) ;
                        }
                        if ( dq < tbest )
                        {
                            tbest = dq ;
                            Xb [0] = X [0] ;
                            Xb [1] = X [1] ;
                        }
                        smallq = FALSE ;
                        if ( t_iter == 0 && ((dq <= CGZERO) ||
                            (-qcost <= 1.e-12*Ck)) )
                        {
#ifdef PASA
                            UnitNewtonStep = TRUE ;
#endif
                            /* if the quadratic objective is tiny, then
                               perform a unit step if possible */
                            if ( dq > CGZERO ) smallq = TRUE ;
                            break ;
                        }
                        else if ( (dq > tbest) && (tbest < CGZERO))
                        {
                            break ;
                        }
                        /* Trial may decrease if f is infinite or nan
                           We continue of decrease rho as long as there is
                           some possibility that the decrease will lead to
                           a better objective value. If the recent absolute
                           objective values are much smaller than the
                           cost improvement (-qcost) based on the quadratic
                           model, stop trying to decrease rho */
                        if ( (status == CG_FUNCTION_NAN_OR_INF) ||
                             (trial < CGONE) || ((dq >= 0.1*qcost)
                              && (-qcost > 1.e-12*Ck)) )
                        {
                            /* continue but for smaller rho */
                            rho *= Parm->RhoDecay ;
                            continue ;
                        }
                        break ;
                    }
#if 0
t = 0 ;
for (i = 0; i <= 100; i++)
{
    t += .01 ;
    cg_copy (d, rhs, n) ;
    cg_scale (d, d, t, n) ;
    trial = CGONE ;
    status = XXCG(evaluate) (CGZERO, &trial, "f", PASA_CG_COM) ;
    printf ("%10.2f %e\n", t, cgcom->f) ;
}
#endif
                    /* Compare the best trust region step to the Newton
                       step and take the best. If the SS direction does not
                       seem to be a good one, then quite SS and switch back
                       to the Newton CG solve. */
                    int quitSS = TRUE ;
                    if ( t_iter == 0 )
                    {
                        if ( sqrt(X[0]*X[0]+X[1]*X[1]) <= .05*rho )
                        {
                        /* if the minimizer is unconstrained and near zero,
                           then abandon the Hessian-based scheme and return
                           to cg_descent (quitSS = TRUE by default);
                           the iteration is essentially gradient descent,
                           so it is more efficient to use cg_descent. Retain
                           NewtonStep = TRUE (t_iter = 0 => d stores the best
                           search direction and cgcom->f = best value) */
                        }
                        else
                        /* if the minimizer of the quadratic not near the
                           point generated by the SS, then quitSS. Otherwise,
                           accept the move and break. */
                        {
                            CGFLOAT tmax = CGZERO ;
                            CGFLOAT smax = CGZERO ;
                            for (i = 0; i < n; i++)
                            {
                                if ( (t = fabs (rhs [i]-d [i]))>tmax ) tmax = t;
                                if ( (s = fabs (d [i])        )>smax ) smax = s;
                            }
                            if ( smax == CGZERO || tmax/smax > 0.5 )
                            {
                                /* symmetric solve did not yield a good search
                                   direction; return to CG (quitSS = TRUE
                                   by default) */
                                if ( PrintLevel )
                                {
                                    printf ("symmetric solver gave bad search "
                                           "direction, return to cg_descent "
                                           "(norm trust d: %e norm d-new: %e\n",
                                            smax, tmax) ;
                                }
                            }
                            else
                            {
                                quitSS = FALSE ;
                            }
                        }
                    }
                        /* make the best move, if t_iter >= 3,
                           increment trustfails */
                    else if ( ((tbest < f) && (trustfails < 2||t_iter<2))
                              || (smallq == TRUE) )
                    {
                        quitSS = FALSE ;
                        if ( PrintLevel ) printf ("  ---- Do Trust\n") ;
                        if ( t_iter >= TrustIterLimit ) trustfails++ ;
                        /* if t_iter >= 1, restore best objective value and
                           best point (overshot and now restore prior best) */
                        cg_scale (D, Z, Xb [0], n) ;
                        cg_daxpy (D, Z1,Xb [1], n) ;
                        if ( !QuadCost ) cgcom->f = f + tbest ;
#ifdef PASA
                        if ( !QuadCost ) pasacom->f = cgcom->f ;
                        if ( Aexists )
                        {
                            pasa_null_project(d, D, NULL, FALSE, pasacom) ;
                        }
#endif
                    }
                    else /* nothing worked, need to return to cg_descent */
                    {
                        NewtonStep = FALSE ;
                        /* negative gradient = search direction */
                        cg_scale (D, gproj, -CGONE, n) ;
#ifdef PASA
                        if ( Aexists )
                        {
                            pasa_null_project(d, D, NULL, FALSE, pasacom) ;
                        }
#endif
                    }
                    if ( quitSS )
                    {
                        trustfails = 0 ;
                        use_hessian = FALSE ;
                        dCG = TRUE ;
                        iter_since_err_update = 0 ;
                        err_mark = err_decay*err ;
                    }
#ifdef TIME
cgcom->loop_sstrust += sopt_timer () - tic ;
#endif
                }
                dphi0  = cg_dot (g, d, n) ;
                if ( PrintLevel )
                {
                    printf ("dphi0 in SS: %e\n", dphi0) ;
                }
                if (dphi0 > CGZERO ) /*ascent direction */
                {
                    /* switch to -gproj */
#ifdef PASA
                    if ( Aexists )
                    {
                        pasa_null_project (d, gproj, NULL, FALSE, pasacom) ;
                    }
#endif
                    NewtonStep = FALSE ;
                    sopt_scale (d, d, -CGONE, n) ;
                    dphi0 = cg_dot (g, d, n) ;
                    if ( dphi0 > CGZERO )
                    { 
                        status = CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION ;
                        XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                        return (status) ;
                    }
                }
            }
        }
        else if ( LBFGS == 1 )
        {
            if ( IterRestart >= nrestart ) /* restart the l-bfgs method */
            {
                if ( PrintLevel >= 1 )
                {
                    printf ("REstart LBFGS\n") ;
                }
                IterRestart = 0 ;
                IterQuad = 0 ;
                mlast = -1 ;
                memk = 0 ;
                scale = (CGFLOAT) 1 ;

                /* copy xnew to x */
                cg_copy (x, xnew, n) ;

                /* for a quadratic objective, evaluate objective
                   and its gradient from scratch */
                if ( QuadCost == TRUE )
                {
                    t = CGZERO ;
                    XXCG(evaluate) (t, &t, "fg", PASA_CG_COM) ;
                }
                else
                {
                    /* copy gnew to g */
                    cg_copy (g, gnew, n) ;
                }

                /* D = -g or -gproj, search direction before final projection */
                cg_scale (D, gproj, -CGONE, n) ;

#ifdef PASA
                /* d is the final search direction */
                pasa_null_project (d, D, NULL, FALSE, pasacom) ;
                /* delete next ?? */
                /* pasacom->cg_bb_est = scale ;*/
#endif

                /* derivative in search direction without penalty term */
                dphi0  = cg_dot (g, d, n) ;
/*printf ("dphi0: %e\n", dphi0 + pasacom->dp) ;*/
            }
            else
            {
                if ( PrintLevel >= 2 )
                {
                    printf ("ordinary LBFGS update\n") ;
                }
                mlast = (mlast+1) % mem ;
                spp = mlast*n ;
                cg_scale (Sk+spp, D, alpha, n) ;

                /* copy xnew to x */
                cg_copy (x, xnew, n) ;

                cg_step (Yk+spp, gnewproj, gproj, -CGONE, n) ;

                /* copy gnew to g */
                cg_copy (g, gnew, n) ;

                if ( Aexists == TRUE )
                {
                    cg_copy (gproj, gnewproj, n) ;
                    cg_copy (gnew, gnewproj, n) ;
                }

                SkYk [mlast] = alpha*(dphi-dphi0) ;
                if (memk < mem)
                {
                    memk++ ;
                }

                /* calculate Hg = H g, saved in gnew */
                mp = mlast ;  /* memk is the number of vectors in the memory */
                for (j = 0; j < memk; j++)
                {
                    mpp = mp*n ;
                    t = cg_dot (Sk+mpp, gnew, n)/SkYk[mp] ;
                    tau [mp] = t ;
                    cg_daxpy (gnew, Yk+mpp, -t, n) ;
                    mp -=  1;
                    if ( mp < 0 )
                    {
                        mp = mem-1 ;
                    }
                }
                ykyk = cg_dot (Yk+spp, Yk+spp, n) ;
                if ( ykyk > CGZERO ) /* compute new approx to inverse Hessian */
                {
                    scale = SkYk[mlast]/ykyk ;/* approximates inverse Hessian */
#ifdef PASA
                    pasacom->cg_bb_est = scale ;
#endif
                }

                cg_scale (gnew, gnew, scale, n) ;

                for (j = 0; j < memk; j++)
                {
                    mp +=  1 ;
                    if ( mp == mem )
                    {
                        mp = 0 ;
                    }
                    mpp = mp*n ;
                    t = cg_dot (Yk+mpp, gnew, n)/SkYk[mp] ;
                    cg_daxpy (gnew, Sk+mpp, tau [mp]-t, n) ;
                }
                cg_scale (D, gnew, -CGONE, n) ;
#ifdef PASA
                /* transform back to x-space */
                pasa_null_project (d, D, NULL, FALSE, pasacom) ;
#endif

                dphi0 = cg_dot (g, d, n) ;
            }
        } /* end of LBFGS */
        else if ( LBFGS == 0 || LBFGS == 2 ) /* search direction by cg formula*/
        {
            /* check to see whether cg should be restated */
            if ( (IterRestart >= nrestart) || ((IterQuad == qrestart)
                                           && (IterQuad != IterRestart)) )
            {
                if ( PrintLevel >= 1 )
                {
                    printf ("Restart CG\n") ;
                }

                IterRestart = 0 ;
                IterQuad = 0 ;
                beta = CGZERO ;
                scale = (CGFLOAT) 1 ;

                /* set x = xnew, g = gnew */
                cg_copy (x, xnew, n) ;

                /* for a quadratic objective, evaluate objective
                   and its gradient from scratch when restart is performed */
                if ( QuadCost == TRUE )
                {
                    t = CGZERO ;
                    XXCG(evaluate) (t, &t, "fg", PASA_CG_COM) ;
                }
                else
                {
                    /* copy gnew to g */
                    cg_copy (g, gnew, n) ;
                }

                /* D is the search direction before the final projection */
                cg_scale (D, gproj, -CGONE, n) ;

#ifdef PASA
                /* delete next ?? */
                /* pasacom->cg_bb_est = scale ;*/
                /* d is the final search direction */
                pasa_null_project (d, D, NULL, FALSE, pasacom) ;
#endif

                /* compute square of 2-norm of d (used in bb formula) */
                dnorm2 = cg_dot (d, d, n) ;

                /* derivative in search direction without penalty term */
                dphi0  = cg_dot (g, d, n) ;
            }
            else /* ordinary cg update without restart */
            {
                if ( PrintLevel >= 2 )
                {
                    printf ("ordinary cg update\n") ;
                }
                /* set x = xnew */
                cg_copy (x, xnew, n) ;

                /* compute: ykPyk  = (newproj - oldproj)*(newproj - oldproj),
                            gkPyk  =  newproj           *(newproj - oldproj)
                   update: oldproj = newproj */
                cg_update_beta (gproj, gnewproj, &gkPyk, &ykPyk, n) ;
                dkyk = dphi - dphi0 ;

#ifdef PASA
                /* For bound constrained problems, we set g = gnew in
                   update_beta since g = gproj and gnew = gnewproj.
                   When linear constraints present, we need to set g = gnew. */
                if ( Aexists == TRUE )
                {
                    cg_copy (g, gnew, n) ;
                }
                pasacom->cg_bb_est = alpha*dnorm2/dkyk ;
#endif


                if ( Parm->AdaptiveTheta )
                {
                    t = 2. - CGONE/(0.1*QuadTrust + CGONE) ;
                }
                else
                {
                    t = Parm->theta ;
                }
                beta = (gkPyk - t*dphi*ykPyk/dkyk)/dkyk ;

                /* faster: initialize dnorm2 = gnorm2 at start, then
                           dnorm2 = gnorm2 + beta**2*dnorm2 - 2.*beta*dphi
                           gnorm2 = ||g_{k+1}||^2
                           dnorm2 = ||d_{k+1}||^2
                           dpi = g_{k+1}' d_k */

                /* lower bound for beta is BetaLower*d_k'g_k/ ||d_k||^2 */
                beta = CGMAX (beta, Parm->BetaLower*dphi0/dnorm2) ;

                /* update search direction D = -gproj + beta*Dold */
                cg_update_d (D, gproj, beta, n) ;

#ifdef PASA
                /* project D into null space to obtain search direction */
                pasa_null_project (d, D, NULL, FALSE, pasacom) ;
#endif

                /* derivative in the new search direction */
                dphi0 = cg_dot (d, g, n) ;

                dnorm2 = cg_dot (d, d, n) ;
            }
        }   /* search direction by cg formula has been computed */
        else /* LBFGS = 3, limited memory CG formula */
        {
            if ( RestartAfterHessian )
            {
                if ( Subspace )  /* Subspace exists */
                {
                    Subspace = FALSE ;
                    FirstFull = FALSE ;
                    InvariantSpace = FALSE ;
                }
                use_memory = TRUE ;
                StartCheck = Stat->iter ;
                Restart = TRUE ;  /* restart in Full Space */
                /* Earlier, we evaluated gnewproj is Aexists and copied
                   either gnewproj or gnew to gproj. Later, in the
                   restart, gnorm2 is needed. */
#ifdef PASA
                    if ( Aexists )
                    {
                        gnorm2 = cg_dot (gproj, gproj, n) ;
                    }
#endif
            }
            else if ( use_memory == FALSE )
            {
                if ( (IterRestart >= nrestart) || ((IterQuad == qrestart)
                     && (IterQuad != IterRestart)) )
                {
                    Restart = TRUE ;
                }
            }
            else /* memory is used */
            {
                if ( Subspace ) /* the iteration is in the subspace */
                {
                    IterSubRestart++ ;
                    if ( PrintLevel >= 2 )
                    {
                        printf ("iteration in subspace, IterSubRestart: %i\n",
                                 IterSubRestart) ;
                    }

                    /* compute projection of g into subspace */
                    mp = SkFstart ;
                    j = nsub - mp ;

                    /* multiply basis vectors by new gradient */
                    cg_matvec (wsub, SkF, gnewproj, nsub, n, 0) ;

                    /* rearrange wsub and store in gsubtemp
                       (elements associated with old vectors should
                        precede elements associated with newer vectors */
                    cg_copy0 (gsubtemp, wsub+mp, j) ;
                    cg_copy0 (gsubtemp+j, wsub, mp) ;

                    /* solve Rk'y = gsubtemp */
                    cg_trisolve (gsubtemp, Rk, mem, nsub, 0) ;
                    gsubnorm2 = cg_dot0 (gsubtemp, gsubtemp, nsub) ;
                    /*  gnorm2 = cg_dot (gtemp, gtemp, n) ; */

                    gnorm2 = cg_dot (gnewproj, gnewproj, n) ;
                    ratio = sqrt(gsubnorm2/gnorm2) ;
                    if ( ratio < CGONE - Parm->eta1  )  /* Exit Subspace */
                    {
                       if ( PrintLevel >= 1 )
                       {
                           printf ("CG iter: %li exit subspace ratio: %e "
                                   "eta1: %e\n",
                                   (LONG) Stat->iter, ratio, Parm->eta1) ;
                       }
                       FirstFull = TRUE ; /* first iteration in full space */
                       Subspace = FALSE ; /* leave the subspace */
                       InvariantSpace = FALSE ;
                       /* check the subspace condition for SubCheck iterations
                          starting from the current iteration (StartCheck) */
                       StartCheck = Stat->iter ;
                    }
                    else
                    {
                       /* Check if a restart should be done in subspace */
                       if ( IterSubRestart == nrestartsub )
                       {
                           Restart = TRUE ;
                       }
                    }
                }
                else  /* in full space */
                {
                    if ( (IterRestart == 1) || FirstFull )
                    {
                        memk = 0 ;
                    }
                    if ( (memk == 1) && InvariantSpace )
                    {
                         memk = 0 ;
                         InvariantSpace = FALSE ;
                    }
                    if ( PrintLevel >= 2 )
                    {
                        printf ("iteration in full space, memk: %i\n", memk) ;
                    }
                    if (memk < mem )
                    {
                        memk_is_mem = FALSE ;
                        SkFstart = 0 ;
                        /* SkF stores basis vector of the form alpha*d. We
                           factor SkF = Zk*Rk where Zk has orthonormal columns
                           and Rk is upper triangular. Zk is not stored;
                           wherever it is needed, we use SkF * inv (Rk) */
                        if (memk == 0)
                        {
                            mlast = 0 ;  /* starting pointer in the memory */
                            memk = 1 ;   /* dimension of current subspace */

                            t = sqrt(dnorm2) ;
                            zeta = alpha*t ;
                            Rk [0] = zeta ;
                            cg_scale (SkF, d, alpha, n) ;

                            Yk [0] = (dphi - dphi0)/t ;
                            gsub [0] = dphi/t ;
                            SkYk [0] = alpha*(dphi-dphi0) ;
                            FirstFull = FALSE ;
                            if ( IterRestart > 1 )
                            {
                               /* Need to save g for later correction of first
                                  column of Yk. Since g does not lie in the
                                  subspace and the first column is dense
                               cg_copy (gkeep, g, n) ; */
                               cg_copy (gkeep, gproj, n);
                               /* Also store dot product of g with the first
                                 direction vector -- this saves a later dot
                                 product when we fix the first column of Yk */
                               stgkeep = dphi0*alpha ;
                               d0isg = FALSE ;
                            }
                            else d0isg = TRUE ;
                        }
                        else
                        {
                            mlast = memk ; /* starting pointer in the memory */
                            memk++ ;       /* total number of Rk in the memory*/
                            mpp = mlast*n ;
                            spp = mlast*mem ;
                            cg_scale (SkF+mpp, d, alpha, n) ;

                            /* check if the alphas are far from 1 */
                            if ( !FastLA || (fabs(alpha-5.05)>4.95) ||
                                            (fabs(alphaold-5.05)>4.95) )
                            {
                                /* multiply basis vectors by new direction */
                                cg_matvec (Rk+spp, SkF, SkF+mpp, mlast,n,0);

                                /* solve Rk'y = wsub to obtain the components of
                                   the new direction vector relative to the
                                   orthonormal basis Z in S = ZR, store in
                                   next column of Rk */
                                cg_trisolve (Rk+spp, Rk, mem, mlast, 0) ;
                            }
                            else /* alphas are close to 1 */
                            {
                                t1 = -alpha ;
                                t2 = beta*alpha/alphaold ;
                                for (j = 0; j < mlast; j++)
                                {
                                    Rk [spp+j] = t1*gsub [j] +t2*Rk [spp-mem+j];
                                }
                            }

                            t = alpha*alpha*dnorm2 ;
                            t1 = cg_dot0 (Rk+spp, Rk+spp, mlast) ;
                            if (t <= t1)
                            {
                                zeta = t*1.e-10 ;
                                Stat->NegDiag = TRUE ;
                            }
                            else zeta = sqrt(t-t1);
                            Rk [spp+mlast] = zeta ;
                            /* t = cg_dot0 (Zk+mlast*n, g, n)*/
                            t = - zeta/alpha ;
                            Yk [spp-mem+mlast] = t ;
                            gsub [mlast] = t ;

                            /* multiply basis vectors by new gradient */
                            cg_matvec (wsub, SkF, gnewproj, mlast, n, 0) ;

                            /* exploit dphi for last multiply */
                            wsub [mlast] = alpha*dphi ;
                            /* solve for new gsub */
                            cg_trisolve (wsub, Rk, mem, memk, 0) ;

                            /* subtract old gsub from new gsub = column of Yk */
                            cg_Yk (Yk+spp, gsub, wsub, NULL, memk) ;

                            SkYk [mlast] = alpha*(dphi-dphi0) ;
                        }
                    }
                    else  /* memk = mem */
                    {
                        memk_is_mem = TRUE ;
                        mlast = mem-1 ;
                        cg_scale (stemp, d, alpha, n) ;

                        /* compute projection of s_k = alpha_k d_k into subspace
                           check if the alphas are far from 1 */
                        if (!FastLA||(fabs(alpha-5.05)>4.95)||
                                     (fabs(alphaold-5.05)>4.95))
                        {
                            mp = SkFstart ;
                            j = mem - mp ;

                            /* multiply basis vectors by sk */
                            cg_matvec (wsub, SkF, stemp, mem, n, 0) ;
                            /* rearrange wsub and store in Re = end col Rk */
                            cg_copy0 (Re, wsub+mp, j) ;
                            cg_copy0 (Re+j, wsub, mp) ;

                            /* solve Rk'y = Re */
                            cg_trisolve (Re, Rk, mem, mem, 0) ;
                        }
                        else /* alphas close to 1 */
                        {
                            t1 = -alpha ;
                            t2 = beta*alpha/alphaold ;
                            for (j = 0; j < mem; j++)
                            {
                                Re [j] = t1*gsub [j] + t2*Re [j-mem] ;
                            }
                        }

                        /* t = 2-norm squared of s_k */
                        t = alpha*alpha*dnorm2 ;
                        /* t1 = 2-norm squared of projection */
                        t1 = cg_dot0 (Re, Re, mem) ;
                        if (t <= t1)
                        {
                            zeta = t*1.e-10 ;
                            Stat->NegDiag = TRUE ;
                        }
                        else zeta = sqrt(t-t1);

                        /* dist from new search direction to prior subspace*/
                        Re [mem] = zeta ;

                        /* projection of prior g on new orthogonal
                           subspace vector */
                        t = -zeta/alpha ; /* t = cg_dot(Zk+mpp, g, n)*/
                        gsub [mem] = t ;
                        Yk [memsq] = t ;  /* also store it in Yk */

                        spp = memsq + 1 ;
                        mp = SkFstart ;
                        j = mem - mp ;

                        /* multiply basis vectors by gnew */
                        cg_matvec (vsub, SkF, gnewproj, mem, n, 0) ;

                        /* rearrange and store in wsub */
                        cg_copy0 (wsub, vsub+mp, j) ;
                        cg_copy0 (wsub+j, vsub, mp) ;

                        /* solve Rk'y = wsub */
                        cg_trisolve (wsub, Rk, mem, mem, 0) ;
                        wsub [mem] = (alpha*dphi
                                     - cg_dot0 (wsub, Re, mem))/zeta;

                        /* add new column to Yk, store new gsub */
                        cg_Yk (Yk+spp, gsub, wsub, NULL, mem+1) ;

                        /* store sk (stemp) at SkF+SkFstart */
                        cg_copy (SkF+SkFstart*n, stemp, n) ;
                        SkFstart++ ;
                        if ( SkFstart == mem ) SkFstart = 0 ;

                        mp = SkFstart ;
                        for (k = 0; k < mem; k++)
                        {
                            spp = (k+1)*mem + k ;
                            t1 = Rk [spp] ;
                            t2 = Rk [spp+1] ;
                            t = sqrt(t1*t1 + t2*t2) ;
                            t1 = t1/t ;
                            t2 = t2/t ;

                            /* update Rk */
                            Rk [k*mem+k] = t ;
                            for (j = (k+2); j <= mem; j++)
                            {
                                spp1 = spp ;
                                spp = j*mem + k ;
                                t3 = Rk [spp] ;
                                t4 = Rk [spp+1] ;
                                Rk [spp1] = t1*t3 + t2*t4 ;
                                Rk [spp+1] = t1*t4 - t2*t3 ;
                            }
                            /* update Yk */
                            if ( k < 2 ) /* mem should be greater than 2 */
                            {
                                /* first 2 rows are dense */
                                spp = k ;
                                for (j = 1; j < mem; j++)
                                {
                                    spp1 = spp ;
                                    spp = j*mem + k ;
                                    t3 = Yk [spp] ;
                                    t4 = Yk [spp+1] ;
                                    Yk [spp1] = t1*t3 + t2*t4 ;
                                    Yk [spp+1] = t1*t4 -t2*t3 ;
                                }
                                spp1 = spp ;
                                spp = mem*mem + 1 + k ;
                                t3 = Yk [spp] ;
                                t4 = Yk [spp+1] ;
                                Yk [spp1] = t1*t3 + t2*t4 ;
                                Yk [spp+1] = t1*t4 -t2*t3 ;
                            }
                            else if ( (k == 2) && (2 < mem-1))
                            {
                                spp = k ;

                                /* col 1 dense since the oldest direction
                                    vector has been dropped */
                                j = 1 ;
                                spp1 = spp ;
                                spp = j*mem + k ;
                                /* single nonzero percolates down the column */
                                t3 = Yk [spp] ;  /* t4 = 0. */
                                Yk [spp1] = t1*t3 ;
                                Yk [spp+1] = -t2*t3 ;
                                /* process rows in Hessenberg part of matrix */
                                for (j = 2; j < mem; j++)
                                {
                                    spp1 = spp ;
                                    spp = j*mem + k ;
                                    t3 = Yk [spp] ;
                                    t4 = Yk [spp+1] ;
                                    Yk [spp1] = t1*t3 + t2*t4 ;
                                    Yk [spp+1] = t1*t4 -t2*t3 ;
                                }
                                spp1 = spp ;
                                spp = mem*mem + 1 + k ;
                                t3 = Yk [spp] ;
                                t4 = Yk [spp+1] ;
                                Yk [spp1] = t1*t3 + t2*t4 ;
                                Yk [spp+1] = t1*t4 -t2*t3 ;
                            }
                            else if ( k < (mem-1) )
                            {
                                spp = k ;

                                /* process first column */
                                j = 1 ;
                                spp1 = spp ;
                                spp = j*mem + k ;
                                t3 = Yk [spp] ;  /* t4 = 0. */
                                Yk [spp1] = t1*t3 ;
                                Yk [spp+1] = -t2*t3 ;

                                /* process rows in Hessenberg part of matrix */
                                j = k-1 ;
                                spp = (j-1)*mem+k ;
                                spp1 = spp ;
                                spp = j*mem + k ;
                                t3 = Yk [spp] ;
                                Yk [spp1] = t1*t3 ; /* t4 = 0. */
                                /* Yk [spp+1] = -t2*t3 ;*/
                                /* Theoretically this element is zero */
                                for (j = k; j < mem; j++)
                                {
                                    spp1 = spp ;
                                    spp = j*mem + k ;
                                    t3 = Yk [spp] ;
                                    t4 = Yk [spp+1] ;
                                    Yk [spp1] = t1*t3 + t2*t4 ;
                                    Yk [spp+1] = t1*t4 -t2*t3 ;
                                }
                                spp1 = spp ;
                                spp = mem*mem + 1 + k ;
                                t3 = Yk [spp] ;
                                t4 = Yk [spp+1] ;
                                Yk [spp1] = t1*t3 + t2*t4 ;
                                Yk [spp+1] = t1*t4 -t2*t3 ;
                            }
                            else /* k = mem-1 */
                            {
                                spp = k ;

                                /* process first column */
                                j = 1 ;
                                spp1 = spp ;
                                spp = j*mem + k ;
                                t3 = Yk [spp] ; /* t4 = 0. */
                                Yk [spp1] = t1*t3 ;

                                /* process rows in Hessenberg part of matrix */
                                j = k-1 ;
                                spp = (j-1)*mem+k ;
                                spp1 = spp ;
                                spp = j*mem + k ;
                                t3 = Yk [spp] ; /* t4 = 0. */
                                Yk [spp1] = t1*t3 ;

                                j = k ;
                                spp1 = spp ;
                                spp = j*mem + k ; /* j=mem-1 */
                                t3 = Yk [spp] ;
                                t4 = Yk [spp+1] ;
                                Yk [spp1] = t1*t3 + t2*t4 ;

                                spp1 = spp ;
                                spp = mem*mem + 1 + k ; /* j=mem */
                                t3 = Yk [spp] ;
                                t4 = Yk [spp+1] ;
                                Yk [spp1] = t1*t3 + t2*t4 ;
                            }
                            /* update g in subspace */
                            if ( k < (mem-1) )
                            {
                                t3 = gsub [k] ;
                                t4 = gsub [k+1] ;
                                gsub [k] = t1*t3 + t2*t4 ;
                                gsub [k+1] = t1*t4 -t2*t3 ;
                            }
                            else /* k = mem-1 */
                            {
                                t3 = gsub [k] ;
                                t4 = gsub [k+1] ;
                                gsub [k] = t1*t3 + t2*t4 ;
                            }
                        }

                        /* update SkYk */
                        for (k = 0; k < mlast; k++) SkYk [k] = SkYk [k+1] ;
                        SkYk [mlast] = alpha*(dphi-dphi0) ;
                    }

                    /* calculate t = ||gsub||/||gtemp||  */
                    gsubnorm2 = cg_dot0 (gsub, gsub, memk) ;

                    gnorm2 = cg_dot (gnewproj, gnewproj, n) ;
                    ratio = sqrt (gsubnorm2/gnorm2) ;
                    if ( ratio > CGONE-Parm->eta2 ) InvariantSpace = TRUE ;

                    /* check to see whether to enter subspace */
                    if ( ((memk > 1) && InvariantSpace) ||
                         ((memk == mem) && (ratio > CGONE-Parm->eta0)) )
                    {
                        Stat->NumSub++ ;
                        if ( PrintLevel >= 1 )
                        {
                            if ( InvariantSpace )
                            {
                                printf ("CG iter: %li invariant space of "
                                        "dimension %i\n",
                                        (LONG) Stat->iter, memk) ;
                            }
                            else
                            {
                                printf ("CG iter: %li using memory with %i "
                                        "vectors in memory\n",
                                         (LONG) Stat->iter, memk) ;
                            }
                        }
                        /* if the first column is dense, we need to correct it
                           now since we do not know the entries until the basis
                           is determined */
                        if ( !d0isg && !memk_is_mem )
                        {
                            wsub [0] = stgkeep ;
                            /* mlast = memk -1 */
                            cg_matvec (wsub+1, SkF+n, gkeep, mlast, n, 0) ;
                            /* solve Rk'y = wsub */
                            cg_trisolve (wsub, Rk, mem, memk, 0) ;
                            /* corrected first column of Yk */
                            Yk [1] -= wsub [1] ;
                            cg_scale0 (Yk+2, wsub+2, -CGONE, memk-2) ;
                        }
                        if ( d0isg && !memk_is_mem ) DenseCol1 = FALSE ;
                        else                         DenseCol1 = TRUE ;

                        Subspace = TRUE ;
                        /* reset subspace skipping to 0 in test
                           for InvariantSpace */
                        SubSkip = 0 ;
                        IterSubRestart = 0 ;
                        nsub = memk ; /* dimension of subspace */
                        nrestartsub = (int) (((CGFLOAT) nsub)*
                                                            Parm->restart_fac) ;
                        mp_begin = mlast ;
                        memk_begin = nsub ;
                        SkFlast = (SkFstart+nsub-1) % mem ;
                        cg_copy0 (gsubtemp, gsub, nsub) ;
                        /* Rk contains the sk for subspace, initialize Sk = Rk*/
                        cg_copy (Sk, Rk, (int) mem*nsub) ;
                    }
                    else
                    {
                       /* deleted 3 spots where checked for ratio < 10 */
                       if ( (IterRestart == nrestart) ||
                          ((IterQuad == qrestart) && (IterQuad != IterRestart)))
                       {
                           Restart = TRUE ;
                       }
                    }
                } /* done checking the full space */
            } /* done using the memory */

            if ( Subspace ) /* compute search direction in subspace */
            {
                Stat->IterSub++ ;
                if ( PrintLevel >= 2 )
                {
                    printf(" Subspace Iteration, IterSub: %ld\n",
                          (LONG) Stat->IterSub);
                }

                /* set x = xnew and g = gnew */
                cg_copy (x, xnew, n) ;
                cg_copy (g, gnew, n) ;
                if ( Aexists )
                {
                    cg_copy(gproj, gnewproj,n) ;
                }

                if ( Restart ) /*restart in subspace*/
                {
                    Restart = FALSE ;
                    RestartAfterHessian = FALSE ;
                    IterRestart = 0 ;
                    IterSubRestart = 0 ;
                    IterQuad = 0 ;
                    mp_begin = -1 ;
                    memk_begin = 0 ;
                    memk = 0 ;
                    scale = (CGFLOAT) 1 ;

                    if ( PrintLevel >= 1 )
                    {
                        printf ("REstart subspace in cg\n") ;
                    }

                    /* search direction d = -Zk gsub, gsub = Zk' g, dsub = -gsub
                       => d =  Zk dsub = SkF (Rk)^{-1} dsub */
                    cg_scale0 (dsub, gsubtemp, -CGONE, nsub) ;
                    cg_copy0 (gsub, gsubtemp, nsub) ;
                    cg_copy0 (vsub, dsub, nsub) ;
                    cg_trisolve (vsub, Rk, mem, nsub, 1) ;
                    /* rearrange and store in wsub */
                    mp = SkFlast ;
                    j = nsub - (mp+1) ;
                    cg_copy0 (wsub, vsub+j, mp+1) ;
                    cg_copy0 (wsub+(mp+1), vsub, j) ;
                    cg_matvec (D, SkF, wsub, nsub, n, 1) ;

                    /* !!! CHECK THIS !!! */
                   /* dnorm2 = pasa_dot(D, D, n) ; */
                    dnorm2 = gsubnorm2 ;

#ifdef PASA
                    /* transform direction back to x-space */
                    pasa_null_project (d, D, NULL, FALSE, pasacom) ;
                    /* delete next ?? */
                    /* pasacom->cg_bb_est = scale ;*/
#endif

                    /* derivative in new search direction without penalty term*/
                    dphi0 = cg_dot (d, g, n) ;
                }
                else  /* continue in subspace without restart */
                {
                    mlast_sub = (mp_begin + IterSubRestart) % mem ;

                    if (IterSubRestart > 0 ) /* not 1st iteration in subspace */
                    {
                        /* add new column to Yk memory,
                           calculate yty, Sk, Yk and SkYk */
                        spp = mlast_sub*mem ;
                        cg_scale0 (Sk+spp, dsub, alpha, nsub) ;
                        /* yty = (gsubtemp-gsub)'(gsubtemp-gsub),
                           set gsub = gsubtemp */
                        cg_Yk (Yk+spp, gsub, gsubtemp, &yty, nsub) ;
                        SkYk [mlast_sub] = alpha*(dphi - dphi0) ;
                    }
                    else
                    {
                        yty = cg_dot0 (Yk+mlast_sub*mem,
                                           Yk+mlast_sub*mem, nsub) ;
                    }
                    if ( yty > CGZERO )
                    {
                        scale = SkYk [mlast_sub]/yty ;
#ifdef PASA
                        pasacom->cg_bb_est = scale ;
#endif
                    }

                    /* calculate gsubtemp = H gsub */
                    mp = mlast_sub ;
                    /* memk = size of the L-BFGS memory in subspace */
                    memk = CGMIN (memk_begin + IterSubRestart, mem) ;
                    l1 = CGMIN (IterSubRestart, memk) ;
                    /* l2 = number of triangular columns in Yk with a zero */
                    l2 = memk - l1 ;
                    /* l1 = number of dense column in Yk (excluding first) */
                    l1++ ;
                    l1 = CGMIN (l1, memk) ;

                    /* process dense columns */
                    for (j = 0; j < l1; j++)
                    {
                        mpp = mp*mem ;
                        t = cg_dot0 (Sk+mpp, gsubtemp, nsub)/SkYk[mp] ;
                        tau [mp] = t ;
                        /* update gsubtemp -= t*Yk+mpp */
                        cg_daxpy0 (gsubtemp, Yk+mpp, -t, nsub) ;
                        mp-- ;
                        if ( mp < 0 ) mp = mem-1 ;
                    }

                    /* process columns from triangular (Hessenberg) matrix */
                    for (j = 1; j < l2; j++)
                    {
                        mpp = mp*mem ;
                        t = cg_dot0 (Sk+mpp, gsubtemp, mp+1)/SkYk[mp] ;
                        tau [mp] = t ;
                        /* update gsubtemp -= t*Yk+mpp */
                        if ( mp == 0 && DenseCol1 )
                        {
                            cg_daxpy0 (gsubtemp, Yk+mpp, -t, nsub) ;
                        }
                        else
                        {
                            cg_daxpy0(gsubtemp, Yk+mpp, -t, CGMIN(mp+2,nsub));
                        }
                        mp-- ;
                        if ( mp < 0 ) mp = mem-1 ;
                    }
                    cg_scale0 (gsubtemp, gsubtemp, scale, nsub) ;

                    /* process columns from triangular (Hessenberg) matrix */
                    for (j = 1; j < l2; j++)
                    {
                        mp++ ;
                        if ( mp == mem ) mp = 0 ;
                        mpp = mp*mem ;
                        if ( mp == 0 && DenseCol1 )
                        {
                            t = cg_dot0 (Yk+mpp, gsubtemp, nsub)/SkYk[mp] ;
                        }
                        else
                        {
                            t = cg_dot0 (Yk+mpp, gsubtemp,
                                                   CGMIN(mp+2,nsub))/SkYk[mp] ;
                        }
                        /* update gsubtemp += (tau[mp]-t)*Sk+mpp */
                        cg_daxpy0 (gsubtemp, Sk+mpp, tau [mp] - t, mp+1) ;
                    }

                    /* process dense columns */
                    for (j = 0; j < l1; j++)
                    {
                        mp++ ;
                        if ( mp == mem ) mp = 0 ;
                        mpp = mp*mem ;
                        t = cg_dot0 (Yk+mpp, gsubtemp, nsub)/SkYk [mp] ;
                        /* update gsubtemp += (tau[mp]-t)*Sk+mpp */
                        cg_daxpy0 (gsubtemp, Sk+mpp, tau [mp] - t, nsub) ;
                    } /* done computing H gsubtemp */

                    /* compute d = Zk dsub = SkF (Rk)^{-1} dsub */
                    cg_scale0 (dsub, gsubtemp, -CGONE, nsub) ;
                    cg_copy0 (vsub, dsub, nsub) ;

                    cg_trisolve (vsub, Rk, mem, nsub, 1) ;

                    /* rearrange and store in wsub */
                    mp = SkFlast ;
                    j = nsub - (mp+1) ;
                    cg_copy0 (wsub, vsub+j, mp+1) ;
                    cg_copy0 (wsub+(mp+1), vsub, j) ;
                    cg_matvec (D, SkF, wsub, nsub, n, 1) ;

                    /* !!! CHECK THIS !!! */
                    /*  dnorm2 = pasa_dot(D, D, n) ; */
                    dnorm2 = cg_dot0 (dsub, dsub, nsub);

#ifdef PASA
                    /* transform direction back to x-space */
                    pasa_null_project (d, D, NULL, FALSE, pasacom) ;
#endif

                    /* derivative in new search direction without penalty term*/
                    dphi0 = cg_dot (d, g, n) ;
                }
            } /* end of subspace search direction */
            else  /* compute the search direction in the full space */
            {
                if ( PrintLevel >= 1 )
                {
                    printf ("compute search direction in full space\n") ;
                }
                /* check to see whether cg should be restated */
                if ( Restart )
                {
                    if ( PrintLevel >= 1 )
                    {
                        printf ("Restart CG/limited memory, use_memory: %i\n",
                                 use_memory) ;
                    }
                    Restart = FALSE ;
                    RestartAfterHessian = FALSE ;
                    IterRestart = 0 ;
                    IterQuad = 0 ;
                    scale = (CGFLOAT) 1 ;

                    /* set x = xnew, g = gnew */
                    cg_copy (x, xnew, n) ;
                    cg_copy (g, gnew, n) ;

                    /*gnewproj and gnorm2 were already computed */
                    if ( use_memory )
                    {
                        if ( Aexists && !RestartAfterHessian )
                        {
                            cg_copy (gproj, gnewproj, n) ;
                        }
                    }
                    else
                    {
                        gnorm2 = cg_dot (gproj, gproj, n) ;
                    }

                    /* D is the search direction before the final projection */
                    cg_scale (D, gproj, -CGONE, n) ;

#ifdef PASA
                    /* transform direction back to x-space */
                    pasa_null_project (d, D, NULL, FALSE, pasacom) ;
                    dnorm2 = cg_dot (d, d, n) ;
#else
                    dnorm2 = gnorm2 ;
#endif
                    /* derivative in search direction without penalty term */
                    dphi0  = cg_dot (g, d, n) ;
                    beta = CGZERO ;
                }
                else if ( !FirstFull ) /* ordinary cg update without restart */
                {
                    if ( PrintLevel >= 2 )
                    {
                        printf ("ordinary cg update in limited memory\n") ;
                    }
                    /* set x = xnew */
                    cg_copy (x, xnew, n) ;

                    if ( use_memory == FALSE )
                    {
                        gnorm2 = cg_dot(gnewproj, gnewproj, n);
                    }

                    /* compute: ykPyk  = (newproj - oldproj)*(newproj - oldproj)
                                gkPyk  =  newproj           *(newproj - oldproj)
                       update: oldproj = newproj */
                    cg_update_beta (gproj, gnewproj, &gkPyk, &ykPyk, n) ;

                    /* For bound constrained problems, we set g = gnew in
                       update_beta since g = gproj and gnew = gnewproj. When
                       linear constraints present, we need to set g = gnew. */
                    if ( Aexists == TRUE )
                    {
                        cg_copy (g, gnew, n) ;
                    }

                    dkyk = dphi - dphi0 ;
                    scale = alpha*dnorm2/dkyk ;

                    if ( Parm->AdaptiveTheta )
                    {
                        t = 2. - CGONE/(0.1*QuadTrust + CGONE) ;
                    }
                    else
                    {
                        t = Parm->theta ;
                    }
                    beta = (gkPyk - t*dphi*ykPyk/dkyk)/dkyk ;

                    /* faster: initialize dnorm2 = gnorm2 at start, then
                       dnorm2 = gnorm2 + beta**2*dnorm2 - 2.*beta*dphi
                       gnorm2 = ||g_{k+1}||^2
                       dpi = g_{k+1}' d_k */

                    /* lower bound for beta is BetaLower*d_k'g_k/ ||d_k||^2 */
                    beta = CGMAX (beta, Parm->BetaLower*dphi0/dnorm2) ;

                    /* update search direction D = -gproj + beta*Dold */
                    cg_update_d (D, gproj, beta, n) ;

#ifdef PASA
                    /* transform direction back to x-space */
                    pasa_null_project (d, D, NULL, FALSE, pasacom) ;
                    pasacom->cg_bb_est = scale ;
#endif
                    dnorm2 = cg_dot(d, d, n) ;
                    /* derivative in new search direction without penalty term*/
                    dphi0 = cg_dot (d, g, n) ;
                }
                else /* FirstFull = TRUE, precondition after leaving subspace*/
                {
                    if ( PrintLevel >= 1 )
                    {
                        printf ("first full iteration after subspace exit\n") ;
                    }
                    /* set x = xtemp */
                    cg_copy (x, xnew, n) ;

                    /* compute: ykPyk  = (newproj - oldproj)*(newproj - oldproj)
                                gkPyk  =  newproj           *(newproj - oldproj)
                       update: oldproj = newproj */
                    cg_update_beta (gproj, gnewproj, &gkPyk, &ykPyk, n) ;

                    /* For bound constrained problems, we set g = gnew in
                       update_beta since g = gproj and gnew = gnewproj. When
                       linear constraints present, we need to set g = gnew. */
                    if ( Aexists == TRUE )
                    {
                        cg_copy (g, gnew, n) ;
                    }

                    mlast_sub = (mp_begin + IterSubRestart) % mem ;
                    /* save Sk */
                    spp = mlast_sub*mem ;
                    cg_scale0 (Sk+spp, dsub, alpha, nsub) ;
                    /* calculate yty, save Yk, set gsub = gsubtemp */

                    cg_Yk (Yk+spp, gsub, gsubtemp, &yty, nsub) ;
                    ytg = cg_dot0  (Yk+spp, gsub, nsub) ;
                    t = alpha*(dphi - dphi0) ;
                    SkYk [mlast_sub] = t ;

                    /* scale = t/ykyk ; */
                    if ( yty > CGZERO )
                    {
                        scale = t/yty ;
                    }

                    /* calculate gsubtemp = H gsub */
                    mp = mlast_sub ;
                    /* memk = size of the L-BFGS memory in subspace */
                    memk = CGMIN (memk_begin + IterSubRestart, mem) ;
                    l1 = CGMIN (IterSubRestart, memk) ;
                    /* l2 = number of triangular columns in Yk with a zero */
                    l2 = memk - l1 ;
                    /* l1 = number of dense column in Yk (excluding first) */
                    l1++ ;
                    l1 = CGMIN (l1, memk) ;

                    /* process dense columns */
                    for (j = 0; j < l1; j++)
                    {
                        mpp = mp*mem ;
                        t = cg_dot0 (Sk+mpp, gsubtemp, nsub)/SkYk[mp] ;
                        tau [mp] = t ;
                        /* update gsubtemp -= t*Yk+mpp */
                        cg_daxpy0 (gsubtemp, Yk+mpp, -t, nsub) ;
                        mp-- ;
                        if ( mp < 0 ) mp = mem-1 ;
                    }
                    /* process columns from triangular (Hessenberg) matrix */
                    for (j = 1; j < l2; j++)
                    {
                        mpp = mp*mem ;
                        t = cg_dot0 (Sk+mpp, gsubtemp, mp+1)/SkYk[mp] ;
                        tau [mp] = t ;
                        /* update gsubtemp -= t*Yk+mpp */
                        if ( mp == 0 && DenseCol1 )
                        {
                            cg_daxpy0 (gsubtemp, Yk+mpp, -t, nsub) ;
                        }
                        else
                        {
                            cg_daxpy0 (gsubtemp, Yk+mpp,-t, CGMIN(mp+2,nsub));
                        }
                        mp-- ;
                        if ( mp < 0 ) mp = mem-1 ;
                    }
                    cg_scale0 (gsubtemp, gsubtemp, scale, nsub) ;

                    /* process columns from triangular (Hessenberg) matrix */
                    for (j = 1; j < l2; j++)
                    {
                        mp++ ;
                        if ( mp == mem ) mp = 0 ;
                        mpp = mp*mem ;
                        if ( mp == 0 && DenseCol1 )
                        {
                            t = cg_dot0 (Yk+mpp, gsubtemp, nsub)/SkYk[mp] ;
                        }
                        else
                        {
                            t = cg_dot0 (Yk+mpp, gsubtemp,
                                                    CGMIN(mp+2,nsub))/SkYk[mp] ;
                        }
                        /* update gsubtemp += (tau[mp]-t)*Sk+mpp */
                        cg_daxpy0 (gsubtemp, Sk+mpp, tau [mp] - t, mp+1) ;
                    }

                    /* process dense columns */
                    for (j = 0; j < l1; j++)
                    {
                        mp++ ;
                        if ( mp == mem ) mp = 0 ;
                        mpp = mp*mem ;
                        t = cg_dot0 (Yk+mpp, gsubtemp, nsub)/SkYk [mp] ;
                        /* update gsubtemp += (tau[mp]-t)*Sk+mpp */
                        cg_daxpy0 (gsubtemp, Sk+mpp, tau [mp] - t, nsub) ;
                    } /* done computing H gsubtemp */

                    /* compute beta */
                    dkyk = dphi - dphi0 ;
#ifdef PASA
                    pasacom->cg_bb_est = alpha*dnorm2/dkyk ;
#endif

                    if ( Parm->AdaptiveTheta )
                    {
                        t = 2. - CGONE/(0.1*QuadTrust + CGONE) ;
                    }
                    else
                    {
                        t = Parm->theta ;
                    }
                    /* Theoretically t1 = ykPyk - yty */
                    t1 = CGMAX(ykPyk-yty, CGZERO) ;
                    if ( ykPyk > CGZERO )
                    {
                        scale = (alpha*dkyk)/ykPyk ; /* = sigma */
                    }
                    beta = scale*((gkPyk - ytg) - t*dphi*t1/dkyk)/dkyk ;
                    /* beta = MAX (beta, Parm->BetaLower*dphi0/dnorm2) ; */
                    beta = CGMAX (beta, Parm->BetaLower*(dphi0*alpha)/dkyk) ;

                    /* compute search direction
                       d = -Zk (H - sigma)ghat - sigma g + beta d
                       Note: d currently contains last 2 terms so only need
                       to add the Zk term. Above gsubtemp = H ghat */

                    /* form vsub = sigma ghat - H ghat = sigma ghat - gsubtemp*/
                    cg_scale0 (vsub, gsubtemp, -CGONE, nsub) ;
                    cg_daxpy0 (vsub, gsub, scale, nsub) ;
                    cg_trisolve (vsub, Rk, mem, nsub, 1) ;

                    /* rearrange vsub and store in wsub */
                    mp = SkFlast ;
                    j = nsub - (mp+1) ;
                    cg_copy0 (wsub, vsub+j, mp+1) ;
                    cg_copy0 (wsub+(mp+1), vsub, j) ;


                    /* save old direction d in gnew */
                    cg_copy (gnew, D, n) ;

                    /* D = Zk (sigma - H)ghat */
                    cg_matvec (D, SkF, wsub, nsub, n, 1) ;

                    /* incorporate the new g and old d terms in new d */
                    cg_daxpy (D, gproj, -scale, n) ;
                    cg_daxpy (D, gnew, beta, n) ;
#ifdef PASA
                    /* transform direction back to x-space */
                    pasa_null_project (d, D, NULL, FALSE, pasacom) ;
#endif
                    dnorm2 = cg_dot(d, d, n) ;
                    /* derivative in new search direction without penalty term*/
                    dphi0 = cg_dot (d, g, n) ;

                }   /* end of preconditioned step */
            }
        }   /* search direction has been computed */

        t = CGZERO ;
#ifdef PASA
        if ( use_penalty == TRUE )
        {
            t = pasacom->dp = cg_dot (d, gpen, n) ;
        }
        err = pasacom->e ;
#else
        err = gnorm ;
#endif

        /* restart when search direction not a descent direction */
        if ( dphi0 + t  >= CGZERO )
        {
            if ( PrintLevel >= 1 )
            {
                printf ("REstart CG due to bad search direction\n") ;
            }

#ifdef PASA
            /* if CG encounters a bad search direction in PASA, then
               return to gradproj since we may have reached the optimum
               point over the active manifold */
            return (PASA_GRAD_PROJ) ;
#endif

            IterRestart = 0 ;
            IterQuad = 0 ;
            mlast = -1 ;
            memk = 0 ;
            beta = CGZERO ;
            Restart = FALSE ;

            /* D is the search direction before the final projection */
            cg_scale (D, gproj, -CGONE, n) ;

#ifdef PASA
            /* d is the final search direction */
            pasa_null_project (d, D, NULL, FALSE, pasacom) ;
#endif
            /* compute square of 2-norm of d (used in bb formula) */
            dnorm2 = cg_dot (d, d, n) ;

            /* derivative in search direction without penalty term */
            dphi0  = cg_dot (g, d, n) ;

            t = CGZERO ;
#ifdef PASA
            if ( use_penalty == TRUE )
            {
                t = pasacom->dp = cg_dot (gpen, d, n) ;
            }
#endif

            if ( dphi0 + t > CGZERO )
            {
                status = CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION ;
                XXCG(wrapup) (status, FALSE, PASA_CG_COM) ;
                return (status) ;
            }
        }

        /* test for slow convergence */
        if ( (f < fbest) || (err < gbest) )
        {
            nslow = 0 ;
            if ( f < fbest )
            {
                fbest = f ;
            }
            if (err < gbest )
            {
                gbest = err ;
            }
        }
        else
        {
            nslow++ ;
        }
        if ( nslow > slowlimit )
        {
            status = CG_NO_COST_OR_GRADIENT_IMPROVEMENT ;
            XXCG(wrapup)(status, TRUE, PASA_CG_COM) ;
            return (status) ;
        }

#ifndef NDEBUG
        if ( Parm->debug )
        {
            if ( f > cgcom->f0 + Parm->debugtol*Ck )
            {
                Stat->newf = f ;
                Stat->oldf = cgcom->f0 ;
                status = CG_DEBUGGER_IS_ON_AND_FUNCTION_VALUE_INCREASES ;
                XXCG(wrapup) (status, TRUE, PASA_CG_COM) ;
                return (status) ;
            }
        }
#endif

    }
    status = XXCG(wrapup) (CG_ITERATIONS_EXCEED_MAXITS, FALSE, PASA_CG_COM) ;
    return (status) ;
}

/* ==========================================================================
   === cg_setup =============================================================
   ==========================================================================
    Generate a pointer to a CGdata structure with default values for all
    parameters, and NULL or EMPTY values for all inputs. Return NULL pointer if
    there is not enough memory.
   ========================================================================== */
CGdata * XXCG(setup) (void)
{
    int status ;
    /* Initialize data structure */
    CGdata *Data ;

    status = CG_OK ;
    /* Allocate memory for parameters and statistics, initialize parameters */
    Data = (CGdata  *) cg_malloc (&status, 1, sizeof (CGdata)) ;
    Data->Parm = (CGparm *) cg_malloc (&status, 1, sizeof(CGparm)) ;
    Data->Stat = (CGstat *) cg_malloc (&status, 1, sizeof(CGstat)) ;
    /* create the cg structure for storing the Hessian of a quadratic */
    Data->H = (SOPT_matrix *) cg_malloc (&status, 1, sizeof (SOPT_matrix)) ;
    if ( status == SOPT_OUT_OF_MEMORY )
    {
        Data = NULL ;
        return (Data) ;
    }

    cg_default (Data->Parm) ;
    sopt_matrix_default (Data->H) ;

    /* set remaining values in the CGdata structure to EMPTY (= -1) or NULL */
    Data->x = NULL ;
    Data->n = EMPTY ;

    /* evaluate the function f at x: value (f, x, n) */
    Data->value = NULL ;

    /* evaluate the gradient at x: grad  (g, x, n) */
    Data->grad = NULL ;

    /* evaluate the function f and gradient g at x: valgrad (f, g, x, n) */
    /* NULL => use value & grad */
    Data->valgrad = NULL ;

    /* If the Hessian-based implementation of cg_descent is employed, then the
       hessian function is for evaluating the Hessian of the objective
       at a given point. See cg_descent.h and sopt.h for the format. */
    Data->hessian = NULL ;

    /* Optional estimate for number of nonzeros in the Hessian.
       EMPTY => determine hnnz by evaluating Hessian at a point. */
    Data->hnnz = EMPTY ;

    /* size n, the linear term for a quadratic or linear objective is c'*x */
    Data->c = NULL ;

   /* If the problem is a QP, then instead of providing routines
       to evaluate the objective function and its gradient, the user
       could simply provide the Hessian of the objective
       and any linear term in the objective (see inputs Data->H
       and data->c above).  If the linear term is not given, then
       it is taken to be zero.  Whenever the objective value or its
       gradient is needed by cg_descent, it will be computed internally
       using the provided Hessian.  See sopt.h for the SOPT_matrix formats. */

    /* For a QP, another option is to provide a routine for computing the
       product between the Hessian and a vector.  hprod (p, x, n) evaluates
       the Hessian times a vector. Here the problem dimension n and x are
       given, while the routine should generate p = H*x where H is the
       Hessian of the objective. */
    Data->hprod = NULL ;

    /* Data->Work can be used for a pointer to a real work space */
    Data->Work = NULL ; /* NULL => let cg_descent allocate real work space. */

    /* If the user does not provide a pointer to x, then an array of size
       n is malloc'd and set to zero.  pasa_terminate will free any allocated
       memory. The x_created pointer is used internally by CG_DESCENT. */
    Data->x_created = NULL ;

    /* If the sparse matrix in H was created by cg_descent, then
       H_created is set to TRUE. */
    Data->H_created = FALSE;

    /* Return pointer to CGdata structure */
    return (Data) ;
}

/* ==========================================================================
   === cg_terminate =========================================================
   ==========================================================================
    Free allocated memory associated with CGdata structure
   ========================================================================== */
void XXCG(terminate)
(
    CGdata **DataHandle
)
{
    CGdata *Data ;
    Data = *DataHandle ;
    /* Free memory allocated in pasa_setup */
    cg_free (Data->Parm) ;
    cg_free (Data->Stat) ;
    if ( Data->x_created != NULL )
    {
        if ( Data->x_created == Data->x )
        {
            cg_free (Data->x) ;
            Data->x = Data->x_created = NULL ;
        }
        else /* user created x */
        {
            cg_free (Data->x_created) ;
            Data->x_created = NULL ;
        }
    }
    if ( Data->H_created )
    {
        cg_free (Data->H->p) ;
        cg_free (Data->H->i) ;
        cg_free (Data->H->x) ;
    }
    cg_free (Data->H) ;
    cg_free (Data) ;
}

int XXCG(wrapup)
(
    int       status,
    int   fastreturn, /* T => return after printing status */
#ifdef PASA
    PASAcom *pasacom,
#endif
    CGcom     *cgcom
)
{
    CGparm *Parm ;
    CGstat *Stat ;

    Stat = cgcom->Stat ;
    Parm = cgcom->Parm ;

    Stat->status = status ;
    Stat->f = cgcom->f ;
    if ( Parm->PrintStatus )
    {
        cg_print_status (cgcom->cgdata) ;
        printf ("\n") ;
        fflush (stdout) ;
    }
    if ( fastreturn ) return (status) ;

    if ( Parm->PrintStat )
    {
        cg_print_stat (cgcom->cgdata) ;
    }
#ifdef TIME
    printf ("time hessian: %19.10f\n", cgcom->loop_hessian) ;
    printf ("time    init: %19.10f\n", cgcom->loop_init) ;
    printf ("time     prp: %19.10f\n", cgcom->loop_prp) ;
    printf ("time cgtrust: %19.10f\n", cgcom->loop_cgtrust) ;
    printf ("time SSmumps: %19.10f\n", cgcom->loop_SSmumps) ;
    printf ("time  factor: %19.10f\n", cgcom->loop_factor) ;
    printf ("time sstrust: %19.10f\n", cgcom->loop_sstrust) ;
#endif

#ifdef PASA
    int const use_napheap = pasacom->use_napheap ;
    int const use_pproj = pasacom->use_pproj ;
    int const PrintLevel = Parm->PrintLevel ;
    int const loExists = pasacom->loExists ;
    int const hiExists = pasacom->hiExists ;
    int const Aexists = pasacom->Aexists ;
    /* the cg error is the local error */
    Stat->err = pasacom->e ;

    /* if the penalty term is used and the objective is nonquadratic,
       then the unpenalized objective is stored in f_orig */
    if ( (cgcom->QuadCost == FALSE) && (pasacom->use_penalty == TRUE) )
    {
        pasacom->f = pasacom->f_orig ;
    }

    /* If status = CG_HITS_BOUNDARY, then the iterate hit the boundary
       of the feasible region.  In this case, we compress the problem,
       update the factorization, and restart either CG or active grad_proj */
    if ( status == CG_HITS_BOUNDARY )
    {
        int blk, blks ;
        PASAINT cols, i, iri, j, k, l, m, ncol, nrow, ncoldel, ni,
                p, pp, q, *Ap, *Anz, *Ai, *ATp, *ATi, *col_start, *ir ;
        PASAFLOAT t, xj, *a, *Ax, *ATx, *b, *bl, *bu ;
        PPprob *Prob ;
        PPwork *Work ;
#ifndef NDEBUG
        int NoRows, *ActiveCols ;
#endif
        int             downdate_ok = TRUE ;
        PASAFLOAT const maxstep     = pasacom->maxstep ;
        PASAFLOAT const maxbndstep  = pasacom->maxbndstep ;
        PASAFLOAT const maxconstep  = pasacom->maxconstep ;
        PASAINT   const maxbndindex = pasacom->maxbndindex ;
        PASAINT   const maxconindex = pasacom->maxconindex ;
        PASAINT                  nf = pasacom->nf ;
        PASAINT                  nc = pasacom->nc ;
        PASAINT              *ifree = pasacom->ifree ;
        PASAINT         *bound_cols = pasacom->bound_cols ;
        PASAFLOAT                *x = pasacom->x ;
        PASAFLOAT                *g = pasacom->g ;
        PASAFLOAT               *lo = pasacom->lo ;
        PASAFLOAT               *hi = pasacom->hi ;

        if ( use_pproj )
        {
            Prob = pasacom->ppcom->Prob ;
            Work = pasacom->ppcom->Work ;
            ATp = Work->ATp ;
            ATi = Work->ATi ;
            ATx = Work->ATx ;
            Ap = Prob->Ap ;
            Anz = Prob->Anz ;
            Ai = Prob->Ai ;
            Ax = Prob->Ax ;
            ni = Prob->ni ;
            ir = Work->ir ;
            b  = pasacom->b ;
            bl = pasacom->bl ;
            bu = pasacom->bu ;
            ASSERT (Work->ncoladd == 0) ;
            ASSERT (Work->ncoldel == 0) ;
            ASSERT (Work->nrowadd == 0) ;
            ASSERT (Work->nrowdel == 0) ;
            ncol = Prob->ncol ;
            nrow = Prob->nrow ;
            ASSERT (nrow == pasacom->nrow) ;
        }
        else if ( use_napheap == TRUE )
        {
            a = pasacom->nap_a ;
        }

        if ( PrintLevel >= 1 )
        {
            printf ("step: %e", maxstep) ;
            if ( pasacom->Bounds )
            {
                if ( maxbndindex != 0 )
                {
                    j = (maxbndindex < 0) ? -(maxbndindex+1) : maxbndindex - 1 ;
                    printf (" maxbndstep: %e maxbndindex: %ld ucol: %ld",
                        maxbndstep, (LONG) maxbndindex, (LONG) ifree [j]) ;
                }
            }
            if ( Aexists )
            {
                if ( maxconindex != 0 )
                {
                    if ( use_pproj )
                    {
                        j = (maxconindex < 0) ? -maxconindex : maxconindex ;
                        printf (" maxconstep: %e maxconindex: %ld row: %ld",
                               maxconstep, (LONG) pasacom->maxconindex,
                                           (LONG) Prob->ineq_row [j]) ;
                    }
                    else /* use_napheap */
                    {
                        printf (" maxconstep: %e maxconindex: %ld row: 0",
                               maxconstep, (LONG) pasacom->maxconindex) ;
                    }
                }
            }
            printf ("\n") ;
        }

#if 0
        /* If the problem has only bound constraints (no inequalities)
           and an expansion step was performed, then remove all the bound
           variables from the problem. Currently, not used. */
        if ( !Aexists && (maxstep > maxbndstep) )
        {
            cols = 0 ;
            for (j = 0; j < nf; j++)
            {
                PASAFLOAT const Xj = x [j] ;
                /* store bound Xj in userx and remove from the problem */
                if ( loExists && (Xj == lo [j]) )
                {
                    k = ifree [j] ;
                    pasacom->userx [k] = Xj ;
                    /* store bound columns in user coordinates */
                    bound_cols [nc] = -(k + 1) ;
                    nc++ ;
                }
                else if ( hiExists && (Xj == hi [j]) )
                {
                    k = ifree [j] ;
                    pasacom->userx [k] = Xj ;
                    /* store bound columns in user coordinates */
                    bound_cols [nc] = k ;
                    nc++ ;
                }
                else /* Xj is free and retained */
                {
                    /* retain associated x, g, lo, hi, and ifree entries */
                    x [cols] = Xj ;
                    g [cols] = g [j] ;
                    if ( loExists ) lo [cols] = lo [j] ;
                    if ( hiExists ) hi [cols] = hi [j] ;
                    ifree [cols] = ifree [j] ;
                    cols++ ;
                }
            }
            pasacom->nf = cols ;
            pasacom->nc = nc ;
            ASSERT (nc+cols == pasacom->ncol) ;
            if ( PrintLevel >= 1 )
            {
                printf ("cg wrapup nc: %ld nf: %ld\n",
                         (LONG) nc, (LONG) cols) ;
            }
        }
        else if ( maxstep < maxconstep )
#endif
        /* maxstep <= maxconstep => a variable reached its bound, remove it */
        if ( maxbndstep <=  maxconstep )
        {
            if ( maxbndindex < 0 )
            {
                j = -(maxbndindex+1) ;
                /* set x [j] = lo [j] to account for rounding errors in x */
                k = ifree [j] ;
                xj = pasacom->userx [k] = lo [j] ;
                bound_cols [nc] = -(k + 1) ;
                if ( use_napheap )
                {
                    PASAFLOAT const aj = a [j] ;
                    pasacom->nap_a2 -= aj*aj ;
                    t = aj*xj ;
                    pasacom->nap_bl -= t ;
                    pasacom->nap_bu -= t ;
                }
            }
            else /* maxbndindex > 0 */
            {
                j = maxbndindex - 1 ;
                /* set x [j] = hi [j] to account for rounding errors in x */
                k = ifree [j] ;
                xj = pasacom->userx [k] = hi [j] ;
                bound_cols [nc] = k ;
                if ( use_napheap )
                {
                    PASAFLOAT const aj = a [j] ;
                    pasacom->nap_a2 -= aj*aj ;
                    t = aj*xj ;
                    pasacom->nap_bl -= t ;
                    pasacom->nap_bu -= t ;
                }
            }
            pasacom->nc++ ;
#ifndef NDEBUG
            if ( use_pproj )
            {
                NoRows = (pasacom->ppcom->Work->RLinkUp [nrow] == nrow) ;
                ActiveCols = pasacom->ppcom->Check->ActiveCols ;
            }
#endif
            /* delete the jth element of x, g, lo, hi, and ifree
               if napheap is used, then delete the jth element of a */
            pasacom->nf-- ; /* a free variable is being bound */
            nf = pasacom->nf ;
            k =  j + 1 ;
            l = nf - j ; /* number of elements to move down */
            sopt_copyi_noblas (ifree+j, ifree+k, l) ;
            sopt_copyx_noblas (    x+j,     x+k, l) ;
            sopt_copyx_noblas (    g+j,     g+k, l) ;
            if ( use_napheap )sopt_copyx_noblas ( a+j,  a+k, l) ;
            if ( loExists   ) sopt_copyx_noblas (lo+j, lo+k, l) ;
            if ( hiExists   ) sopt_copyx_noblas (hi+j, hi+k, l) ;
#ifndef NDEBUG
            if ( use_pproj && !NoRows )
            {
                sopt_copyi_noblas (ActiveCols+j, ActiveCols+k, l) ;
            }
#endif
            ASSERT (pasacom->nc + nf == pasacom->ncol) ;
            if ( use_napheap )
            {
                if ( pasacom->nap_a2 <= .01*pasacom->nap_a2save )
                {
                    pasacom->nap_a2save = pasacom->nap_a2 = pasa_dot (a, a, nf);
                }
            }
            if ( PrintLevel >= 1 )
            {
                printf ("cg wrapup nc: %ld nf: %ld\n",
                         (LONG) pasacom->nc, (LONG) nf) ;
            }

#ifndef NOPPROJ
            /* if A exists, compress the matrix and update the factorization*/
            if ( use_pproj )
            {
                /* update the factorization if the column has active rows
                   and the column has not yet been removed */
                if ( Anz [j] && pasacom->update_factor )
                {
#ifndef NDEBUG
                    /* prevents error message in pproj check routines */
                    Work->ib [j] = 1 ;
#endif
                    ncoldel = 1 ;
                    i = ncol - ncoldel ;
                    Work->ColmodList [i] = j ;
                    Work->ColmodFlag [j] = i ;
                    Work->ncoldel = ncoldel ;
                    Work->Annz -= Anz [j] ;
                    pasacom->ppcom->Parm = pasacom->pprojparm ;
                    downdate_ok = pproj_modcol (pasacom->ppcom, 0, 0, -1, NULL,
                                  Ap+j, Anz+j, NULL, NULL, 1) ;
                    if ( !downdate_ok )
                    {
                        if ( PrintLevel )
                        {
                            printf ("CG: failure in update of factorization "
                                    "after deleting a column\n") ;
                        }
                        /* the matrix should be refactored from scratch*/
                        Work->fac = FALSE ;
                    }
                }
                /* After removing column j, shift elements of ib down */
                sopt_copy_int (Work->ib+j, Work->ib+k, l) ;

                /* adjust the right side of the inequalities to account
                   for the bound variable */
                if ( xj != CGZERO )
                {
                    q = Ap [j+1] ;
                    for (p = Ap [j]; p < q; p++)
                    {
                        i = Ai [p] ;
                        iri = ir [i] ;
                        t = Ax [p]*xj ;
                        if ( iri == 0 ) /* equality constraint */
                        {
                            b [i] -= t ;
                        }
                        else /* inequality constraint */
                        {
                            if ( iri < 0 )
                            {
                                iri = -iri ;
                            }
                            else if ( iri > ni )
                            {
                                iri -= ni ;
                            }
                            bl [iri] -= t ;
                            bu [iri] -= t ;
                        }
                    }
                }

                /* When columns are deleted, the column start of each block is
                   reduced by the number of proceeding columns that were
                   deleted.  Find the first block with column start > j.
                   The following code terminates since col_start [blks] =
                   number of columns in matrix. */
                col_start = Work->col_start ;
                blks = Work->blks ;
                for (blk = 1; blk <= blks; blk++)
                {
                    if ( col_start [blk] > j )
                    {
                        break ;
                    }
                }
                l = col_start [blk] ;
                m = Ap [j] ;
                p = Ap [j+1] ;
                cols = j ;
                for (k = j+1; k < ncol; k++)
                {
                    if ( k == l )
                    {
                        col_start [blk] = cols ;
                        blk++ ;
                        while ( col_start [blk] == l )
                        {
                            col_start [blk] = cols ;
                            blk++ ;
                        }
                        l = col_start [blk] ;
                    }
                    q = Ap [k+1] ;
                    Anz [cols] = Anz [k] ;

                    /* save column k in matrix */
                    Ap [k] = Ap [cols] + q - p ;
                    for (; p < q; p++)
                    {
                        Ai [m] = Ai [p] ;
                        Ax [m] = Ax [p] ;
                        m++ ;
                    }
                    cols = k ;
                }
                Prob->ncol = cols ;
                Work->A->ncol = cols ;
                Work->A->nzmax = Ap [cols] ;
                col_start [blk] = cols ;
                while ( blk < blks )
                {
                    blk++ ;
                    col_start [blk] = cols ;
                }

                /* Update AT matrix. This can be done much more
                   efficiently if AT uses ATnz */
                pp = 0 ;
                p = 0 ;
                for (i = 0; i < nrow; i++)
                {
                    q = ATp [i+1] ;
                    for (; p < q; p++)
                    {
                        k = ATi [p] ;
                        if ( k < j ) /* keep the element unchanged */
                        {
                            ATi [pp] = k ;
                            ATx [pp] = ATx [p] ;
                            pp++ ;
                        }
                        else if ( k > j )  /* decrease element by 1 */
                        {
                            ATi [pp] = k - 1 ;
                            ATx [pp] = ATx [p] ;
                            pp++ ;
                        }
                        /* if k = j, then omit the element */
                    }
                    ATp [i+1] = pp ;
                }
                pasa_copyi (Work->AFTp, ATp, nrow+1) ;
#ifndef NDEBUG
                pasa_checkA (pasacom, TRUE) ;
#endif
            }
#endif
        }
        else /* an inequality has become active */
        {
#ifndef NOPPROJ
            int  *sol_to_blk ;
            CGINT nr, nrowadd, p1, p2, p2m1, q1, q2, row,
                 *bound_rows, *ineq_row, *RLinkDn, *RLinkUp ;

            bound_rows = pasacom->bound_rows ;

            if ( use_napheap )
            {
                /* maxbndindex = -1 if at lower bound, +1 if at upper bound */
                bound_rows [0] = pasacom->nap_constraint = 2 ;
                if ( maxconindex < 0 ) /* lower bound active */
                {
                    pasacom->nap_bu = pasacom->nap_bl ;
                }
                else                  /* upper bound active */
                {
                    pasacom->nap_bl = pasacom->nap_bu ;
                }
                pasacom->nr = 1 ;
            }
            else /* use_pproj == TRUE */
            {
#ifndef NDEBUG
                pasa_checkA (pasacom, TRUE) ;
#endif
                ineq_row = Prob->ineq_row ;
                nr = pasacom->nr ;
                if ( maxconindex < 0 ) /* lower bound active */
                {
                    k = -maxconindex ;
                    row = ineq_row [k] ;
                    b [row] = bl [k] ;
                    bound_rows [nr] = -(row+1) ;
                }
                else                  /* upper bound active */
                {
                    k = maxconindex ;
                    row = ineq_row [k] ;
                    b [row] = bu [k] ;
                    bound_rows [nr] = row + 1 ;
                }
                ir [row] = 0 ;
                if ( pasacom->update_factor )
                {
#ifndef NDEBUG
                    int *ActiveRows = pasacom->ppcom->Check->ActiveRows ;
                    if ( ActiveRows [row] == 1 )
                    {
                        printf ("row %ld is currently active, cannot re-add "
                                "in cg_descent\n", (LONG) row) ;
                        pasa_error (-1, __FILE__, __LINE__, "stop") ;
                    }
                    ActiveRows [row] = 1 ;
#endif
                    Work->nactive++ ;
                    Work->ATnz += ATp [row+1] - ATp [row] ;
                    nrowadd = 1 ;
                    i = nrow - nrowadd ;
                    Work->RowmodList [i]   = row ;
                    Work->RowmodFlag [row] = i ;
                    Work->nrowadd = nrowadd ;

                    /* add row to linked lists */
                    RLinkDn = Work->RLinkDn ;
                    RLinkUp = Work->RLinkUp ;

                    /* this can be done more efficiently if a significant number
                       of rows are active by searching ir above and below the
                       new active row */
                    i = nrow ;
                    while ( (i = RLinkUp [i]) < row )

                    ASSERT (i != row) ;
                    j = RLinkDn [i] ;
                    RLinkDn [i] = row ;
                    RLinkDn [row] = j ;
                    RLinkUp [j] = row ;
                    RLinkUp [row] = i ;

                    /* Update A to put the new active row in the active part of
                       A. Recall that the active rows in column must be in
                       increasing order. */
                    q = ATp [row+1] ;
                    for (p = ATp [row]; p < q; p++)
                    {
                        j = ATi [p] ;
                        p1 = Ap [j] ;
                        p2 = q1 = p1 + Anz [j] ;
                        q2 = Ap [j+1] ;
                        /* find row in the inactive part of matrix */
                        for (; p2 < q2; p2++)
                        {
                            if ( Ai [p2] == row )
                            {
                                /* move 1st inactive row to location of
                                   row in column */
                                Ai [p2] = Ai [q1] ;
                                t = Ax [p2] ;
                                Ax [p2] = Ax [q1] ;
                                break ;
                            }
                        }
                        ASSERT (p2 < q2) ;
                        /* find where to insert row in active part of column */
                        p2 = q1 ;
                        for (; p2 > p1; )
                        {
                            p2m1 = p2 - 1 ;
                            if ( Ai [p2m1] < row ) /* row goes at location p2 */
                            {
                                break ;
                            }
                            Ai [p2] = Ai [p2m1] ;
                            Ax [p2] = Ax [p2m1] ;
                            p2 = p2m1 ;
                        }
                        Ai [p2] = row ;
                        Ax [p2] = t ;
                        Anz [j]++ ;
                    }
                    pasacom->ppcom->Parm = pasacom->pprojparm ;
                    downdate_ok = pproj_modrow (pasacom->ppcom, 0, TRUE, FALSE,
                              +1, NULL, NULL, NULL, NULL) ;
                    if ( !downdate_ok )
                    {
                        if ( PrintLevel )
                        {
                            printf ("CG: failure in update of factorization "
                                    "after adding a row\n") ;
                        }
                        /* the matrix should be refactored from scratch*/
                        Work->fac = FALSE ;
                    }
                }
                sol_to_blk = Work->sol_to_blk ;
                for (i = 1; i < k; i++)
                {
                    row = ineq_row [i] ;
                    if ( ir [row] > ni ) ir [row]-- ; /* ni decreased by 1 */
                }
                /* The bounds as well as the block numbers associated with each
                   inequality need to be shifted down. */
                PASAINT km = k ;
                for (k = k+1; k <= ni; k++)
                {
                    bl [km] = bl [k] ;
                    bu [km] = bu [k] ;
                    sol_to_blk [km] = sol_to_blk [k] ;
                    row = ineq_row [k] ;
                    /* subtract 2 for inactive rows since ni and row decreased*/
                    if ( ir [row] > ni ) ir [row] -= 2 ;
                    /* adjust by 1 for active rows since row changed by 1 */
                    else if ( ir [row] > 0 ) ir [row]-- ;
                    else                     ir [row]++ ;
                    ineq_row [km] = row ;
                    km = k ;
                }
                ineq_row [ni] = nrow ;
                sol_to_blk [km] = Work->blks ;
                Prob->ni-- ;    /* one fewer inequality to check */
                pasacom->nr++ ; /* one more bound inequality */
            }
#endif
        }

#if 0
        if ( downdate_ok ) /* check the error if the factorization valid */
        {
            pasa_null_project (NULL, pasacom->g, NULL, TRUE, pasacom) ;
            /* it might save some time to only perform this call to
               checktol when pasacom->e <= pasacom->testtol (similar to
               what is done in the main program). */
            if ( PrintLevel >= 1 )
            {
                printf ("Local error: %e testtol: %e Global error: %e\n",
                        pasacom->e, pasacom->testtol, pasacom->E) ;
            }

            status = pasa_cg_checktol (pasacom->x, pasacom, cgcom) ;
            if ( (status == PASA_ERROR_TOLERANCE_SATISFIED) ||
                 (status == PASA_GRAD_PROJ) )
            {
                return (status) ;
            }
        }
#endif
        /* multiply pasacom->switchfactor by pasaparm->switchdecay, but do
           not let switchfactor drop beneath switchlower */
        PASAFLOAT switchlower = pasacom->pasaparm->switchlower ;
        PASAFLOAT switchdecay = pasacom->pasaparm->switchdecay ;
        t = pasacom->switchfactor ;
        /* if t = switchlower, then do nothing */
        if ( (t > switchlower) && (maxstep > PASAZERO) )
        {
            t *= switchdecay ;
            if ( t < switchlower )
            {
                t = switchlower/pasacom->switchfactor ;
                pasacom->testtol *= t ;
                pasacom->switchfactor = switchlower ;
            }
            else
            {
                pasacom->switchfactor = t ;
                pasacom->testtol *= switchdecay ;
            }
            /* do not let testtol drop below grad_tol */
            pasacom->testtol = PASAMAX (pasacom->testtol, pasacom->grad_tol) ;
        }
        /* if the downdate failed, the matrix needs to be refactored; we go
           to the activeGP to do this */
        if ( (pasacom->pasaparm->use_activeGP == TRUE) || !downdate_ok )
        {
            return (PASA_ACTIVE_GP) ;
        }
        return (PASA_CG_DESCENT) ; /* restart CG */
        /* this ends the ifdef PASA: either terminate, return to grad_proj,
                                     return to active grad_proj, or restart CG*/
    }
    /* possible values for status: PASA_ACTIVE_GP, PASA_GRAD_PROJ,
                  PASA_ERROR_TOLERANCE_SATISFIED or an error flag */
#else
    /* this starts the pure CG routine, remove any regularization from f */
    if ( (cgcom->QuadCost == TRUE) && (Parm->QPshift > CGZERO) )
    {
        CGFLOAT t = cg_dot (cgcom->x, cgcom->x, cgcom->n) ;
        Stat->f -= 0.5*t*t*Parm->QPshift ;
    }
    /* free any memory malloc'd by cg */
    if ( (cgcom->cgdata->Work == NULL) && (cgcom->work_created != NULL) )
    {
        cg_free (cgcom->work_created) ;
    }
    if ( cgcom->deriv_mode > 1 )
    {
#ifdef USE_MUMPS
        DMUMPS_STRUC_C *id = cgcom->id ;
        if ( id != NULL )
        {
            id->job = JOB_END ;
            dmumps_c (id) ;
            cg_free (id) ;
        }
#endif
#ifdef USE_HARWELL
        cholmod_common *cmm = cgcom->cmm ;
        cholmod_sparse *Achol = cgcom->Achol ;
        LBL_Data *lbl = cgcom->lbl ;
        FILE *lblfile = cgcom->lblfile ;
        CGINT *Cparent = cgcom->Cparent ;
        CGINT *Cmember = cgcom->Cmember ;
        if ( cmm != NULL )
        {
            CHOLMOD (finish) (cmm) ;
            cg_free (cmm) ;
        }

        if ( Achol != NULL )
        {
            if ( Achol->p != NULL )
            {
                cg_free (Achol->p) ;
                cg_free (Achol->i) ;
            }
            cg_free (Achol) ;
            if ( Cparent != NULL )
            {
                cg_free (Cparent) ;
                cg_free (Cmember) ;
            }
        }

        if ( lbl != NULL )
        {
            /* did not use malloc'd space in irn and jcn */
            lbl->irn = lbl->ma57->irn = NULL ;
            lbl->jcn = lbl->ma57->jcn = NULL ;
            LBL_Finalize (cgcom->lbl) ;
        }
        if ( lblfile != NULL ) fclose (lblfile) ;
#endif
        CGFLOAT *basisWork = cgcom->basisWork ;
        CGINT *basisiWork = cgcom->basisiWork ;
        SOPT_Cmatrix *C = cgcom->C ;
        cg_free (cgcom->speed) ;
        if ( basisWork != NULL )
        {
            cg_free (basisWork) ;
            cg_free (basisiWork) ;
        }
        /* free work space for CG compressed matrix C */
        if ( C->p == C->malloc_p )
        {
            cg_free (C->p) ;
            cg_free (C->i) ;
            cg_free (C->x) ;
        }
        if ( C->rows == C->malloc_rows )
        {
            cg_free (C->rows) ;
            cg_free (C->cols) ;
            cg_free (C->vals) ;
        }
        if ( C->map   != NULL ) cg_free (C->map) ;
        if ( C->Smap  != NULL ) cg_free (C->Smap) ;
        if ( C->col   != NULL ) cg_free (C->col) ;
        if ( C->Hwork != NULL ) cg_free (C->Hwork) ;
        if ( C->rhs   != NULL ) cg_free (C->rhs) ;
        cg_free (C) ;

        /* if problem not QuadCost, then return user's cgdata->H in an
           initialized state */
        if ( !cgcom->QuadCost )
        {
            sopt_matrix_default (cgcom->cgdata->H) ;
        }
    }
#endif
    return (status) ;
}

int XXCG(evaluate)
(
    CGFLOAT alpha_good, /* a value of alpha for which function is finite */
    CGFLOAT     *Alpha, /* stepsize along the search direction */
    char         *what, /* fg = eval func & grad, g = grad only,f = func only */
#ifdef PASA
    PASAcom *pasacom,
#endif
    CGcom       *Com
)
{
    int const QuadCost = Com->QuadCost ;
#ifdef PASA
    int status ;
    /* pasa is used */
    status = pasa_evaluate (alpha_good, Alpha, pasacom, what) ;
    if ( QuadCost )
    {
        if ( *Alpha = -2 ) Com->dq = pasacom->dq ;
    }
    else /* not a quadratic */
    {
        if ( !strcmp (what, "f") )
        {
            Com->f = pasacom->f ;
        }
        else if ( !strcmp (what, "g") )
        {
            Com->df = pasacom->df ;
        }
        else /* fg */
        {
            Com->f = pasacom->f ;
            Com->df = pasacom->df ;
        }
    }
    if ( status == PASA_OK )
    {
        return (CG_CONTINUE) ;
    }
    else
    {
        return (CG_FUNCTION_NAN_OR_INF) ;
    }
#else
    /* pasa is NOT being used */
    int i ;
    CGINT n ;
    CGFLOAT alpha, df, t, *d, *f, *g, *x, *gnew, *xnew ;
    CGparm *Parm ;
    CGstat *Stat ;

    Parm = Com->Parm ;
    alpha = *Alpha ;
    n = Com->n ;
    x = Com->x ;
    g = Com->g ;
    f = &(Com->f) ;
    Stat = Com->Stat ;

    /* initial evaluation of objective and gradient */
    if ( alpha == CGZERO )
    {
        Stat->ngrad++ ;
        if ( QuadCost ) /* quadratic objective */
        {
            if ( Com->hprod_format )
            {
                cg_builtin_hprod (g, x, n, Com->H->p, Com->H->i, Com->H->x) ;
            }
            else Com->hprod (g, x, n) ; /* g = Q*x */

            /* regularize the Hessian if QPshift != 0 */
            if ( Com->QPshift )    /* g = g + QPshift*x */
            {
                cg_step (g, g, x, Com->QPshift, n) ;
            }

            *f = 0.5*cg_dot (g, x, n) ;      /* f = .5*x'Qx */
            if ( Com->c != NULL )
            {
                *f += cg_dot(x, Com->c, n) ; /* f = .5*x'*Qx + c'*x */
                /* store gradient in userg */
                cg_step (g, g, Com->c, CGONE, n) ; /* g = Qx+c */
            }
        }
        else
        {
            Stat->nfunc++ ;
            if ( Com->valgrad != NULL )
            {
                Com->valgrad (f, g, x, n) ;
            }
            else
            {
                Com->grad  (g, x, n) ;
                Com->value (f, x, n) ;
            }
        }
        return (CG_CONTINUE) ;
    }

    d = Com->d ;
    /* alpha = -1 or -2: evaluate Hessian times direction for a quadratic
       alpha =       -2: also evaluate cost change at x + d */
    if ( alpha < CGZERO )
    {
        Stat->ngrad++ ;
        if ( Com->hprod_format )
        {
            cg_builtin_hprod (Com->Qd, d, n, Com->H->p, Com->H->i, Com->H->x) ;
        }
        else Com->hprod (Com->Qd, d, n) ;
        /* add the regularization term if QPshift != 0 */
        if ( Com->QPshift )
        {
            cg_step (Com->Qd, Com->Qd, d, Com->QPshift, n) ;
        }
        if ( alpha == -2 ) /* also evaluate cost change */
        {
            Com->dq = cg_dot (Com->g, d, n) + 0.5*cg_dot (Com->Qd, d, n) ;
        }
        return (CG_CONTINUE) ;
    }

    /* take a step along the search direction */
    gnew = Com->gnew ;
    xnew = Com->xnew ;
    cg_step (xnew, x, d, alpha, n) ;
    if ( !strcmp (what, "fg") )     /* compute function and gradient */
    {
        Stat->nfunc++ ;
        Stat->ngrad++ ;
        if ( Com->valgrad != NULL )
        {
            Com->valgrad (f, gnew, xnew, n) ;
        }
        else
        {
            Com->grad  (gnew, xnew, n) ;
            Com->value (f, xnew, n) ;
        }
        df = Com->df = cg_dot (gnew, d, n) ;
        if ( (*f == *f) && (*f < CGINF) && (*f > -CGINF) &&
             (df == df) && (df < CGINF) && (df > -CGINF) )
        {
            return (CG_CONTINUE) ;
        }
        t = Parm->cg_infdecay ;
        for (i = 0; i < Parm->cg_ninf_tries; i++)
        {
            Stat->nfunc++ ;
            Stat->ngrad++ ;
            alpha = alpha_good + t*(alpha - alpha_good) ;
            t *= Parm->cg_infdecay_rate ;
            cg_step (xnew, x, d, alpha, n) ;
            if ( Com->valgrad != NULL )
            {
                Com->valgrad (f, gnew, xnew, n) ;
            }
            else
            {
                Com->grad  (gnew, xnew, n) ;
                Com->value (f, xnew, n) ;
            }
            df = Com->df = cg_dot (gnew, d, n) ;
            if ( (*f == *f) && (*f < CGINF) && (*f > -CGINF) &&
                 (df == df) && (df < CGINF) && (df > -CGINF) )
            {
                break ;
            }
        }
    }
    else if ( !strcmp (what, "f") ) /* compute function value */
    {
        Stat->nfunc++ ;
        Com->value (f, xnew, n) ;
        if ( (*f == *f) && (*f < CGINF) && (*f > -CGINF) )
        {
            return (CG_CONTINUE) ;
        }
        t = Parm->cg_infdecay ;
        for (i = 0; i < Parm->cg_ninf_tries; i++)
        {
            Stat->nfunc++ ;
            alpha = alpha_good + t*(alpha - alpha_good) ;
            t *= Parm->cg_infdecay_rate ;
            cg_step (xnew, x, d, alpha, n) ;
            Com->value (f, xnew, n) ;
            if ( (*f == *f) && (*f < CGINF) && (*f > -CGINF) )
            {
                break ;
            }
        }
    }
    else
    {
        Stat->ngrad++ ;
        Com->grad (gnew, xnew, n) ;
        df = Com->df = cg_dot (gnew, d, n) ;
        if ( (df == df) && (df < CGINF) && (df > -CGINF) )
        {
            return (CG_CONTINUE) ;
        }
        t = Parm->cg_infdecay ;
        for (i = 0; i < Parm->cg_ninf_tries; i++)
        {
            Stat->ngrad++ ;
            alpha = alpha_good + t*(alpha - alpha_good) ;
            t *= Parm->cg_infdecay_rate ;
            cg_step (xnew, x, d, alpha, n) ;
            Com->grad (gnew, xnew, n) ;
            df = Com->df = cg_dot (gnew, d, n) ;
            if ( (df == df) && (df < CGINF) && (df > -CGINF) )
            {
                break ;
            }
        }
    }
    if ( i == Parm->cg_ninf_tries )
    {
        return (CG_FUNCTION_NAN_OR_INF) ;
    }
    else
    {
        *Alpha = alpha ;
        return (CG_CONTINUE) ;
    }
#endif
}

/* =========================================================================
   ==== cg_Wolfe ===========================================================
   =========================================================================
   Check whether the Wolfe or the approximate Wolfe conditions are satisfied
   Return:
       CG_WOLFE_OK     (Wolfe or approximate Wolfe conditions satisfied)
       CG_WOLFE_NOT_OK (Wolfe condition does not hold)
   ========================================================================= */
int XXCG(Wolfe)
(
    CGFLOAT  alpha, /* stepsize */
    CGFLOAT      f, /* function value associated with stepsize alpha */
    CGFLOAT   dphi, /* derivative value associated with stepsize alpha */
    CGcom   *cgcom
)
{
    if ( dphi >= cgcom->wolfe_lo )
    {
        /* test original Wolfe conditions */
        if ( f - cgcom->f0 <= alpha*cgcom->wolfe_hi )
        {
            if ( cgcom->Parm->PrintLevel >= 2 )
            {
                printf ("Wolfe conditions hold\n") ;
/*              printf ("wolfe f: %25.15e f0: %25.15e df: %25.15e\n",
                         f, cgcom->f0, dphi) ;*/
            }
            return (CG_WOLFE_OK) ;
        }
        /* test approximate Wolfe conditions */
        else if ( cgcom->approxstep )
        {
/*          if ( cgcom->Parm->PrintLevel >= 2 )
            {
                printf ("f:    %e fpert:    %e ", f, cgcom->fpert) ;
                if ( f > cgcom->fpert ) printf ("(fail)\n") ;
                else                  printf ("(OK)\n") ;
                printf ("dphi: %e hi bound: %e ", dphi, cgcom->awolfe_hi) ;
                if ( dphi > cgcom->awolfe_hi ) printf ("(fail)\n") ;
                else                         printf ("(OK)\n") ;
            }*/
            if ( (f <= cgcom->fpert) && (dphi <= cgcom->awolfe_hi) )
            {
                if ( cgcom->Parm->PrintLevel >= 2 )
                {
                    printf ("Approximate Wolfe conditions hold\n") ;
/*                  printf ("f: %25.15e fpert: %25.15e dphi: %25.15e awolf_hi: "
                       "%25.15e\n", f, cgcom->fpert, dphi, cgcom->awolfe_hi) ;*/
                }
                return (CG_WOLFE_OK) ;
            }
        }
    }
    return (CG_WOLFE_NOT_OK) ;
}

/* =========================================================================
   ==== cg_line ============================================================
   =========================================================================
   Approximate Wolfe line search routine
   Return:
       CG_WOLFE_OK
       CG_WOLFE_CONDITIONS_NOT_SATISFIED
       CG_SLOPE_ALWAYS_NEGATIVE
       CG_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS
   ========================================================================= */
int XXCG(line)
(
    int       repeat, /* TRUE => Wolfe search failed, retry using approxstep */
#ifdef PASA
    PASAcom *pasacom,
#endif
    CGcom     *cgcom
)
{
    int approxstep, AvoidFeval, iter, ngrow, PrintLevel, qa, qa0, qb, qb0,
        status, toggle ;
    CGFLOAT alpha, a, a1, a2, b, bmin, da, db, d0, d1, d2, df, f, fa, fb,
           a0, b0, da0, db0, fa0, ExpandSafe, fb0, maxstep, rho, t, width ;
    char *s1, *s2, *fmt1, *fmt2, *fmt3 ;
    CGparm *Parm ;

    approxstep = cgcom->approxstep ;
    Parm = cgcom->Parm ;
    PrintLevel = Parm->PrintLevel ;
    if ( PrintLevel >= 2 )
    {
        if ( approxstep )
        {
            printf ("Approximate Wolfe line search, AvoidFeval: %i "
                    "repeat: %i QuadOK: %i\n",
                     cgcom->AvoidFeval, repeat, cgcom->QuadOK) ;
            printf ("=============================\n") ;
        }
        else
        {
            printf ("Wolfe line search, AvoidFeval: %i\n"
                    "repeat: %i QuadOK: %i\n",
                     cgcom->AvoidFeval, repeat, cgcom->QuadOK) ;
            printf ("=================\n") ;
        }
    }
    maxstep = cgcom->maxstep ;
    status = CG_CONTINUE ;

    /* We check the Wolfe conditions at alpha if it was computed
       by a quadratically accurate formula. However, below we can change
       this decision if the value of the derivative implies that the
       Wolfe conditions cannot be satisfied.

       As the cg_descent converges, the function values typically
       approach a constant value. When this happens, the cubic interpolation
       step in the line search loses its accuracy and it is better to
       use a secant step based on the derivative of the objective function.
       AvoidFeval is TRUE when the function values are close together.  */

    AvoidFeval = cgcom->AvoidFeval ;

    /* repeat = FALSE => implies that this is the initial attempt to
       satisfy the Wolfe conditions */
    if ( repeat == FALSE )
    {
        alpha = cgcom->alpha ;
        /* It is advantageous to have both the function value and
           derivative except when QuadOK = F and AvoidFeval = T.
           In this case, we will not check the Wolfe conditions, and
           since the function values have converged, the line search
           is based on the derivative. */
        if ( (cgcom->QuadOK == FALSE) && (AvoidFeval == TRUE) )
        {
            a0 = alpha ;
            /* if the value of the objective at alpha exists, we
               nonetheless make note of it */
            if ( (cgcom->f_at >= CGZERO) && (alpha == cgcom->f_at) )
            {
                qb = TRUE ;
            }
            else
            {
                qb = FALSE ;
            }
            /* evaluate derivative if not present */
            if ( (cgcom->df_at < CGZERO) || (alpha != cgcom->df_at) )
            {
                status = XXCG(evaluate) (CGZERO, &alpha, "g", PASA_CG_COM) ;
                if ( status == CG_FUNCTION_NAN_OR_INF )
                {
                    return (status) ;
                }
                /* note that the value of alpha returned by the evaluation
                   routine may be different from the input alpha due to
                   nan's or inf's */
                if ( qb && (alpha != a0 ) ) /* evaluate f and g at same point */
                {
                    status = XXCG(evaluate) (CGZERO, &alpha, "fg", PASA_CG_COM);
                    if ( status == CG_FUNCTION_NAN_OR_INF )
                    {
                        return (status) ;
                    }
                }
            }
            /* even though function values may have converged since
               AvoidFeval is TRUE, we will evaluate the function value
               if the derivative at alpha is much smaller than the
               derivative at the starting point. */
            if ( !qb && cgcom->df < Parm->BigDfactor*cgcom->df0 )
            {
                qb = TRUE ;
                a0 = alpha ; /* the value of alpha where g evaluated */
                status = XXCG(evaluate) (CGZERO, &alpha, "f", PASA_CG_COM) ;
                if ( status == CG_FUNCTION_NAN_OR_INF )
                {
                    return (status) ;
                }
                if ( alpha != a0 ) /* evaluate f and g at same point */
                {
                    status = XXCG(evaluate) (CGZERO, &alpha, "fg", PASA_CG_COM);
                    if ( status == CG_FUNCTION_NAN_OR_INF )
                    {
                        return (status) ;
                    }
                }
            }
        }
        else /* we want to have both the objective value and derivative */
        {
            qb = TRUE ;
            a0 = alpha ;
            /* if function has already been computed, only need derivative */
            if ( (cgcom->f_at >= CGZERO) && (alpha == cgcom->f_at) )
            {
                status = XXCG(evaluate) (CGZERO, &alpha, "g", PASA_CG_COM) ;
            }
            /* if derivative has already been computed, only need value */
            else if ( (cgcom->df_at >= CGZERO) && (alpha == cgcom->df_at) )
            {
                status = XXCG(evaluate) (CGZERO, &alpha, "f", PASA_CG_COM) ;
            }
            /* otherwise both value and derivative are needed */
            else
            {
                a0 = -CGONE ;
                status = XXCG(evaluate) (CGZERO, &alpha, "fg", PASA_CG_COM) ;
            }
            /* since the evaluation routine decreases alpha when an nan
               or inf is encountered, it is possible in the first 2 cases
               above that the function and derivative were computed at
               different points */
            if ( (a0 >= CGZERO) && (a0 != alpha) )
            {
                status = XXCG(evaluate) (CGZERO, &alpha, "fg", PASA_CG_COM) ;
            }
            /* if the function was nan or inf everywhere it was evaluated,
               then return with error */
            if ( status == CG_FUNCTION_NAN_OR_INF )
            {
                return (status) ;
            }
        }

        if ( qb == TRUE )
        {
            fb = cgcom->f ;
            /* make any necessary adjustments for a Wolfe line search */
            if ( !approxstep )
            {
                fb -= alpha*cgcom->wolfe_hi ;
            }

            /* if the objective value has converged, then it is better to
               use gradient to satisfy stopping conditions */
            if (fabs (cgcom->f-cgcom->f0) <= Parm->CostConverge*fabs (cgcom->f))
            {
                AvoidFeval = TRUE ;
            }
            else
            {
                AvoidFeval = FALSE ;
            }
        }
    }
    else /* this is second attempt to satisfy Wolfe conditions */
    {
        cgcom->df = cgcom->savedf ;
        qb = cgcom->saveqb ;
        alpha = cgcom->savealpha ;
        if ( qb == TRUE )
        {
            fb = cgcom->f = cgcom->savefb ;
        }
    }
    b = alpha ;

    if ( approxstep )
    {
        db = cgcom->df ;
        d0 = da = cgcom->df0 ;
    }
    else
    {

/*printf ("df: %e\n", cgcom->df) ;*/
/*printf ("wolfe_hi: %e\n", cgcom->wolfe_hi) ;*/
        db = cgcom->df - cgcom->wolfe_hi ;
        d0 = da = cgcom->df0 - cgcom->wolfe_hi ;

        /* save data in case we repeat the line search using approxstep */
        cgcom->savedf = cgcom->df ;
        cgcom->savealpha = alpha ;
        cgcom->saveqb = qb ;
        if ( qb == TRUE )
        {
            cgcom->savefb = cgcom->f ;
        }
    }

    a = CGZERO ;
    a1 = CGZERO ;
    d1 = d0 ;
    fa = cgcom->f0 ; /* alpha = 0 so no adjustment for Wolfe step */
    if ( PrintLevel >= 1 )
    {
        fmt1 = "%9s %2s a: %13.6e b: %13.6e fa: %13.6e fb: %13.6e "
               "da: %13.6e db: %13.6e\n" ;
        fmt2 = "%9s %2s a: %13.6e b: %13.6e fa: %13.6e fb:  x.xxxxxxxxxx "
               "da: %13.6e db: %13.6e\n" ;
        fmt3 = "%9s %2s a: %13.6e b: %13.6e fa:  x.xxxxxxxxxx "
               "fb: %13.6e da: %13.6e db: %13.6e\n" ;
        if ( cgcom->QuadOK ) s2 = "OK" ;
        else               s2 = "" ;
        if ( qb ) printf (fmt1, "start    ", s2, a, b, fa, fb, da, db);
        else      printf (fmt2, "start    ", s2, a, b, fa, da, db) ;
    }

    /* if a quadratic interpolation step performed, then f was evaluated at b,
       check Wolfe conditions */
    if ( (cgcom->QuadOK == TRUE) && (cgcom->f <= cgcom->f0) )
    {
        status = XXCG(Wolfe) (b, cgcom->f, cgcom->df, cgcom) ;
        if ( status == CG_WOLFE_OK )
        {
#ifdef PASA
            if ( alpha == maxstep )
            {
                return (CG_HITS_BOUNDARY) ;
            }
#endif
            return (status) ; /* Wolfe conditions hold */
        }
    }

    /* For a nice function, the Wolfe line search might terminate above
       in the quadratic interpolation step.  If it does not terminate there
       and a full Wolfe line search is performed, then set cg->Wolfe to TRUE */
    if ( !approxstep )
    {
        cgcom->Wolfe = TRUE ;
    }

    /*Find initial interval [a,b] such that
      da <= 0, db >= 0, fa <= fpert = [(f0 + eps*fabs (f0)) or (f0 + eps)] */
    rho = Parm->rho ;
    ExpandSafe = Parm->ExpandSafe ;
    ngrow = 1 ;
    qa = TRUE ;
    while ( db < CGZERO ) /* need db >= 0 */
    {
        /* evaluate function at b if not yet evaluated there */
        if ( AvoidFeval == FALSE )
        {
            if ( !qb )
            {
                b0 = b ;
                status = XXCG(evaluate) (a, &b, "f", PASA_CG_COM) ;
                /* Note that the initial value of is reduced in the evaluation
                   routine if nan or infinite function value encountered. */
                if ( status == CG_FUNCTION_NAN_OR_INF )
                {
                    return (status) ;
                }
                if ( b0 != b )
                {
                    status = XXCG(evaluate) (a, &b, "fg", PASA_CG_COM) ;
                    if ( status == CG_FUNCTION_NAN_OR_INF )
                    {
                        return (status) ;
                    }
                    db = cgcom->df ;
                }

                /* since the function was evaluated, we again check
                   on its reliability */
                if ( fabs (cgcom->f-cgcom->f0) <=
                                            Parm->CostConverge*fabs (cgcom->f) )
                {
                    AvoidFeval = TRUE ;
                }
                else
                {
                    AvoidFeval = FALSE ;
                }

                if ( approxstep )
                {
                    fb = cgcom->f ;
                    db = cgcom->df ;
                }
                else
                {
                    fb = cgcom->f - b*cgcom->wolfe_hi ;
                    db = cgcom->df - cgcom->wolfe_hi ;
                }
                qb = TRUE ;
            }
            if ( fb > cgcom->fpert ) /* contract interval [a, b] */
            {
                status = XXCG(contract) (&a, &fa, &da, &b, &fb,&db,PASA_CG_COM);
                if ( status == CG_WOLFE_OK )
                {
#ifdef PASA
                    if ( cgcom->alpha == maxstep )
                    {
                        /* cg hits boundary, restart cg_descent*/
                        return (CG_HITS_BOUNDARY) ;
                    }
#endif
                    return (status) ; /* Wolfe conditions hold */
                }
                else if ( status == CG_INTERVAL_OK ) /* db >= 0 */
                {
                    /* we now have fa <= fpert, da <= 0, and db >= 0,
                       go to the line search below */
                    break ;
                }
                else if ( (status == CG_EXCESSIVE_UPDATING_OF_PERT_EPS) ||
                          (status == CG_FUNCTION_NAN_OR_INF) )
                {
                    return (status) ;
                }
                /* else new fpert generated so that fb <= fpert */
            }
        }

        /* fb <= fpert and db < 0 */
#ifdef PASA
        if ( b == maxstep ) /* b cannot grow further in pasa */
        {
            /* cg hits boundary, restart cg_descent*/
            cgcom->alpha = b ;
            return (CG_HITS_BOUNDARY) ; /* restart cg_descent */
        }
#endif

        /* expansion phase */
        ngrow++ ;
        if ( ngrow > Parm->maxsteps )
        {
            return (CG_SLOPE_ALWAYS_NEGATIVE) ;
        }
        /* update interval (a replaced by b) */
        a = b ;
        if ( qb == TRUE )
        {
            fa = fb ;
            qa = TRUE ;
        }
        else
        {
            qa = FALSE ;
        }
        da = db ;
        /* store old values of a and corresponding derivative */
        d2 = d1 ;
        d1 = da ;
        a2 = a1 ;
        a1 = a ;

        bmin = rho*b ;
        if ( (ngrow == 2) || (ngrow == 3) || (ngrow == 6) )
        {
            if ( d1 > d2 )
            {
                /* Use secant to estimate the spot where derivative vanishes
                   since secant moves us to the right. Based on the
                   value of the derivative at a0 = 0, a1 and a2, we can
                   obtain an estimate for the third derivative. If the
                   estimate is positive, then the second derivative of the
                   first derivative is positive, and the secant method takes
                   us past the root of the first derivative. Note that when
                   ngrow = 2, a0 = a1 = 0, so we cannot estimate the second
                   derivative. */
                if ( (ngrow == 2) || ((d1-d2)/(a1-a2) >= (d2-d0)/a2) )
                {
                    b = a1 - (a1-a2)*(d1/(d1-d2)) ;
                }
                else
                {
                    /* The estimate of the second derivative is negative,
                       so the secant method is short of the the true root.
                       SecantAmp is a parameter that increases the step. */
                    {
                        b = a1 - Parm->SecantAmp*(a1-a2)*(d1/(d1-d2)) ;
                    }
                }
                /* safeguard growth */
                t = ExpandSafe*a1 ;
                if ( b > t )
                {
                    b = t ;
#if 0
                    /* rho smaller than ExpandSafe, set rho = ExpandSafe */
                    if ( rho < ExpandSafe )
                    {
                        rho = ExpandSafe ;
                    }
                    ExpandSafe *= Parm->RhoGrow ;
#endif
                }
            }
            else /* secant method makes no sense, must be far from the root,
                    increase rho and hence bmin above */
            {
                rho *= Parm->RhoGrow ;
            }
        }
        else
        {
            rho *= Parm->RhoGrow ;
        }
        b = CGMAX (bmin, b) ;
        if ( b > maxstep )
        {
            b = maxstep ;
        }
        status = XXCG(evaluate) (a, &b, "fg", PASA_CG_COM) ;
        if ( status == CG_FUNCTION_NAN_OR_INF )
        {
            return (status) ;
        }
        cgcom->alpha = b ;
        qb = TRUE ;
        /* since the function was evaluated, we again check on its reliability*/
        if ( fabs (cgcom->f-cgcom->f0) <= Parm->CostConverge*fabs (cgcom->f) )
        {
            AvoidFeval = TRUE ;
        }
        else
        {
            AvoidFeval = FALSE ;
        }
        if ( approxstep )
        {
            fb = cgcom->f ;
            db = cgcom->df ;
        }
        else
        {
            db = cgcom->df - cgcom->wolfe_hi ;
            fb = cgcom->f - b*cgcom->wolfe_hi ;
        }
        if ( PrintLevel >= 2 )
        {
            if ( cgcom->QuadOK ) s2 = "OK" ;
            else                 s2 = "" ;
            if ( qa == TRUE )
            {
                printf (fmt1, "expand   ", s2, a, b, fa, fb, da, db) ;
            }
            else
            {
                printf (fmt3, "expand   ", s2, a, b, fb, da, db) ;
            }
        }
    }

    /* We now have fa <= fpert, da <= 0, db >= 0; hence, the iteration is
       trapped in [a, b] and will not hit the boundary where a new
       constraint becomes active */
    toggle = 0 ;
    width = b - a ;
    for (iter = 0; iter < Parm->maxsteps; iter++)
    {
        /* determine the next iterate */
        if ( (toggle == 0) || ((toggle == 2) && ((b-a) <= width)) )
        {
            cgcom->QuadOK = TRUE ;
            if ( cgcom->UseCubic && qa && qb && (AvoidFeval == FALSE) )
            {
                s1 = "cubic 0  " ;
                alpha = XXCG(cubic) (a, fa, da, b, fb, db) ;
                if ( (alpha <= a) || (alpha >= b) ) /* use secant method */
                {
                    s1 = "secant 0 " ;
                    if      ( -da  < db ) alpha = a - (a-b)*(da/(da-db)) ;
                    else if (  da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                    else                  alpha = -1. ;
                }
            }
            else
            {
                s1 = "secant 0 " ;
                if      ( -da  < db ) alpha = a - (a-b)*(da/(da-db)) ;
                else if (  da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                else                  alpha = -1. ;
            }
            width = Parm->stepdecay*(b - a) ;
        }
        else if ( toggle == 1 ) /* another variation of cubic interpolation */
        {
            cgcom->QuadOK = TRUE ;
            if ( cgcom->UseCubic && (AvoidFeval == FALSE) )
            {
                s1 = "cubic 1  " ;
                if ( (alpha == a) && (a-a0 < b-a) && qa0 )
                {
                    /* a is most recent iterate and is closer to a0 than to b */
                    alpha = XXCG(cubic) (a0, fa0, da0, a, fa, da) ;
                }
                else if ( alpha == a && qb )
                {
                    /* a is most recent iterate and is closer to b than to a0 */
                    alpha = XXCG(cubic) (a, fa, da, b, fb, db) ;
                }
                else if ( alpha == b ) /* b is most recent iterate */
                {
                    /* check if b is closer to b0 than to a */
                    if ( (b0 - b < b - a) && qb0 )
                    {
                        alpha = XXCG(cubic) (b, fb, db, b0, fb0, db0) ;
                    }
                    else if ( qb && qa )      /* b is closer to a than to b0 */
                    {
                        alpha = XXCG(cubic) (a, fa, da, b, fb, db) ;
                    }
                    else
                    {
                        alpha = -1. ;
                    }
                }
                else
                {
                    alpha = -1. ;
                }

                /* if alpha no good, use cubic between a and b */
                if ( (alpha <= a) || (alpha >= b) )
                {
                    if ( qa && qb  )
                    {
                        alpha = XXCG(cubic) (a, fa, da, b, fb, db) ;
                    }
                    else
                    {
                        alpha = -1. ;
                    }
                }

                /* if alpha still no good, use secant method between a and b */
                if ( alpha < CGZERO )
                {
                    s1 = "secant 1 " ;
                    if      ( -da  < db ) alpha = a - (a-b)*(da/(da-db)) ;
                    else if (  da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                    else                  alpha = -1. ;
                }
            }
            else /* ( use secant ) */
            {
                s1 = "secant 1b" ;
                if ( (alpha == a) && (da > da0) )/* use a0 if possible*/
                {
                    alpha = a - (a-a0)*(da/(da-da0)) ;
                }
                else if ( db < db0 )/* try b0 */
                {
                    alpha = b - (b-b0)*(db/(db-db0)) ;
                }
                else /* secant based on a and b */
                {
                    if      ( -da  < db ) alpha = a - (a-b)*(da/(da-db)) ;
                    else if (  da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                    else                  alpha = -1. ;
                }

                if ( (alpha <= a) || (alpha >= b) )
                {
                    if      ( -da  < db ) alpha = a - (a-b)*(da/(da-db)) ;
                    else if (  da != db ) alpha = b - (a-b)*(db/(da-db)) ;
                    else                  alpha = -1. ;
                }
            }
        }
        else
        {
            alpha = .5*(a+b) ; /* use bisection if b-a decays slowly */
            s1 = "bisection" ;
            cgcom->QuadOK = FALSE ;
        }

        if ( (alpha <= a) || (alpha >= b) )
        {
            alpha = .5*(a+b) ;
            s1 = "bisection" ;
            if ( (alpha == a) || (alpha == b) )
            {
                return (CG_WOLFE_CONDITIONS_NOT_SATISFIED) ;
            }
            cgcom->QuadOK = FALSE ; /* bisection was used */
        }

        if ( toggle == 0 ) /* save values for next iteration */
        {
            a0 = a ;
            b0 = b ;
            da0 = da ;
            db0 = db ;
            if ( qa )
            {
                fa0 = fa ;
                qa0 = TRUE ;
            }
            else
            {
                qa0 = FALSE ;
            }
            if ( qb )
            {
                fb0 = fb ;
                qb0 = TRUE ;
            }
            else
            {
                qb0 = FALSE ;
            }
        }

        toggle++ ;
        if ( toggle > 2 ) toggle = 0 ;

        status = XXCG(evaluate) (a, &alpha, "fg", PASA_CG_COM) ;
        if ( status == CG_FUNCTION_NAN_OR_INF )
        {
            return (status) ;
        }

        /* since the function was evaluated, we again check
           on its reliability */
        if ( fabs (cgcom->f-cgcom->f0) <= Parm->CostConverge*fabs (cgcom->f) )
        {
            AvoidFeval = TRUE ;
        }
        else
        {
            AvoidFeval = FALSE ;
        }
        f = cgcom->f ;
        df = cgcom->df ;
        if ( cgcom->QuadOK )
        {
            status = XXCG(Wolfe) (alpha, f, df, cgcom) ;
            if ( status == CG_WOLFE_OK )
            {
                cgcom->alpha = alpha ;
                if ( PrintLevel >= 2 )
                {
                    printf ("             a: %13.6e f: %13.6e df: %13.6e %1s\n",
                             alpha, f, df, s1) ;
                }
                return (status) ;
            }
        }
        if ( !approxstep )
        {
            f -= alpha*cgcom->wolfe_hi ;
            df -= cgcom->wolfe_hi ;
        }
        if ( df >= CGZERO )
        {
            b = alpha ;
            fb = f ;
            db = df ;
            qb = TRUE ;
        }
        else if ( f <= cgcom->fpert )
        {
            a = alpha ;
            fa = f ;
            da = df ;
        }
        else /* df < 0 and f > fpert try to contract interval [a, alpha] */
        {
            CGFLOAT B, fB, dB ;
            B = alpha ;
            fB = f ;
            dB = df ;
            qb = TRUE ;
            /* contract interval [a, alpha] */
            status = XXCG(contract) (&a, &fa, &da, &B, &fB, &dB, PASA_CG_COM) ;
            if ( (status == CG_WOLFE_OK) ||
                 (status == CG_EXCESSIVE_UPDATING_OF_PERT_EPS) ||
                 (status == CG_FUNCTION_NAN_OR_INF) )
            {
                return (status) ;
            }
            if ( status == CG_NEW_PERT )
            {
                a = alpha ;
                fa = f ;
                da = df ;
            }
            else /* interval OK, returned B has positive derivative */
            {
                toggle = 0 ;
                b = B ;
                fb = fB ;
                db = dB ;
            }
        }
        if ( PrintLevel >= 2 )
        {
            if ( cgcom->QuadOK ) s2 = "OK" ;
            else                 s2 = "" ;
            if ( qa && !qb )     printf (fmt2, s1, s2, a, b, fa, da, db) ;
            else if ( qa && qb ) printf (fmt1, s1, s2, a, b, fa, fb, da, db) ;
            else /* !qa && qb */ printf (fmt3, s1, s2, a, b, fb, da, db) ;
        }
    }
    return (CG_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS) ;
}

/* =========================================================================
   ==== cg_contract ========================================================
   =========================================================================
   The input for this routine is an interval [a, b] with the property that
   fa <= fpert, da <= 0, db <= 0, and fb >= fpert. The returned status is
   Return:
       CG_WOLFE_OK            (Wolfe conditions are satisfied)
       CG_INTERVAL_OK         (a subinterval [a, b] is generated with the
                               property that fa <= fpert, da <= 0 and db >= 0)
       CG_NEW_PERT            (db < 0 and a new fpert is generated so that
                               fb < fpert -- hence, b can grow in the line
                               search)
       CG_EXCESSIVE_UPDATING_OF_PERT_EPS
       CG_FUNCTION_NAN_OR_INF

   NOTE: The input arguments are unchanged when status = CG_NEW_PERT since
         we need to move the original b to the right
   ========================================================================= */
int XXCG(contract)
(
    CGFLOAT       *A, /* left side of bracketing interval */
    CGFLOAT      *fA, /* function value at a */
    CGFLOAT      *dA, /* derivative at a */
    CGFLOAT       *B, /* right side of bracketing interval */
    CGFLOAT      *fB, /* function value at b */
    CGFLOAT      *dB, /* derivative at b */
#ifdef PASA
    PASAcom *pasacom,
#endif
    CGcom     *cgcom
)
{
    int approxstep, iter, PrintLevel, toggle, status ;
    CGFLOAT a, alpha, b, old, da, db, df, dold, f0, f, fa, fb, fold, t, width ;
    char *s ;
    CGparm *Parm ;

    approxstep = cgcom->approxstep ;
    Parm = cgcom->Parm ;
    PrintLevel = Parm->PrintLevel ;
    a = *A ;
    fa = *fA ;
    da = *dA ;
    b = *B ;
    fb = *fB ;
    db = *dB ;
    f0 = cgcom->f0 ;
    t = Parm->eps_grow*(fb-f0) ;
    toggle = 0 ;
    width = CGZERO ;
    for (iter = 0; iter < Parm->ncontract; iter++)
    {
        if ( (toggle == 0) || ((toggle == 2) && ((b-a) <= width)) )
        {
            /* cubic based on bracketing interval */
            alpha = XXCG(cubic) (a, fa, da, b, fb, db) ;
            toggle = 0 ;
            width = Parm->stepdecay*(b-a) ;
            cgcom->QuadOK = TRUE ;
        }
        else if ( toggle == 1 )
        {
            cgcom->QuadOK = TRUE ;
            /* cubic based on most recent iterate and smallest value */
            if ( (old < a) && (a - old < b - a) )
            {
                /* a is most recent iterate and a is closer to old than to b */
                alpha = XXCG(cubic) (a, fa, da, old, fold, dold) ;
            }
            else if ( (b < old) && (old - b < b - a) )
            {
                /* b is most recent iterate and b is closer to old than to a */
                alpha = XXCG(cubic) (b, fb, db, old, fold, dold) ;
            }
            else
            {
                alpha = XXCG(cubic) (a, fa, da, b, fb, db) ;
            }
        }
        else
        {
            alpha = .5*(a+b) ; /* use bisection if b-a decays slowly */
            cgcom->QuadOK = FALSE ;
        }

        if ( (alpha <= a) || (alpha >= b) )
        {
            alpha = .5*(a+b) ;
            cgcom->QuadOK = FALSE ; /* bisection was used */
        }

        toggle++ ;
        if ( toggle > 2 ) toggle = 0 ;

        status = XXCG(evaluate) (a, &alpha, "fg", PASA_CG_COM) ;
        if ( status == CG_FUNCTION_NAN_OR_INF )
        {
            return (status) ;
        }

        f = cgcom->f ;
        df = cgcom->df ;

        if ( cgcom->QuadOK )
        {
            status = XXCG(Wolfe) (alpha, f, df, cgcom) ;
            if ( status == CG_WOLFE_OK )
            {
                cgcom->alpha = alpha ;
                return (status) ;
            }
        }
        if ( !approxstep )
        {
            f -= alpha*cgcom->wolfe_hi ;
            df -= cgcom->wolfe_hi ;
        }
        if ( df >= CGZERO ) /* done */
        {
            *B = alpha ;
            *fB = f ;
            *dB = df ;
            *A = a ;
            *fA = fa ;
            *dA = da ;
            return (CG_INTERVAL_OK) ;
        }
        if ( f <= cgcom->fpert ) /* update a using alpha */
        {
            old = a ;
            a = alpha ;
            fold = fa ;
            fa = f ;
            dold = da ;
            da = df ;
        }
        else                     /* update b using alpha */
        {
            old = b ;
            b = alpha ;
            fold = fb ;
            fb = f ;
            dold = db ;
            db = df ;
        }
        if ( PrintLevel >= 2 )
        {
            if ( cgcom->QuadOK ) s = "OK" ;
            else                 s = "" ;
            printf ("contract  %2s a: %13.6e b: %13.6e fa: %13.6e fb: "
                    "%13.6e da: %13.6e db: %13.6e\n", s, a, b, fa, fb, da, db) ;
        }
    }

    /* cg_contract was unable to either satisfy the Wolfe conditions or
       obtain db >= 0. Return to the line search with a larger value for
       fpert. Check to see if the Wolfe conditions are satisfied. */
    f0 = cgcom->f0 ;
    if ( cgcom->PertRule == 1 )
    {
        if ( f0 != CGZERO )
        {
            cgcom->pert_eps = t/fabs (f0) ;
            cgcom->fpert = f0 + t ;
        }
        else
        {
            cgcom->fpert = 2.*fabs (*fB) ;
        }
    }
    else /* PertRule = 0 */
    {
        cgcom->pert_eps = t ;
        cgcom->fpert = f0 + t ;
    }
    /* count the number of times that pert_eps was modified */
    cgcom->neps++ ;

    /* check to see if the Wolfe line search conditions are now satisfied at
       the last iterate alpha */
    status = XXCG(Wolfe) (alpha, f, df, cgcom) ;
    if ( status == CG_WOLFE_OK )
    {
        cgcom->alpha = alpha ;
        return (status) ;
    }

    /* check to see if the Wolfe line search conditions are now satisfied at
       the right side of the original interval */
    status = XXCG(Wolfe) (*B, *fB, *dB, cgcom) ;
    if ( status == CG_WOLFE_OK )
    {
        cgcom->alpha = *B ;
        cgcom->f = *fB ;
        cgcom->df = *dB ;
        return (status) ;
    }

    /* If Wolfe conditions do not hold, then we have *fB > fpert, do not
       modify the returned arguments of cg_contract. */
    if ( cgcom->neps >= Parm->neps )
    {
        return (CG_EXCESSIVE_UPDATING_OF_PERT_EPS) ;
    }
    else
    {
        return (CG_NEW_PERT) ;
    }
}

/* =========================================================================
   ==== cg_cubic ===========================================================
   =========================================================================
   Compute the minimizer of a Hermite cubic. If the computed minimizer
   outside [a, b], return -1 (it is assumed that a >= 0).
   ========================================================================= */
CGFLOAT XXCG(cubic)
(
    CGFLOAT  a,
    CGFLOAT fa, /* function value at a */
    CGFLOAT da, /* derivative at a */
    CGFLOAT  b,
    CGFLOAT fb, /* function value at b */
    CGFLOAT db  /* derivative at b */
)
{
    CGFLOAT c, d1, d2, delta, t, v, w ;
    delta = b - a ;
    if ( delta == CGZERO ) return (a) ;
    v = da + db - 3.*(fb-fa)/delta ;
    t = v*v - da*db ;
    if ( t < CGZERO ) /* complex roots, use secant method */
    {
         if ( fabs (da) < fabs (db) ) c = a - (a-b)*(da/(da-db)) ;
         else if ( da != db )         c = b - (a-b)*(db/(da-db)) ;
         else                         c = -1 ;
         return (c) ;
    }

    if ( delta > CGZERO ) w = sqrt(t) ;
    else                  w =-sqrt(t) ;
    d1 = da + v - w ;
    d2 = db + v + w ;
    if ( (d1 == CGZERO) && (d2 == CGZERO) ) return (-1.) ;
    if ( fabs (d1) >= fabs (d2) ) c = a + delta*da/d1 ;
    else                          c = b - delta*db/d2 ;
    return (c) ;
}

void XXCG(builtin_hprod)
(
    CGFLOAT       *Hd, /* Hd = H*d */
    CGFLOAT        *d,
    CGINT   const   n, /* number of entries in d, H is n by n */
    CGINT   const *Hp, /* column pointers for Hessian */
    CGINT   const *Hi, /* row indices for Hessian */
    CGFLOAT const *Hx  /* nonzero numerical entries in Hessian */
)
{
    CGINT j, p ;
    cg_initx (Hd, CGZERO, n) ;
    p = 0 ;
    for (j = 0; j < n; j++)
    {
        CGINT const q = Hp [j+1] ;
        CGFLOAT const t = d [j] ;
        if ( t ) /* if t != 0 */
        {
            for (; p < q; p++)
            {
                Hd [Hi [p]] += t*Hx [p] ;
            }
        }
        else p = q ;
    }
}

#ifdef PASA
/* =========================================================================
   ==== pasa_cg_maxstep ====================================================
   =========================================================================
   Determine the maximum feasible step in the search direction.
   ========================================================================= */
CGFLOAT XXCG(maxstep)
(
    PASAcom *pasacom,
    CGcom     *cgcom
)
{
    CGINT boundindex, constraintindex, i, j, n, ni, p, q, row,
         *ineq_row, *ATp, *ATi ;
    CGFLOAT maxstep, s, t, *Adk, *Axk, *ATx, *bl, *bu, *d, *lo, *hi, *x ;
    PPprob *Prob ;
    PPwork *Work ;

/*printf ("check maxstep\n") ;*/
    int const loExists = pasacom->loExists ;
    int const hiExists = pasacom->hiExists ;
    int const use_napheap = pasacom->use_napheap ;
    int const use_pproj = pasacom->use_pproj ;
    int const use_penalty = pasacom->use_penalty ;
    int const Aexists = pasacom->Aexists ;
    n = pasacom->nf ;
    lo = pasacom->lo ;
    hi = pasacom->hi ;
    d = cgcom->d ;
    x = cgcom->x ;
    boundindex = constraintindex = 0 ;
    maxstep = PASAINF ;
    if ( pasacom->Bounds )
    {
        PASAFLOAT dmax = PASAZERO ;
        for (j = 0; j < n; j++)
        {
            PASAFLOAT const dj = d [j] ;
            if ( hiExists )
            {
                if ( dj > CGZERO )
                {
                    PASAFLOAT const step = (hi [j]-x [j])/dj ;
                    if ( step <= maxstep )
                    {
                        /* when a bound becomes immediately active, make active
                           the bound corresponding to the largest dj */
                        if ( step <= PASAZERO )
                        {
                            maxstep = PASAZERO ;
                            if ( dj > dmax )
                            {
                                dmax = dj ;
                                boundindex = j + 1 ;
                            }
                        }
                        else
                        {
                            maxstep = step ;
                            boundindex = j + 1 ;
                        }
                    }
                }
                else if ( loExists && (dj < PASAZERO) )
                {
                    PASAFLOAT const step = (lo [j]-x [j])/dj ;
                    if ( step <= maxstep )
                    {
                        /* when a bound becomes immediately active, make active
                           the bound corresponding to the largest dj */
                        if ( step <= PASAZERO )
                        {
                            maxstep = PASAZERO ;
                            if ( -dj > dmax )
                            {
                                dmax = -dj ;
                                boundindex = -(j + 1) ;
                            }
                        }
                        else
                        {
                            maxstep = step ;
                            boundindex = -(j + 1) ;
                        }
                    }
                }
            }
            /* else loExists since bounds exist */
            else if ( dj  < CGZERO )
            {
                PASAFLOAT const step = (lo [j]-x [j])/dj ;
                if ( step <= maxstep )
                {
                    /* when a bound becomes immediately active, make active
                       the bound corresponding to the largest dj */
                    if ( step <= PASAZERO )
                    {
                        maxstep = PASAZERO ;
                        if ( -dj > dmax )
                        {
                            dmax = -dj ;
                            boundindex = -(j + 1) ;
                        }
                    }
                    else
                    {
                        maxstep = step ;
                        boundindex = -(j + 1) ;
                    }
                }
            }
        }
    }
    pasacom->maxbndstep = maxstep ;
    pasacom->maxbndindex = boundindex ;

    /* Next, find the maxstep based on the inequality constraints.
       If the penalty technique is employed, we also need to compute
       the products Ad and A'Ad. If we are not using the penalty
       technique, then we only need the product between the inactive
       rows of A and d. */
    if ( Aexists )
    {
        maxstep = PASAINF ;
        if ( use_pproj )
        {
            bl = pasacom->bl ;
            bu = pasacom->bu ;
            Work = pasacom->ppcom->Work ;
            ATp = Work->ATp ;
            ATi = Work->ATi ;
            ATx = Work->ATx ;
            Axk = pasacom->Axk ;
            Adk = pasacom->Adk ;

            if ( use_penalty == TRUE )
            {
                PPINT l, nrow, *ir ;
                PPFLOAT *AAd ;
                AAd = cgcom->AAd ;
                pasa_initx (AAd, PASAZERO, n) ;
                ir = Work->ir ;
                pasacom->Ad2 = CGZERO ;
                nrow = pasacom->nrow ;
                p = 0 ;
                l = 0 ;
                for (i = 0; i < nrow; i++)
                {
                    s = CGZERO ;
                    q = ATp [i+1] ;
                    for (; p < q; p++)
                    {
                        s += ATx [p]*d [ATi [p]] ;
                    }
                    if ( ir [i] == 0 )
                    {
                        pasacom->Ad2 += s*s ;
                        p = ATp [i] ;
                        for (; p < q; p++)
                        {
                            AAd [ATi [p]] += ATx [p]*s ;
                        }
                    }
                    else
                    {
                        l++ ;
                        Adk [l] = s ;
                        if ( s > CGZERO )
                        {
                            t = (bu [l] - Axk [l])/s ;
                            if ( t < maxstep )
                            {
                                maxstep = t ;
                                constraintindex = l ;
                            }
                        }
                        else if ( s < CGZERO )
                        {
                            t = (bl [l] - Axk [l])/s ;
                            if ( t < maxstep )
                            {
                                maxstep = t ;
                                constraintindex = -l ;
                            }
                        }
                        else
                        {
                            continue ;
                        }
                    }
                }
                pasacom->Ad2 *= pasacom->penalty ;
            }
            else /* focus on the inactive inequalities */
            {
                Prob = pasacom->ppcom->Prob ;
                ni = Prob->ni ;
                if ( ni > 0 )
                {
                    ineq_row = Prob->ineq_row ;
                    for (i = 1; i <= ni; i++)
                    {
                        row = ineq_row [i] ;
                        s = CGZERO ;
                        q = ATp [row+1] ;
                        for (p = ATp [row]; p < q; p++)
                        {
                            s += ATx [p]*d [ATi [p]] ;
                        }
                        Adk [i] = s ;
                        if ( s > CGZERO )
                        {
                            t = (bu [i] - Axk [i])/s ;
                            if ( t < maxstep )
                            {
                                maxstep = t ;
                                constraintindex = i ;
                            }
                        }
                        else if ( s < CGZERO )
                        {
                            t = (bl [i] - Axk [i])/s ;
                            if ( t < maxstep )
                            {
                                maxstep = t ;
                                constraintindex = -i ;
                            }
                        }
                        else
                        {
                            continue ;
                        }
                    }
                }
            }
        }
        else if ( use_napheap == TRUE )
        {
            CGFLOAT *a ;
            maxstep = PASAINF ;
            a = pasacom->nap_a ;
            s = cg_dot (a, d, n) ;
            if ( use_penalty == TRUE ) /* => linear constraint is active */
            {
                pasacom->Ad2 = s*s*pasacom->penalty ;
                cg_scale (cgcom->AAd, a, s, n) ;
            }
            else if ( pasacom->nap_constraint == 0 ) /* inactive linear const */
            {
                pasacom->Adk [1] = s ;
                /* determine the maxstep */
                if ( s > CGZERO )
                {
                    t = (pasacom->nap_bu - pasacom->Axk [1])/s ;
                    if ( t < maxstep )
                    {
                        maxstep = t ;
                        constraintindex = 1 ;
                    }
                }
                else if ( s < CGZERO )
                {
                    t = (pasacom->nap_bl - pasacom->Axk [1])/s ;
                    if ( t < maxstep )
                    {
                        maxstep = t ;
                        constraintindex = -1 ;
                    }
                }
            }
        }
        maxstep = PASAMAX (maxstep, PASAZERO) ;
        pasacom->maxconstep = maxstep ;
        maxstep = PASAMIN (maxstep, pasacom->maxbndstep) ;
    }
    pasacom->maxconindex = constraintindex ;
    cgcom->maxstep = pasacom->maxstep = maxstep ;

    if ( cgcom->Parm->PrintLevel >= 1 )
    {
        if ( pasacom->Bounds )
        {
            printf ("maxbndstep: %e maxbndindex: %ld ",
                     pasacom->maxbndstep, (LONG) pasacom->maxbndindex) ;
        }
        if ( Aexists )
        {
            printf ("maxconstep: %e maxconindex: %ld\n",
                 pasacom->maxconstep, (LONG) pasacom->maxconindex) ;
        }
        else printf ("\n") ;
    }
    return (maxstep) ;
}

#if 0
/* =========================================================================
   ==== cg_checktol ========================================================
   =========================================================================
   This routine is entered when the nominal stopping criterion
   pasacom->e <= switchfactor*pasacom->E is satisfied where pasacom->E is the
   global error bound before entering CG. The goal is to decide what
   is done next. There are four possibilities:

   1. Terminate with status = PASA_ERROR_TOLERANCE_SATISFIED if the global
      error is sufficiently small.
   2. Terminate CG and return to grad_proj where bound variables are
      freed.
   3. If use_penalty is TRUE, then we could update the multiplier and
      penalty term and continue CG but with a restart.
   4. Continue the current iteration but with an updated (smaller)
      pasacom->E value.

   Recall that the optimization problem is

   min f(x) subject to lo <= x <= hi, bl <= Ax <= bu.

   We reformulate this as

   min f(x) subject to lo <= x <= hi, bl <= b <= bu, Ax - b = 0.

   Let L(x, y) = f(x) + y'(Ax - b) be the Lagrangian.  The first-order
   optimality conditions are that x and b should be feasible and that

       L'_j (x, y)  = 0 if lo_j < x_j < hi_j,
       L'_j (x, y) <= 0 if        x_j = hi_j,
       L'_j (x, y) >= 0 if lo_j = x_j,
           y_i      = 0 if bl_i < (Ax)_i < bu_i,
           y_i     <= 0 if bl_i = (Ax)_i
           y_i     >= 0 if        (Ax)_i = bu_i

   In null_project (previously executed) we determined y to minimize
   ||nabla f(x) - B'y||_2 and then stored in Com->e the quantity
   ||nabla f(x) - B'y||_inf (the absolute largest component), where B
   is the submatrix of A corresponding to the active constraints at x.
   Thus Com->e stores the local error in the manifold associated with
   the active constraints at x. To obtain an (over) estimate for the
   global error at x, set the components of y corresponding to the
   inactive inequalities to 0.  Let maxsignerr denote the maximum
   violation of the sign constraints in the first-order
   optimality conditions, and let maxeqerr denote the maximum amount
   that the current iterate violates the polyhedral constraint.
   Our decision about what to do next is based on the following rules:

   A. If MAX (pasacom->e, maxsignerr, maxeqerr) <= pasacom->grad_tol,
      then we terminate with status = PASA_ERROR_TOLERANCE_SATISFIED.
   B. If pasacom->e <= switchfactor*maxsignerr, then we branch to grad_proj.
      The rationale is that we wish to determine the active constraints
      at the solution and reduce maxsignerr.
   C. If pasacom->e < switchfactor*maxeqerr and use_penalty is FALSE,
      then we branch to grad_proj.  The rationale is that there is no
      mechanism in CG for reducing the equation error when use_penalty
      is FALSE. On the other hand, grad_proj will project onto the
      feasible set and remove the error associated with constraint violation.
   D. If pasacom->e < switchfactor*maxeqerr and use_penalty is TRUE,
      then we update the multiplier and penalty term, and restart CG.
      The rationale is that we wish to reduce the violation in the
      constraint by utilizing the multiplier structure of the objective.
      The stopping condition is pasacom->e <=
      switchfactor*MAX (maxsignerr, maxeqerr, old pasacom->e).
   E. If pasacom->e >= switchfactor*maxeqerr, then continue
      CG with the new stopping condition pasacom->e <=
      switchfactor*MAX (maxsignerr, maxeqerr, old pasacom->e).

   When Aexists = FALSE and there are only bound constraints, maxeqerr = 0.
   ========================================================================= */
int XXCG(checktol) /* return:
                               PASA_ERROR_TOLERANCE_SATISFIED
                               PASA_GRAD_PROJ
                               CG_RESTART
                               CG_CONTINUE */
(
    CGFLOAT       *x, /* current iterate */
    PASAcom *pasacom,
    CGcom     *cgcom
)
{
    int Aexists, PrintLevel, relerr ;
    CGINT j, k, nc, nrow, *ATi, *ATp, *bound_cols, *RLinkUp ;
    CGFLOAT E, Lambda, maxeqerr, maxsignerr, penalty, residual, switchfactor,
            testtol_old, *ATx, *userg, *lambda, *temp ;
    PASAparm *Parm ;      /* parameters for PASA */
    PPwork *Work ;
    CGstat *Stat ;

    Parm = pasacom->pasaparm ; /* parameters */
    const int use_penalty = pasacom->use_penalty ;
    const int use_pproj = pasacom->use_pproj ;
    const int use_napheap = pasacom->use_napheap ;
    penalty = Parm->penalty ;
    switchfactor = pasacom->switchfactor ;
    PrintLevel = Parm->PrintLevel ;
    Aexists = pasacom->Aexists ;
    Stat = cgcom->Stat ;

    if ( pasacom->QP )
    {
        userg = pasacom->gtot ;
    }
    else
    {
        userg = pasacom->userg ;
    }

    /* we use the max violation of the first order optimality conditions
       relative to multiplier signs to determine whether to decrease testtol */
    maxsignerr = CGZERO ;
    maxeqerr = CGZERO ;
    nc = pasacom->nc ;
    bound_cols = pasacom->bound_cols ;
    nrow = pasacom->nrow ;
    if ( use_pproj )
    {
        Work = pasacom->ppcom->Work ;
        RLinkUp = Work->RLinkUp ;
        /* the problem is effectively bound constrained if there are no
           active equations */
        if ( RLinkUp [nrow] == nrow )
        {
            Aexists = FALSE ;
        }
    }
    else if ( use_napheap == TRUE )
    {
        if ( (pasacom->nap_constraint == 0) || (pasacom->nap_a2 == PASAZERO) )
        {
            Aexists = FALSE ;
        }
    }
    if ( Aexists == FALSE ) /* only bound constraints present */
    {
        for (k = 0; k < nc; k++)
        {
            j = bound_cols [k] ;
            if ( j < 0 ) /* at lower bound */
            {
                j = -(j+1) ;
                maxsignerr = CGMAX (maxsignerr, -userg [j]) ;
            }
            else
            {
                maxsignerr = CGMAX (maxsignerr, userg [j]) ;
            }
        }
    }
#ifndef NOPPROJ
    else if ( use_pproj ) /* linear equations/inequalities present*/
    {
        CGINT i, l, m, n, nr, p, q,
             *Ap, *Anz, *Ai, *bound_rows, *invperm, *ir, *RLinkDn ;
        CGFLOAT r, s, t, absAx, *Ax, *Axk, *b ;
        PPprob *Prob ;

        Prob = pasacom->ppcom->Prob ;

        /* If use_penalty is FALSE, then we can use lambda stored in
           pasacom->lambda for the multiplier. This lambda was computed
           during the last call to pasa_null_project.  If use_penalty is TRUE
           and we have not started the first iteration, then we can use
           the lambda computed at the top of CG (lambda_pen). Otherwise,
           we must update the lambda computed in null_project to take into
           account the fact that the complete multiplier is the sum
           of the multiplier computed in null_project, the multiplier
           lambda_pen obtained at the start of CG, and the penalty term
           p*(b-Bx). Instead of making these adjustments, we compute
           the multiplier from scatch here. Note though that time would
           be saved if the multiplier estimate was obtained by the update
           process. */
        lambda = pasacom->lambda ;
        if ( use_penalty == TRUE )
        {
            if ( Stat->iter == cgcom->FirstIter )
            {
                lambda = pasacom->lambda_pen ;
            }
            else /* compute lambda from scratch */
            {
                Ap = Prob->Ap ;
                Anz = Prob->Anz ;
                Ai = Prob->Ai ;
                Ax = Prob->Ax ;
                pasa_initx (lambda, CGZERO, nrow) ;

                /* use code from start of CG to estimate multiplier
                   NOTE: we could estimate the multiplier using the first-order
                         multiplier method update */
                n = pasacom->nf ;
                for (j = 0; j < n; j++)
                {
                    t = pasacom->g [j] ;
                    if ( t != PASAZERO )
                    {
                        k = Ap [j] ;
                        l = k + Anz [j] ;
                        for (; k < l; k++)
                        {
                            lambda [Ai [k]] += t*Ax [k] ;
                        }
                    }
                }
                RLinkDn = Work->RLinkDn ;
                pproj_lsol (Work->L, lambda, RLinkUp [nrow], nrow, RLinkUp) ;
                k = RLinkUp [nrow] ;
                /* momentarily set the initial RLinkDn to -1, this simplifies
                   indexing in dltsolve */
                RLinkDn [k] = -1 ;
                l = RLinkDn [nrow] ;
                pproj_dltsol (Work->L, lambda, lambda, l, k, RLinkDn) ;
                RLinkDn [k] = nrow ; /* restore RLinkDn */
            }
        }
        /* lambda is the negative of the multiplier described in the comments */

        Ap = pasacom->Copy->Ap ;
        Ax = pasacom->Copy->Ax ;
        Ai = pasacom->Copy->Ai ;
        invperm = pasacom->invperm ;
        ir  = Work->ir ;
        for (k = 0; k < nc; k++)
        {
            j = bound_cols [k] ;
            if ( j < 0 ) /* at lower bound */
            {
                m = -(j+1) ;      /* user indexing */
                l = invperm [m] ; /* convert from user index to pproj index */
            }
            else
            {
                m = j ;
                l = invperm [j] ;
            }
            t = CGZERO ;
            q = Ap [l+1] ;
            for (p = Ap [l]; p < q; p++)
            {
                t -= lambda [Ai [p]]*Ax [p] ;
            }
            t += userg [m] ;
            if ( j < 0 ) /* variable at lower bound */
            {
                maxsignerr = CGMAX (maxsignerr, -t) ;
                if ( (PrintLevel >= 2) && (t < CGZERO) )
                {
                    printf ("lower bound %ld wrong sign: %e\n", (LONG) m, -t) ;
                }
            }
            else
            {
                maxsignerr = CGMAX (maxsignerr, t) ;
                if ( (PrintLevel >= 2) && (t > CGZERO) )
                {
                    printf ("upper bound %ld wrong sign: %e\n", (LONG) m, t) ;
                }
            }
        }
        nr = pasacom->nr ;
        bound_rows = pasacom->bound_rows ;
        for (k = 0; k < nr; k++)
        {
            j = bound_rows [k] ;
            if ( j < 0 ) /* at lower bound */
            {
                m = -(j+1) ;      /* row index */
                maxsignerr = CGMAX (maxsignerr, -lambda [m]) ;
                if ( (PrintLevel >= 2) && (lambda [m] < CGZERO) )
                {
                    printf ("lower ineqn %ld wrong sign: %e\n",
                           (LONG) m, -lambda [m]) ;
                }
            }
            else
            {
                m = j - 1 ;
                maxsignerr = CGMAX (maxsignerr, lambda [m]) ;
                if ( (PrintLevel >= 2) && (lambda [m] > CGZERO) )
                {
                    printf ("upper ineqn %ld wrong sign: %e\n",
                           (LONG) m, lambda [m]) ;
                }
            }
        }


        relerr = TRUE ;
        if ( pasacom->pprojparm->stop_condition == 1 ) /* use absolute error */
        {
            relerr = FALSE ;
        }

        /* maxeqerr is computed at top CG, otherwise compute it here */
        if ( Stat->iter > cgcom->FirstIter )
        {
            absAx = CGZERO ;   /* stores max_i sum_j |A_{ij}x_j| */
            ATp = Work->ATp ;
            ATi = Work->ATi ;
            ATx = Work->ATx ;
            temp = Work->arrayd ;
            pasa_initx (temp, PASAZERO, nrow) ;
            Axk = pasacom->Axk ;
            b = pasacom->b ;
            p = 0 ;
            l = 0 ;
            for (i = 0; i < nrow; i++)
            {
                r = s = CGZERO ;
                q = ATp [i+1] ;
                for (; p < q; p++)
                {
                    t = ATx [p]*x [ATi [p]] ;
                    s += t ;
                    r += fabs (t) ;
                }
                if ( r > absAx )
                {
                    absAx = r ;
                }
                if ( ir [i] == 0 )
                {
                    temp [i] = t = b [i] - s ;
                    if ( fabs (t) > maxeqerr )
                    {
                        maxeqerr = fabs (t) ;
                    }
                }
                else
                {
                    l++ ;
                    Axk [l] = s ;
                }
            }
            ASSERT (l == pasacom->ppcom->Prob->ni) ;
        }
        else /* extract absAx from pproj */
        {
            absAx = Work->absAx ;
            maxeqerr = cgcom->maxeqerr ;
        }
        if ( relerr == TRUE )
        {
            if ( absAx > CGZERO )
            {
                maxeqerr /= absAx ;
            }
        }
    }
#endif
    else if ( use_napheap == TRUE )
    {
        /* If use_penalty is FALSE, then we can use lambda stored in
           pasacom->lambda for the multiplier. This lambda was computed
           during the last call to pasa_null_project.  If use_penalty is TRUE
           and we have not started the first iteration, then we can use
           the lambda computed at the top of CG (lambda_pen). Otherwise,
           we must update the lambda computed in null_project to take into
           account the fact that the complete multiplier is the sum
           of the multiplier computed in null_project, the multiplier
           lambda_pen obtained at the start of CG, and the penalty term
           p*(b-Bx). Instead of making these adjustments, we compute
           the multiplier from scatch here. Note though that time would
           be saved if the multiplier estimate was obtained by the update
           process. */
        PASAINT m ;
        PASAFLOAT absAx, s, t, *a ;
        a = pasacom->nap_a ;
        Lambda = *(pasacom->lambda) ;
        if ( use_penalty == TRUE )
        {
            if ( Stat->iter == cgcom->FirstIter )
            {
                Lambda = *(pasacom->lambda_pen) ;
            }
            else /* compute lambda from scratch */
            {
                /* multiplier estimate = a'g/a'a */
                Lambda = pasa_dot (a, pasacom->g, pasacom->nf)/pasacom->nap_a2 ;
            }
        }
        for (k = 0; k < nc; k++)
        {
            j = bound_cols [k] ;
            if ( j < 0 ) /* at lower bound */
            {
                m = -(j+1) ;      /* user indexing */
            }
            else
            {
                m = j ;
            }
            t = userg [m] - Lambda*pasacom->nap_auser [m] ;
            if ( j < 0 ) /* variable at lower bound */
            {
                maxsignerr = CGMAX (maxsignerr, -t) ;
                if ( (PrintLevel >= 2) && (t < CGZERO) )
                {
                    printf("lower bound %ld (user index) wrong sign: %e\n",
                          (LONG) m, -t) ;
                }

            }
            else
            {
                maxsignerr = CGMAX (maxsignerr, t) ;
                if ( (PrintLevel >= 2) && (t > PASAZERO) )
                {
                    printf("upper bound %ld (user index) wrong sign: %e\n",
                          (LONG) m, t) ;
                }

            }
        }
        if ( pasacom->nr > 0 )
        {
            if ( pasacom->bound_rows [0] < 0 ) /* at lower bound */
            {
                maxsignerr = PASAMAX (maxsignerr, -Lambda) ;
                if ( (PrintLevel >= 2) && (Lambda < PASAZERO) )
                {
                    printf ("lower ineqn %i wrong sign: %e\n", 0, -Lambda) ;
                }
            }
            else                           /* at upper bound */
            {
                maxsignerr = PASAMAX (maxsignerr, Lambda) ;
                if ( (PrintLevel >= 2) && (Lambda > PASAZERO) )
                {
                    printf ("upper ineqn %i wrong sign: %e\n", 0, Lambda) ;
                }
            }
        }

        relerr = TRUE ;
        absAx = CGZERO ;
        s = CGZERO ;
        for (j = 0; j < pasacom->nf; j++)
        {
            t = a [j]*x [j] ;
            s += t ;
            absAx += fabs (t) ;
        }
        residual = pasacom->nap_bl - s ;
        maxeqerr = fabs (residual) ;
        if ( relerr == TRUE )
        {
            if ( absAx > CGZERO )
            {
                maxeqerr /= absAx ;
            }
        }
    }

    /* TEST A: If MAX (pasacom->e, maxsignerr, maxeqerr) <= pasacom->grad_tol,
       then we terminate with status = PASA_ERROR_TOLERANCE_SATISFIED. */
    if ( PrintLevel >= 1 )
    {
        printf ("maxsignerr: %e maxeqerr: %e\n", maxsignerr, maxeqerr) ;
    }
    E = CGMAX (maxeqerr, maxsignerr) ;
    E = CGMAX (E, pasacom->e) ;
    testtol_old = pasacom->testtol ;
    if ( E < pasacom->E )
    {
        pasacom->E = E ;
        pasacom->testtol = E*switchfactor ;
        pasacom->testtol = CGMAX (pasacom->testtol, pasacom->grad_tol) ;
    }
    if ( PrintLevel >= 1 )
    {
        printf ("cg error bound E: %e Global error: %e testtol: %e\n",
                 E, pasacom->E, pasacom->testtol) ;
    }
    if ( pasacom->E <= pasacom->grad_tol )
    {
        if ( PrintLevel >= 1 )
        {
            printf ("checktol: A. Error tolerance satisfied\n") ;
        }
        return (PASA_ERROR_TOLERANCE_SATISFIED) ;
    }

    /* TEST B: If pasacom->e <= switchfactor*maxsignerr, then we branch
       to grad_proj.  The rationale is that we wish to determine the active
       constraints at the solution and reduce maxsignerr. */
    if ( pasacom->e <= pasacom->testtol )
    {
        if ( PrintLevel >= 1 )
        {
            printf ("checktol: B. Branch to grad_proj\n") ;
        }
        return (PASA_GRAD_PROJ) ;
    }

    /* TEST C: If pasacom->e < switchfactor*maxeqerr and use_penalty is FALSE,
       then we branch to grad_proj.  The rationale is that there is no
       mechanism in CG for reducing the equation error when use_penalty
       is FALSE. On the other hand, grad_proj will project onto the
       feasible set and remove the error associated with constraint violation.*/
    if ( (pasacom->e < switchfactor*maxeqerr) && (use_penalty == FALSE) )
    {
        if ( PrintLevel >= 1 )
        {
            printf ("checktol: C. Branch to grad_proj\n") ;
        }
        return (PASA_GRAD_PROJ) ;
    }

    /* TEST D: If pasacom->e <= switchfactor*maxeqerr and use_penalty is TRUE,
       then we update the multiplier and penalty term, and restart CG.
       The rationale is that we wish to reduce the violation in the
       constraint by utilizing the multiplier structure of the objective.
       The stopping condition is pasacom->e <=
       switchfactor*MAX (maxsignerr, maxeqerr, pasacom->e). */

    /* only perform test D in the middle of CG, not at the start */
    if ( Stat->iter > cgcom->FirstIter )
    {
        if ( (pasacom->e < switchfactor*maxeqerr) && (use_penalty == TRUE) )
        {
            CGINT i, p, q ;
            CGFLOAT c, fp, s, t, *gpen ;
            gpen = cgcom->gpen ;
            c = 0.5*penalty ;
            if ( use_pproj )
            {
                fp = CGZERO ;
                pasa_initx (gpen, PASAZERO, pasacom->nf) ;
                i = nrow ;
                while ( (i = RLinkUp [i]) < nrow )
                {
                    t = temp [i] ;
                    /* fp = penalty term in the objective */
                    fp += t*(c*t + lambda [i]) ;
                    s = -(penalty*t + lambda [i]) ;
                    q = ATp [i+1] ;
                    for (p = ATp [i]; p < q; p++)
                    {
                        gpen [ATi [p]] += ATx [p]*s ;/* gradient of penalty */
                    }
                }
            }
            else if ( use_napheap == TRUE )
            {
                fp = residual*(c*residual + Lambda) ;
                s = -(penalty*residual + Lambda) ;
                pasa_scale (gpen, pasacom->nap_a, s, pasacom->nf) ;
            }
            pasacom->fp = fp ;

            /* since the penalty term has changed, recompute objective value */
            pasacom->f = pasacom->f_orig + fp ;

            /* if the testtol did not change, then make it smaller */
            if ( testtol_old == pasacom->testtol )
            {
                pasacom->testtol *= switchfactor ;
                /* do not let testtol drop below grad_tol */
                pasacom->testtol = PASAMAX(pasacom->testtol, pasacom->grad_tol);
            }

            if ( PrintLevel >= 1 )
            {
                printf ("checktol: D. Restart CG with tolerance %e\n",
                         pasacom->testtol) ;
            }
            return (CG_RESTART) ;
        }
    }

    /* TEST E: If pasacom->e > switchfactor*maxeqerr, then continue CG.
       The stopping condition is pasacom->e <=
       switchfactor*MAX (maxsignerr, maxeqerr, pasacom->e). */
    if ( PrintLevel >= 1 )
    {
        printf ("checktol: E. Continue CG with new tolerance %e "
                "(old tolerance: %e\n", pasacom->testtol, testtol_old) ;
    }
    return (CG_CONTINUE) ;
}
#endif
#endif
#ifndef PASA
/* ==========================================================================
   === cg_basis =============================================================
   ==========================================================================
   Compute orthonormal basis for columns of A, store the basis in Z
   The algorithm is based on a QR factorization of A, with column permuations.
   Columns are treated as linearly dependent if the distance from the
   normalized column to the space spanned by the prior vectors is <= QRcutoff.
   ========================================================================== */
CGINT cg_basis
(
    CGFLOAT       *Z, /* dense matrix of orthonormal basis vectors */
    CGFLOAT       *T, /* upper triangular matrix stored by columns */
    CGFLOAT       *A, /* A*P*T = Z dense matrix, size nrow by rank */
    CGINT      *perm, /* column # of A for column i of Z is perm [i] */
    int   PrintLevel,
    CGINT       nrow,
    CGINT       ncol,
    CGFLOAT QRcutoff, /* QRcutoff for discarding dependent columns */
    CGFLOAT    *Work  /* work space for Householder vectors, copy of A,
                           column norms of A, column norms of partially
                           factorized A */
)
{
    CGINT i, j, k, m ;
    CGFLOAT s, t, *Aj, *Ak, *Hj, *Hk, *Tk, *Zk ;
    CGFLOAT *work = Work ;
    CGFLOAT mindiag = CGINF ;
    /* when the max norm of the column squared exceeds colmax_recompute,
       recompute the square norms of the columns */
    CGINT const minrowcol = CGMIN (nrow, ncol) ;
    /* normalize columns of A to be unit vectors */
    CGFLOAT        *H = work ;  work += (minrowcol*(2*nrow - minrowcol+1))/2 ;
    CGFLOAT    *Anorm = work ;  work += ncol ;
    CGFLOAT  *colnorm = work ;  work += ncol ;

    Hj = H ;
    Aj = A ;
    for (j = 0; j < ncol; j++)
    {
        perm [j] = j ;
        t = cg_dot (Aj, Aj, nrow) ;
        if ( t != CGZERO )
        {
            CGFLOAT const u = CGONE/sqrt (t) ;
            Anorm [j] = u ;
            colnorm [j] = CGONE ;
            cg_scale (Aj, Aj, u, nrow) ;
        }
        else
        {
            Anorm [j] = CGZERO ;
            colnorm [j] = CGZERO ;
        }
        Aj += nrow ;
    }
    CGFLOAT colmax = CGONE ;
    CGFLOAT colmax_recompute = 1.e-5 ;
    k = 0 ;
    /* replac A by A * diag (Anorm) */
    for (j = 0; j < minrowcol; j++)
    {
        /* recompute column norms squared if colmax <= colmax_recompute */
        if ( (colmax <= colmax_recompute) && (j > 0) )
        {
            colmax = -1 ; /* reevaluate colmax */
            Aj = A+(j*nrow) ;
            for (i = j; i < minrowcol; i++)
            {
                t = cg_dot (Aj+j, Aj+j, nrow-j) ;
                colnorm [i] = t ;
                if ( t > colmax )
                {
                    colmax = t ;
                    k = i ;
                }
                Aj += nrow ;
            }
            colmax_recompute = 1.e-5 * colmax ;
        }
        if ( sqrt(colmax) <= QRcutoff ) /* done */
        {
            if ( PrintLevel )
            {
                printf ("break at col %ld due to small colmax: %e cutoff: %e\n",
                         (LONG) j, sqrt (colmax), QRcutoff) ;
            }
            for (i = j; i < ncol; i++) perm [i] = EMPTY ;
            break ;
        }
        /* swap columns k and j where k has the max norm if colnorm [j] <
           0.01*colmax -- if possible, try to use the provided order */
        Aj = A+nrow*j ; /* start of column j */
        if ( (k != j) && (colnorm [j] < .01*colmax) )
        {
            colnorm [k] = colnorm [j] ;
            colnorm [j] = colmax ;
            i = perm [k] ;
            perm [k] = perm [j] ;
            perm [j] = i ;
            t = Anorm [k] ;
            Anorm [k] = Anorm [j] ;
            Anorm [j] = t ;
            Ak = A+nrow*k ; /* start of column k */
            /* swap columns j and k */
            t = CGZERO ;
            for (i = 0; i < j; i++)
            {
                s = Ak [i] ;
                Ak [i] = Aj [i] ;
                Aj [i] = s ;
            }
            for (; i < nrow; i++)
            {
                s = Ak [i] ;
                Ak [i] = Aj [i] ;
                Aj [i] = s ;
                t += s*s ;
            }
        }
        else /* k = j */
        {
            t = cg_dot (Aj+j, Aj+j, nrow-j) ;
        }
        /* compute Householder matrix for column j */
        if ( t == CGZERO )
        {
            if ( PrintLevel )
            {
                printf ("remainder of columns are zero: %ld\n", (LONG) j) ;
            }
            for (i = j; i < ncol; i++) perm [i] = EMPTY ;
            break ;
        }
        s = sqrt(t) ;
        if ( s < mindiag ) mindiag = s ;
        t = 1/sqrt(2*s*(s+fabs(Aj [j]))) ;
        if ( Aj [j] < CGZERO )
        {
            Hj [0] = (Aj [j] - s)*t ;
            Aj [j] = s ;
        }
        else
        {
            Hj [0] = (Aj [j] + s)*t ;
            Aj [j] = -s ;
        }
        cg_scale (Hj+1, Aj+j+1, t, nrow-j-1) ;
        /* apply Household to A and update column norms */
        colmax = -1 ;
        for (i = j+1; i < minrowcol; i++)
        {
            Aj += nrow ;
            t = 2*cg_dot (Hj, Aj+j, nrow-j) ;
            cg_daxpy (Aj+j, Hj, -t, nrow-j) ;
            colnorm [i] -= Aj [j]*Aj [j] ;
            if ( colnorm [i] > colmax )
            {
                colmax = colnorm [i] ;
                k = i ;
            }
        }
        Hj += nrow-j ;
    }
    if ( PrintLevel )
    {
        printf ("min diagonal of R: %e\n", mindiag) ;
    }
    /* determine number of linearly independent vectors */
    m = j ; /* returned at end of code */
    /* Let P denote the permutation matrix corresponding to the column
       swapping during the QR factorization. We have computed factorization

               A*D*P = H1*H2 ... Hm * [R above 0]

       or equivalently,

               A*P(P*D*P)inv(R) = first m columns of H1*H2 ... Hm

       The first m columns of the Householder product correspond to the basis.
       The product (P*D*P)*inv(R) is returned in T, while the permutation is
       is in output perm.  */

    Tk = T;
    for (k = 0; k < m; k++)
    {
        /* compute column k of P*D*P*inv (R) and store in T */
        Ak = A+k*nrow ;
        t = CGONE/Ak [k] ;
        Tk [k] = t*Anorm [k] ; /* multiply by Dk */
        for (j = k-1; j >= 0; j--) Tk [j] = -t*Ak [j] ;
        Ak -= nrow ;
        for (j = k-1; j >= 0; j--)
        {
            t = Tk [j]/Ak [j] ;
            Tk [j] = t*Anorm [j] ; /* multiply by Dj */
            for (i = j-1; i >= 0; i--) Tk [i] -= t*Ak [i] ;
            Ak -= nrow ;
        }
        Tk += k+1 ;
    }

    /* evaluate the m columns of Z, the basis vectors */
    Hk = H ;
    Zk = Z ;
    for (k = 0; k < m; k++)
    {
        Hj = Hk ;
        t = 2*Hj [0] ;
        cg_initx (Zk, CGZERO, k) ;
        cg_scale (Zk+k, Hj, -t, nrow-k) ;
        Zk [k] += CGONE ;
        for (j = k-1; j >= 0; j--)
        {
            Hj -= (nrow - j) ;
            t = 2*cg_dot (Zk+j, Hj, nrow-j) ;
            cg_daxpy (Zk+j, Hj, -t, nrow-j) ;
        }
        Zk += nrow ;
        Hk += nrow-k ;
    }
    return (m) ;
}
#endif

/* =========================================================================
   ==== pasa_cg_PHPz =======================================================
   =========================================================================
   Compute the product P*H*P*z
   ========================================================================= */
CGFLOAT * XXCG(PHPz)
(
    CGFLOAT       *z,
    CGFLOAT      *p1,
    CGFLOAT      *p2,
    int      Aexists,
    int      triples,
    SOPT_Cmatrix  *C
#ifdef PASA
    , PASAcom *pasacom
#endif
)
{
    CGFLOAT *pz ;
    pz = z ;
    /* operate using P */
    if ( Aexists )
    {
#ifdef PASA
        pasa_null_project (p1, z, NULL, FALSE, pasacom);
#endif
        pz = p1 ;
    }
    /* multiply pz by H and store in p2 */
    if ( triples )
    {
        cg_initx (p2, CGZERO, C->n) ;
        CGINT   *rows = C->rows ;
        CGINT   *cols = C->cols ;
        CGFLOAT *vals = C->vals ;
        CGFLOAT *tp2 = p2-1 ;
        CGFLOAT *tz = z-1 ;
        CGINT i, row, col ;
        for (i = 0; i < C->nnz; i++)
        {
            row = rows [i] ;
            col = cols [i] ;
            tp2 [row] += tz [col]*vals [i] ;
            if ( row != col )
            {
                tp2 [col] += tz [row]*vals [i] ;
            }
        }
    }
    else cg_builtin_hprod (p2, pz, C->n, C->p, C->i, C->x) ;
    pz = p2 ;
#ifdef PASA
    if ( Aexists )
    {
        /* compute final projection in P*H*P*z */
        pasa_null_project (p1, p2, NULL, FALSE, pasacom) ;
        pz = p1 ;
    }
#endif
    return (pz) ;
}

int XXCG(debug_C_in_CG)
(
    SOPT_Cmatrix *C,
    SOPT_matrix  *H
)
{
    /* check that each element of the Hessian coincides with an
       element of the compressed Hessian. First copy C to CDebug */
    SOPTINT i, j, k, l, p, p1 ;
    SOPTFLOAT t ;
    int status = 0 ;
    int const fortran = H->fortran ;
    CGINT n = C->n ;
    CGINT Cnnz = C->p [n] ;
    int const compress = (C->col == NULL) ? FALSE : TRUE ;
    CGFLOAT *CxDebug = sopt_malloc (&status, Cnnz, sizeof (SOPTFLOAT)) ;
    sopt_copyx (CxDebug, C->x, Cnnz) ;
    if ( status ) return (sopt_convert_error (LCG, status)) ;
#ifdef PASA
    CGINT *Ccol = C->col ;
#endif
    if ( C->Hformat == 3 )
    {
        for (k = 0; k < H->nnz; k++)
        {
            if ( fortran )
            {
                i = H->rows [k] - 1 ;
                j = H->cols [k] - 1 ;
            }
            else
            {
                i = H->rows [k] ;
                j = H->cols [k] ;
            }
            if ( compress )
            {
#ifdef PASA
                if ( (Ccol [i] == EMPTY) || (Ccol [j] == EMPTY)) continue;
                i = Ccol [i] ;
                j = Ccol [j] ;
#endif
            }
            t = H->vals [k] ;
            CGINT const q1 = C->p [j+1] ;
            int found = FALSE ;
            for (p = C->p [j]; p < q1; p++)
            {
                if ( C->i [p] == i )
                {
                    found = TRUE ;
                    if ( CxDebug [p] != t )
                    {
                        printf ("row: %li col: %li error %e in val: %e\n",
                                (LONG) i, (LONG) j, C->x[p]-t, t) ;
                        sopt_error (-1, __FILE__, __LINE__, "stop") ;
                    }
                    else
                    {
                        CxDebug [p] = CGZERO ;
                    }
                }
            }
            if ( found == FALSE )
            {
                printf ("row: %li col: %li missing in C\n", (LONG) i, (LONG) j);
                sopt_error (-1, __FILE__, __LINE__, "stop") ;
            }
            if ( i == j ) continue ;
            CGINT const q2 = C->p [i+1] ;
            found = FALSE ;
            for (p = C->p [i]; p < q2; p++)
            {
                if ( C->i [p] == j )
                {
                    found = TRUE ;
                    if ( CxDebug [p] != t )
                    {
                        printf ("Row: %li Col: %li error %e in val: %e\n",
                                (LONG) i, (LONG) j, C->x[p]-t, t);
                        sopt_error (-1, __FILE__, __LINE__, "stop");
                    }
                    else
                    {
                        CxDebug [p] = CGZERO ;
                    }
                }
            }
            if ( found == FALSE )
            {
                printf ("row: %li col: %li missing in C\n", (LONG) i, (LONG) j);
                sopt_error (-1, __FILE__, __LINE__, "stop") ;
            }
        }
        for (j = 0; j < n; j++)
        {
            CGINT q = C->p [j+1] ;
            for (p = C->p [j]; p < q; p++)
            {
                if ( CxDebug [p] != CGZERO )
                {
                    printf ("row: %li col: %li val: %e not in triples\n",
                             (LONG) C->i [p], (LONG) j, C->x [p]) ;
                    sopt_error (-1, __FILE__, __LINE__, "stop") ;
                }
            }
        }
    }
    else if ( C->Hformat == 0 )
    {
        for (j = 0; j < H->ncol; j++)
        {
            l = j ;
#ifdef PASA
            if ( (l = Ccol [j]) == EMPTY ) continue ;
#endif
            CGINT const q = H->p [j+1] ;
            k = 0 ;
            for (p = H->p [j]; p < q; p++)
            {
                i = H->i [p] ;
#ifdef PASA
                if ( (i = Ccol [i]) == EMPTY ) continue ;
#endif
                /* look for i in column l of C */
                CGINT const q1 = C->p [l+1] ;
                for (p1 = C->p [l]; p1 < q1; p1++)
                {
                    if ( i == C->i [p1] ) break ;
                }
                if ( p1 < q1 ) /* found it */
                {
                    if ( C->x [p1] != H->x [p] ) /* values do not match */
                    {
                        printf ("in row %li and column %li of C, the value "
                                "in Cx = %e does not match value %e in Hx\n",
                                (LONG) i, (LONG) l, C->x [p1], H->x [p]) ;
                    }
                    else
                    {
                        if ( CxDebug [p1] < CGINF ) CxDebug [p1] = CGINF ;
                        else
                        {
                            printf ("row index %li repeats in column %li "
                                    "of C\n", (LONG) i, (LONG) l) ;
                            sopt_error (-1, __FILE__, __LINE__, "stop") ;
                        }
                    }
                }
                else /* error since index i not found in column l of C */
                {
                    printf ("Row: %li Col: %li missing from C in H->C\n",
                            (LONG) i, (LONG) l) ;
                    sopt_error (-1, __FILE__, __LINE__, "stop");
                }
                k++ ;
            }
            if ( C->p [l+1] - C->p [l] != k )
            {
                printf ("C->p [%li] - C->p [%li] != %li in _CG\n",
                        (LONG) l+1, (LONG) l, (LONG) k) ;
                sopt_error (-1, __FILE__, __LINE__, "stop");
            }
        }
        /* check for sorted cols */
        p = 0 ;
        for (j = 0; j < n; j++)
        {
            CGINT const q = C->p [j+1] ;

            if ( p == q ) continue ;
            for (p++; p < q; p++)
            {
                if ( C->i [p-1] >= C->i [p] )
                {
                    printf ("column %li in C not sorted %li >= %li at %li\n",
                       (LONG) j, (LONG) C->i [p-1], (LONG) C->i [p], (LONG) p) ;
                    sopt_error (-1, __FILE__, __LINE__, "stop");
                }
            }
        }
        /* check that CxDebug is completely infinity */
        for (i = 0; i < Cnnz; i++)
        {
            if ( CxDebug [i] < CGINF )
            {
                printf ("element %li in C->x not in H\n", (LONG) i) ;
                sopt_error (-1, __FILE__, __LINE__, "stop");
            }
        }
    }
    sopt_free (CxDebug) ;
    return (status) ;
}

int XXCG(debug_C_in_SS)
(
    SOPT_Cmatrix *C,
    SOPT_matrix  *H
)
{
    CGINT i, j, k, l ;
    /* check the compressed Hessian */
    int status = 0 ;
    int const fortran = H->fortran ;
#ifdef PASA
    CGINT I, J ;
    int const compress = (C->col == NULL) ? FALSE : TRUE ;
    CGINT *Ccol = C->col ;
#endif
    if ( C->Hformat == 3 )
    {
        l = 0 ;
        for (k = 0; k < H->nnz; k++)
        {
            i = H->rows [k] ;
            j = H->cols [k] ;
#ifdef PASA
            if ( compress )
            {
                I = i ;
                J = j ;
                if ( fortran ) /* if H based on Fortran, lower indices */
                {
                    I-- ;
                    J-- ;
                }
                if ( ((I = Ccol [I]) < 0) || ((J = Ccol [J]) < 0) ) continue;
                if ( fortran )
                {
                    i = I + 1 ;
                    j = J + 1 ;
                }
            }
#endif
            if ( i > j ) continue ;
            if ( !fortran )
            {
                i++ ;
                j++ ;
            }
            if ( C->rows [l] != i )
            {
                printf ("%li != %li comp. row wrong at k: %li l: %li\n",
                         (LONG) C->rows [l], (LONG) i, (LONG) k, (LONG) l) ;
                sopt_error (-1, __FILE__, __LINE__, "stop") ;
            }
            if ( C->cols [l] != j )
            {
                printf ("%li != %li comp. col wrong at k: %li l: %li\n",
                         (LONG) C->cols [l], (LONG) j, (LONG) k, (LONG) l) ;
                sopt_error (-1, __FILE__, __LINE__, "stop") ;
            }
            if ( H->vals [k] != C->vals [l] )
            {
                printf ("%e != %e comp. val wrong at k: %li l: %li\n",
                         C->vals [l], H->vals [k], (LONG) k, (LONG) l) ;
                sopt_error (-1, __FILE__, __LINE__, "stop") ;
            }
            l++ ;
        }
    }
    return (status) ;
}

#if defined (USE_MUMPS) || defined (USE_HARWELL)
int XXCG(checkSS)
(
#ifdef USE_MUMPS
    DMUMPS_STRUC_C *id,
#elif USE_HARWELL
    LBL_Data      *lbl,
#endif
    CGFLOAT      *vals,
    CGFLOAT       *rhs,
    CGFLOAT  *rhsDebug,
    CGINT const    nnz
)
{
    CGINT i, j, k, n ;
    CGFLOAT s, t ;
    int status = 0 ;
#ifdef USE_MUMPS
    n = id->n ;
    CGINT   *rows = id->irn ;
    CGINT   *cols = id->jcn ;
#elif USE_HARWELL
    n = lbl->ma57->n ;
    CGINT *rows = lbl->irn ;
    CGINT *cols = lbl->jcn ;
#endif
/*sopt_printxMATLAB (rhsDebug, n, "rhs") ;
printf ("A = [\n") ;*/
    CGFLOAT *absAx = sopt_malloc (&status, n, sizeof (CGFLOAT)) ;
    if ( status ) return (CG_OUT_OF_MEMORY) ;
    sopt_initx (absAx, CGZERO, n) ;
    for (k = 0; k < nnz; k++)
    {
        i = rows [k] - 1 ;
        j = cols [k] - 1 ;
/*printf ("%i %i %25.15e\n", i+1, j+1, vals [k]) ;*/
        t = rhs [j]*vals [k] ;
        absAx [i] += fabs (t) ;
        rhsDebug [i] -= t ;
        if ( i != j  )
        {
/*printf ("%i %i %25.15e\n", j+1, i+1, vals [k]) ;*/
            t = rhs [i]*vals [k] ;
            absAx [j] += fabs (t) ;
            rhsDebug [j] -= t ;
        }
    }
/*printf ("] ;\n") ;*/
    s = sopt_sup_normx (absAx, n) ;
    if ( s != CGZERO )
    {
        t = sopt_sup_normx (rhsDebug, n)/s ;
    }
    else
    {
        t = sopt_sup_normx (rhsDebug, n) ;
    }
    printf ("Relative Error in symmetric solve: %e max absAx: %e\n\n", t, s) ;
    printf ("nnz in Hessian: %li\n", (LONG) nnz) ;
    sopt_free (rhsDebug) ;
    sopt_free (absAx) ;
    return (status) ;
}
#endif

#ifdef USE_HARWELL
void XXCG(harwell_fact_data)
(
    LBL_Data *lbl
)
{
    printf ("Number of entries in factors:           %li\n",
            (LONG) lbl->ma57->info [13]) ;
    printf ("Storage for real data in factors:       %li\n",
            (LONG) lbl->ma57->info [14]) ;
    printf ("Min length of fact:                     %li\n",
            (LONG) lbl->ma57->info [16]) ;
    printf ("Min length of ifact:                    %li\n",
            (LONG) lbl->ma57->info [17]) ;
    printf ("Min length of fact without compression :%li\n",
            (LONG) lbl->ma57->info [18]) ;
    printf ("Min length of ifact without compression:%li\n",
            (LONG) lbl->ma57->info [19]) ;
    printf ("Order of largest frontal matrix:        %li\n",
            (LONG) lbl->ma57->info [20]) ;
    printf ("Number of 2x2 numerical pivots:         %li\n",
            (LONG) lbl->ma57->info [21]) ;
    printf ("Number of fully-summed variables:       %li\n",
            (LONG) lbl->ma57->info [22]) ;
    printf ("Number of negative eigenvalues:         %li\n",
            (LONG) lbl->ma57->info [23]) ;
    printf ("Rank of the factorization:              %li\n",
            (LONG) lbl->ma57->info [24]) ;
    printf ("Floating point flops for elimin:        %e\n",
            lbl->ma57->rinfo [3]) ;
}
#endif
