/* Contents:

    1. SSMboundary - sequential subspace method for ||x|| = r
    2. SSMinterior - SSM for ||x|| < r, solve Ax = -b
    3. SSMminresP  - minimum residual algorithm applied to SQP system
    4. SSMminres   - minimum residual algorithm applied to Ax + b = 0 */

#include "SSM.h"

/* ==========================================================================
   === SSMboundary ==========================================================
   ==========================================================================
    Starting with a guess for the solution of the following problem,
    use SSM to refine the solution estimate until the error tolerance
    is fulfilled:

         Minimize x'Ax + 2x'b subject to x'x = r^2

    For the SQP step in SSM we need to solve a linear system of the form

 (SQP)  P(A + mu I)P y = -P(b + (A+ muI)x)

    where x = prior iterate, P = I - X*X', X = x/||x||
    The updated iterate and multiplier are

            x_new = x + y
            mu_new = rho (x_new), rho (x) = - (b+Ax)'x/||x||^2.
                            
    We use preconditioned MINRES to solve (SQP). Since we
    know that A + mu I should be positive semidefinite at the optimal
    mu, we choose mu large enough to ensure that A + mu I
    is positive definite. As a consequence, diag (A) + mu I is positive.
    To apply SSOR preconditioning, we need the diagonal of P(A + mu I)P
    positive. Since A + mu I is positive definite, x'P(A + mu I)Px = 0
    if and only if x is a multiple of w.  The diagonal of P(A + mu I)P
    has a zero if and only if ei'P(A + mu I)P*ei = 0 for some column
    ei of the identity. Hence, diag (P(A + mu I)P) is positive if and
    only if w is not a column of the identity */

int SSMboundary /* return 0 (error tolerance satisfied)
                                -1 (min residual convergence failure in SQP)
                                -2 (SSM failed to converge)
                                -3 (error decay in SSM too slow)
                                -5 (failure of QR diagonalization)
                                -6 (insufficient space in QR method)
                                -8 (solution in interior ||x|| < r) */
(
    SSMFLOAT    *x, /* estimated solution to ball problem, a'x = 0 */
    SSMcom    *Com  /* SSMcom structure */
)
{
    int ssm_limit, done, status, BndLessThan, PrintLevel ;
    SSMINT j, n ;
    SSMFLOAT s, t, fj, vj, ball_tol, mu, mulow, rad, prior_err,
          *eig_res, *eig_guess, *xj, *r0, *Ax, *Av, *b, *V, *v ;
    SSMParm *Parm ; /* parameter structure */
    SSMProblem *PB ; /* problem structures  */

    PB = Com->PB ;
    n = PB->n ;
    rad = PB->rad ;
    b = PB->b ;
    Parm = Com->Parm ;
    BndLessThan = Parm->BndLessThan ;
    PrintLevel = Parm->PrintLevel ;
    ball_tol = Com->tol ;
    ssm_limit = Com->ssm_limit ;
    xj = Com->wx1 ;
    eig_res = Com->wx2 ;
    eig_guess = Com->wx3 ;
    V = Com->V ;
    v = Com->v ;
    Av = Com->Av ;
    Ax = Com->Ax ;
    r0 = Com->r0 ;

    while (Com->error > ball_tol)
    {
        Restart:
        prior_err = Com->error ;
        Com->ssm_its++ ;
        if ( Com->ssm_its > ssm_limit )
        {
            if ( PrintLevel >= 1 )
            {
                printf ("SSM refinement did not converge within %i \n"
                        "iterations as required by Parm->ssm_its_fac %e\n",
                         (int) ssm_limit, Parm->ssm_its_fac) ;
            }
            return (-2) ;
        }
        /* convert x to unit vector to facilitate projection computation */
        t = SSMONE/rad ;
        for (j = 0; j < n; j++) V [j] = t*x [j] ;

        /* A + mulow I should be positive definite */
        mulow = SSMMIN (Com->eig_error - Com->emin, -Parm->eig_lower) ;
        if ( mulow < -PB->Dmin ) mulow = -PB->Dmin + Parm->eps*PB->Amax ;
        mu = Com->mu ;
        if ( mulow > mu ) mu = mulow + (mulow-mu) ;/* increase mu if necessary*/

        if ( PrintLevel < 0 )
        {
            printf ("before minresP\n");
            fflush (stdout) ;
        }
        done = SSMminresP (xj, r0, x, Ax, b, V, mu, rad, PB, Com) ;
        if ( PrintLevel < 0 )
        {
            printf ("before subspace\n");
            fflush (stdout) ;
        }
        /* solve subspace problem, optimize over x, v, r0, xj */
        status = SSMsubspace (x, rad, Ax, v, r0, xj, NULL, 4, SSMFALSE, Com) ;
        if ( PrintLevel < 0 )
        {
            printf ("after subspace \n");
            fflush (stdout) ;
        }
        if ( PrintLevel >= 2 )
        {
            printf ("\nSSM iteration: %i ball_err: %e ball_tol: %e\n",
                     (int) Com->ssm_its, Com->error, ball_tol) ;
            printf ("emin: %e mu: %e eig_err: %e\n",
                     Com->emin, Com->mu, Com->eig_error) ;
        }
        if ( status ) return (status) ; /* failure of QR diagonalization */
        mu = Com->mu ;
        if ( done == 1 )
        {
            /* check if constraint inactive */
            if ( !Com->Active && (SSMMAX(-rad*Com->mu, SSMZERO) > ball_tol) )
                return (-8) ;
            /* otherwise constraint is active. Check to see if there
               was excessive growth of error in the evaluation of
               ball_err in minresP */
            if ( Com->error > ball_tol ) goto Restart ;
            return (0) ;
        }
        if ( done < 0 ) return (done) ;


        /* exit if constraint inactive */
        if ( BndLessThan && !Com->Active) return (-8) ;

        mu = Com->mu ;
        mulow = SSMMIN (Com->eig_error - Com->emin, -Parm->eig_lower) ;
        if ( mulow < -PB->Dmin ) mulow = -PB->Dmin + Parm->eps*PB->Amax ;
        /* check if eigenvalue should be refined when mulow > mu */
        if ( mulow > mu )
        {
            if ( PrintLevel >= 2 )
            {
                printf ("mulow %e > mu %e eig_err: %e ball_err: %e\n",
                      mulow, mu, Com->eig_error, Com->error) ;
            }
            mu = mulow + (mulow - mu) ;
            /* refine eigenvalue estimate when eigen error > x error */
            if ( Com->eig_error > Parm->eig_refine*Com->error/rad )
            {
                if ( PrintLevel >= 2 ) printf ("refine eigenpair\n") ;
                /* be sure that A + mu I is positive definite */
                mu += ((SSMFLOAT) 100)*Parm->eps*PB->Amax ;

                if ( Parm->IPM == SSMTRUE ) /* apply inverse power method */
                {
                    /* prepare for inverse power method iteration applied to
                       (A + mu I)x = v, current estimate for eigenvector.
                       Initial guess is x = s v where s is chosen to
                       minimize the residual ||(A + mu I)(s v) - v||.
                       Defining f = (A + mu I)v, we have s = f'v/f'f.
                       We now compute s, the residual (A + mu I)(s v) - v,
                       and sv the starting guess for the inverse power scheme */
                    s = SSMZERO ;
                    t = SSMZERO ;
                    for (j = 0; j < n; j++)
                    {
                        vj = v [j] ;
                        fj = Av [j] + mu*vj ;
                        s += fj*vj ;
                        t += fj*fj ;
                    }
                    s = s/t ;
                    t = mu*s - SSMONE ;
                    for (j = 0; j < n; j++)
                    {
                        vj = v [j] ;
                        /* starting guess */
                        eig_guess [j] = s*vj ;
                        /* = residual (A + mu I)(s v) - v */
                        eig_res [j] = s*Av [j] + t*vj ;
                    }
                    /* MINRES applied to (A + mu I)x = v */
                    if ( PrintLevel < 0 )
                    {
                        printf ("before minres\n");
                        fflush (stdout) ;
                    }
                    SSMminres (xj, eig_res, eig_guess, 1.e20, mu,
                               SSMTRUE, PB, Com);
                    if ( PrintLevel < 0 )
                    {
                        printf ("after minres\n");
                        fflush (stdout) ;
                    }
                    /* minimize over subspace spanned by previous vectors and
                       xj = solution - v */
                }
                else         /* use SQP method to improve eigenvalue estimate */
                {
                    /* The smallest eigenvalue is the solution to the
                       problem: min x'Ax subject to ||x|| = 1. The
                       multiplier associated with the constraint is the
                       eigenvalue. We apply one iteration of the SQP
                       method to the SQP system using the current
                       eigenvalue and eigenvector estimates as the
                       linearization points. We solve the SQP system
                       using SSMminresP (MINRES with projection) */
               
                    s = SSMZERO ;
                    for (j = 0; j < n; j++)
                    {
                        t = v [j] ;
                        s += t*t ;
                    }
                    s = sqrt (s) ;
                    if ( s > SSMZERO ) s = SSMONE/s ;
                    for (j = 0; j < n; j++) v [j] *= s ;
                    for (j = 0; j < n; j++) Av [j] *= s ;
                    if ( PrintLevel < 0 )
                    {
                        printf ("before eig minresP\n");
                        fflush (stdout) ;
                    }
                    SSMminresP (xj, Av, v, Av, NULL, v, mu, SSMONE, PB, Com);
                    if ( PrintLevel < 0 )
                    {
                        printf ("after eig minresP\n");
                        fflush (stdout) ;
                    }
                }
                /* add the new eigenvector estimate to the subspace and
                   reoptimize */
                if ( PrintLevel < 0 )
                {
                    printf ("before subspace\n");
                    fflush (stdout) ;
                }
                status = SSMsubspace (x, rad, Ax, v, r0, xj, xj, 5, SSMFALSE,
                                      Com) ;
                if ( PrintLevel < 0 )
                {
                    printf ("after subspace\n");
                    fflush (stdout) ;
                }
                if ( PrintLevel >= 2 )
                {
                    printf ("\nSSM iteration: %i ball_err: %e ball_tol: %e\n",
                             (int) Com->ssm_its, Com->error, ball_tol) ;
                    printf ("emin: %e mu: %e eig_err: %e\n",
                             Com->emin, Com->mu, Com->eig_error) ;
                }
                if ( status ) return (status) ; /* QR diagonalization failure */
            }
        }
        if ( Com->error <= ball_tol )
        {
            /* check if constraint inactive */
            if ( !Com->Active && (SSMMAX(-rad*mu, SSMZERO) > ball_tol) )
                return (-8) ;
            /* otherwise constraint is active and error tolerance satisfied */
            return (0) ;
        }

        /* return if convergence is slow */
        if ( Com->error > prior_err*Parm->ssm_decay ) return (-3) ;
    }
    return (0) ;
}

/* ==========================================================================
   === SSMinterior ==========================================================
   ==========================================================================
    Starting with an interior guess for the solution of the following problem,
    apply precondition MINRES to Ax = -b.

         Minimize x'Ax + 2x'b subject to x'x <= r^2

    If the iterate exceeds the norm constraint, then switch to
    SSMboundary. Otherwise, continue to apply MINRES iterations
    until convergence is achieved. If the convergence criterion is not
    satisfied, then we return to the Lanczos process and try to compute
    a better starting guess */

int SSMinterior /* return 0 (error tolerance satisfied)
                                -1 (minimum residual algorithm did not converge)
                                -5 (failure of QR diagonalization)
                                -6 (insufficient space in QR method)
                                -9 (solution on boundary ||x|| = r) */
(
    SSMFLOAT    *x, /* estimated solution to sphere constrained problem */
    SSMcom    *Com  /* SSMcom structure */
)
{
    int done, status, PrintLevel ;
    SSMINT j, n ;
    SSMFLOAT rad, t, normx, ball_tol, *xj, *y ;
    SSMParm *Parm ; /* parameter structure */
    SSMProblem *PB ;

    PB = Com->PB ;
    n = PB->n ;
    rad = PB->rad ;
    Parm = Com->Parm ;
    PrintLevel = Parm->PrintLevel ;
    ball_tol = Com->tol ;
    xj = Com->wx1 ;
    y = Com->MINRES ;

    if ( PB->Dmin == SSMZERO )
    {
        t = Parm->eps*PB->Amax ;
        for (j = 0; j < n; j++)
        {
            if ( PB->D [j] == SSMZERO ) PB->D [j] = t ;
        }
        PB->Dmin = t ;
    }
    Com->mu = SSMZERO ;

    done = SSMminres (xj, Com->r0, x, rad, SSMZERO, SSMFALSE, PB, Com) ;
    normx = Com->normx ;
    if ( done == 0  )               /* error tolerance achieved */
    {
        if ( normx <= rad )
        {
            for (j = 0; j < n; j++) x [j] = y [j] ;
            return (0) ;
        }
        else if ( normx - rad <= ball_tol )
        {
            t = rad/normx ;
            for (j = 0; j < n; j++) x [j] = t*y [j] ;
            return (0) ;
        }
        /* else treat as boundary solution */
        Com->error = SSMMAX (Com->error, normx - rad) ;
    }
    if ( done == -1 ) return (-1) ; /* MINRES convergence failure */

    /* solution lies on boundary - to ensure convergence, minimize objective
       over subspace, then return to boundary routine */
    normx = SSMZERO ;
    for (j = 0; j < n; j++)
    {
        t = x [j] ;
        normx += t*t ;
    }
    normx = sqrt (normx) ;
    status = SSMsubspace (x, normx, Com->Ax, Com->v, Com->r0, xj,
                          NULL, 4, SSMFALSE,Com);
    if ( PrintLevel >= 1 ) printf ("interior solution reaches boundary\n") ;

    if ( status ) return (status) ;  /* QR method convergence failure */

    /* normalize x if constraint not active */
    if ( !Com->Active )
    {
        normx = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = x [j] ;
            normx += t*t ;
        }
        normx = sqrt (normx) ;
        if ( normx > SSMZERO )
        {
            t = rad/normx ;
            for (j = 0; j < n; j++) x [j] *= t ;
        }
    }
    return (-9) ;                             /* solution on boundary */
}
 
/* ==========================================================================
   === SSMminresP (minimum residual algorithm with projection) =============
   ==========================================================================
    Apply MINRES to the SQP system:

       (A+muI)y + (x) nu = -[(A+muI)x + b]
         x'y             = 0

    Define X = x/||x|| and P = I - XX'. Make the change of variables
    y = Pz. Hence, the SQP system is equivalent to solving
   
            P(A + mu I)Pz = -P r0 where r0 = (A+muI)x + b.

    We solve this system using Algorithm 3 in the paper
    W. W. Hager, Iterative methods for nearly singular linear systems,
    SIAM Journal on Scientific Computing, 22 (2000), pp. 747-766.
   ==========================================================================*/

int SSMminresP /* return 1 if convergence tolerance met,
                                 0 if decay tolerance met
                                -1 for min residual convergence failure
                                   (too many iterations) */
(
    SSMFLOAT   *xj, /* computed SQP iterate */
    SSMFLOAT   *r0, /* b + Ax, cost gradient at starting point */
    SSMFLOAT    *x, /* solution estimate */
    SSMFLOAT   *Ax, /* A*x */
    SSMFLOAT    *b,
    SSMFLOAT    *X, /* x/||x|| */
    SSMFLOAT    mu, /* multiplier for the constraint */
    SSMFLOAT   rad, /* radius r */
    SSMProblem *PB, /* problem specification */
    SSMcom    *Com  /* pointer to SSMcom structure */
)
{
    int it, it_limit, PrintLevel ;
    SSMINT j, n ;
    SSMFLOAT ball_tol, ball_err, sqp_err, s, t, tj,
          u, uj, uj1, ej, c, z1, nu, dj,
          *Ay, *r, *aj, *vj, *vj1, *wj, *pr, *rj1, *rj2, *zj1, *zj2,
          *d, *D, *ds, *p, *q, *y, Qj1 [4], Qj2 [2] ;
    SSMParm *Parm ;
    Parm = Com->Parm ;
    PrintLevel = Parm->PrintLevel ;
    n = PB->n ;
    D = PB->D ;
    ball_tol = Com->tol ;
    y = Com->MINRES ;
    Ay = y+(1*n) ;
    r  = y+(2*n) ;
    aj = y+(3*n) ;
    vj = y+(4*n) ;
    vj1= y+(5*n) ;
    wj = y+(6*n) ;
    pr = y+(7*n) ;
    rj1= y+(8*n) ;
    rj2= y+(9*n) ;
    zj1= y+(10*n) ;
    zj2= y+(11*n) ;

    d = Com->SSOR ; /* d  = diag (A + mu*I) */
    ds= d+(1*n) ;
    p = d+(2*n) ;
    q = d+(3*n) ; /* q = (A + mu*I)w */

    s = SSMZERO ;
    for (j = 0; j < n; j++) s += X [j]*r0 [j] ;
    for (j = 0; j < n; j++) r [j] = s*X [j] - r0 [j] ;
    t = SSMZERO ;
    for (j = 0; j < n; j++) t += X [j]*Ax [j] ;
    t /= rad ;
    u = SSMONE/rad ;
    for (j = 0; j < n; j++)
    {
        s = u*Ax [j] ;
        p [j]  = s - X [j]*t ;
        q [j]  = s + X [j]*mu ;
    }
    /* compute diagonal of SSOR matrix */
    s = SSMINF ;
    for (j = 0; j < n; j++)
    {
        t = mu + D [j] - (q [j] + p [j])*X [j] ;
        if ( t < s ) s = t ;
        d [j] = t ;
    }
#ifndef NDEBUG
    if ( s <= SSMZERO )
    {
        DPRINT (("Negative diagonal in SSOR: %e\n", s)) ;
        DPRINT (("this only occurs in unusual situations\n"
                "(for example, the eigenvalue estimate is poor)\n")) ;
/*      if ( s < SSMZERO ) SOPTERROR ("stop!") ;*/
    }
#endif

    if ( s <= SSMZERO ) /* try to make diagonal positive */
    {
        DPRINT (("Warning: negative diagonal %e in SSOR, try to recover\n", s));
        /* if mu is increased by t, the diagonal increases by t(1-Xj^2) */
        u = SSMZERO ;
        c = (SSMFLOAT) 100 * Parm->eps * PB->Amax ; /* min diagonal element */
        for (j = 0; j < n; j++)
        {
            if ( d [j] <= c )
            {
                t = SSMONE - X [j]*X [j] ;
                if ( t > SSMZERO )
                {
                    t = (c-d [j])/t ;
                    if ( t > u ) u = t ;
                }
                else
                {
                    u = -SSMONE ;
                    break ;
                }
            }
        }
        if ( u < SSMZERO ) /* the matrix cannot be fixed, force d > 0 */
        {
            DPRINT (("Warning: Xj/||X|| > 1??, X [%i] = %e\n", j, X [j]));
            s = c - s ;
            for (j = 0; j < n; j++) d [j] += s ; /* make diagonal >= c */
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                t = X [j] ;
                q [j] += t*u ;
                d [j] += u*(SSMONE-t*t) ;
            }
        }
    }
    for (j = 0; j < n; j++) ds [j] = sqrt (d [j]) ;

    SSMSSORmultP (vj, r, wj, X, aj, mu, (int) 1, PB, Com) ;
    uj = SSMZERO ;
    ej = SSMZERO ;
    for (j = 0; j < n; j++)
    {
        ej += vj [j]*vj [j] ;
        xj [j] = SSMZERO ;
        Ay [j] = SSMZERO ;
        rj1 [j] = SSMZERO ;
        rj2 [j] = SSMZERO ;
        zj1 [j] = SSMZERO ;
        zj2 [j] = SSMZERO ;
        vj1 [j] = SSMZERO ;
    }
    ej = sqrt (ej) ;
    for (j = 0; j < n; j++) vj [j] /= ej ;
    Qj1 [0] = SSMONE ;
    Qj1 [1] = SSMZERO ;
    Qj1 [2] = SSMZERO ;
    Qj1 [3] = SSMONE ;
    Qj2 [0] = SSMZERO ;
    Qj2 [1] = SSMZERO ;
    /* do not perform more iterations than was done in Lanczos process */
    if ( PB->DT->nalloc == Parm->Lanczos_bnd ) it_limit = Com->minres_limit ;
    else               it_limit = SSMMIN (Com->minres_limit, 2*PB->DT->nalloc) ;
    for (it = 0; it < it_limit; it++)
    {
        Com->minres_its++ ;
        SSMSSORmultP (pr, vj, wj, X, aj, mu, (int) 0, PB, Com) ;
        dj = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            dj += pr [j]*vj [j] ;
        }
        uj1 = uj ;
        uj = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = pr [j] - vj [j]*dj - vj1 [j]*uj1 ;
            pr [j] = t ;
            uj += t*t ;
        }
        uj = sqrt (uj) ;
        for (j = 0; j < n; j++)
        {
            vj1 [j] = vj [j] ;
            vj [j] = pr [j]/uj ;
        }
        tj = Qj2 [0]*uj1 ;
        t  = Qj2 [1]*uj1 ;
        uj1 = Qj1 [0]*t + Qj1 [2]*dj ;
        dj  = Qj1 [1]*t + Qj1 [3]*dj ;
        Qj2 [0] = Qj1 [2] ;
        Qj2 [1] = Qj1 [3] ;
        t = sqrt (dj*dj + uj*uj) ;
        if ( t != SSMZERO )
        {
            c = dj/t ;
            s = uj/t ;
        }
        else
        {
            c = SSMONE ;
            s = SSMZERO ;
        }
        Qj1 [0] = c ;
        Qj1 [1] =-s ;
        Qj1 [2] = s ;
        Qj1 [3] = c ;
        z1 = ej*c ;
        ej = -s*ej ;
        dj = c*dj + s*uj ;
        u = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = (wj [j] - zj1 [j]*uj1 - zj2 [j]*tj)/dj ;
            zj2 [j] = zj1 [j] ;
            zj1 [j] = t ;
            s = xj [j] + t*z1 ;
            xj [j] = s ;
            t = x [j] + s ;
            y [j] = t ;
            u += t*t ;
            t = (aj [j] - rj1 [j]*uj1 - rj2 [j]*tj)/dj ;
            rj2 [j] = rj1 [j] ;
            rj1 [j] = t ;
            Ay [j] += t*z1 ;
        }

        u = rad/sqrt (u) ;
        /* y = x_old + xj, the updated SQP iteration, u is the scale
           factor to ensure that y has the correct norm */
        for (j = 0; j < n; j++)
        {
            y [j] *= u ;
            wj [j] = u*(Ax [j] + Ay [j]) ;
        }
        t = SSMZERO ;
        if ( b == NULL )
        {
            for (j = 0; j < n; j++) t += wj [j]*y[j] ;
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                s = b [j] + wj [j] ; /* wj = Ay */
                wj [j] = s ;
                t += s*y[j] ;
            }
        }
        nu = -t/(rad*rad) ;
        Com->mu = nu ;

        /* ball_err = ||b + Ay + nu*y||, nu chosen to minimize residual */
        ball_err = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = y [j]*nu + wj [j] ;
            ball_err += t*t ;
        }
        ball_err = sqrt (ball_err) ;
        Com->error = ball_err ;
        /* include factor .5 in ball_tol since ball_err polluted by error */
        if ( ball_err <= .5*ball_tol ) return (1) ;
        s = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = r0 [j] + Ay [j] + xj [j]*mu ;
            /* r = b + A(x+xj) + mu*xj = r0 + (A + mu I)xj, xj orthogonal to x*/
            r [j] = t ;
            s += t*X [j] ;
        }

        /* the error in the the SQP linear system is sqp_err = ||Pr||,
           which we compute below.  The error along x can be set to zero
           by an appropriate choice of nu in the SQP system
           (A + mu I)xj + nu x = -r0  */
        sqp_err = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = r [j] - X [j]*s ;
            sqp_err += t*t ;
        }
        sqp_err = sqrt (sqp_err) ;
        if ( PrintLevel >= 3 )
        {
            printf ("minres it: %i err: %e ball_err: %e ball_tol: %e\n",
                     (int) it, sqp_err, ball_err, ball_tol) ;
        }
        if ( sqp_err <= Parm->sqp_decay*ball_err ) return (0) ;
    }
    if ( PrintLevel >= 2 )
    {
        printf ("after %i minres iterations, convergence not achieved\n",
                 (int) it) ;
    }
    return (-1) ;
}
/* ==========================================================================
   === SSMminres ===========================================================
   ==========================================================================
    Solve the unconstrained optimization problem

         Minimize x'(A + mu I)x + 2x'b

   This differs from SSMminresP in that there is no projection.
   If this is routine called from the interior point routine,
   then mu is zero. If this is called by the boundary routine,
   then we are refining the estimate for the eigenvector associated with the
   smallest eigenvalue. In the refinement process, we apply one iteration
   of the inverse power method. When called from the interior point
   routine, convergence is achieved either when the requested error
   tolerance has been satisfied, or when the solution norm exceeds
   the constraint radius. For the inverse power method, convergence is
   achieved when the relative residual for the eigenvector estimate
   decreases by the factor eig_decay.
   ==========================================================================*/

int SSMminres /* return 0 if ||Ax + b|| <= ball_tol
                               -1 for convergence failure (too many iterations)
                               -2 if ||x|| > r */
(
    SSMFLOAT   *xj, /* computed SQP iterate */
    SSMFLOAT   *r0, /* b + Ax, cost gradient at starting point */
    SSMFLOAT    *x, /* solution estimate */
    SSMFLOAT   rad, /* radius of sphere */
    SSMFLOAT    mu, /* safeguard, (A + mu I) positive definite */
    int        IPM, /* TRUE (inverse power method), FALSE (interior point) */
    SSMProblem *PB, /* problem specification */
    SSMcom    *Com  /* pointer to SSMcom structure */
)
{
    int it, it_limit, PrintLevel ;
    SSMINT j, n ;
    SSMFLOAT ball_tol, ball_err, eig_err, eig_tol, normx, s, t, tj,
          uj, uj1, ej, c, z1, dj,
          *Ay, *r, *aj, *vj, *vj1, *wj, *pr, *rj1, *rj2, *zj1, *zj2,
          *d, *D, *ds, *y, Qj1 [4], Qj2 [2] ;
    SSMParm *Parm ;
    Parm = Com->Parm ;
    eig_tol = Parm->eig_decay ;
    PrintLevel = Parm->PrintLevel ;
    n = PB->n ;
    D = PB->D ;
    ball_tol = Com->tol ;
    y = Com->MINRES ;
    Ay = y+(1*n) ;
    r  = y+(2*n) ;
    aj = y+(3*n) ;
    vj = y+(4*n) ;
    vj1= y+(5*n) ;
    wj = y+(6*n) ;
    pr = y+(7*n) ;
    rj1= y+(8*n) ;
    rj2= y+(9*n) ;
    zj1= y+(10*n) ;
    zj2= y+(11*n) ;

    d = Com->SSOR ; /* d  = diag (A + mu*I) */
    ds= d+(1*n) ;

#ifndef NDEBUG
    s = SSMINF ;
    for (j = 0; j < n; j++) s = SSMMIN (s, D [j] + mu) ;
    if ( s <= SSMZERO )
    {
        DPRINT (("Nonpositive diagonal in SSOR: %e\n", s)) ;
        SOPTERROR ("stop!") ;
    }
#endif

    for (j = 0; j < n; j++) r [j] = -(r0 [j] + mu*x [j]) ;
    for (j = 0; j < n; j++)
    {
        t = D [j] + mu ;
        ds [j] = sqrt (t) ;
        d [j] = t ;
    }

    SSMSSORmult (vj, r, wj, aj, d, ds, mu, SSMTRUE, PB, Com) ;
    uj = SSMZERO ;
    ej = SSMZERO ;
    for (j = 0; j < n; j++)
    {
        ej += vj [j]*vj [j] ;
        xj [j] = SSMZERO ;
        Ay [j] = SSMZERO ;
        rj1 [j] = SSMZERO ;
        rj2 [j] = SSMZERO ;
        zj1 [j] = SSMZERO ;
        zj2 [j] = SSMZERO ;
        vj1 [j] = SSMZERO ;
    }
    ej = sqrt (ej) ;
    for (j = 0; j < n; j++) vj [j] /= ej ;
    Qj1 [0] = SSMONE ;
    Qj1 [1] = SSMZERO ;
    Qj1 [2] = SSMZERO ;
    Qj1 [3] = SSMONE ;
    Qj2 [0] = SSMZERO ;
    Qj2 [1] = SSMZERO ;
    it_limit = Com->minres_limit ;
    for (it = 0; it < it_limit; it++)
    {
        Com->minres_its++ ;
        SSMSSORmult (pr, vj, wj, aj, d, ds, mu, SSMFALSE, PB, Com) ;
        dj = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            dj += pr [j]*vj [j] ;
        }
        uj1 = uj ;
        uj = SSMZERO ;
        for (j = 0; j < n; j++)
        {
            t = pr [j] - vj [j]*dj - vj1 [j]*uj1 ;
            pr [j] = t ;
            uj += t*t ;
        }
        uj = sqrt (uj) ;
        for (j = 0; j < n; j++)
        {
            vj1 [j] = vj [j] ;
            vj [j] = pr [j]/uj ;
        }
        tj = Qj2 [0]*uj1 ;
        t  = Qj2 [1]*uj1 ;
        uj1 = Qj1 [0]*t + Qj1 [2]*dj ;
        dj  = Qj1 [1]*t + Qj1 [3]*dj ;
        Qj2 [0] = Qj1 [2] ;
        Qj2 [1] = Qj1 [3] ;
        t = sqrt (dj*dj + uj*uj) ;
        if ( t != SSMZERO )
        {
            c = dj/t ;
            s = uj/t ;
        }
        else
        {
            c = SSMONE ;
            s = SSMZERO ;
        }
        Qj1 [0] = c ;
        Qj1 [1] =-s ;
        Qj1 [2] = s ;
        Qj1 [3] = c ;
        z1 = ej*c ;
        ej = -s*ej ;
        dj = c*dj + s*uj ;
        for (j = 0; j < n; j++)
        {
            t = (wj [j] - zj1 [j]*uj1 - zj2 [j]*tj)/dj ;
            zj2 [j] = zj1 [j] ;
            zj1 [j] = t ;
            s = xj [j] + t*z1 ;
            xj [j] = s ;
            t = x [j] + s ;
            y [j] = t ;
            t = (aj [j] - rj1 [j]*uj1 - rj2 [j]*tj)/dj ;
            rj2 [j] = rj1 [j] ;
            rj1 [j] = t ;
            Ay [j] += t*z1 ;
        }
        /* compute ball_err = ||b + A(x+xj)|| = ||r0 + Axj|| if b != NULL
                            = ||A(x+xj)||                    otherwise
                   sqp_err  = ||r0 + mu(x+xj)||              if b != NULL
                            = ||A(x+xj) + mu(x+xj)||         otherwise */

        /* termination conditions for Ax = b and for the inverse power method
           are slightly different. For Ax = b, we continue until the residual
           satisfies the convergence criterion. For the inverse power
           method, we require that the residual norm is less than eig_tol */

        if ( IPM == SSMTRUE )
        {
            eig_err = SSMZERO ;
            normx = SSMZERO ;
            for (j = 0; j < n; j++)
            {
                s = y [j] ;
                t = r0 [j] + Ay [j] + mu*s ;   /* Ay = Axj */
                eig_err += t*t ;
                normx += s*s ;
            }
            normx = sqrt (normx) ;
            eig_err = sqrt (eig_err)/SSMMAX (SSMONE, normx) ;
            if ( PrintLevel >= 3 )
            {
                printf ("EIGRES it: %i eig_res: %e eig_tol: %e\n",
                     (int) it, eig_err, eig_tol) ;
            }
            if ( (normx > rad) || (eig_err <= eig_tol) ) return (0) ;
            /* return when constraint is violated
               or the relative residual <= eig_tol */
        }
        else
        {
            ball_err = SSMZERO ;
            normx = SSMZERO ;
            for (j = 0; j < n; j++)
            {
                t = r0 [j] + Ay [j] ;   /* Ay = Axj */
                ball_err += t*t ;
                s = y [j] ;
                normx += s*s ;
            }
            ball_err = sqrt (ball_err) ;
            normx = sqrt (normx) ;
            Com->normx = normx ;
            Com->error = ball_err ;
            if ( PrintLevel >= 3 )
            {
                printf ("MINRES it: %i ball_err: %e ball_tol: %e\n",
                         (int) it, ball_err, ball_tol) ;
            }
            /* include factor .5 since ball_err polluted by errors */
            if ( ball_err <= .5*ball_tol ) return (0) ;
            /* return when constraint is sufficiently violated */
            if ( normx > rad*Parm->radius_flex ) return (-2) ;
        }
    }
    if ( PrintLevel >= 2 )
    {
        printf ("after %i iterations, MINRES has not converged\n", it) ;
    }
    return (-1) ;
}
