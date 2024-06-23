/* Contents:

    1. SSMmult        - compute (A+D)x by rows
    2. SSMGivensMult  - compute x'G_1 G_2 ..., G_i = ith rotation in
                        diagonalization of a tridiagonal matrix
    3. SSMtGivensMult - compute G_1 G_2 ... x
    4. SSMDenseMult   - multiply dense matrix and vector 
    5. SSMtDenseMult  - multiply dense matrix transpose and vector
    6. SSMSSORmultP   - SSOR preconditioning operation with projection
    7. SSMSSORmult    - SSOR preconditioning operation */


#include "SSM.h"
/* ========================================================================== */
/* === SSMmult ============================================================== */
/* ========================================================================== */

/* Evaluate p = (A+D)x */

void SSMmult
(
    SSMFLOAT  *p, /* output vector of size n */
    SSMFLOAT  *x, /* input vector of size n */
    SSMFLOAT *Ax, /* numerical values in A excluding diagonal */
    SSMFLOAT  *D, /* diagonal of matrix */
    SSMINT   *Ai, /* row indices for each column of A */
    SSMINT   *Ap, /* Ap [j] = start of column j */
    SSMINT     n, /* dimension of matrix */
    SSMcom  *Com
)
{
    SSMINT j, k, l ;
    SSMFLOAT t ;
    k = 0 ;
    for (j = 0; j < n; j++)
    {
        l = Ap [j+1] ;
        t = D [j]*x [j] ;
        for (; k < l; k++)
        {
            t += Ax [k]*x [Ai [k]] ;
        }
        p [j] = t ;
    }
    Com->mults++ ;
}

/* ========================================================================== */
/* === SSMGivensMult ======================================================== */
/* ========================================================================== */

/* Multiply row vector x on the left by the Givens rotations generated
   during the QR algorithm.  The rotations are applied in the same order
   that they were generated during the QR algorithm. */

void SSMGivensMult
(
    SSMFLOAT   *x,    /* the vector to which the rotations are applied */
    SSMDiag   *DT     /* diagonalization structure for tridiagonal matrix */
)
{   

    SSMINT ngivens, nscale ;
    int pass, npass ;
    SSMINT i, j, k, f, n, *gk, *gf, *gi ;
    short *gz, *Gz ;
    SSMFLOAT xxt, xxs, *gx, *gy, *gs, *Gx, *Gy, *Gs ;

    /* ---------------------------------------------------------------------- */
    /* Read in the needed arrays */
    /* ---------------------------------------------------------------------- */
    
    /* input */
    npass = DT->its ;   /* number of QR iterations */
    if ( npass == 0 ) return ; /* diagonal matrix */
    n = DT->n ;
    gx = DT->gx ;       /* fast Givens x factor */
    gy = DT->gy ;       /* fast Givens y factor */
    gz = DT->gz ;       /* rotation type 1 or 2 */
    gs = DT->gs ;       /* scaling factors */
    gi = DT->gi ;       /* row indices for scaling factors */
    gk = DT->gk ;       /* index of last diagonal element in QR iteration */
    gf = DT->gf ;       /* index of first diagonal element in QR iteration */

    ngivens = 0 ;
    nscale = 0 ;
    for (pass = 0; pass < npass; pass++)
    {
        k = gk [pass] ;
        f = gf [pass] ;

        /* apply Givens rotation */
        Gx = gx+(ngivens-f) ;
        Gy = gy+(ngivens-f) ;
        Gz = gz+(ngivens-f) ;
        nscale++ ;
        i = abs (gi [nscale]) - 1 ;
        for (j = gf [pass]; j < k; j++)
        {
            if ( i == j )
            {
                if ( gi [nscale] > 0 )
                {
                    x [j] *= gs [nscale] ;
                    nscale++ ;
                    i = abs (gi [nscale]) - 1 ;
                    if ( i == j )
                    {
                        x [j+1] *= gs [nscale] ;
                        nscale++ ;
                        i = abs (gi [nscale]) - 1 ;
                    }
                }
                else
                {
                    x [j+1] *= gs [nscale] ;
                    nscale++ ;
                    i = abs (gi [nscale]) - 1 ;
                }
            }
            xxt = x [j] ;
            xxs = x [j+1] ;
            if ( Gz [j] == 1 )
            {
                x [j]   = xxt + Gx [j] * xxs ;
                x [j+1] = xxs - Gy [j] * xxt ;
            }
            else
            {
                x [j]   = Gx [j] * xxt + xxs ;
                x [j+1] = Gy [j] * xxs - xxt ;
            }
        }
        ngivens += (k-f) ;
    }

    /* final scaling */
    nscale++ ;
    Gs = gs+nscale ;
    for (i = 0; i < n ; i++) x [i] *= Gs [i] ;
}

/* ==========================================================================
   === SSMtGivensMult =======================================================
   ==========================================================================

   Multiply a column vector x on the right by the Givens rotations generated
   during the QR algorithm.  The rotations are applied in the reverse order
   that they were generated during the QR algorithm. */

void SSMtGivensMult
(
    SSMFLOAT    *x,    /* the vector to which the rotations are applied */
    SSMDiag    *DT     /* diagonalization structure for tridiagonal matrix */
)
{
    SSMINT jj, nscale, *gj ;
    int npass ;
    SSMINT i, j, k, kk, kf, n, *gi, *gk, *gf ;
    short *gz, *Gz ;
    SSMFLOAT *gx, *gy, *gs, *Gx, *Gy, *Gs ;

    /* input */
    npass = DT->its ;   /* number of QR iterations */
    if ( npass == 0 ) return ; /* diagonal matrix */
    gx = DT->gx ;       /* fast Givens x factor */
    gy = DT->gy ;       /* fast Givens y factor */
    gz = DT->gz ;       /* rotation type 1 or 2 */
    gs = DT->gs ;       /* scaling factors */
    gi = DT->gi ;       /* indexed by QR pass, points into scale factors */
    gj = DT->gj ;       /* indexed by QR pass, points into Givens factors */
    gk = DT->gk ;       /* index of last diagonal element in QR iteration */
    gf = DT->gf ;       /* index of first diagonal element in QR iteration */
    nscale = DT->nscale ; /* number of scaling operations */
    n = DT->n ;

    Gs = gs+nscale ;
    for (j = 0 ; j < n ; j++)
    {
        x [j] *= Gs [j] ;
    }
    nscale-- ;

    for (k = npass-1 ; k >= 0 ; k--)
    {
        jj = gj [k] ;
        kk = gk [k] -1 ;
        kf = gf [k] ;
        Gx = gx+(jj-kf) ;
        Gy = gy+(jj-kf) ;
        Gz = gz+(jj-kf) ;
        nscale-- ;
        i = abs (gi [nscale]) - 1 ;

        for (j = kk; j >= kf; j--)
        {
            double xt = x [j] ;
            double xs = x [j+1] ;

            if (Gz [j] == 1)
            {
                x [j]   = xt - Gy [j] * xs ;
    	        x [j+1] = xs + Gx [j] * xt ;
            }
            else
            {
                x [j]   = Gx [j] * xt - xs ;
                x [j+1] = Gy [j] * xs + xt ;
            }
            /* scale */
            if ( j == i )
            {
                if ( gi [nscale] < 0 )
                {
                    x [j+1] *= gs [nscale] ;
                    nscale-- ;
                    i = abs (gi [nscale]) - 1 ;
                    if ( i == j )
                    {
                        x [j] *= gs [nscale] ;
                        nscale-- ;
                        i = abs (gi [nscale]) - 1 ;
                    }
                }
                else
                {
                    x [j] *= gs [nscale] ;
                    nscale-- ;
                    i = abs (gi [nscale]) - 1 ;
                }
            }
        }
    }
}

/* ==========================================================================
   === SSMDenseMult =========================================================
   ==========================================================================

   Compute y = Vx where V is m by n */

void SSMDenseMult
(
    SSMFLOAT *y,     /* m by 1 product, output */
    SSMFLOAT *x,     /* n by 1 given vector */
    SSMFLOAT *V,     /* dense m by n matrix */
    SSMINT    m,     /* number of rows */
    SSMINT    n      /* number of columns */
)
{
    SSMINT i, j ;
    SSMFLOAT *Vp, t ;

    if ( n < 1 ) return ;
    t = x [0] ;
    for (i = 0; i < m; i++) y [i] = V [i]*t ;
    Vp = V+m ;
    for (j = 1; j < n; j++)
    {
        t = x [j] ;
        for (i = 0; i < m; i++) y [i] += Vp [i]*t ;
        Vp += m ;
    }
}

/* ==========================================================================
   === SSMtDenseMult ==========================================================
   ==========================================================================

   Compute y' = x'V where V is m by n */

void SSMtDenseMult
(
    SSMFLOAT *y,     /* n by 1 product, output */
    SSMFLOAT *x,     /* m by 1 given vector */
    SSMFLOAT *V,     /* dense m by n matrix */
    SSMINT    m,     /* number of rows */
    SSMINT    n      /* number of columns */
)
{
    SSMINT i, j ;
    SSMFLOAT *Vp, t ;

    Vp = V ;
    for (j = 0; j < n; j++)
    {
        t = 0. ;
        for (i = 0; i < m; i++)
        {
            t += Vp [i]*x [i] ;
        }
        y [j] = t ;
        Vp += m ;
    }
}

/* ==========================================================================
   === SSMSSORmultP =========================================================
   ==========================================================================

    Multiply b by the matrix associated with SSOR
    preconditioning for a linear system with matrix P(A + mu I)P
    where P projects a vector into the space orthogonal to x.
    See the documentation for the SSMSSORmult algorithm where we
    give the SSOR preconditioner for a matrix A. In order to handle
    the projection, we use the strategy explained on page 203 of the
    following paper (see Algorithms 2 and 3):
    W. W. Hager, Minimizing a quadratic over a sphere, SIAM Journal on
    Optimization, 12 (2001), pp. 188-208. Below d = diag (A) + mu,
    s = sqrt(d), p and q are vectors associated with the
    projection, and w = x/||x||. */

void SSMSSORmultP
(
    SSMFLOAT    *y,  /* the resulting vector */
    SSMFLOAT    *b,  /* vector to be multiplied by SSOR matrix */
    SSMFLOAT   *wj,  /* the first half of the SSOR multiplication operation */
    SSMFLOAT    *w,  /* w = x/||x|| */
    SSMFLOAT   *aj,  /* aj = Awj, product of A with wj */
    SSMFLOAT    mu,  /* multiplier */
    int    startup,  /* = 1 for starting multiplication, 0 otherwise */
    SSMProblem *PB,  /* problem specification */
    SSMcom    *Com   /* SSMcom structure */
)
{
    SSMINT k, l, *Ap, *Ap1, *Ai, *Au ;
    SSMINT j, n, n1 ;
    SSMFLOAT r, t, u, xj, yj, *Ax, *d, *p, *q, *s ;
    n  = PB->n ;  /* dimension of A */
    Ax = PB->x ;  /* numerical values in A (input) */
    Ai = PB->i ;  /* row indices for each column of A */
    Ap = PB->p ;  /* Ap [j] = start of column j */
    Au = PB->u ;  /* Au [j] = location right after last nonzero above diagonal*/
    Ap1= Ap+1 ;
    d = Com->SSOR ; /* d  = diag (C), C = (A + mu*I) */
    s = d+(1*n) ; /* s = sqrt (d) */
    p = d+(2*n) ; /* p = A*w - (w'*Aw)w */
    q = d+(3*n) ; /* q = (A + mu*I)w */
    n1 = n - 1 ;
/*  s = 0 ;
    t = 0 ;
    for i = n : -1 : 2
        y(i) = ( y(i) + w(i)*s + p(i)*t ) / d(i) ;
        s = s + q(i)*y(i) ;
        t = t + w(i)*y(i) ;
%       y(1:i-1) = y(1:i-1) - A(1:i-1,i)*y(i) ;
        y = y - U (:,i)*y(i) ; U = upper triangle of A
    end
    y(1) = ( y(1) + w(1)*s + p(1)*t ) / d(1) ; */

    if ( startup )
    {
        for (j = 0; j < n; j++) y [j] = b [j] ;
        goto StartUp ;
    }

    for (j = 0; j < n; j++) y [j] = s [j]*b [j] ;
    u = SSMZERO ;
    t = SSMZERO ;
    for (j = n1; ; j--)
    {
        xj = w [j] ;
        yj = (y [j] + xj*u + q [j]*t)/d [j] ;
        y [j] = yj ;
        if ( j <= 0 ) break ;

        t += (xj   * yj) ;
        u += p [j] * yj ;

        /* y = y - U (:,j) * yj */
        l = Au [j] ;
        for (k = Ap [j]; k < l; k++)
        {
            y [Ai [k]] -= Ax [k]*yj ;
        }
    }
    Com->mults += SSMHALF ;

    r = SSMZERO ;
    for (j = 0; j < n; j++) r += y [j]*w [j] ;
    /* note: aj and wj below are returned arrays */
    for (j = 0; j < n; j++)
    {
        t = y [j] - r*w [j] ;
        y [j] = t ;
        wj [j] = t ;
    }

    SSMmult (aj, wj, Ax, PB->D, Ai, Ap, n, Com) ; /* aj = A*wj */
    r = SSMZERO ;
    for (j = 0; j < n; j++)
    {
        t = aj [j] + y [j]*mu ;
        y [j] = t ;
        r += w [j]*t ;
    }
    for (j = 0; j < n; j++) y [j] -= r*w [j] ;

/*  s = 0 ;
    t = 0 ;
    for i = 1:n-1
        y(i) = ( y(i) + w(i)*s + q(i)*t) / d(i) ;
        s = s + p(i)*y(i) ;
        t = t + w(i)*y(i) ;
%       y(i+1:n) = y(i+1:n) - A(i+1:n,i)*y(i) ;
        y = y - L (: ,i)*y(i) ;
    end
    y(n) = ( y(n) + w(n)*s + q(n)*t ) / d(n) ; */

    StartUp: /* compute sqrt (d)*lower triangular SSOR stuff */
    u = SSMZERO ;
    t = SSMZERO ;
    for (j = 0; ; j++)
    {
        xj = w [j] ;

        yj = (y [j] + xj*u + p [j]*t) / d [j] ;
        y [j] = s [j]*yj ;
        if ( j >= n1 ) break ;

        t += xj    * yj ;
        u += q [j] * yj ;

        /* y = y - L (:,j) * yj */
        l = Ap1 [j] ;
        for (k = Au [j]; k < l ; k++)
        {
            y [Ai [k]] -= Ax [k]*yj ;
        }
    }
    Com->mults += SSMHALF ;
}

/* ==========================================================================
   === SSMSSORmult ==========================================================
   ==========================================================================

    Multiply b by the matrix associated with SSOR
    preconditioning for a linear system with matrix A.
    This code is based on formula (6.3) in the following paper:
 ***W. W. Hager, Iterative methods for nearly singular systems,
    SIAM Journal on Scientific Computing, 22 (2000), pp. 747-766.
    If this formula is combined with (6.2), we see that in the
    minimal residual algorithm, each iteration involves multiplication
    by the matrix

                    - sqrt(D) inv(L) A inv(L)' sqrt(D)

    Here D is the diagonal of A and L is the lower triangular
    matrix formed by elements of A on the diagonal and below the diagonal.
    In the code below we ignore the leading "-" sign. To correct for this,
    we add the correction term to x in the calling routine SSMSSOR
    (instead of subtract the term). See Algorithm 3 in the paper ***. */

void SSMSSORmult
(
    SSMFLOAT    *y,  /* the resulting vector */
    SSMFLOAT    *b,  /* vector to be multiplied by SSOR matrix */
    SSMFLOAT   *wj,  /* the first half of the SSOR multiplication operation */
    SSMFLOAT   *aj,  /* aj = Awj, product of A with wj */
    SSMFLOAT    *d,  /* diagonal of matrix */
    SSMFLOAT    *s,  /* square root of d */
    SSMFLOAT    mu,  /* diagonal safeguard */
    int    startup,  /* = TRUE (starting multiplication), FALSE (otherwise) */
    SSMProblem *PB,  /* problem specification */
    SSMcom    *Com   /* SSMcom structure */
)
{
    SSMINT k, l, *Ai, *Ap, *Ap1, *Au ;
    SSMINT j, n, n1 ;
    SSMFLOAT yj, *Ax ;
    n  = PB->n ;  /* dimension of A */
    Ax = PB->x ;  /* numerical values in A (input) */
    Ai = PB->i ;  /* row indices for each column of A */
    Ap = PB->p ;  /* Ap [j] = start of column j */
    Au = PB->u ;  /* Au [j] = location right after last nonzero above diagonal*/
    Ap1= Ap+1 ;
    n1 = n - 1 ;
/*  s = 0 ;
    t = 0 ;
    for i = n : -1 : 2
        y(i) = y(i) / d(i) ;
        y = y - U (:,i)*y(i) ; U = upper triangle of A
    end
    y(1) = y(1) / d(1) ; */

    if ( startup == SSMTRUE )
    {
        for (j = 0; j < n; j++) y [j] = b [j] ;
        goto StartUp ;
    }

    /* multiply by sqrt (D) */
    for (j = 0; j < n; j++) y [j] = s [j]*b [j] ;
    
    /* multiply by inv(L)' */
    for (j = n1; ; j--)
    {
        yj = y [j]/d [j] ;
        y [j] = yj ;
        if ( j <= 0 ) break ;

        /* y = y - U (:,j) * yj */
        l = Au [j] ;
        for (k = Ap [j]; k < l; k++)
        {
            y [Ai [k]] -= Ax [k]*yj ;
        }
    }
    Com->mults += SSMHALF ;

    /* note: aj and wj below are returned arrays */
    for (j = 0; j < n; j++) wj [j] = y [j] ;

    /* multiply by A, aj = A*wj */
    SSMmult (aj, wj, Ax, PB->D, Ai, Ap, n, Com) ;
    if ( mu == SSMZERO ) for (j = 0; j < n; j++) y [j] = aj [j] ;
    else                 for (j = 0; j < n; j++) y [j] = aj [j] + y [j]*mu ;

/*  for i = 1:n-1
        y(i) = y(i) / d(i) ;
        y = y - L (: ,i)*y(i) ;
    end
    y(n) = y(n) / d(n) ; */

    StartUp: /* compute sqrt (d)*lower triangular SSOR stuff */
    /* multiply by sqrt(D) inv(L) */
    for (j = 0; ; j++)
    {
        yj = y [j] / d [j] ;
        y [j] = s [j]*yj ;            /* y_j multiplied by sqrt(D_j) */
        if ( j >= n1 ) break ;

        /* y = y - L (:,j) * yj */
        l = Ap1 [j] ;
        for (k = Au [j]; k < l ; k++)
        {
            y [Ai [k]] -= Ax [k]*yj ;
        }
    }
    Com->mults += SSMHALF ;
}
