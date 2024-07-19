#include "cg_descent.h"

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Utility Routines
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* =========================================================================
   ==== cg_matvec ==========================================================
   =========================================================================
   Compute y = A*x or A'*x where A is a dense rectangular matrix
   ========================================================================= */
void cg_matvec
(
    CGFLOAT *y, /* product vector */
    CGFLOAT *A, /* dense matrix */
    CGFLOAT *x, /* input vector */
    CGINT    n, /* number of columns of A */
    CGINT    m, /* number of rows of A */
    int      w  /* T => y = A*x, F => y = A'*x */
)
{
/* if the blas have not been installed, then hand code the product */
#ifdef NOBLAS
    CGINT j, l ;
    l = 0 ;
    if ( w )
    {
        cg_scale0 (y, A, x [0], (int) m) ;
        for (j = 1; j < n; j++)
        {
            l += m ;
            cg_daxpy0 (y, A+l, x [j], (int) m) ;
        }
    }
    else
    {
        for (j = 0; j < n; j++)
        {
            y [j] = cg_dot0 (A+l, x, (int) m) ;
            l += m ;
        }
    }
#else

/* if the blas have been installed, then possibly call gdemv */
    if ( !w && (w || (m >= MATVEC_START)) )
    {
        BLAS_INT int_one = 1 ;
        BLAS_INT M, N ;
        CGFLOAT float_one, float_zero ;
        M = (BLAS_INT) m ;
        N = (BLAS_INT) n ;
        float_zero = (CGFLOAT) 0 ;
        float_one = (CGFLOAT) 1 ;
        /* only use transpose mult with blas
          SOPT_DGEMV ("n", &M, &N, one, A, &M, x, blas_one, zero, y,blas_one);*/
        SOPT_DGEMV ("t", &M, &N, &float_one , A, &M, x, &int_one, &float_zero,
                   y, &int_one) ;
    }
    else
    {
        CGINT j, l ;
        l = 0 ;
        if ( w )
        {
            cg_scale (y, A, x [0], m) ;
            for (j = 1; j < n; j++)
            {
                l += m ;
                cg_daxpy (y, A+l, x [j], m) ;
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                y [j] = cg_dot0 (A+l, x, (int) m) ;
                l += m ;
            }
        }
    }
#endif

    return ;
}

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       End of routines that could use the BLAS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* =========================================================================
   ==== cg_update_2 ========================================================
   =========================================================================
   Set gold = gnew (if not equal), compute 2-norm^2 of gnew, and optionally
      set d = -gnew
   ========================================================================= */
CGFLOAT cg_update_2
(
    CGFLOAT *gold, /* old g */
    CGFLOAT *gnew, /* new g */
    CGFLOAT    *d, /* d */
    CGINT       n  /* length of vectors */
)
{
    CGINT i, n5 ;
    CGFLOAT s, t ;
    t = CGZERO ;
    n5 = n % 5 ;

    if ( gold != NULL )
    {
        for (i = 0; i < n5; i++)
        {
            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
        }
        for (; i < n; )
        {
            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            gold [i] = s ;
            d [i] = -s ;
            i++ ;
        }
    }
    else
    {
        for (i = 0; i < n5; i++)
        {
            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
        }
        for (; i < n; )
        {
            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;

            s = gnew [i] ;
            t += s*s ;
            d [i] = -s ;
            i++ ;
        }
    }
    return (t) ;
}

/* =========================================================================
   ==== cg_update_beta =====================================================
   =========================================================================
   compute: ykPyk  = (newproj - oldproj)*(newproj - oldproj),
            gkPyk  =  newproj           *(newproj - oldproj)
   update: oldproj = newproj
   ========================================================================= */
void cg_update_beta
(
    CGFLOAT *oldproj,
    CGFLOAT *newproj,
    CGFLOAT   *GkPyk,
    CGFLOAT   *YkPyk,
    CGINT          n  /* length of vectors */
)
{
    CGINT i, n5 ;
    CGFLOAT gkPyk, ykPyk ;
    gkPyk = CGZERO ;
    ykPyk = CGZERO ;
    n5 = n % 5 ;

    for (i = 0; i < n5; i++)
    {
        CGFLOAT t, p ;

        t = newproj [i] ;
        p = t - oldproj [i] ;
        oldproj [i] = t ;

        gkPyk +=  t*p ;
        ykPyk +=  p*p ;
    }
    for (; i < n; )
    {
        CGFLOAT t, p ;

        t = newproj [i] ;
        p = t - oldproj [i] ;
        oldproj [i] = t ;

        gkPyk +=  t*p ;
        ykPyk +=  p*p ;
        i++ ;
        /* 1 ------------------------- */

        t = newproj [i] ;
        p = t - oldproj [i] ;
        oldproj [i] = t ;

        gkPyk +=  t*p ;
        ykPyk +=  p*p ;
        i++ ;
        /* 2 ------------------------- */

        t = newproj [i] ;
        p = t - oldproj [i] ;
        oldproj [i] = t ;

        gkPyk +=  t*p ;
        ykPyk +=  p*p ;
        i++ ;
        /* 3 ------------------------- */

        t = newproj [i] ;
        p = t - oldproj [i] ;
        oldproj [i] = t ;

        gkPyk +=  t*p ;
        ykPyk +=  p*p ;
        i++ ;
        /* 4 ------------------------- */

        t = newproj [i] ;
        p = t - oldproj [i] ;
        oldproj [i] = t ;

        gkPyk +=  t*p ;
        ykPyk +=  p*p ;
        i++ ;
        /* 5 --------------------------- */
    }
    *GkPyk = gkPyk ;
    *YkPyk = ykPyk ;
    return ;
}

/* =========================================================================
   ==== cg_update_d ========================================================
   =========================================================================
   Set d = -gproj + beta*d
   ========================================================================= */
void cg_update_d
(
    CGFLOAT      *d,
    CGFLOAT  *gproj,
    CGFLOAT    beta,
    CGINT         n  /* length of vectors */
)
{
    CGINT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++)
    {
        d [i] = -gproj [i] + beta*d [i] ;
    }
    for (; i < n; )
    {
        d [i] = -gproj [i] + beta*d [i] ;
        i++ ;

        d [i] = -gproj [i] + beta*d [i] ;
        i++ ;

        d [i] = -gproj [i] + beta*d [i] ;
        i++ ;

        d [i] = -gproj [i] + beta*d [i] ;
        i++ ;

        d [i] = -gproj [i] + beta*d [i] ;
        i++ ;
    }
}

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Start of limited memory CG routines  (+ matvec above)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* =========================================================================
   ==== cg_trisolve ========================================================
   =========================================================================
   Solve Rx = y or R'x = y where R is a dense upper triangular matrix
   ========================================================================= */
void cg_trisolve
(
    CGFLOAT *x, /* right side on input, solution on output */
    CGFLOAT *R, /* dense matrix */
    int      m, /* leading dimension of R */
    int      n, /* dimension of triangular system */
    int      w  /* T => Rx = y, F => R'x = y */
)
{
    int i, l ;
    if ( w )
    {
        l = m*n ;
        for (i = n; i > 0; )
        {
            i-- ;
            l -= (m-i) ;
            x [i] /= R [l] ;
            l -= i ;
            cg_daxpy0 (x, R+l, -x [i], i) ;
        }
    }
    else
    {
        l = 0 ;
        for (i = 0; i < n; i++)
        {
            x [i] = (x [i] - cg_dot0 (x, R+l, i))/R [l+i] ;
            l += m ;
        }
    }

/* equivalent to:
    BLAS_INT M, N ;
    M = (BLAS_INT) m ;
    N = (BLAS_INT) n ;
    if ( w ) SOPT_DTRSV ("u", "n", "n", &N, R, &M, x, blas_one) ;
    else     SOPT_DTRSV ("u", "t", "n", &N, R, &M, x, blas_one) ; */

    return ;
}

/* =========================================================================
   ==== cg_Yk ==============================================================
   =========================================================================
   Compute y = gnew - gold, set gold = gnew, compute y'y
   ========================================================================= */
void cg_Yk
(
    CGFLOAT    *y, /*output vector */
    CGFLOAT *gold, /* initial vector */
    CGFLOAT *gnew, /* search direction */
    CGFLOAT  *yty, /* y'y */
    CGINT       n  /* length of the vectors */
)
{
    CGINT n5, i ;
    CGFLOAT s, t ;
    n5 = n % 5 ;
    if ( (y != NULL) && (yty == NULL) )
    {
        for (i = 0; i < n5; i++)
        {
            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
        }
        for (; i < n; )
        {
            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;

            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;

            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;

            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;

            y [i] = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            i++ ;
        }
    }
    else if ( (y == NULL) && (yty != NULL) )
    {
        s = CGZERO ;
        for (i = 0; i < n5; i++)
        {
            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
        }
        for (; i < n; )
        {
            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            s += t*t ;
            i++ ;
        }
        *yty = s ;
    }
    else
    {
        s = CGZERO ;
        for (i = 0; i < n5; i++)
        {
            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
        }
        for (; i < n; )
        {
            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;

            t = gnew [i] - gold [i] ;
            gold [i] = gnew [i] ;
            y [i] = t ;
            s += t*t ;
            i++ ;
        }
        *yty = s ;
    }
    return ;
}

/* =========================================================================
   ==== cg_update_inf2 =====================================================
   =========================================================================
   Set gold = gnew, compute inf-norm of gnew & 2-norm of gnew, set d = -gnew
   ========================================================================= */
CGFLOAT cg_update_inf2
(
    CGFLOAT   *gold, /* old g */
    CGFLOAT   *gnew, /* new g */
    CGFLOAT      *d, /* d */
    CGFLOAT *gnorm2, /* 2-norm of g */
    CGINT         n  /* length of vectors */
)
{
    CGINT i, n5 ;
    CGFLOAT gnorm, s, t ;
    gnorm = CGZERO ;
    s = CGZERO ;
    n5 = n % 5 ;

    for (i = 0; i < n5; i++)
    {
        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ; 
        gold [i] = t ;
        d [i] = -t ;
    }
    for (; i < n; )
    {
        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ; 
        gold [i] = t ;
        d [i] = -t ;
        i++ ;

        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ; 
        gold [i] = t ;
        d [i] = -t ;
        i++ ;

        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
        i++ ;

        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
        i++ ;

        t = gnew [i] ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
        s += t*t ;
        gold [i] = t ;
        d [i] = -t ;
        i++ ;
    }
    *gnorm2 = s ;
    return (gnorm) ;
}

#if 0
CGFLOAT cg_daxpy_dot
(
    CGFLOAT      *x, /* new x */
    CGFLOAT      *d, /* direction */
    CGFLOAT const s, /* scalar */
    CGINT         n  /* length of vectors */
)
{
    CGINT i, n5 ;
    CGFLOAT t ;
    n5 = n % 5 ;

    t = CGZERO ;
    for (i = 0; i < n5; i++)
    {
        CGFLOAT const di = d [i] ;
        x [i] += s*di ;
        t += di*di ;
    }
    for (; i < n; )
    {
        CGFLOAT const d1 = d [i] ;
        x [i] += s*d1 ;
        t += d1*d1 ;
        i++ ;

        CGFLOAT const d2 = d [i] ;
        x [i] += s*d2 ;
        t += d2*d2 ;
        i++ ;

        CGFLOAT const d3 = d [i] ;
        x [i] += s*d3 ;
        t += d3*d3 ;
        i++ ;

        CGFLOAT const d4 = d [i] ;
        x [i] += s*d4 ;
        t += d4*d4 ;
        i++ ;

        CGFLOAT const d5 = d [i] ;
        x [i] += s*d5 ;
        t += d5*d5 ;
        i++ ;
    }
    return (t) ;
}
#endif

/* The rationale for the speed computations is as follows:

   Assume that the error decays like c*a^k where k is the
   iteration number.  Two different methods have different
   a and c, however, when deciding if method 1 is faster
   than method 2, the "c" does not matter since we consider
   the two methods starting from the same error when k = 0.
   That is c = err and k = 0 when the methods start.

   The formula for a is

   a = [error decay factor over the time window]^(1/ #its)

   where the error decay factor is error at the end of
   the time window divided by the error at the start of
   the time window.

   If a_i is the error decay factor for method i, then the
   number of iterations of method 2 needed to achieve the
   same error reduction as one iteration of method 1 is the
   solution k to a2^k = a1. That is, k2 = log(a1)/log(a2).
   The time to perform these iterations is

       (time per iteration for method 2)*k2

    After rearrangement, method 1 is faster than method 2 if

       (time per iteration for method 1)/ log (a1) >
       (time per iteration for method 2)/ log (a2)

    Note the this formula simplifies since the number of
    iterations in a time window appears in both the numerator
    and the denominator so they cancel. In summary, method 1
    is faster than method 2 if
Thank you for reviewing the original submission of the paper below. We have now received the revised manuscript. I am hoping you can take a look at the revision and see if your previous comments were addressed. Thank you for your help.

Best regards and best wishes for 2023,

William W. Hager
       (time of method 1)/log (decay factor for 1) >
       (time of method 2)/log (decay factor for 2) */

void cg_update_speed
(
    CGFLOAT     *speed, /* speed of the 3 methods */
    int         Method, /* PureCG, NewtonCG, NewtonSS */
    CGFLOAT     minerr, /* min error */
    CGFLOAT     maxerr, /* max error */
    CGFLOAT        tic, /* start time */
    CGFLOAT        toc  /* end   time */
)
{
    if ( minerr == maxerr ) return ;
    speed [Method] = (toc - tic)/log (minerr/maxerr) ;
    return ;
}



/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       End of limited memory CG utility routines
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
