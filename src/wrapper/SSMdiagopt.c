/* Contents:

    1. SSMdiagopt      -  main routine to optimize diagonal QP with ||x|| = r.
    2. SSMdiagF        -  evaluate F (see below)
    3. SSMdiagF1       -  evaluate F'
   ==========================================================================
   === SSMdiagopt ===========================================================
   ==========================================================================
    Find a global minimizer for the diagonalized QP

    (DQP)     minimize sum d_i x_i^2 + 2 f_i x_i
              subject to ||x|| <= r or ||x|| = r

    If the constraint is ||x|| <= r, then we first check to see if the
    constraint is active. The constraint is inactive if:

        1. d_i >= 0 for all i
        2. f_i =  0 if d_i = 0
        3. sum { (f_i/d_i)^2: d_i > 0 } <= r

    If the constraint is inactive, we store the solution and exit.
    Otherwise, we search for values of mu for which the KKT conditions hold:

        (d_i+mu)x_i + f_i = 0  or x_i = -f_i/(d_i+mu)

    and ||x|| = r. In other words,
                          _        _
                         |  f_i^2   |
    (*)     F (mu) = sum | ---------| - r^2 = 0 .
                         |(d_i+mu)^2|
                          -        -
    Note that mu >= -dmin = - min (d).  We define d_shift = d - dmin
    and we replace d by d_shift. Instead of mu >= -dmin, the
    constraint is mu >= 0.  We start with upper and
    lower bounds for the root of (*) and perform a Newton step
    on the left side and a secant step on the right side. The
    Newton step should produce an iterate on the left of the root
    and the secant iterate should yield a point on the right of
    the root. The initial bounds on the root are obtained as follows:
    Neglecting those i for which d_shift_i > 0 yields

              mu >= sqrt (sum (f_i^2: d_shift_i = 0)) / r

    Replacing all d_shift_i by 0 gives

              mu <= ||f||/ r

    The algorithm for computing the root mu is as follows in MATLAB notation:

    dmin = min (d) ;
    dshift = d - dmin ;            % dshift >= 0
    dzero = find (dshift <= 0) ;   % dzero = list of all (now zero) entries
    dpositive = find (d >  0) ;    % dpositive = list of all pos. entries

    x = zeros (n, 1) ;
    if ( isempty(dpositive) )      % all elements of d are identical
        normf = norm(f) ;
        if ( norm(f) > 0 )         % f is nonzero and d is a multiple of 1
            x = -(r/normf)*f ;
        else
            x(1) = r ;             % all elements of f are zero, x can be
        end                        % any norm r vector
        A = (normf/r) - emin ;
        return
    end

    x (dpositive) = -f (dpositive) ./ e (dpositive) ;
    normx = norm (x) ;
    f0 = norm (f (dzero)) ;        % part of f corresponding to dshift = 0
    f2 = f.^2 ;                    % square of f

    if ( (f0 == 0) & (normx <= r) )% mu = -dmin satisfies KKT conditions
                                   % x for dpositive as above, rest of x
                                   % chosen to satisfy norm constraint
        x (dzero (1)) = sqrt (r^2 - normx^2) ;
        mu = -dmin ;

    else                           % mu > -dmin

        A = f0 / r ;               % A <= root (root lower bound) corresponds
                                   % to discarding dpositive terms in (*)
        if ( A == 0 ) % shift A positive to avoid pole in (*) -- see code
        Fl = sum ((f./(A+d)).^2) - r^2 ; %Fl = F (A)

        if ( Fl < 0 )              % F not positive => mu = A

            x = -f./(d+A) ;
            mu = A - dmin ;

        else
            B = norm (f) / r ;    % right of root (upper bound) corresponds
                                   % replacing dshift by zero
            Fr = sum ((f./(B+d)).^2) - r^2 ; % F2 = F (B)
            while (abs ((B - A) / A) > tol)

                st = .5*Fl / sum (f2./((A+d).^3));% Newton step
                A = A + st ;                      % left side of root
                Fl = sum ((f./(A+d)).^2) - r^2 ;
                B = A - Fl*(B-A)/(Fr-Fl) ;     % secant step right side
                Fr = sum ((f./(B+d)).^2) - r^2 ;
                if ( Fl < 0 )
                    break       % impossible, mu = A
                elseif ( Fr > 0 )
                    mu = B ;    % impossible, mu = B
                    break
                end
            end
            x = -f ./ (A + d) ;
            mu = A - emin ;
        end
    end
   ========================================================================== */
#include "SSM.h"

void SSMdiagopt
(
    SSMFLOAT     *x,  /* n-by-1 solution vector (output) */
    SSMFLOAT      r,  /* radius of sphere */
    int BndLessThan,  /* TRUE means ||x|| <= r, FALSE means ||x|| = r */
    SSMDiagOpt  *DO,  /* diagonal optimization structure */
    SSMcom     *Com
)
{
    SSMINT i, imin, n ;
    SSMFLOAT dmin, dmax, fmax, f0, Fl, Fr, Ft, A, B, Aold, Flold, normf, normx2,
          rr, s, t, fi, tol, width, deps, feps, *d, *f, *f2, *dshift ;

    rr = r*r ;                        /* rr = radius squared */
    n = DO->n ;
    d = DO->d ;
    /* find dmin */
    dmin = SSMINF ;
    dmax = SSMZERO ;
    for (i = 0; i < n; i++)
    {
        t = d [i] ;
        if ( t < dmin )
        {
            dmin = t ;
            imin = i ;
        }
        if ( fabs (t) > dmax ) dmax = fabs (t) ;
    }
    DO->dmin = dmin ;
    DO->imin = imin ;

    fmax = SSMZERO ;
    f = DO->f ;
    for (i = 0; i < n; i++) if ( fabs (f [i]) > fmax ) fmax = fabs (f [i]) ;

    /* Treat components of d near dmin as equal to dmin, the corresponding
       components of f that are essentially zero are set to zero.
       Also, compute dshift = d - dmin, f2 = f^2, normf = ||f||,
       f0 = ||{f [i]: d [i] = dmin}|| */
    t = Com->Parm->diag_eps*n ;
    deps = dmax*t ;
    feps = fmax*t ;

    /* if problem is essentially positive semidefinite,
       check if optimal solution satisfies ||x|| <= r */
    if ( (dmin >= -deps) && BndLessThan )
    {
        s = SSMZERO ;
        for (i = 0; i < n; i++)
        {
            if ( fabs (d [i]) > deps )
            {
                t = -f [i]/d [i] ;
                s += t*t ;
                x [i] = t ;
            }
            else
            {
                /* if d_i essentially vanishes but not f_i, treat as ||x|| = r*/
                if ( fabs (f [i]) > feps )
                {
                    if ( d [i] > SSMZERO )
                    {
                        t = -f [i]/d [i] ;
                        s += t*t ;
                        x [i] = t ;
                    }
                }
                /* if d_i and f_i both vanish to within tolerance, set x_i = 0*/
                else x [i] = SSMZERO ; 
            }
        }
        if ( (s <= rr) && (i == n) )
        {
            DO->mu = SSMZERO ;
            return ;
        }
    }

    /* ||x|| = r */
    normx2 = SSMZERO ;
    f0 = SSMZERO ;
    normf = SSMZERO ;
    f2 = DO->f2 ;
    dshift = DO->dshift ;
    for (i = 0; i < n; i++)
    {
        fi = f [i] ;
        s = fi*fi ;
        f2 [i] = s ;
        normf += s ;
        t = d [i] - dmin ;
        dshift [i] = t ;
        if ( t == SSMZERO ) f0 += s ;
        else { t = fi/t ; normx2 += t*t ; }
    }
    f0 = sqrt (f0) ;        /* norm of f for indices with dshift_i = 0*/
    normf = sqrt (normf) ;  /* ||f|| */

    /* trivial case, x [0] = +- r */
    if ( n == 1 )
    {
        if ( f [0] > SSMZERO ) x [0] = -r ;
        else                   x [0] =  r ;
        DO->mu = -(d [0] + f [0]/x [0]) ;
        return ;
    }

    /* Compute the global minimizer */

    tol = DO->tol ;       /* relative accuracy of solution */

    /*------------------------------------------------------------------------*/
    /* solve the problem when dshift is nonzero */
    /*------------------------------------------------------------------------*/

    B = normf/r ;                     /* B is upper bound for mu */
    /*PRINTF (("f0: %e B: %e r: %e\n", f0, B, r) ; */
    if ( f0 == SSMZERO )              /* search for positive lower bound A */
    {
        i = 0 ;
        if ( normx2 > rr )            /* mu > 0 */
        {
            A = B ;                   /* start right of root */
            for (i = 0; i < 20; i++)
            {
                A = A*Com->Parm->shrink ;/* multiply upper bound by shrink */
                Fl = SSMdiagF (A, f2, dshift, rr, n) ;
                if ( Fl > SSMZERO ) break ; /* F (A) > 0 => A is left of root */
                B = A ;
            }
        }
        if ( (i == 20) || (normx2 <= rr) ) /* mu = 0 */
        {
            for (i = 0; i < n; i++)
            {
                if ( dshift [i] > SSMZERO )
                {
                    x [i] = -f [i]/dshift [i] ;
                }
                else
                {
                    x [i] = SSMZERO ;
                }
            }
            x [imin] = sqrt (rr - normx2) ;
            DO->mu = -dmin ;
            return ;
        }
    }
    else                       /* f0 > 0 */
    {
        A = f0/r ;             /* A stores lower bound for mu */
        Fl = SSMdiagF (A, f2, dshift, rr, n) ;
#ifndef NDEBUG
        if ( Fl < -1.e-10 )
        {
            DPRINT (("Fl: %e < 0 in SSMdiagopt when it should be "
                    "positive at location 1\n", Fl)) ;
        }
#endif
        if ( Fl < SSMZERO ) goto Exit ; /* Fl should be > 0, Fl < 0 => mu = A */
    }

#ifndef NDEBUG
    if ( A <= 0 )
    {
        DPRINT (("A = %e has wrong sign in SSMdiagopt\n", A)) ;
        SOPTERROR ("stop") ;
    }
    if ( A > B )
    {
        DPRINT (("A = %e > B = %e in SSMdiagopt\n", A, B)) ;
        SOPTERROR ("stop") ;
    }
#endif

    /*------------------------------------------------------------------------*/
    /* the root is bracketed, 0 < A < mu < B */
    /*------------------------------------------------------------------------*/
    Fr = SSMdiagF (B, f2, dshift, rr, n) ;
    if ( Fr > SSMZERO )        /* Fr should be < 0, Fr > 0 => mu = B */
    {
        A = B ;
        goto Exit ;
    }
#ifndef NDEBUG
    if ( Fr > tol )
    {
        DPRINT (("Fr = %e has wrong sign in SSMdiagopt, location 2\n", Fr)) ;
        SOPTERROR ("stop") ;
    }
#endif

    /*------------------------------------------------------------------------*/

    width = ((SSMFLOAT) 2)*fabs (B-A) ;
/*  printf ("A: %30.15e B: %30.15e\n", A, B) ;*/
    while ( SSMTRUE )
    {
        Aold = A ;
        A -= Fl/SSMdiagF1 (A, f2, dshift, n) ;   /* Newton step */
        Flold = Fl ;
        Fl = SSMdiagF (A, f2, dshift, rr, n) ;   /* function value at A */
/*      printf ("A: %30.15e B: %30.15e\n", A, B) ; */
/*      printf ("Fl: %30.15e Fr: %30.15e\n", Fl, Fr) ; */
        if ( ((B-A)/A <= tol) || (B-A <= tol*tol) ) break ;
        if ( Fl < SSMZERO ) /* with perfect precision, this would not happen */
        {
            B = A ;
            Fr = Fl ;
            A = Aold ;
            Fl = Flold ;
        }
#ifndef NDEBUG
        if ( A > B )
        {
            DPRINT (("Newton step = %e > B = %e in SSMdiagopt\n", A, B)) ;
            SOPTERROR ("stop") ;
        }
/*      if ( Fr == Fl )
        {
            DPRINT (("Fr = Fl = %e in SSMdiagopt\n", Fl) ;
            SOPTERROR ("stop") ;
        }*/
        if ( Fl < -.1 )
        {
            DPRINT (("Fl = %e has wrong sign in SSMdiagopt\n", Fl)) ;
            SOPTERROR ("stop") ;
        }
#endif
        if ( Fr == Fl ) break ;
        B = A - Fl*(B-A)/(Fr-Fl) ;             /* secant step */
        Fr = SSMdiagF (B, f2, dshift, rr, n) ; /* function value at B */
#ifndef NDEBUG
        if ( Fr >= .1 )
        {
            DPRINT (("Fr = %e has wrong sign in SSMdiagopt, B = %e\n", Fr, B)) ;
        }
#endif
        if ( (Fr >= SSMZERO) || ((B-A)/A <= tol) )
        {
            A = B ;                            /* B = best estimate of mu */
            break ;
        }
        if ( B - A > width ) /* bisection step when slow decay of width */
        {
            t = SSMHALF*(A+B) ;
            if ( (t <= A) || (t >= B) ) break ;
            Ft = SSMdiagF (t, f2, dshift, rr, n) ;/* function value at t */
            if ( Ft > SSMZERO )
            {
                A = t ;
                Fl = Ft ;
            }
            else
            {
                B = t ;
                Fr = Ft ;
            }
        }
        width *= SSMHALF ;
    }

    /* best estimate for mu is stored in A */
    Exit:
    for (i = 0; i < n; i++) x [i] = -f [i]/(A + dshift [i]) ;
    DO->mu = A - dmin ;

#ifndef NDEBUG
    A = DO->mu ;
    t = SSMZERO ;
    for (i = 0; i < n; i++) t += x [i]*x [i] ;
    t = SSMMAX (sqrt (t)-r, SSMZERO) ;
    if ( t/r > Com->Parm->check_kkt )
    {
        DPRINT (("diagopt solution violates norm constraint: %e > r = %e\n",
                  t, r)) ;
        SOPTERROR ("stop") ;
    }
    s = t = SSMZERO ;
    for (i = 0; i < n; i++)
    {
        s = SSMMAX (s, fabs (d [i])) ;
        t = SSMMAX (t, fabs ((d [i]+A)*x [i] + f [i])) ;
    }
    s = SSMMAX (s, fabs (A)) ;
    if ( s > SSMZERO )
    {
        if ( t > s*r*Com->Parm->check_kkt )
        {
            DPRINT (("global solution does not satisfy KKT, err = %e\n", t)) ;
            SOPTERROR ("stop") ;
        }
    }
#endif
}

/* ==========================================================================
   Evaluate the function  (*):
                          _        _
                         |  f_i^2   |
            F (mu) = sum |----------| - r^2
                         |(d_i+mu)^2|
                          -        -
   ========================================================================== */
SSMFLOAT SSMdiagF
(
    SSMFLOAT    mu,    /* the multiplier */
    SSMFLOAT   *f2,    /* f_i^2 */
    SSMFLOAT    *d,
    SSMFLOAT    rr,    /* radius of sphere squared */
    SSMINT       n     /* dimension */
)
{
    SSMINT i ;
    SSMFLOAT F, t ;
    F = SSMZERO ;
    for (i = 0; i < n; i++)
    {
        t = mu + d [i] ;
        F += f2 [i]/(t*t) ;
    }
    F -= rr ;
    return (F) ;
}

/* ==========================================================================
   Evaluate the function F':
                               _        _
                              |   f_i^2  |
            F' (mu) = - 2 sum | ---------| .
                              |(d_i+mu)^3|
                               -        -
   ========================================================================== */
SSMFLOAT SSMdiagF1
(
    SSMFLOAT   mu,    /* the multiplier */
    SSMFLOAT  *f2,    /* f_i^2 */
    SSMFLOAT   *d,
    SSMINT      n     /* dimension */
)
{
    SSMINT i ;
    SSMFLOAT F, t ;
    F = 0 ;
    for (i = 0; i < n; i++)
    {
        t = mu + d [i] ;
        F += f2 [i]/(t*t*t) ;
    }
    F = -2.*F ;
    return (F) ;
}
