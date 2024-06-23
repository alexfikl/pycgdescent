/* Contents:

   1. SSMtridiag    - reduce A to tridiagonal matrix by orthog. similarity
   2. SSMtriHouse   - full Householder reduction to tridiagonal form
   3. SSMtriLanczos - partial Lanczos tridiagonalization

   ==========================================================================
   === SSMtridiag ===========================================================
   ==========================================================================

    Given a sparse symmetric matrix A, contruct a matrix V with
    orthonormal columns such that V'AV is tridiagonal:
                           _
                          | d_1  u_1   0    0    .  .  .
                          |
            -V'AV = T =   | u_1  d_2   u_2   0   .  .  .
                          |
                          |  0   u_2   d_3   0   .  .  .
                          |
                          |  .    .     .    .   .  .  .

    The "-" sign is needed since the matrix is -(A+D).
    If the dimension of A is <= House_limit, then the Householder
    process is used to obtain a square orthogonal matrix V. Otherwise,
    the Lanczos process is employed until either a small u_i is
    encountered or max_its iterations are performed, which ever
    occurs first. The matrix P is an orthogonal projection into the space
    orthogonal to a given vector "a". If "a" is NULL, then
    we assume that it is the vector of all ones. If the numerical
    entries of A are null, then we assume that all the elements are one.
    In the Lanczos process, the starting column of V is a vector of
    the form Px where x is chosen randomly.  This choice ensures that
    the columns of V are all contained in the range of P.
   ========================================================================== */
#include "SSM.h"

int SSMtridiag /* return 0 (process completes)
                               -7 (starting Lanczos vector vanishes) */
(
    SSMFLOAT    *x, /* current solution estimate, ignored in Householder */
    int         it, /* iteration number */
    SSMLanczos *LH,
    SSMProblem *PB,
    SSMcom    *Com
)
{
    SSMINT j, k, l, n, n2, *Ai, *Ap, *Ap1 ;
    int status ;
    SSMFLOAT t, nnz, rand_max, *b, *Ax, *V, *D, *start_vector ;
    SSMParm *Parm ;

    Parm = Com->Parm ;
    n = PB->n ;
    b = PB->b ;
    Ax = PB->x ;
    Ai = PB->i ;
    Ap = PB->p ;
    D  = PB->D ;
    Ap1= Ap+1 ;
    V = LH->V ;
    /* use Householder */
    nnz = (SSMFLOAT) Ap [n] ;
    t = (SSMFLOAT) n ;
    if ( (n <= LH->House_limit)
         || ((t*t*t < 20*nnz*PB->DT->nalloc) && (PB->DT->nalloc == n)) )
    {
        /* create dense matrix */
        n2 = n*n ;
        for (k = 0; k < n2; k++) V [k] = SSMZERO ;
        k = 0 ;
        for (j = 0; j < n; j++)
        {
             l = Ap1 [j] ;
             for (; k < l; k++) V [Ai [k]] = Ax [k] ;
             V [j] = D [j] ;
             V = V+n ;
        }

        SSMtriHouse (LH->V, LH->d, LH->u, Com->wx1, n) ;
        LH->ncols = n ;
        status = 0 ;
        /* set bound on number of qr its */
    }
    else
    {
        /* Compute Lanczos starting point.  For starting iteration,
           starting guess is
           1. user guess is x if x not NULL
           2. b if b is nonzero
           3. random otherwise */
        start_vector = LH->V ;
        if ( it == 0 )
        {
            if ( x == NULL )
            {
                for (j = 0; j < n; j++)
                {
                    if ( b [j] != SSMZERO ) break ;
                }
                if ( j < n ) start_vector = b ;
                else
                {
                    rand_max = (SSMFLOAT) RAND_MAX ;
                    for (j = 0; j < n; j++)
                    {
                        start_vector [j] =
                                     (((SSMFLOAT) rand ())/rand_max)-SSMHALF ;
                    }
                }
            }
            else start_vector = x ;
        }
        /* after the first iteration, start_vector is the projected
           KKT residual which is stored in Com->MINRES+(8*n) */
        else start_vector = Com->MINRES+(8*n) ;
        status = SSMtriLanczos (start_vector, LH, PB, Com) ;
        if ( Parm->PrintLevel >= 1 )
        {
            printf ("number Lanczos: %i upper limit: %i\n",
                (int) LH->ncols, (int) LH->max_its) ;
        }
    }
    return (status) ;
}

/* ==========================================================================
   === SSMtriHouse ==========================================================
   ==========================================================================
  
    Reduce A to tridiagonal form by Householder orthogonal similarity
    transformations:
                                 _
                                | d_1  u_1   0    0    .  .  .
                                |
             V'AV = T =         | u_1  d_2   u_2   0   .  .  .
                                |
                                |  0   u_2   d_3   0   .  .  .
                                |
                                |  .    .     .    .   .  .  .

   where V is computed as a product of Householder matrices. A is
   stored as a dense symmetric matrix, and V is stored as a dense matrix.
   The output matrix V overwrites the input matrix A.

n = size (A, 1) ;
d = zeros(n,1) ;
u = zeros(n,1) ;

if ( n == 1 ) % {

    d(1) = A(1,1) ;
    V(1) = 1 ;

elseif ( n == 2 ) % } {

    d(1) = A(1,1) ;
    d(2) = A(2,2) ;
    u(1) = A(2,1) ;
    V = zeros(2) ;
    V(1,1) = 1 ;
    V(2,2) = 2 ;

else % } {

    W = zeros(n,n) ; Householder vectors (stored in A)
    for j = 1:n-2 % {

        %-----------------------------------------------------------------------
        % compute v
        %-----------------------------------------------------------------------

        % accesses only the strictly lower triangular part of A:
        v = house (A(:,j),j+1) ;
        % v (1:j) is zero

        %-----------------------------------------------------------------------
        % compute column j+1 of W
        %-----------------------------------------------------------------------

        % column 1 of W is zero.  W is strictly lower triangular
        W(:,j+1) = v ;

        %-----------------------------------------------------------------------
        % compute x
        %-----------------------------------------------------------------------

        % since v (j) is zero, we can skip column j of A
        % also note that A (1:(j-1),(j+1):n) is zero
        % accesses just lower part of A to compute this
        x = zeros (n,1) ;
        x (j:n) = A (j:n,(j+1):n) * v ((j+1):n) ;

        %-----------------------------------------------------------------------
        % compute a
        %-----------------------------------------------------------------------

        % x (1:(j-1)) is zero and v (1:j) is zero
        a = .5 * ( (v ((j+1):n))' * x ((j+1):n)) ;

        %-----------------------------------------------------------------------
        % compute b
        %-----------------------------------------------------------------------

        b = zeros (n,1) ;
        b (j:n) = a * v (j:n) - x (j:n) ;

        %-----------------------------------------------------------------------
        % update A
        %-----------------------------------------------------------------------

        % just update the lower triangular part of A
        A(j:n,j:n) = A(j:n,j:n) + v(j:n)*b(j:n)' + b(j:n)*v(j:n)' ;

        %-----------------------------------------------------------------------
        % save diagonal and off-diagonal
        %-----------------------------------------------------------------------

        d(j) = A(j,j) ;
        u(j) = A(j+1,j) ;

    end % }
    d(n-1) = A(n-1,n-1) ;
    u(n-1) = A(n,n-1) ;
    d(n) = A(n,n) ;

    %---------------------------------------------------------------------------
    % project W to V
    %---------------------------------------------------------------------------

    V = zeros(n,n) ;
    for j = 1:n % {
        V(:,j) = qj (W, j) ;
    end % }

end % }

   ========================================================================== */

void SSMtriHouse
(
    SSMFLOAT *A, /* n-by-n matrix dense symmetric matrix (input)
                               dense orthogonal matrix (output) */
    SSMFLOAT *d, /* n-by-1 vector (output) */
    SSMFLOAT *u, /* n-by-1 vector (output) */
    SSMFLOAT *x, /* n-by-1 vector (workspace) */
    SSMINT    n
)
{
    SSMFLOAT a, hj, y, s, t, *Ak, *Aj ;
    SSMINT i, j, jp1, k ;

    if (n == 1)
    {

        d [0] = A [0] ;
        u [0] = 0 ;
        A [0] = 1 ;
        return ;

    }
    else if (n == 2)
    {

        d [0] = A [0] ;
        d [1] = A [3] ;
        u [0] = A [1] ;
        A [0] = 1 ;
        A [1] = 0 ;
        A [2] = 0 ;
        A [3] = 1 ;
        return ;

    }

     /* start Householder reduction for matrix of dimension >= 3 */
     /* initialize current column of A */
    Aj = A ;
    for (j = 0 ; j < n-2 ; j++)
    {
        d [j] = Aj [j] ;
        s = SSMZERO ;
        jp1 = j+1 ;
        for (i = jp1; i < n; i++) s += Aj [i]*Aj [i] ;
        if ( s == SSMZERO )
        {
            u [j] = SSMZERO ;
            Aj += n ;
            continue ;   /* the Householder vector is 0 as is subcolumn of A */
        }

        /* scaling factor for Household matrix */
        s = sqrt (s) ;
        hj = Aj [jp1] ;
        t = 1/sqrt (s*(s + fabs (hj))) ;
        if ( hj >= 0 )
        {
            hj += s ;    /* hj is (j+1)st element of Household matrix */
            y = -s ;     /* y = u_j = (j+1)st element of (Householder * Aj) */
        }
        else
        {
            hj -= s ;
            y = s ;
        }
        u [j] = y ;

        /* ----------------------------------------------------------- */
        /* overwrite column j of A with Householder vector */
        /* ----------------------------------------------------------- */

        Aj [jp1] = hj ;
        for (i = jp1; i < n; i++) Aj [i] *= t ;
    
        /* x ((j+1):n) = A ((j+1):n,(j+1):n) * house ((j+1):n) */
        for (i = jp1; i < n; i++) x [i] = SSMZERO ;
        Ak = Aj ;
        for (k = jp1 ; k < n ; k++)
        {
    	    /* Ak is column k of A */
            Ak += n ;
            /* diagonal */
            {
                SSMFLOAT vk = Aj [k] ;
                SSMFLOAT xk = Ak [k] * vk ;
    
                /* dot product with row k of A and saxpy with column k of A */
                for (i = k+1 ; i < n ; i++)
                {
                    SSMFLOAT aki = Ak [i] ;
                    xk += aki * Aj [i] ;
                    x [i] += aki * vk ;
                }
                x [k] += xk ;
            }
        }

        /* -------------------------------------------------------------- */
        /* a = .5 * ( (v ((j+1):n))' * x ((j+1):n)) ; */
        /* -------------------------------------------------------------- */

        a = 0 ;
        for (i = jp1; i < n; i++) a += Aj [i] * x [i] ;
        a /= 2 ;

        /* ----------------------------------------------------- */
        /* b (j:n) = a * v (j:n) - x (j:n) ; write b back into x */
        /* ----------------------------------------------------- */

        for (i = jp1; i < n; i++)
        {
            x [i] = a*Aj [i] - x [i] ;
        }

        /* ----------------------------------------------------------------- */
        /* A(j+1:n,j+1:n) = A(j+1:n,j+1:n)
                            + v(j+1:n)*b(j+1:n)' + b(j+1:n)*v(j+1:n)' */
        /* ----------------------------------------------------------------- */

        /* note that x holds b */
        /* update only the lower triangular part of A */
        Ak = Aj ;
        for (k = jp1 ; k < n ; k++)
        {
            Ak += n ;
            {
                SSMFLOAT vk = Aj [k] ;
                SSMFLOAT xk = x [k] ;
                Ak [k] += 2.*(vk * xk) ;
                for (i = k+1; i < n; i++)
                {
                    Ak [i] += (Aj [i] * xk) + (x [i] * vk) ;
                }
            }
        }
        Aj += n ;
    }

    /* ------------------------------------------------------------------ */
    /* d(n-1) = A(n-1,n-1) ; u(n-1) = A(n,n-1) ; u(n) = 0 ; d(n) = A(n,n) */
    /* ------------------------------------------------------------------ */

    d [n-2] = Aj [n-2] ;
    u [n-2] = Aj [n-1] ;
    Aj += n ;           /* Aj = last column of matrix */
    d [n-1] = Aj [n-1] ;
    u [n-1] = 0 ;

    /* ------------------------------------------------------------------ */
    /* compute columns of V and store them in A */
    /* start with last column and work to first column */
    /* ------------------------------------------------------------------ */

    /* last column is special */
    Ak = Aj-(n+n) ;    /* first column with a Householder vector */
    t = Ak [n-1] ;
    Aj [n-2] =   - t*Ak [n-2] ;
    Aj [n-1] = 1 - t*Ak [n-1] ;
    for (j = n-4; j >= 0; j--)
    {
        /* compute the new element in Aj */
        Ak -= n ;
        t = SSMZERO ;
        for (i = j+2; i < n; i++)
        {
            t += Aj [i]*Ak [i] ;
        }
        Aj [j+1] = -Ak [j+1]*t ;

        /* update the remaining elements in Aj */
        for (i = j+2; i < n; i++)
        {
            Aj [i] -= Ak [i]*t ;
        }
    }
    Aj [0] = 0 ;

    /* compute next to last through 2nd column of V */
    for (j = n-2; j > 0; j--)
    {
        Aj -= n ;
        /* startup, column of identity times preceding Householder */
        Ak = Aj - n ;
        t = Ak [j] ;
        for (i = j; i < n; i++) Aj [i] = -Ak [i]*t ;
        Aj [j]++ ;
        for (k = j-2; k >= 0; k--)
        {
            /* compute the new element in Aj */
            Ak -= n ;
            t = SSMZERO ;
            for (i = k+2; i < n; i++)
            {
                t += Aj [i]*Ak [i] ;
            }
            Aj [k+1] = -Ak [k+1]*t ;

            /* update the remaining elements in Aj */
            for (i = k+2; i < n; i++)
            {
                Aj [i] -= Ak [i]*t ;
            }
        }
        Aj [0] = SSMZERO ;
    }
    /* store the first columns of V = first column of identity */
    Aj -= n ;
    for (i = 0; i < n; i++) Aj [i] = 0 ;
    Aj [0] = 1 ;
}

/* ========================================================================== */
/* === SSMtriLanczos ======================================================== */
/* ========================================================================== 
    Partially reduce P(A+D)P to tridiagonal form by the Lanczos process.
    The maximum number of Lanczos iterations is Com->Lanczos_max_its.  The
    Lanczos process is terminated when there is a small
    off-diagonal element (|u_j| <= Amax*Parm->Lanczos_tol, where Amax is the
    absolute largest element in A + D).
                           _
                          | d_1  u_1   0    0    .  .  .
                          |
            -V'AV = T =   | u_1  d_2   u_2   0   .  .  .
                          |
                          |  0   u_2   d_3   0   .  .  .
                          |
                          |  .    .     .    .   .  .  .

   ========================================================================== */

int SSMtriLanczos /* return 0 (process completes)
                                  -7 (starting Lanczos vector vanishes) */
(
    SSMFLOAT *start_vector, /* starting point */
    SSMLanczos         *LH, /* Lanczos structure */
    SSMProblem         *PB, /* Problem specification */
    SSMcom            *Com  /* SSMcom structure */
)
{
    SSMINT *Ap ;
    SSMINT i, j, n ;
    SSMINT *Ai ;
    SSMFLOAT s, t, uj, *Ax, *d, *u, *r, *Vj, *Vjprior ;

    /* ---------------------------------------------------------------------- */
    /* Read in the needed arrays */
    /* ---------------------------------------------------------------------- */

    /* output */
    Vj= LH->V ;    /* dense SSMFLOAT matrix of orthonormal vectors */
    d = LH->d ;    /* SSMFLOAT array containing diagonal of tridiagonal matrix*/
    u = LH->u ;    /* SSMFLOAT array containing subdiagonal tridiagonal matrix*/

    /* problem specification */
    n  = PB->n ;  /* problem dimension */
    Ax = PB->x ;  /* numerical values for edge weights */
    Ai = PB->i ;  /* adjacent vertices for each node */
    Ap = PB->p ;  /* points into Ax or Ai */

    /* work arrays */
    r = Com->wx1 ;

/* Lanczos iteration (orthogonal vectors stored in q_1, q_2, ... )

    q_0 = 0, r = q_1 (starting vector)
    for j = 1:n
        u_{j-1} = ||r||                  (u_0 is discarded)
        q_j = r/u_{j-1}                  (normalize vector)
        d_j = q_j'Aq_j                   (diagonal element)
        r = (A-d_jI)q_j - u_{j-1}q_{j-1} (subtract projection on prior vectors)
    end
*/

    s = SSMZERO ;
    for (i = 0; i < n; i++)
    {
        t = start_vector [i] ;
        s += t*t ;
    }
    if ( s == SSMZERO )
    {
        if ( Com->Parm->PrintLevel > 0 )
        {
            printf ("starting vector in SSMLanczos vanishes!!\n") ;
            LH->ncols = 0 ;
            return (-7) ;
        }
    }
    s = sqrt (s) ;

    /* normalize starting vector */
    t = SSMONE/s ;
    for (i = 0; i < n; i++) Vj [i] = t*start_vector [i] ;

    SSMmult (r, Vj, Ax, PB->D, Ai, Ap, n, Com) ; /* startup: compute r = A*Vj */

    /* compute d [0] */
    s = SSMZERO ;
    for (i = 0; i < n; i++) s += r [i]*Vj [i] ;
    d [0] = s ;
    for (i = 0; i < n; i++) r [i] -= s* Vj [i] ; /* 1st iteration r computed */

    /* the main Lanczos iteration */
    j = 0 ;
    while ( 1 )
    {
        s = 0;
        for (i = 0; i < n; i++) s += r [i]*r [i] ;
        uj = sqrt (s) ;
        u [j] = uj ;
        if ( uj < LH->tol ) break ;
        Vjprior = Vj ;
        Vj += n ;
        for (i = 0; i < n; i++) Vj [i] = r[i]/uj ;

        SSMmult (r, Vj, Ax, PB->D, Ai, Ap, n, Com) ; /* compute q = Ar */

        /* compute d_j */
        s = SSMZERO ;
        for (i = 0; i < n; i++) s += r [i]*Vj [i] ;
        j++ ;
        d [j] = s ;
        if ( j+1 == LH->max_its ) break ; /* max number of iteration */

        /* compute new residual r */

        for (i = 0; i < n; i++) r [i] -= s*Vj [i] + uj*Vjprior [i] ;

    }
    u [j] = SSMZERO ;
    LH->ncols = j+1 ;                        /* number of columns in V */
    return (0) ;

/*  PRINTF (("Lanczos cols: %i\n", j+1)) ;
    SSMINT k ;
    Vj = V ;
    for (k = 0; k < j+1; k++)
    {
        PRINTF (("V (:, %i) = [\n", k+1)) ;
        for (i = 0; i < n; i++) PRINTF (("%25.15e\n", Vj [i])) ;
        PRINTF (("] ;\n")) ;
        Vj += n ;
    }
    PRINTF (("A = [\n"));
    for (k = 0; k < n; k++) for (i = Ap [k]; i < Ap [k+1]; i++)
        PRINTF (("%i %i 1\n", Ai [i]+1, k+1)) ; */

}
