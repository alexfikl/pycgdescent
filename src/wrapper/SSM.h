#ifndef _SSM_H_
#define _SSM_H_

#include "sopt.h"

/* ========================================================================== */
/* === macros =============================================================== */
/* ========================================================================== */

#define SSMFLOAT SOPTFLOAT
#define SSMINT   SOPTINT
#define SSMFALSE SuiteOPTfalse
#define SSMTRUE SuiteOPTtrue

#define SSMONE ((SSMFLOAT) 1)
#define SSMZERO ((SSMFLOAT) 0)
#define SSMHALF ((SSMFLOAT) .5)

#define SSMMIN SOPTMIN
#define SSMMAX SOPTMAX

#define SSMINF inf

#define ssm_free sopt_free
#define ssm_malloc sopt_malloc

#define SSM_STATUS_OK                                             (0)
#define SSM_TOO_MANY_ITERATIONS_IN_MINRES                       (401)
#define SSM_VIOLATES_SPHERE_CONSTRAINT                          (402)
#define SSM_RESTARTS_EXCEEDS_LIMIT                              (403)
#define SSM_DIMENSION_NONPOSITIVE                               (404)
#define SSM_QR_DIAGONALIZATION_FAILS                            (405)
#define SSM_NOT_ENOUGH_SPACE_FOR_QR_DIAGONALIZATION             (406)
#define SSM_STARTING_LANCZOS_VECTOR_VANISHES                    (407)
#define SSM_SOLUTION_IN_INTERIOR_OF_SPHERE                      (408)
#define SSM_SOLUTION_ON_BOUNDARY_OF_SPHERE                      (409)
#define SSM_OUT_OF_MEMORY                                       (410)
#define SSM_START_MESSAGES                                      (400)
#define SSM_END_MESSAGES                                        (499)

/* debug printf uses DPRINT */
#ifndef NDEBUG
#define DPRINT(e) printf e
#else
#define DPRINT(e)
#endif

/* -------------------------------------------------------------------------- */
/* SSM version information */
/* -------------------------------------------------------------------------- */
#define SSM_DATE "May 26, 2023"
#define SSM_MAIN_VERSION 2
#define SSM_SUB_VERSION 0
#define SSM_SUBSUB_VERSION 0

/* ========================================================================== */
/* === structures =========================================================== */
/* ========================================================================== */

typedef struct SSMParm      /* parameters */
{
    unsigned int seed ;     /* seed to start srand; do not use srand if zero */

    SSMFLOAT         eps ; /* machine epsilon for FLOAT  */
    int     Lanczos_dim1 ; /* Lanczos iterations L:           */
    int     Lanczos_dim2 ; /*   if      ( n <= dim1 ) L = n-1 */
    int       Lanczos_L1 ; /*   else if ( n <= dim2 ) L = L1  */
    int       Lanczos_L2 ; /*   else             L = L1 + L2*(log10 (n/dim2)) */
    int      Lanczos_bnd ; /* upper bound on number of Lanczos iterations */
    SSMFLOAT Lanczos_tol ; /* stop Lanczos when |u_j|/max |A| <= tol */
    int      House_limit ; /* use Householder when n <= limit */
    int           qr_its ; /* max number QR algorithm iterations*/
    SSMFLOAT  qr_tol_rel ; /* relative eigenvalue convergence tol*/
    SSMFLOAT  qr_tol_abs ; /* absolute eigenvalue convergence tol*/
    SSMFLOAT    qr_lower ; /* lower bound for qr scaling factor */
    SSMFLOAT diag_optim_tol ; /* accuracy of multiplier for diagonal matrix */
    SSMFLOAT    diag_eps ; /* bound on error in diagonal of tridiag matrix */
    SSMFLOAT   check_tol ; /* tolerance used when checking QP structure */
    SSMFLOAT   check_kkt ; /* tolerance used when checking KKT error */
    SSMFLOAT  eig_refine ; /* factor determines when to refine eig in SSM */
    SSMFLOAT radius_flex ; /* flex factor for radius in interior routines */
    SSMFLOAT   sqp_decay ; /* factor determines when to stop SQP iteration */
    SSMFLOAT   eig_decay ; /* factor determines when to stop inverse iteration*/
    SSMFLOAT   eig_lower ; /* lower bound smallest eigenvalue (-INF default) */
    int              IPM ; /* TRUE (eigenvalue estimate using inverse power
                              method), FALSE (... using SQP method) */
    SSMFLOAT   ssm_decay ; /* required error decay factor in SSM, else restart*/
    SSMFLOAT      shrink ; /* mu multiplied by shrink in diagopt */
    SSMFLOAT minres_its_fac ; /* max number minres iterations in SSM is fac*n */
    SSMFLOAT ssm_its_fac ; /* max number SSM iterations is fac*n */
    SSMFLOAT grow_Lanczos ; /* factor to grow number of Lanczos iterations */
    int         nrestart ; /* number of Lanczos restarts attempted */
    int       PrintLevel ; /* Level 0  = no print, ... , Level 2 = max print*/
    int       PrintParms ; /* TRUE means print parameter values */
    int       PrintFinal ; /* TRUE means print status and statistics */
    int      BndLessThan ; /* TRUE means ||x|| <= r, FALSE means ||x|| = r */
} SSMParm ;

typedef struct ssm_stat_struct /* statistics returned to user */
{
    SSMINT        mults ; /* number of matrix/vector multiplications */
    int         ssm_its ; /* number of SSM iterations */
    int      minres_its ; /* number of minimum residual iterations */
    int     restart_its ; /* number of times Lanczos restarted */
    SSMFLOAT      error ; /* 2-norm of KKT error */
} ssm_stat ;

typedef struct SSMDiag     /* diagonalization of tridiagonal matrix */
{
    SSMINT             n ; /* dimension of matrix */
    SSMFLOAT        tol1 ; /* stop when u_{k-1}^2 <= qr_tol1*|d_k| */
    SSMFLOAT        tol2 ; /* stop when |u_{k-1}| <= qr_tol2 */
    int          max_its ; /* max QR iterations (default: 4n) */
    int              its ; /* actual number of QR iterations */
    SSMINT        nscale ; /* number of scaling operations */
    SSMINT       ngivens ; /* number of Givens rotation */
    SSMFLOAT          *e ; /* eigenvalues of tridiagonal matrix, size nalloc */
    SSMINT        nalloc ; /* space allocated for diagonal, must be >= n */
    SSMFLOAT    qr_lower ; /* lower bound for scaling in qr method */
    SSMFLOAT         *gx ; /* first Givens factor, size=max_its*nalloc */
    SSMFLOAT         *gy ; /* second Givens factors, same size as gx */
    short            *gz ; /* 1 or 2, type of Givens transform, size of gx */
    SSMFLOAT         *gs ; /* fast Givens scale factors, same size as gx */
    SSMINT           *gi ; /* row indices associated with scale factors */
    SSMINT           *gj ; /* indexed by QR pass, points into Givens factors*/
    SSMINT           *gk ; /* index of last diagonal element in QR iteration */
    SSMINT           *gf ; /* index of first diagonal element in QR iteration */
} SSMDiag ;

typedef struct SSMLanczos  /* Lanczos/Householder tridiagonalization structure*/
{
    int          max_its ; /* upper bound on the number of Lanczos iterations */
    SSMFLOAT         tol ; /* termination tolerance for u_j in Lanczos iter */
    int      House_limit ; /* Householder tridiag if n <= House_limit */
    SSMFLOAT          *V ; /* orthonormal vectors, V'AV = T is tridiagonal */
    int            ncols ; /* number of cols in the Lanczos orthogonal matrix*/
    int            nrows ; /* number of rows in the Lanczos orthogonal matrix*/
    SSMFLOAT          *d ; /* diagonal of tridiagonal T */
    SSMFLOAT          *u ; /* upper diagonal of tridiagonal T */
} SSMLanczos ;

typedef struct SSMDiagOpt  /* optimization of diagonal quadratic */
{
    SSMINT             n ; /* dimension */
    SSMFLOAT          *d ; /* diagonal matrix */
    SSMFLOAT          *f ; /* linear term in cost function */
    SSMFLOAT         *f2 ; /* square of linear term in cost function */
    SSMFLOAT     *dshift ; /* d - min (d) */
    SSMFLOAT        dmin ; /* dmin = min (d) */
    SSMFLOAT        imin ; /* index of minimal component of d */
    SSMFLOAT         tol ; /* relative error tolerance for optimal multiplier */
    SSMFLOAT          mu ; /* multiplier associated with minimizer */
} SSMDiagOpt ;

typedef struct SSMProblem  /* problem specification */
{
    SSMINT             n ; /* problem dimension */
    SSMFLOAT          *x ; /* edge weights */
    SSMINT            *i ; /* adjacent vertices for each node */
    SSMINT            *u ; /* location last super diag in col (used in SSOR) */
    SSMINT            *p ; /* points into x and i, p [0], ..., p [n] */
    SSMFLOAT          *D ; /* diagonal of objective matrix */
    SSMFLOAT          *b ; /* linear term in objective function */
    SSMFLOAT         rad ; /* radius of sphere */
    SSMFLOAT        Amax ; /* absolute largest element in matrix */
    SSMFLOAT        Dmin ; /* smallest diagonal element */
    SSMDiag          *DT ; /* diagonalization structure for tridiagonal matrix*/
    SSMDiagOpt       *DO ; /* diagonal optimization structure */
    SSMLanczos       *LH ; /* Lanczos/Householder structure */
} SSMProblem ;

typedef struct SSMcom
{
    SSMINT          *wi1 ; /* work array of size n */
    SSMINT          *wi2 ; /* work array of size n */
    SSMFLOAT        *wx1 ; /* work array of size n, wx1 */
    SSMFLOAT        *wx2 ; /* work array of size n, wx2 = wx1+n */
    SSMFLOAT        *wx3 ; /* work array of size n, wx3 = wx2+2n */
    SSMFLOAT       *SSOR ; /* work array of size 4n for SSOR d, s, p, q */
    SSMFLOAT     *MINRES ; /* work array of size 12n needed for MINRES */
    SSMFLOAT          *W ; /* work array of size 5n for Householder vectors */
    SSMFLOAT          *V ; /* work array of size 5n for orthogonal vectors */
    SSMFLOAT        *VAV ; /* work array of size 25 for storing V'AV */
    SSMFLOAT          *v ; /* current eigenvector estimate */
    SSMFLOAT         *Av ; /* A times v */
    SSMFLOAT         *Ax ; /* A times x */
    SSMFLOAT         *r0 ; /* Ax + b */
    int           Active ; /* TRUE (||x|| = r), FALSE (||x|| < r) */
    SSMFLOAT        emin ; /* current minimum eigenvalue estimate */
    SSMFLOAT          mu ; /* current multiplier estimate */
    SSMFLOAT       normx ; /* norm of x */
    SSMFLOAT         tol ; /* KKT error tolerance for solution */
    SSMProblem       *PB ; /* problem specification */
    SSMProblem  *PBdense ; /* structure for dense subproblems */
    SSMParm        *Parm ; /* parameters */
    SSMINT         mults ; /* number of matrix/vector products */
    int          ssm_its ; /* number of SSM iterations */
    int        ssm_limit ; /* limit on number of SSM iterations */
    int       minres_its ; /* number of minimum residual iterations */
    int     minres_limit ; /* limit on minres iterations */
    SSMFLOAT       error ; /* 2-norm of KKT error */
    SSMFLOAT   eig_error ; /* relative 2-norm of eigenvector error */

/* check info */
    SSMFLOAT    cost_old ; /* cost previous time objective function checked */
    SSMFLOAT    emin_old ; /* prior eigenvalue estimate */
} SSMcom ;

/* prototypes */

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
    SSMFLOAT    *x, /* solution (output) */
    ssm_stat *Stat, /* NULL means do not return statistics */

/* input: */
    SSMINT        n, /* size of x, constraint lo <= a'x <= hi */
    SSMFLOAT    *Ax, /* numerical entries in matrix, packed array */
    SSMINT      *Ai, /* row indices for each column */
    SSMINT      *Ap, /* points into Ax and Ai, Ap [j] = start of column j */
    SSMFLOAT     *b, /* linear term in objective function */
    SSMFLOAT      r, /* radius of ball */
    SSMFLOAT    tol, /* KKT error tolerance for solution */
    SSMFLOAT *guess, /* starting, NULL means no guess given */
    SSMParm  *UParm/* user parameters, NULL means use default parameters */
) ;

void SSMdefault
(
    SSMParm  *Parm /* pointer to parameter structure */
) ;

void SSM_print_status
(
    int status  /* status from SSM */
) ;

void SSMprint_parms
(
    SSMParm *Parm /* SSMparm structure to be printed */
) ;

int SSMallocate
(
    SSMcom    *Com, /* pointer to SSMcom structure */
    SSMProblem *PB, /* problem specification */
    SSMParm  *Parm, /* parameters, needed to determine allocation */
    SSMINT       n  /* problem dimension */
) ;

int SSMreallocate
(
    SSMFLOAT    *x, /* current estimate of solution */
    SSMcom    *Com  /* SSMcom structure */
) ;

void SSMdestroy
(
    SSMcom *Com  /* SSMcom structure to free */
) ;

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
) ;

void SSMorth
(
    SSMFLOAT  *w, /* k-th orthonormal vector */
    SSMFLOAT  *W, /* packed matrix storing Householder vectors */
    SSMFLOAT  *x, /* k-th new vector */
    SSMINT     k,
    SSMINT     n  /* dimension of x */
) ;

int SSMdiag /* return: 0 if convergence tolerance satisfied
                              -5 (algorithm did not converge)
                              -6 (insufficient space allocated for Givens) */
(
    SSMFLOAT *din, /* diagonal of tridiagonal matrix */
    SSMFLOAT *uin, /* superdiagonal of tridiagonal matrix, use uin (0:n-2) */
    SSMINT      n, /* dimension of tridiagonal matrix */
    SSMDiag   *DT, /* structure for storing diagonalization */
    SSMINT   *wi1, /* work array of size n */
    SSMFLOAT *wx1, /* work array of size n */
    SSMFLOAT *wx2, /* work array of size n */
    SSMFLOAT *wx3  /* work array of size n */
) ;

int SSMballdense /*return 0 (error tolerance satisfied)
                                 -5 (failure of QR diagonalization)
                                 -6 (insufficient space in QR diagonalization)*/
(
    SSMFLOAT     *x, /* n-by-1 solution vector (output) */
    SSMProblem  *PB, /* problem structure associated with A */
    SSMcom     *Com /* SSMcom structure */
) ;

#ifndef NDEBUG
int SSMcheckKKT
(
    SSMFLOAT    *x,  /* n-by-1 solution vector (output) */
    SSMFLOAT     r,  /* radius of sphere */
    SSMProblem *PB,  /* problem specification */
    SSMcom    *Com   /* pointer to SSMcom structure */
) ;
#endif

SSMFLOAT SSMdiagF
(
    SSMFLOAT    mu,    /* the multiplier */
    SSMFLOAT   *f2,    /* f_i^2 */
    SSMFLOAT    *d,
    SSMFLOAT    rr,    /* radius of sphere squared */
    SSMINT       n     /* dimension */
) ;

SSMFLOAT SSMdiagF1
(
    SSMFLOAT    mu,    /* the multiplier */
    SSMFLOAT   *f2,    /* f_i^2 */
    SSMFLOAT    *d,
    SSMINT       n     /* dimension */
) ;

void SSMmineig
(
    SSMFLOAT *emin, /* estimated smallest eigenvalue */
    SSMFLOAT    *v, /* associated eigenvector */
    SSMFLOAT    *w, /* work array of size n */
    SSMProblem *PB  /* problem specification */
) ;

int SSMtridiag /* return 0 (process completes)
                               -7 (starting Lanczos vector vanishes) */
(
    SSMFLOAT    *x, /* current solution estimate, ignored in Householder */
    int         it, /* iteration number */
    SSMLanczos *LH,
    SSMProblem *PB,
    SSMcom    *Com
) ;

void SSMtriHouse
(
    SSMFLOAT *A,      /* n-by-n matrix dense symmetric matrix (input)
                                     dense orthogonal matrix (output) */
    SSMFLOAT *d,      /* n-by-1 vector (output) */
    SSMFLOAT *u,      /* n-by-1 vector (output) */
    SSMFLOAT *x,      /* n-by-1 vector (workspace) */
    SSMINT    n
) ;

int SSMtriLanczos /* return 0 (process completes)
                                  -7 (starting Lanczos vector vanishes) */
(
    SSMFLOAT    *x, /* starting point, NULL = random starting point */
    SSMLanczos *LH, /* Lanczos structure */
    SSMProblem *PB, /* Problem specification */
    SSMcom    *Com
) ;

void SSMdiagopt
(
    SSMFLOAT     *x,  /* n-by-1 solution vector (output) */
    SSMFLOAT      r,  /* radius of sphere */
    int BndLessThan,  /* TRUE means ||x|| <= r, FALSE means ||x|| = r */
    SSMDiagOpt  *DO,  /* diagonal optimization structure */
    SSMcom     *Com
) ;

void SSMmult
(
    SSMFLOAT  *p,    /* output vector of size n */
    SSMFLOAT  *x,    /* input vector of size n */
    SSMFLOAT *Ax,    /* numerical values in A excluding diagonal */
    SSMFLOAT  *D,    /* diagonal of matrix */
    SSMINT   *Ai,    /* row indices for each column of A */
    SSMINT   *Ap,    /* Ap [j] = start of column j */
    SSMINT     n,    /* dimension of matrix */
    SSMcom  *Com
) ;

void SSMGivensMult
(
    SSMFLOAT   *x,  /* the vector to which the rotations are applied */
    SSMDiag   *DT   /* diagonalization structure */
) ;

void SSMtGivensMult
(
    SSMFLOAT   *x,  /* the vector to which the rotations are applied */
    SSMDiag   *DT   /* diagonalization structure */
) ;

void SSMDenseMult
(
    SSMFLOAT *y,     /* m by 1 product, output */
    SSMFLOAT *x,     /* n by 1 given vector */
    SSMFLOAT *V,     /* dense m by n matrix */
    SSMINT    m,     /* number of rows */
    SSMINT    n      /* number of columns */
) ;

void SSMtDenseMult
(
    SSMFLOAT *y,     /* n by 1 product, output */
    SSMFLOAT *x,     /* m by 1 given vector */
    SSMFLOAT *V,     /* dense m by n matrix */
    SSMINT    m,     /* number of rows */
    SSMINT    n      /* number of columns */
) ;

void SSMSSORmultP
(
    SSMFLOAT    *y,  /* the resulting vector */
    SSMFLOAT    *b,  /* vector to be multiplied by SSOR matrix */
    SSMFLOAT   *wj,  /* the first half of the SSOR multiplication operation */
    SSMFLOAT    *X,  /* normalized x */
    SSMFLOAT   *aj,  /* aj = Awj, product of A with wj */
    SSMFLOAT    mu,  /* multiplier */
    int    startup,  /* = 1 for starting multiplication, 0 otherwise */
    SSMProblem *PB,  /* problem specification */
    SSMcom    *Com   /* SSMcom structure */
) ;

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
) ;

int SSMboundary /* return 0 (error tolerance satisfied)
                                -1 (min residual convergence failure in SQP)
                                -2 (SSM failed to converge)
                                -3 (error decay in SSM too slow)
                                -5 (failure of QR diagonalization)
                                -6 (solution in interior ||x|| < r) */
(
    SSMFLOAT    *x, /* estimated solution to ball problem, a'x = 0 */
    SSMcom    *Com  /* SSMcom structure */
) ;

int SSMinterior /* return 0 (error tolerance satisfied)
                                -1 (minimum residual algorithm did not converge)
                                -5 (failure of QR diagonalization)
                                -7 (solution on boundary ||x|| = r) */
(
    SSMFLOAT    *x, /* estimated solution to ball problem, a'x = 0 */
    SSMcom    *Com  /* SSMcom structure */
) ;

int SSMminresP /* return 1 if convergence tolerance met,
                                 0 if decay tolerance met
                                -1 for min residual convergence failure
                                   (too many iterations) */
(
    SSMFLOAT   *xj, /* computed SQP iterate */
    SSMFLOAT   *r0, /* b + Ax, cost gradient at starting point */
    SSMFLOAT    *x, /* solution estimate */
    SSMFLOAT   *Ax, /* A*x */
    SSMFLOAT    *b, /* NULL if b vanishes */
    SSMFLOAT    *X, /* x/||x|| */
    SSMFLOAT    mu, /* multiplier for the constraint */
    SSMFLOAT   rad, /* radius r */
    SSMProblem *PB, /* problem specification */
    SSMcom    *Com  /* pointer to SSMcom structure */
) ;

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
) ;

int SSMsubspace /* return: 0 (error tolerance satisfied)
                                 -5 (failure of QR diagonalization) */
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
) ;
#endif
