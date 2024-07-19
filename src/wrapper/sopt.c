#include "sopt.h"
#define SOPTZERO ((SOPTFLOAT) 0)
#define SOPTONE  ((SOPTFLOAT) 1)

/* =========================================================================
   ============================== sopt_free ================================
   ========================================================================= */
void sopt_free
(
    void * p
)
{

#ifdef MATLAB_MEX_FILE
    mxFree (p) ;
#else
    free (p) ;
#endif

}

/* =========================================================================
   ============================== sopt_malloc ==============================
   ========================================================================= */
void * sopt_malloc
(
    int *status,
    SOPTINT   n,
    int    size
)
{
    void *p ;
    if ( n > 0 )
    {
#ifdef MATLAB_MEX_FILE
        p = mxMalloc ((LONG) n * size) ;
#else
        p = malloc ((LONG) n * size) ;
#endif
    }
    else
    {
        return (NULL) ;
    }
    if      ( p == NULL )                     *status = SOPT_OUT_OF_MEMORY ;
    else if ( *status != SOPT_OUT_OF_MEMORY ) *status = 0 ;
    return (p) ;
}

/* =========================================================================
   ============================== sopt_lmalloc =============================
   ========================================================================= */
void * sopt_lmalloc
(
    int *status,
    LONG      n,
    int    size
)
{
    void *p ;
    if ( n > 0 )
    {
#ifdef MATLAB_MEX_FILE
        p = mxMalloc (n * size) ;
#else
        p = malloc (n * size) ;
#endif
    }
    else
    {
        return (NULL) ;
    }
    if      ( p == NULL )                     *status = SOPT_OUT_OF_MEMORY ;
    else if ( *status != SOPT_OUT_OF_MEMORY ) *status = 0 ;
    return (p) ;
}

void sopt_matrix_default
(
    SOPT_matrix *A
)
{
    A->by_rows = NULL ;
    A->by_cols = NULL ;
    A->rows      = NULL ;
    A->cols      = NULL ;
    A->vals      = NULL ;
    A->nnz       = EMPTY ;
    A->sym       = TRUE ;
    A->fortran   = TRUE ;
    A->p         = NULL ;
    A->i         = NULL ;
    A->x         = NULL ;
}

void sopt_Cmatrix_default
(
    SOPT_Cmatrix *C
)
{
    C->n             = EMPTY ;
    C->nnzfull       = EMPTY ;
    C->nnzsym        = EMPTY ;
    C->rows          = NULL ;
    C->cols          = NULL ;
    C->vals          = NULL ;
    C->nnz           = EMPTY ;
    C->p             = NULL ;
    C->i             = NULL ;
    C->x             = NULL ;
    C->col           = NULL ;
    C->map           = NULL ;
    C->Smap          = NULL ;
    C->SparsityFixed = EMPTY ;
    C->Hformat       = EMPTY ;
    C->malloc_p      = NULL ;
    C->malloc_rows   = NULL ;
    C->Hwork         = NULL ;
    C->rhs           = NULL ;
}

/* ========================================================================== */
/* ====== sopt_error ======================================================== */
/* ========================================================================== */
/* when -g compiler option is used, prints line number of error
   ========================================================================== */
void sopt_error
(
    int status,
    const char *file,
    int line,
    const char *message
)
{
    if (status < 0)
    {
        printf ("file: %s line: %d status: %d %s\n",
                 file, line, status, message) ;
        fflush (stdout) ;
#ifdef MATLAB_MEX_FILE
        mexErrMsgTxt (message) ;
#else
        ASSERT (0) ;
        abort ( ) ;
#endif
    }
}

/* ========================================================================== */
/* ====== sopt_timer ======================================================== */
/* ========================================================================== */
/* Returns current walltime in seconds.  */
/* ========================================================================== */

#ifdef SUITEOPT_TIMER
#include <time.h>

double sopt_timer ( void )
{
    struct timespec t ;
    clock_gettime (CLOCK_MONOTONIC, &t) ;
    return ( (double) t.tv_sec + 1.e-9*((double) t.tv_nsec) ) ;
}

#else

double sopt_timer ( void )
{
    return ( (double) 0 ) ;
}
#endif

/* ==========================================================================
   === sopt_print_TF ========================================================
   ==========================================================================
    Print TRUE if TF is TRUE, otherwise print FALSE
   ========================================================================== */
void sopt_print_TF 
(
    int TF /* TRUE or FALSE (or EMPTY if not True or False) */
)
{
    if ( TF == SuiteOPTtrue )
    {
        printf ("TRUE\n") ;
    }
    else if ( TF == SuiteOPTfalse )
    {
        printf ("FALSE\n") ;
    }
    else
    {
        printf ("EMPTY\n") ;
    }
}

/* ==========================================================================
   === sopt_printA ==========================================================
   ==========================================================================
    Print a sparse matrix
   ========================================================================== */
void sopt_printA
(
    SOPTINT  ncol, /* number of cols in A */
    SOPTINT   *Ap, /* size ncol+1, column pointers */
    SOPTINT  *Anz, /* if NULL, A is packed; otherwise gives # nonzeros in cols*/
    SOPTINT   *Ai, /* size Ap [ncol], row indices for A */
    SOPTFLOAT *Ax, /* size Ap [ncol], numerical entries of A */
    char    *what  /* name of the matrix */
)
{
    SOPTINT j, p, q ;
    printf ("%s =\n", what) ;
    if ( Anz == NULL ) /* matrix is packed */
    {
        p = 0 ;
        for (j = 0; j < ncol; j++)
        {
            q = Ap [j+1] ;
            for (; p < q; p++)
            {
                printf ("%ld %ld %e\n", (LONG) Ai [p], (LONG) j, Ax [p]) ;
            }
        } 
    }
    else /* Anz gives the number of nonzeros in each column of A */
    {
        for (j = 0; j < ncol; j++)
        {
            p = Ap [j] ;
            q = Ap [j] + Anz [j] ;
            for (; p < q; p++)
            {
                printf ("%ld %ld %e\n", (LONG) Ai [p], (LONG) j, Ax [p]) ;
            }
        }
    }
}

/* ==========================================================================
   === sopt_printAMATLAB ====================================================
   ==========================================================================
    Print text that can be fed to MATLAB to construct a sparse matrix.
    The sparse matrix is constructed from triples and rows/columns start at 1.
   ========================================================================== */
void sopt_printAMATLAB
(
    SOPTINT   const ncol, /* number of cols in A */
    SOPTINT   const  *Ap, /* size ncol+1, column pointers */
    SOPTINT   const  *Ai, /* size Ap [ncol], row indices for A */
    SOPTINT   const *Anz, /* if NULL, A packed; otherwise # nonzeros in cols */
    SOPTFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    char           *what  /* name of the matrix */
)
{
    SOPTINT j, p, q ;
    printf ("Temp = [\n") ;
    for (j = 0; j < ncol; j++)
    {
        if ( Anz == NULL ) q = Ap [j+1] ;
        else               q = Ap [j] + Anz [j] ;
        {
            for (p = Ap [j]; p < q; p++)
            {
                if ( Ax == NULL )
                {
                    printf ("%ld %ld 1\n", (LONG) Ai [p]+1, (LONG) j+1) ;
                }
                else
                {
                    printf ("%ld %ld %25.15e\n",
                        (LONG) Ai [p]+1, (LONG) j+1, Ax [p]) ;
                }
            }
        }
    }
    printf ("] ;\n") ;
    printf ("%s = sparse (Temp(:, 1), Temp(:, 2), Temp(:, 3)) ;\n", what) ;
}

/* ==========================================================================
   === sopt_printx ==========================================================
   ==========================================================================
    Print a float vector.
   ========================================================================== */
void sopt_printx
(
    SOPTFLOAT const *x, /* numerical entries in the vector */
    SOPTINT   const  n, /* dimension of the vector */
    char          *what  /* name of the vector */
)
{
    SOPTINT i ;
    printf ("%s =\n", what) ;
    for (i = 0; i < n; i++)
    {
        printf ("%ld %25.15e\n", (LONG) i, x [i]) ;
    }
}

/* ==========================================================================
   === sopt_printxMATLAB ====================================================
   ==========================================================================
    Print text that can be fed to MATLAB to construct a dense vector.
   ========================================================================== */
void sopt_printxMATLAB
(
    SOPTFLOAT *x, /* numerical entries in the vector */
    SOPTINT    n, /* dimension of the vector */
    char   *what  /* name of the vector */
)
{
    SOPTINT j ;
    printf ("%s = [\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%26.16e\n", x [j]) ;
    }
    printf ("] ;\n") ;
}

/* ==========================================================================
   === sopt_printi ==========================================================
   ==========================================================================
    Print a SOPTINT vector.
   ========================================================================== */
void sopt_printi
(
    SOPTINT *i, /* array of SOPTINT's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of array */
)
{
    SOPTINT j ;
    printf ("%s =\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%ld %ld\n", (LONG) j, (LONG) i [j]) ;
    }
}

/* ==========================================================================
   === sopt_print_int =======================================================
   ==========================================================================
    Print an int vector.
   ========================================================================== */
void sopt_print_int 
(
    int     *i, /* array of int's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of array */
)
{
    SOPTINT j ;
    printf ("%s =\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%ld %i\n", (LONG) j, i [j]) ;
    }
}

/* ==========================================================================
   === sopt_print_iMATLAB ====================================================
   ==========================================================================
    Print text that can be fed to MATLAB to construct a dense int vector.
   ========================================================================== */
void sopt_printiMATLAB
(
    SOPTINT *i, /* array of SOPTINT's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of the array */
)
{
    SOPTINT j ;
    printf ("%s = [\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%ld\n", (LONG) i [j]) ;
    }
    printf ("] ;\n") ;
}

/* =========================================================================
   ===================== sopt_add ==========================================
   =========================================================================
   add scalar s to each component of array x
   ========================================================================= */

void sopt_add
(
    SOPTFLOAT       *x,  /* array to which s is added */
    SOPTFLOAT const  s,  /* scalar */
    SOPTINT   const  n   /* length of x */
)
{
    SOPTINT j, n5 ;
    SOPTFLOAT *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] += s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
    }
}

/* =========================================================================
   ===================== sopt_addi =========================================
   =========================================================================
   add scalar s to each component of array x
   ========================================================================= */

void sopt_addi
(
    SOPTINT       *x,  /* array to which s is added */
    SOPTINT const  s,  /* scalar */
    SOPTINT const  n   /* length of x */
)
{
    SOPTINT j, n5 ;
    SOPTINT *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] += s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
        *(xj++) += s ;
    }
}

/* =========================================================================
   ===================== sopt_add2i ========================================
   =========================================================================
   add scalar s to each component of array y and stored in x
   ========================================================================= */

void sopt_add2i
(
    SOPTINT        *x,  /* array to which y is stored  */
    SOPTINT  const *y,  /* array to which s is added */
    SOPTINT  const  s,  /* scalar */
    SOPTINT  const  n   /* length of x */
)
{
    SOPTINT j, n5 ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s + y [j] ;
    for (; j < n; j += 5)
    {
        x [j]   = s + y [j] ;
        x [j+1] = s + y [j+1] ;
        x [j+2] = s + y [j+2] ;
        x [j+3] = s + y [j+3] ;
        x [j+4] = s + y [j+4] ;
    }
}

/* =========================================================================
   ===================== sopt_scale ========================================
   =========================================================================
   Scale a SOPTFLOAT array x = s*y
   ========================================================================= */

void sopt_scale
(
    SOPTFLOAT       *x,  /* scaled array */
    SOPTFLOAT const *y,  /* array used for the scaling */
    SOPTFLOAT const  s,  /* scale */
    SOPTINT   const  n   /* length of x */
)
{
    if ( x == y )
    {
        SOPTFLOAT *xj ;
        if ( s == SOPTONE ) return ; /* nothing to do */
#ifndef NOBLAS
        if ( n >= DSCAL_START )
        {
            SOPTFLOAT S = s ;
            BLAS_INT int_one = 1 ;
            BLAS_INT N = (BLAS_INT) n ;
            SOPT_DSCAL (&N, &S, x, &int_one) ;
        }
        else
#endif
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++) x [j] *= s ;
            xj = x+j ;
            for (; j < n; j += 5)
            {
                *(xj++) *= s ;
                *(xj++) *= s ;
                *(xj++) *= s ;
                *(xj++) *= s ;
                *(xj++) *= s ;
            }
        }
    }
    else /* x != y */
    {
        if ( s == SOPTONE ) /* this is a copy */
        {
            sopt_copyx (x, y, n) ;
        }
        else if ( s == -SOPTONE ) /* s = -1 */
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++) x [j] = -y [j] ;
            for (; j < n; j += 5)
            {
                x [j]   = -y [j] ;
                x [j+1] = -y [j+1] ;
                x [j+2] = -y [j+2] ;
                x [j+3] = -y [j+3] ;
                x [j+4] = -y [j+4] ;
            }
        }
        else
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++) x [j] = s*y [j] ;
            for (; j < n; j += 5)
            {
                x [j]   = s*y [j] ;
                x [j+1] = s*y [j+1] ;
                x [j+2] = s*y [j+2] ;
                x [j+3] = s*y [j+3] ;
                x [j+4] = s*y [j+4] ;
            }
        }
    }
}

/* =========================================================================
   ===================== sopt_scale_max ====================================
   =========================================================================
   Scale a SOPTFLOAT array x = s*y and also evaluate ||y||_sup
   ========================================================================= */

SOPTFLOAT sopt_scale_max
(
    SOPTFLOAT       *x,  /* scaled array */
    SOPTFLOAT const *y,  /* array used for the scaling */
    SOPTFLOAT const  s,  /* scale */
    SOPTINT   const  n   /* length of x */
)
{
    SOPTFLOAT ymax, *xj ;
    ymax = SOPTZERO ;
    if ( x == y )
    {
        if ( s == SOPTONE ) /* nothing to do except compute the max */
        {
            ymax = sopt_sup_normx (y, n) ;
        }
#ifndef NOBLAS
        else if ( n >= DSCAL_START )
        {
            ymax = sopt_sup_normx (y, n) ;
            SOPTFLOAT S = s ;
            BLAS_INT int_one = 1 ;
            BLAS_INT N = (BLAS_INT) n ;
            SOPT_DSCAL (&N, &S, x, &int_one) ;
        }
#endif
        else
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++)
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] *= s ;
            }
            xj = x+j ;
            for (; j < n; j += 5)
            {
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
                if ( ymax < fabs(*xj) ) ymax = fabs(*xj) ;
                *(xj++) *= s ;
            }
        }
    }
    else /* x != y */
    {
        if ( s == SOPTONE ) /* this is a copy */
        {
            ymax = sopt_sup_normx (y, n) ;
            sopt_copyx (x, y, n) ;
        }
        else if ( s == -SOPTONE ) /* s = -1 */
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++)
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ;
            }
            for (; j < n; )
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = -y [j] ; j++ ;
            }
        }
        else
        {
            SOPTINT j, n5 ;
            n5 = n % 5 ;
            for (j = 0; j < n5; j++)
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ;
            }
            for (; j < n; )
            {
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
                if ( ymax < fabs(y [j]) ) ymax = fabs(y [j]) ;
                x [j] = s*y [j] ; j++ ;
            }
        }
    }
    return (ymax) ;
}

/* =========================================================================
   === sopt_step ============================================================
   =========================================================================
    Set xnew = x + alpha*d
   ========================================================================= */
void sopt_step
(
    SOPTFLOAT       *xnew, /* updated x vector */
    SOPTFLOAT const    *x, /* current x */
    SOPTFLOAT const    *d, /* search direction */
    SOPTFLOAT const alpha, /* stepsize */
    SOPTINT   const     n  /* dimension */
)
{
    SOPTINT j, n5 ;

    if ( x == NULL ) return ;
    n5 = n % 5 ;     /* n5 = n mod 5 */
    /* check if step size equals 1 */
    if ( alpha == SOPTONE )
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + d [j] ;
        }
        if ( xnew == x )
        {
            for (; j < n; )
            {
                 xnew [j] += d [j] ; j++ ;
                 xnew [j] += d [j] ; j++ ;
                 xnew [j] += d [j] ; j++ ;
                 xnew [j] += d [j] ; j++ ;
                 xnew [j] += d [j] ; j++ ;
            }
        }
        else
        {
            for (; j < n; )
            {
                 xnew [j] = x [j] + d [j] ; j++ ;
                 xnew [j] = x [j] + d [j] ; j++ ;
                 xnew [j] = x [j] + d [j] ; j++ ;
                 xnew [j] = x [j] + d [j] ; j++ ;
                 xnew [j] = x [j] + d [j] ; j++ ;
            }
        }
    }
    else if ( alpha == -SOPTONE )
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] - d [j] ;
        }
        for (; j < n; )
        {
            xnew [j] = x [j] - d [j] ; j++ ;
            xnew [j] = x [j] - d [j] ; j++ ;
            xnew [j] = x [j] - d [j] ; j++ ;
            xnew [j] = x [j] - d [j] ; j++ ;
            xnew [j] = x [j] - d [j] ; j++ ;
        }
    }
    /* else step size is not 1 */
    else
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + alpha*d [j] ;
        }
        for (; j < n; )
        {
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
            xnew [j] = x [j] + alpha*d [j] ; j++ ;
        }
    }
}

/* =========================================================================
   === sopt_step_max =======================================================
   =========================================================================
    Set xnew = x + alpha*d, return dmax = ||d||_sup
   ========================================================================= */
SOPTFLOAT sopt_step_max
(
    SOPTFLOAT       *xnew, /* updated x vector */
    SOPTFLOAT const    *x, /* current x */
    SOPTFLOAT const    *d, /* search direction */
    SOPTFLOAT const alpha, /* stepsize */
    SOPTINT   const     n  /* dimension */
)
{
    SOPTFLOAT dmax ;
    SOPTINT j, n5 ;

    dmax = SOPTZERO ;
    n5 = n % 5 ;     /* n5 = n mod 5 */
    /* check if step size equals 1 */
    if ( alpha == SOPTONE )
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
        }
        if ( xnew == x )
        {
            for (; j < n; )
            {
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] += d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
            }
        }
        else
        {
            for (; j < n; )
            {
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
                xnew [j] = x [j] + d [j] ;
                if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
                j++ ;
            }
        }
    }
    else if ( alpha == -SOPTONE )
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
        }
        for (; j < n; )
        {
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] - d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
        }
    }
    /* else step size is not 1 */
    else
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
        }
        for (; j < n; )
        {
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
            xnew [j] = x [j] + alpha*d [j] ;
            if ( dmax < fabs(d [j]) ) dmax = fabs(d [j]) ;
            j++ ;
        }
    }
    return (dmax) ;
}

/* =========================================================================
   ==== sopt_daxpy =========================================================
   =========================================================================
   Compute x = x + s * d
   ========================================================================= */
void sopt_daxpy
(
    SOPTFLOAT       *x, /* input and output vector */
    SOPTFLOAT const *d, /* direction vector */
    SOPTFLOAT const  s, /* stepsize */
    SOPTINT   const  n  /* length of the vectors */
)
{
#ifndef NOBLAS
    if ( n >= DAXPY_START )
    {
        SOPTFLOAT S = s ;
        BLAS_INT int_one = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        SOPT_DAXPY (&N, &S, (SOPTFLOAT *) d, &int_one, x, &int_one) ;
    }
    else
#endif
    {
        SOPTINT i, n5 ;
        n5 = n % 5 ;
        if ( s == -SOPTONE)
        {
            for (i = 0; i < n5; i++) x [i] -= d [i] ;
            for (; i < n; i += 5)
            {
                x [i]   -= d [i] ;
                x [i+1] -= d [i+1] ;
                x [i+2] -= d [i+2] ;
                x [i+3] -= d [i+3] ;
                x [i+4] -= d [i+4] ;
            }
        }
        else if ( s == SOPTONE )
        {
            for (i = 0; i < n5; i++) x [i] += d [i] ;
            for (; i < n; i += 5)
            {
                x [i]   += d [i] ;
                x [i+1] += d [i+1] ;
                x [i+2] += d [i+2] ;
                x [i+3] += d [i+3] ;
                x [i+4] += d [i+4] ;
            }
        }
        else
        {
            for (i = 0; i < n5; i++) x [i] += s*d [i] ;
            for (; i < n; i += 5)
            {
                x [i]   += s*d [i] ;
                x [i+1] += s*d [i+1] ;
                x [i+2] += s*d [i+2] ;
                x [i+3] += s*d [i+3] ;
                x [i+4] += s*d [i+4] ;
            }
        }
    }
}

/* =========================================================================
   === sopt_copyx ==========================================================
   =========================================================================
   Copy SOPTFLOAT vector y into vector x
   ========================================================================= */
void sopt_copyx
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    SOPTINT   const  n  /* length of vectors */
)
{
    if ( (y == x) || (y == NULL) ) return ;
#ifndef NOBLAS
    if ( n >= DCOPY_START )
    {
        BLAS_INT int_one  = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        SOPT_DCOPY (&N, (SOPTFLOAT *) y, &int_one, x, &int_one) ;
    }
    else
#endif
    {
        SOPTINT i, n5 ;
        n5 = n % 5 ;
        for (i = 0; i < n5; i++) x [i] = y [i] ;
        for (; i < n; i += 5)
        {
            x [i]   = y [i] ;
            x [i+1] = y [i+1] ;
            x [i+2] = y [i+2] ;
            x [i+3] = y [i+3] ;
            x [i+4] = y [i+4] ;
        }
    }
}

/* =========================================================================
   === sopt_lcopyx ==========================================================
   =========================================================================
   Copy SOPTFLOAT vector y into vector x
   ========================================================================= */
void sopt_lcopyx
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    LONG      const  n  /* length of vectors */
)
{
    if ( (y == x) || (y == NULL) ) return ;
#ifndef NOBLAS
    if ( n >= DCOPY_START )
    {
        BLAS_INT int_one  = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        SOPT_DCOPY (&N, (SOPTFLOAT *) y, &int_one, x, &int_one) ;
    }
    else
#endif
    {
        LONG i, n5 ;
        n5 = n % 5 ;
        for (i = 0; i < n5; i++) x [i] = y [i] ;
        for (; i < n; i += 5)
        {
            x [i]   = y [i] ;
            x [i+1] = y [i+1] ;
            x [i+2] = y [i+2] ;
            x [i+3] = y [i+3] ;
            x [i+4] = y [i+4] ;
        }
    }
}

/* =========================================================================
   === sopt_copyx_noblas ===================================================
   =========================================================================
   Copy SOPTFLOAT vector y into vector x
   NOTE: when y is a part of x, cannot use BLAS
   ========================================================================= */
void sopt_copyx_noblas
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    SOPTINT   const  n  /* length of vectors */
)
{
    if ( (y == x) || (y == NULL) || (n == 0) ) return ;
    SOPTINT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) x [i] = y [i] ;
    for (; i < n; i += 5)
    {
        x [i]   = y [i] ;
        x [i+1] = y [i+1] ;
        x [i+2] = y [i+2] ;
        x [i+3] = y [i+3] ;
        x [i+4] = y [i+4] ;
    }
}

/* =========================================================================
   === sopt_copyi_noblas ===================================================
   =========================================================================
   Copy SOPTINT vector y into vector x
   NOTE: when y is a part of x, cannot use BLAS
   ========================================================================= */
void sopt_copyi_noblas
(
    SOPTINT       *x, /* output of copy */
    SOPTINT const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
)
{
    if ( (y == x) || (y == NULL) || (n == 0) ) return ;
    SOPTINT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) x [i] = y [i] ;
    for (; i < n; i += 5)
    {
        x [i]   = y [i] ;
        x [i+1] = y [i+1] ;
        x [i+2] = y [i+2] ;
        x [i+3] = y [i+3] ;
        x [i+4] = y [i+4] ;
    }
}

/* =========================================================================
   ======================== sopt_copyi =====================================
   =========================================================================
   Copy SOPTINT vector y into vector x
   ========================================================================= */
void sopt_copyi
(
    SOPTINT       *x, /* output of copy */
    SOPTINT const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
)
{
    SOPTINT i, n5 ;
    if ( (x == NULL) || (y == x) ) return ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) x [i] = y [i] ;
    for (; i < n; i += 5)
    {
        x [i]   = y [i] ;
        x [i+1] = y [i+1] ;
        x [i+2] = y [i+2] ;
        x [i+3] = y [i+3] ;
        x [i+4] = y [i+4] ;
    }
}

/* =========================================================================
   ======================== sopt_copy_int ==================================
   =========================================================================
   Copy int vector y into vector x
   ========================================================================= */
void sopt_copy_int
(
    int           *x, /* output of copy */
    int     const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
)
{
    SOPTINT i, n5 ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) x [i] = y [i] ;
    for (; i < n; i += 5)
    {
        x [i]   = y [i] ;
        x [i+1] = y [i+1] ;
        x [i+2] = y [i+2] ;
        x [i+3] = y [i+3] ;
        x [i+4] = y [i+4] ;
    }
}

/*  =========================================================================
    ==== sopt_dot ===========================================================
    =========================================================================
    Compute dot product of x and y
    ========================================================================= */
SOPTFLOAT sopt_dot
(
    SOPTFLOAT const *x, /* first vector */
    SOPTFLOAT const *y, /* second vector */
    SOPTINT   const  n  /* length of vectors */
)
{
#ifndef NOBLAS
    if ( n >= DDOT_START )
    {
        BLAS_INT int_one = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        return (SOPT_DDOT (&N, (SOPTFLOAT *) x, &int_one,
                               (SOPTFLOAT *) y, &int_one)) ;
    }
    else
#endif
    {
        SOPTINT i, n5 ;
        SOPTFLOAT t = SOPTZERO ;
        if ( n <= 0 ) return (t) ;
        n5 = n % 5 ;
        for (i = 0; i < n5; i++) t += x [i]*y [i] ;
        for (; i < n; i += 5)
        {
            t += x [i]*y [i] + x [i+1]*y [i+1] + x [i+2]*y [i+2] 
                             + x [i+3]*y [i+3] + x [i+4]*y [i+4] ;
        }
        return (t) ;
    }
}

/* =========================================================================
   ===================== sopt_initx ========================================
   =========================================================================
   Initialize a SOPTFLOAT array
   ========================================================================= */
void sopt_initx
(
    SOPTFLOAT      *x,  /* array to be initialized */
    SOPTFLOAT const s,  /* scalar */
    SOPTINT   const n   /* length of x */
)
{
    SOPTINT j, n5 ;
    SOPTFLOAT *xj ;
    if ( n == 0 ) return ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

/* =========================================================================
   ===================== sopt_initi ========================================
   =========================================================================
   Initialize a SOPTINT array
   ========================================================================= */
void sopt_initi
(
    SOPTINT      *x,  /* array to be initialized */
    SOPTINT const s,  /* scalar */
    SOPTINT const n   /* length of x */
)
{
    SOPTINT j, n5 ;
    SOPTINT *xj ;
    if ( n == 0 ) return ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

/* =========================================================================
   ===================== sopt_init_int =====================================
   =========================================================================
   Initialize an int array
   ========================================================================= */
void sopt_init_int
(
    int          *x,  /* array to be initialized */
    int     const s,  /* scalar */
    SOPTINT const n   /* length of x */
)
{
    SOPTINT j, n5 ;
    int     *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

/* =========================================================================
   ==== sopt_sup_normx =====================================================
   =========================================================================
   Compute sup norm of a SOPTFLOAT vector
   ========================================================================= */
SOPTFLOAT sopt_sup_normx
(
    SOPTFLOAT const *x, /* vector */
    SOPTINT   const  n  /* length of vector */
)
{
#ifndef NOBLAS
    if ( n >= IDAMAX_START )
    {
        BLAS_INT int_one = 1 ;
        BLAS_INT N = (BLAS_INT) n ;
        BLAS_INT i = SOPT_IDAMAX (&N, (SOPTFLOAT *) x, &int_one) ;
        return (fabs (x [i-1])) ; /* adjust for fortran indexing */
    }
    else
#endif
    {
        SOPTINT i, n5 ;
        SOPTFLOAT t ;
        t = SOPTZERO ;
        n5 = n % 5 ;

        for (i = 0; i < n5; i++) if ( t < fabs (x [i]) ) t = fabs (x [i]) ;
        for (; i < n; i += 5)
        {
            if ( t < fabs (x [i]  ) ) t = fabs (x [i]  ) ;
            if ( t < fabs (x [i+1]) ) t = fabs (x [i+1]) ;
            if ( t < fabs (x [i+2]) ) t = fabs (x [i+2]) ;
            if ( t < fabs (x [i+3]) ) t = fabs (x [i+3]) ;
            if ( t < fabs (x [i+4]) ) t = fabs (x [i+4]) ;
        }
        return (t) ;
    }
}

/* =========================================================================
   === sopt_supi ===========================================================
   =========================================================================
   Return largest entry of an SOPTINT vector
   ========================================================================= */
SOPTINT sopt_supi
(
    SOPTINT const *x, /* vector */
    SOPTINT const  n  /* length of vector */
)
{
    SOPTINT xsup ;
    SOPTINT j, n5 ;
    n5 = n % 5 ;             /* n5 = n mod 5 */
    xsup = -SuiteOPTinfint ; /* initializing xsup */
    for (j = 0; j < n5; j++)
    {
        if ( xsup < x [j] ) xsup = x [j] ;
    }
    for (; j < n; j += 5)
    {
        if ( xsup < x [j]   ) xsup = x [j] ;
        if ( xsup < x [j+1] ) xsup = x [j+1] ;
        if ( xsup < x [j+2] ) xsup = x [j+2] ;
        if ( xsup < x [j+3] ) xsup = x [j+3] ;
        if ( xsup < x [j+4] ) xsup = x [j+4] ;
    }
    return (xsup) ;
}

/* ========================================================================== */
/* === sopt_transpose ======================================================= */
/* ========================================================================== */
/*    Transpose a sparse matrix: B = A' (when A is symmetric, Bp = Ap so Bp   */
/*                                       does not need to be computed         */
/* ========================================================================== */
int sopt_transpose
(
    SOPTINT          *Bp, /* size nrow+1, column pointers (output) = Ap if sym*/
    SOPTINT          *Bi, /* size Ap [ncol], row indices of B (output) */
    SOPTFLOAT        *Bx, /* size Ap [ncol], numerical entries of B (output) */
    SOPTINT          *Ap, /* size ncol+1, column pointers */
    SOPTINT   const  *Ai, /* size Ap [ncol], row indices for A */
    SOPTFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    SOPTINT   const nrow, /* number of rows in A */
    SOPTINT   const ncol, /* number of cols in A */
    int       const  sym, /* TRUE if A is symmetric */
    SOPTINT        *work  /* if not NULL, then workspace is provided */
)
{
    SOPTINT i, j, p, *W ;

    int status = 0 ;
    if ( work == NULL )
    {
        W = (SOPTINT *) sopt_malloc (&status, nrow, sizeof (SOPTINT)) ;
        if ( status ) return (status) ;
    }
    else
    {
        W = work ;
    }
    if ( sym )
    {
        sopt_copyi (W, Ap, ncol) ;
    }
    else /* not symmetric */
    {
        /* ================================================================== */
        /* === compute row counts of A ====================================== */
        /* ================================================================== */

        sopt_initi (W, (SOPTINT) 0, nrow) ;
        p = 0 ;
        for (j = 1; j <= ncol; j++)
        {
            SOPTINT const q = Ap [j] ;
            for (; p < q; p++)
            {
                W [Ai [p]]++ ;
            }
        }

        /* ================================================================== */
        /* === compute column pointers of B from the row counts ============= */
        /* ================================================================== */

        Bp [0] = 0 ;
        for (i = 0; i < nrow; i++)
        {
            Bp [i+1] = Bp [i] + W [i] ;
            W [i] = Bp [i] ;
        }
    }

    /* ====================================================================== */
    /* === B = A' =========================================================== */
    /* ====================================================================== */

    p = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
        SOPTINT const q = Ap [j+1] ;
        for (; p < q; p++)
        {
            SOPTINT const pp = W [Ai [p]]++ ;
            Bi [pp] = j ;
            Bx [pp] = Ax [p] ;
        }
    }
    /* if symmetric, we assume that user sets Bp = Ap */
    if ( work == NULL ) sopt_free (W) ; /* nonsymmetric and W malloc'd */
    return (status) ;
}


/* ========================================================================== */
/* === sopt_sort_cols ======================================================= */
/* ========================================================================== */
/*    Perform transposes to sort row indices in each column of matrix         */
/*    If the matrix is symmetric, then one transpose sorts the columns and    */
/*    the transpose matrix AT, which equals the original matrix by symmetry,  */
/*    contains sorted columns. If the matrix is not symmetric, then two       */
/*    transposes are needed to sort the columns. The output matrix A has the  */
/*    sorted columns.                                                         */
/* ========================================================================== */
int sopt_sort_cols /* returned integer:
                         0 if conversion successful
                         SOPT_OUT_OF_MEMORY otherwise */
(
    SOPTINT         *Ap, /* column pointers */
    SOPTINT         *Ai, /* row indices */
    SOPTFLOAT       *Ax, /* numerical values */
    SOPTINT        *ATp, /* row pointers for transpose */
    SOPTINT        *ATi, /* column indices for transpose */
    SOPTFLOAT      *ATx, /* numerical values for transpose */
    SOPTINT        nrow, /* number of rows */
    SOPTINT        ncol, /* number of cols */
    int const       sym  /* TRUE if A is symmetric matrix */
)
{
    int status, sI, sF, mallocAT ;
    status = 0 ;
    mallocAT = FALSE ;
    sI = sizeof (SOPTINT) ;
    sF = sizeof (SOPTFLOAT) ;
    if ( ATi == NULL )
    {
        mallocAT = TRUE ;
        ATi = (SOPTINT *)   sopt_malloc (&status, Ap [ncol], sI) ;
        ATx = (SOPTFLOAT *) sopt_malloc (&status, Ap [ncol], sF);
    }
    if ( !sym ) /* matrix not symmetric */
    {
        if ( mallocAT == TRUE && (ATp != NULL) )
        {
            sopt_free (ATp) ;
            ATp = NULL ;
        }
        if ( ATp == NULL )
        {
            ATp = (SOPTINT *) sopt_malloc (&status, (nrow+1), sI) ;
        }
    }
    if ( status ) return (status) ;
    /* if the matrix is symmetric, then a single transpose sorts the columns */
    status = sopt_transpose (ATp, ATi, ATx, Ap, Ai, Ax, nrow, ncol, sym, NULL) ;
    if ( status || sym ) return (status) ;

    /* for a nonsymmetric matrix, a second transpose sorts the columns */
    status = sopt_transpose (Ap, Ai, Ax, ATp, ATi, ATx, ncol,nrow,sym,NULL);
    if ( mallocAT )
    {
        sopt_free (ATp) ;
        sopt_free (ATi) ;
        sopt_free (ATx) ;
    }
    return (status) ;
}

/* =========================================================================
   === sopt_convert_to_sparse ==============================================
   =========================================================================
   An input Amatrix structure is processed to build the sparse matrix
   component of the structure if it does not yet exist. If either the
   triples arrays or a dense array is provided in the Amatrix structure,
   then they are converted to the sparse matrix structure.
   ========================================================================= */
int sopt_convert_to_sparse /* returned integer:
                                   0 if conversion successful
                                   OUT_OF_MEMORY
                                   ERROR_IN_INPUT_MATRIX */
(
    SOPT_matrix       *A, /* input and output matrix */
    SOPTINT         nrow, /* number of rows in A */
    SOPTINT         ncol, /* number of columns in A */
    int const order_cols, /* TRUE if row indices in each column of the sparse
                             matrix format should be put in increasing order */
    int const   no_zeros, /* if an element of Tx is zero, then delete it */
    int         location  /* location (PASA, LCG, PPROJ) where routine
                             was invoked */
)
{
    int status ;
    status = 0 ;
    A->nrow = nrow ;
    A->ncol = ncol ;
    if ( A->p != NULL )
    {
        return (SOPT_EXISTING_MATRIX_SPARSE) ;
    }
    else if ( (A->rows != NULL) && (A->cols != NULL) && (A->vals != NULL) )
    {
        status = sopt_convert_triple_to_sparse (A, order_cols, no_zeros) ;
    }
    else if ( (A->by_rows != NULL) || (A->by_cols != NULL) )
    {
        status = sopt_convert_dense_to_sparse (A) ;
    }
    else return (SOPT_NO_MATRIX_PROVIDED) ;
    /* If an error occurred, convert the error message based on the
       location where the error occurred. */
    if ( status > 0 ) return (sopt_convert_error (location, status)) ;
    return (SOPT_SPARSE_MATRIX_CREATED) ;
}

int sopt_convert_error
(
    int  location,
    int  status
)
{
    if ( location == LSOPT ) return (status) ;

    if ( location == LPASA )
    {
        if (      status == SOPT_OUT_OF_MEMORY )
                  status  = PASA_OUT_OF_MEMORY ;
        else if ( status == SOPT_ERROR_IN_INPUT_MATRIX )
                  status  = PASA_ERROR_IN_INPUT_MATRIX ;
        else if ( status == SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE )
                  status  = PASA_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE ;
        else if ( status == SOPT_TRIPLES_FORMAT_ERROR )
                  status  = PASA_TRIPLES_FORMAT_ERROR ;
        else if ( status == SOPT_HESSIAN_NOT_COMPUTED )
                  status  = PASA_HESSIAN_NOT_COMPUTED ;
    }
    else if ( location == LCG )
    {
        if (      status == SOPT_OUT_OF_MEMORY )
                  status  =   CG_OUT_OF_MEMORY ;
        else if ( status == SOPT_ERROR_IN_INPUT_MATRIX )
                  status  =   CG_ERROR_IN_INPUT_MATRIX ;
        else if ( status == SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE )
                  status  =   CG_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE ;
        else if ( status == SOPT_TRIPLES_FORMAT_ERROR )
                  status  =   CG_TRIPLES_FORMAT_ERROR ;
        else if ( status == SOPT_HESSIAN_NOT_COMPUTED )
                  status  =   CG_HESSIAN_NOT_COMPUTED ;
    }
    else if ( location == LPPROJ )
    {
        if (      status ==  SOPT_OUT_OF_MEMORY )
                  status  = PPROJ_OUT_OF_MEMORY ;
        else if ( status ==  SOPT_ERROR_IN_INPUT_MATRIX )
                  status  = PPROJ_ERROR_IN_INPUT_MATRIX ;
        else if ( status ==  SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE )
                  status  = PPROJ_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE ;
        
    }
    else if ( location == LCUTE )
    {
        if (      status == SOPT_OUT_OF_MEMORY )
                  status  = CUTE_OUT_OF_MEMORY ;
        else if ( status == SOPT_ERROR_IN_INPUT_MATRIX )
                  status  = CUTE_ERROR_IN_INPUT_MATRIX ;
        else if ( status == SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE )
                  status  = CUTE_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE ;
    }
    else if ( location == LSSM )
    {
        if (      status == SOPT_OUT_OF_MEMORY )
                  status  = SSM_OUT_OF_MEMORY ;
    }
    return (status) ;
}
/* =========================================================================
   === sopt_convert_dense_to_sparse ========================================
   =========================================================================
   convert a dense matrix into a sparse matrix
   Sparse matrices are stored using 3 arrays:

        Ai - row indices of the nonzeros elements. Indices for each
             column should be in increasing order.
        Ax - the nonzeros numerical entries corresponding to elements of Ai
        Ap - the array of column pointers of size ncol + 1.  Ap [j] is
             location of first nonzero in Ax associated with column j and
             Ap [ncol] is the number of nonzeros in the matrix
   ========================================================================= */
int sopt_convert_dense_to_sparse /* returned integer:
                                  0 if conversion successful
                                  SOPT_OUT_OF_MEMORY
                                  SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE
                                  SOPT_ERROR_IN_INPUT_MATRIX */
(
    SOPT_matrix *A
)
{
    int by_rows ;
    SOPTINT i, j, k, nnz, *ap, *ap_shift, *ai ;
    SOPTFLOAT *ax, *dense ;
    SOPTINT nrow = A->nrow ;
    SOPTINT ncol = A->ncol ;
    if ( (nrow <= 0 ) || (ncol <= 0) )
    {
        printf ("\nWhen converting a dense matrix to a sparse matrix, it\n"
                "was found that either the number of rows or the number\n"
                "of columns were not provided.\n") ;
        return (SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE) ;
    }
    int status = 0 ;
    int sI = sizeof (SOPTINT) ;
    int sF = sizeof (SOPTFLOAT) ;
    
    /* count the number of nonzeros in each column */
    A->p = ap = (SOPTINT *) sopt_malloc (&status, (ncol+1), sI) ;
    if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;
    sopt_initi (ap, (SOPTINT) 0, ncol+1) ;
    if ( A->by_rows == NULL )
    {
        dense = A->by_cols ;
        by_rows = FALSE ;
    }
    else
    {
        dense = A->by_rows ;
        by_rows = TRUE ;
    }
    k = 0 ;
    if ( by_rows )
    {
        for (i = 0; i < nrow; i++)
        {
            SOPTINT const l = k + ncol ;
            ap_shift = ap-k ;
            for (; k < l; k++)
            {
                if ( dense [k] != SOPTZERO )
                {
                    ap_shift [k]++ ;
                }
            }
        }
        /* ap now contains the number of nonzeros in each column of A
           change ap to pointers to the start of each column */
        k = 0 ;
        for (j = 0; j < ncol; j++)
        {
            SOPTINT const l = k + ap [j] ;
            ap [j] = k ;
            k = l ;
        }
        ap [ncol] = k ;
        if ( k == 0 )
        {
            sopt_free (A->p) ;
            A->p = A->i = NULL ; A->x = NULL ;
            return (status) ; /* the matrix is all zero */
        }
        A->i = ai = (SOPTINT *)   sopt_malloc (&status, k, sI) ;
        A->x = ax = (SOPTFLOAT *) sopt_malloc (&status, k, sF) ;
        if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;
        k = 0 ;
        for (i = 0; i < nrow; i++)
        {
            SOPTINT const l = k + ncol ;
            ap_shift = ap-k ;
            for (; k < l; k++)
            {
                if ( dense [k] )
                {
                    j = ap_shift [k]++ ;
                    ai [j] = i ;
                    ax [j] = dense [k] ;
                }
            }
        }
        for (j = ncol-1; j > 0; j--)
        {
            ap [j] = ap [j-1] ;
        }
        ap [0] = 0 ;
    }
    else /* matrix is stored by columns */
    {
        /* When a dense matrix stored by columns is converted to a sparse
           matrix by the program sopt_convert_to_sparse, the leading
           dimension must be >= the number of rows in the matrix. */
        ap [0] = 0 ;
        nnz = 0 ;
        for (j = 0; j < ncol; j++)
        {
            SOPTINT const l = k + nrow ;
            for (; k < l; k++)
            {
                if ( dense [k] ) nnz++ ;
            }
            ap [j+1] = nnz ;
        }
        A->i = ai = (SOPTINT *)   sopt_malloc (&status, nnz, sI) ;
        A->x = ax = (SOPTFLOAT *) sopt_malloc (&status, nnz, sF) ;
        if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;
        nnz = k = 0 ;
        for (j = 0; j < ncol; j++)
        {
            SOPTINT const k0 = k ;
            SOPTINT const l = k + nrow ;
            for (; k < l; k++)
            {
                if ( dense [k] )
                {
                    ai [nnz] = k - k0 ;
                    ax [nnz] = dense [k] ;
                    nnz++ ;
                }
            }
        }
    }
    return (status) ;
}
/* =========================================================================
   === sopt_convert_triple_to_sparse =======================================
   =========================================================================
   convert triples Ti (row indices), Tj (column indices), Tx (numerical values)
   to              Ai (row indices), Ap (column pointer), Ax (numerical values)

   If either the nrow or the ncol argument are equal to EMPTY, then their
   values are set to 1 + max(Ti) or 1 + max(Tj) respectively.

   Sparse matrices are stored using 3 arrays:

        Ai - row indices of the nonzeros elements. Indices for each
             column should be in increasing order.
        Ax - the nonzeros numerical entries corresponding to elements of Ai
        Ap - the array of column pointers of size ncol + 1.  Ap [j] is
             location of first nonzero in Ax associated with column j and
             Ap [ncol] is the number of nonzeros in the matrix
   ========================================================================= */
int sopt_convert_triple_to_sparse /* returned integer:
                                   0 if conversion successful
                                   SOPT_OUT_OF_MEMORY
                                   SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE
                                   SOPT_ERROR_IN_INPUT_MATRIX */
(
    SOPT_matrix *A, /* matrix structure where sparse matrix is stored */
    int order_cols, /* = TRUE => order columns in increasing order */
    int   no_zeros  /* = TRUE => delete all zeros, only store nonzero entries */
)
{
    SOPTINT imax, i, j, k, knnz, l, *ai, *ap, *Bp, *Ci, *Cp, *cols, *rows ;
    SOPTFLOAT *ax, *Cx, *vals ;

    /* set default status */
    int status = 0 ;
    int const sI = sizeof (SOPTINT) ;
    int const sF = sizeof (SOPTFLOAT) ;
    /* sym = TRUE  => matrix is symmetric and only elements on main diagonal and
                      and one side are given
             FALSE => all nonzero matrix elements are given*/
    int const sym = A->sym ;
    /* TRUE  => Fortran indexing is used where first row and column are 1
       FALSE => first row and column are 0 */
    int const fortran = A->fortran ;
         
    SOPTINT nrow = A->nrow ;
    SOPTINT ncol = A->ncol ;
    SOPTINT nnz  = A->nnz ;

    /* immediately return if no matrix is input */
    if ( (nrow == 0) || (ncol == 0) || (nnz == 0) )
    {
        return (status) ;
    }

    if ( nnz < 0 )
    {
        printf ("\nWhen converting a matrix in triples format to a sparse\n"
                "matrix, it was found that the number of nonzeros nnz for\n"
                "the triples was negative.\n") ;
        return (SOPT_ERROR_IN_INPUT_MATRIX) ;
    }

    A->p = ap = (SOPTINT *) sopt_malloc (&status, ncol+1, sI) ;
    if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;

    rows = A->rows ;
    cols = A->cols ;
    vals = A->vals ;

    /* count number of nonzeros in each column */
    imax = -1 ;
    sopt_initi (ap, (SOPTINT) 0, ncol) ;
    for (l = 0; l < nnz; l++)
    {
        if ( fortran ) j = cols [l] - 1 ; /* convert Fortran indexing to C */
        else           j = cols [l] ;
    
        if ( j < 0 )
        {
            if ( fortran )
            {
                printf ("\nWhen converting a matrix from triples format to "
                        "sparse format, column index %li equals %li.\n"
                        "Since the indexing parameter for the matrix was "
                        "specified as fortran, the smallest column is 1.\n",
                        (LONG) l,  (LONG) cols [l]) ;
            }
            else
            {
                printf ("\nWhen converting a matrix from triples format to "
                        "sparse format, column index %li equals %li < 0.\n",
                        (LONG) l,  (LONG) cols [l]) ;
            }
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
        else if ( j >= ncol )
        {
            printf ("\nWhen converting a matrix from triples format to sparse\n"
                    "format, a column index %li exceeded the number of\n"
                    "columns %li of the input matrix.\n",
                    (LONG) j, (LONG) ncol) ;
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
        if ( fortran ) i = rows [l] - 1 ;
        else           i = rows [l] ;
        if ( i < 0 )
        {
            if ( fortran )
            {
                printf ("\nWhen converting a matrix from triples format to "
                        "sparse format, row index %li equals %li.\n"
                        "Since the indexing parameter for the matrix was "
                        "specified as fortran, the smallest column is 1.\n",
                        (LONG) l,  (LONG) rows [l]) ;
            }
            else
            {
                printf ("\nWhen converting a matrix from triples format to "
                        "sparse format, row index %li equals %li < 0.\n",
                        (LONG) l,  (LONG) rows [l]) ;
            }
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
        if ( i > imax ) imax = i ;
        if ( sym && (i != j) )
        {
            /* complete matrix by symmetry */
            ap [i]++ ;
            if ( j > imax ) imax = j ;
        }
        if ( (vals [l] == SOPTZERO) && no_zeros )
        {
            /* remove element in column i if matrix symmetric */
            if ( sym && (i != j) ) ap [i]-- ;
        }
        else ap [j]++ ; /* there is another element in column j */
    }
    imax++ ;
    if ( (nrow > 0) && (nrow < imax) )
    {
        printf ("When converting a matrix from triples format to\n"
                "sparse matrix format, it was found that the number of rows\n"
                "(%ld) specified for the matrix was strictly less than the\n"
                "number of rows (%ld) in the input triples.\n\n",
                (LONG) nrow, (LONG) imax) ;
        return (SOPT_ERROR_IN_INPUT_MATRIX) ;
    }
    else if ( nrow <= 0 ) /* if nrow not given, then use the found number */
    {
        nrow = A->nrow = imax ;
    }

    /* ap now contains the number of nonzeros in each column of A;
       change ap to pointers to start of each column */
    k = 0 ;
    for (j = 0; j < ncol; j++)
    {
        l = k + ap [j] ;
    /* flag column containing matrix elements: ap [j] <= 0 means that
       column j contains matrix elements, but none have been put there yet */
        if ( l > k ) ap [j] = -(k+1) ;
        else         ap [j] =  k ;
        k = l ;
    }
    knnz = ap [ncol] = k ;
    if ( k == 0 )
    {
        sopt_free (ap) ;
        return (status) ; /* the matrix is all zero */
    }

    /* matrix   symmetric => build in B
       matrix unsymmetric => build in A */
    SOPTINT   *Bi = sopt_malloc (&status, knnz, sI) ;
    SOPTFLOAT *Bx = sopt_malloc (&status, knnz, sF) ;
    if ( sym )
    {
        Cp = Bp = ap ;
        Ci = Bi ;
        Cx = Bx ;
    }
    else /* unsymmetric */
    {
        A->i = ai = (SOPTINT *)   sopt_malloc (&status, knnz, sI) ;
        A->x = ax = (SOPTFLOAT *) sopt_malloc (&status, knnz, sF) ;
        Bp = sopt_malloc (&status, 1 + A->nrow, sI) ;
        Cp = ap ;
        Ci = ai ;
        Cx = ax ;
    }
    if ( status == SOPT_OUT_OF_MEMORY ) return (status) ;

    int col_out_of_order = 0 ; /* changes to 1 when column found with rows
                                  out-of-order */
    /* store entries in Ci and Cx */
    for (k = 0; k < nnz; k++)
    {
        SOPTFLOAT const Tx = vals [k] ;
        if ( !Tx && no_zeros ) continue ; /* skip zero matrix elements */
        if ( fortran )
        {
            /* convert Fortran to C indexing */
            i = rows [k] - 1 ;
            j = cols [k] - 1 ;
        }
        else
        {
            i = rows [k] ;
            j = cols [k] ;
        }
        l = Cp [j] ;
        if ( l < 0 ) /* no elements have been put in column so far */
        {
            l = -(l+1) ;
            Ci [l] = i ;
            Cp [j] = l + 1 ;
        }
        else
        {
            Cp [j]++ ;
            Ci [l] = i ;
            if ( order_cols )
            {
                if ( i <= Ci [l-1] )
                {
                    col_out_of_order = 1 ;
                }
            }
        }
        Cx [l] = Tx ;
        /* if symmetric matrix & full matrix output, swap i and j if i != j */
        if ( sym && (i != j) )
        {
            l = Cp [i] ;
            if ( l < 0 ) /* no elements have been put in column so far */
            {
                l = -(l+1) ;
                Ci [l] = j ;
                Cp [i] = l + 1 ;
            }
            else
            {
                Cp [i]++ ;
                Ci [l] = j ;
                if ( order_cols )
                {
                    if ( j <= Ci [l-1] )
                    {
                        col_out_of_order = 1 ;
                    }
                }
            }
            Cx [l] = Tx ;
        }
    }

    /* shift ap array up */
    for (j = ncol-1; j > 0; j--)
    {
        Cp [j] = Cp [j-1] ;
    }
    Cp [0] = 0 ;

    /* the matrix is now stored in C = B for symmetric matrix and
                                   C = A for unsymmetrix matrix */

    /* use transposes to sort the columns so that row indices
       are in increasing order when the columns are not yet sorted */
    if ( order_cols && col_out_of_order )
    {
        if ( sym )
        {
            A->i = (SOPTINT *)   sopt_malloc (&status, knnz, sI) ;
            A->x = (SOPTFLOAT *) sopt_malloc (&status, knnz, sF) ;
            status = sopt_sort_cols (Bp, Bi, Bx, A->p, A->i, A->x,
                                     ncol, ncol, sym) ;
        }
        else /* unsymmetric matrix, perform double transpose */
        {
            status = sopt_sort_cols (ap, ai, ax, Bp, Bi, Bx,
                                     nrow, ncol, sym);
        }
    }
    else if ( order_cols && sym )
    {
        /* Symmetric matrix stored in B has ordered column,
           copy B back to A. If the matrix is unsymmetric, then there
           is nothing to do since unsymmtric matrix in A has sorted cols. */
        A->x = Bx ;
        A->i = Bi ;
        return (status) ;
    }
    if ( !sym ) sopt_free (Bp) ;
    sopt_free (Bi) ;
    sopt_free (Bx) ;
    return (status) ;
}

/* ==========================================================================
   === sopt_convert_H_to_C_in_CG ============================================
   ==========================================================================
   Convert a symmetric matrix stored in H into a compressed matrix
   corresponding to the keep set of rows and columns for C.
   For fixed sparsity structure in the CG phase, we just need to
   malloc Cp, Ci, Cx once based on the size of H.  If there are
   no bounds or if cg_descent is stand alone, no malloc needed, just copy
   the pointers. */
int sopt_convert_H_to_C_in_CG /* returned integer:
                              0 if conversion successful
                              SOPT_OUT_OF_MEMORY
                              SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE
                              SOPT_ERROR_IN_INPUT_MATRIX */
(
    SOPT_Cmatrix         *C, /* compressed matrix */
    SOPT_matrix          *H, /* original matrix */
    int const       initial, /* TRUE if the initial iteration in CG */
    int const      location  /* = LPASA or LCG */
)
{
    SOPTINT i, j, k, l, m, p, nnz,
            *Bi, *Ccol, *Cmap, *order, *Fp, *Hp, *Hi, *Hrows, *Hcols, *temp ;
    SOPTFLOAT t, *dense, *Hx, *Hvals ;
    LONG li, le ;

    int const            sI = sizeof (SOPTINT) ;
    int const            sF = sizeof (SOPTFLOAT) ;
    SOPTINT const      ncol = H->ncol ;
    int const       fortran = H->fortran ;
    int const           sym = H->sym ;
    SOPTINT const         n = C->n ;
    int const        compress = (C->col == NULL) ? FALSE : TRUE ;
    int             Hformat = C->Hformat ;
    int const SparsityFixed = C->SparsityFixed ;
    int              status = 0 ;

    if ( Hformat == EMPTY )
    {
        if      ( H->p != NULL )       Hformat = 0 ;
        else if ( H->rows != NULL )    Hformat = 3 ;
        else if ( H->by_rows != NULL ) Hformat = 1 ;
        else if ( H->by_cols != NULL ) Hformat = 2 ;
        else
        {
            status = SOPT_HESSIAN_NOT_COMPUTED ;
        }
        C->Hformat = Hformat ;
    }

    if ( Hformat == 0 )
    {
        Hp = H->p ;
        Hi = H->i ;
        Hx = H->x ;
    }
    else if ( Hformat == 3 )
    {
        Hrows = H->rows ;
        Hcols = H->cols ;
        Hvals = H->vals ;
        /* If triples format is used & sym = TRUE, then only the elements on the
           main diagonal and on one side of the diagonal should be given. If
           fortran = TRUE, then indices should be in Fortran format where the
           first row and column are one. */
    }
    else if ( Hformat == 1 ) dense = H->by_rows ;
    else                     dense = H->by_cols ;

    if ( !Hformat && !compress ) /* CG, HFormat = 0 and no compression */
    {
        C->p = H->p ;
        C->i = H->i ;
        C->x = H->x ;
        return (status) ;
    }

    /* malloc space for Cp, Ci, Cx if necessary */
    if ( (C->p == NULL) || !SparsityFixed )
    {
        if      ( Hformat == 0 ) nnz = Hp [ncol] ;
        else if ( Hformat == 3 ) nnz = 2*H->nnz ;
        else /* dense matrix */
        {
            nnz = 0 ;
            li = 0 ;
            le = (LONG) ncol * (LONG) ncol ;
            nnz = 0 ;
            for (li = 0; li < le; li++)
            {
                if ( dense [li] ) nnz++ ; /* if matrix element nonzero, nnz++ */
            }
        }
        C->nnzfull = nnz ;
        if ( (C->p != NULL) && !SparsityFixed )
        {
            sopt_free (C->p) ;
            sopt_free (C->i) ;
            sopt_free (C->x) ;
            C->p = NULL ;
        }
        if ( C->p == NULL )
        {
            C->p = sopt_malloc (&status, ncol+1, sI) ;
            C->i = sopt_malloc (&status, nnz, sI) ;
            C->x = sopt_malloc (&status, nnz, sF) ;
            if ( status ) return (sopt_convert_error (location, status)) ;
            C->malloc_p = C->p ;
        }
    }
    SOPTINT   *Cp = C->p ;
    SOPTINT   *Ci = C->i ;
    SOPTFLOAT *Cx = C->x ;

    /* malloc space for Ccol and Cmap if necessary
       Ccol [i] = column in C associated with column i in H (EMPTY if i does
                  not appear in C
       Cmap [i] = index in Hx associated with Cx [i] */
    if ( (initial && (compress||(!compress && (C->map==NULL))))||!SparsityFixed)
    {
        if ( (C->map == NULL) && SparsityFixed )
        {
            if ( Hformat == 3 )
            {
                /* space for C->map & 2 other arrays, plus column point array */
                C->map = sopt_malloc (&status, 3*C->nnzfull + ncol + 1, sI) ;
            }
            else
            {
                C->map = sopt_malloc (&status, C->nnzfull, sI) ;
            }
            if ( status ) return (sopt_convert_error (location, status)) ;
        }

        Cmap = C->map ;
        Ccol = C->col ;
        order= C->order ;
        SOPTINT *Cpplus = Cp+1 ;
        Cp [0] = 0 ;
        if ( Hformat == 0 ) /* !compress and Hformat = 0 handled earlier */
        {
            /* bounds exist */
            for (j = 0; j < ncol; j++)
            {
                if ( (k = Ccol [j]) < 0 ) continue ;
                m = 0 ;
                SOPTINT const q = Hp [j+1] ;
                for (p = Hp [j]; p < q; p++)
                {
                    if ( Ccol [Hi [p]] >= 0 ) m++ ;
                }
                Cpplus [k] = m ;
            }
            l = 0 ;
            for (k = 1; k <= n; k++)
            {
                l += Cp [k] ;
                Cp [k] = l ;
            }
            /* Scatter the columns into the rows, starting with
               smallest column. Note that we also take the transpose
               of C (which equals C due to symmetry) in order to
               achieve sorted columns. */
            for (i = 0; i < n; i++)
            {
                j = order [i] ;  /* col in H corresponding to col i in C */
                SOPTINT const q = Hp [j+1] ;
                for (p = Hp [j]; p < q; p++)
                {
                    if ( (k = Ccol [Hi [p]]) >=  0 )
                    {
                        SOPTINT const qq = Cp [k]++ ;
                        Ci [qq] = i ;
                        if ( SparsityFixed ) Cmap [qq] = p ;
                        Cx [qq] = Hx [p] ;
                    }
                }
            }
            for (j = n-1; j > 0; j--)
            {
                Cp [j] = Cp [j-1] ;
            }
            Cp [0] = 0 ;
        }
        else if ( Hformat == 3 )
        {

            sopt_initi (Cp, (SOPTINT) 0, n+1) ;
            SOPTINT const knnz = H->nnz ;
            for (k = 0; k < knnz; k++)
            {
                if ( fortran )
                {
                    i = Hrows [k] - 1 ;
                    j = Hcols [k] - 1 ;
                }
                else
                {
                    i = Hrows [k] ;
                    j = Hcols [k] ;
                }
                if ( compress )
                {
                    if ( ((i = Ccol [i]) < 0)||((j = Ccol [j]) < 0) ) continue;
                }
                Cpplus [j]++ ;
                if ( sym && (i != j) ) Cpplus [i]++ ;
            }
            /* convert the row counts to row starts */
            k = 0 ;
            for (j = 1; j <= n; j++)
            {
                k += Cp [j] ;
                Cp [j] = k ;
            }
            /* scatter to B; transpose back to C to achieve sorted columns */
            if ( SparsityFixed )
            {
                Bi = Cmap+C->nnzfull ;
                temp = Bi+k ;
                Fp = temp+k ;
            }
            else /* Sparsity NOT fixed */
            {
                Bi = sopt_malloc (&status, k, sI) ;
                temp = sopt_malloc (&status, k, sI) ;
                Fp = sopt_malloc (&status, n+1, sI) ;
                if ( status ) return (sopt_convert_error (location, status)) ;
            }
            for (k = 0; k < knnz; k++)
            {
                if ( fortran )
                {
                    i = Hrows [k] - 1 ;
                    j = Hcols [k] - 1 ;
                }
                else
                {
                    i = Hrows [k] ;
                    j = Hcols [k] ;
                }
                if ( compress)
                {
                    if ( ((i = Ccol [i]) < 0)||((j = Ccol [j]) < 0) ) continue;
                }
                if ( i == j )
                {
                    Bi [Cp [i]] = i ;
                    temp [Cp [i]++] = k ;
                }
                else
                {
                    if ( sym )
                    {
                        Bi [Cp [i]] = j ;
                        Bi [Cp [j]] = i ;
                        temp [Cp [i]++] = temp [Cp [j]++] = k ;
                    }
                    else
                    {
                        Bi [Cp [j]] = i ;
                        temp [Cp [j]++] = k ;
                    }
                }
            }
            /* shift Cp array up */
            for (j = n-1; j > 0; j--)
            {
                Cp [j] = Cp [j-1] ;
            }
            Cp [0] = 0 ;

            /* take transpose to obtain C with sorted columns */
            p = 0 ;
            sopt_copyi (Fp, Cp, n+1) ;
            for (j = 0; j < n; j++)
            {
                SOPTINT const q = Cp [j+1] ;
                for (; p < q; p++)
                {
                    k = Fp [Bi [p]]++ ;
                    Ci [k] = j ;
                    l = temp [p] ;
                    Cx [k] = Hvals [l] ;
                    if ( SparsityFixed ) Cmap [k] = l ;
                }
            }

            if ( !SparsityFixed )
            {
                sopt_free (Bi) ;
                sopt_free (temp) ;
                sopt_free (Fp) ;
            }
        }
        else /* dense H */
        {
            nnz = l = k = 0 ;
            Cp [0] = 0 ; 
            nnz = 0 ;
            for (j = 0; j < ncol; j++)
            {
                if ( compress )
                {
                    if ( Ccol [j] >= 0 )
                    {
                        for (i = 0; i < ncol; i++) 
                        {
                            if ( (l = Ccol [i]) >= 0 )
                            {
                                /* m = i+k = location in dense matrix */
                                if ( t = dense [m = (i+k)] )
                                {
                                    Cx [nnz] = t ;
                                    if ( SparsityFixed ) Cmap [nnz] = m ;
                                    Ci [nnz] = l ;
                                    nnz++ ;
                                }
                            }
                        }
                        k += ncol ;
                        l++ ;       /* column number in compressed matrix */
                        Cp [l] = nnz ;
                    }
                }
                else /* no bounds */
                {
                    for (i = 0; i < ncol; i++)
                    {
                        if ( dense [l = (i+k)] )
                        {
                            Cx [nnz] = dense [l] ;
                            if ( SparsityFixed ) Cmap [nnz] = l ;
                            Ci [nnz] = i ;
                            nnz++ ;
                        }
                    }
                    Cp [j+1] = nnz ;
                    k += ncol ;
                }
            }
        } /* end of dense H */
        return (status) ;
    }

    /* evaluate Cx */
    if      ( Hformat == 0 ) Hx = H->x ;
    else if ( Hformat == 3 ) Hx = H->vals ;
    else if ( Hformat == 1 ) Hx = H->by_rows ;
    else                     Hx = H->by_cols ;
    SOPTINT const q = Cp [n] ;
    Cmap = C->map ;
    for (p = 0; p < q; p++)
    {
        Cx [p] = Hx [Cmap [p]] ;
    }
    return (status) ;
}

/* ==========================================================================
   === sopt_convert_H_to_C_in_SS ============================================
   ==========================================================================
   Convert a symmetric matrix stored in H into a compressed matrix C
   corresponding to the keep set of rows and columns for C. The output
   matrix is in triples format.  */
int sopt_convert_H_to_C_in_SS /* returned integer:
                              0 if conversion successful
                              SOPT_OUT_OF_MEMORY */
(
    SOPT_Cmatrix       *C, /* compressed matrix */
    SOPT_matrix        *H, /* original matrix */
    int const     initial, /* TRUE if the initial iteration in SS */
    int const        Cdim, /* max C dim = ncol + nrow */
    int const        Annz, /* max nnz needed for A in SS */
    int const    location  /* = LPASA or LCG */
  
)
{
    SOPTINT i, j, k, l, m, nnz, p,
            *Ccol, *Smap, *Hp, *Hi, *Hrows, *Hcols ;
    SOPTFLOAT t, *dense, *Hx, *Hvals ;
    LONG li, le ;
    int const            sI = sizeof (SOPTINT) ;
    int const            sF = sizeof (SOPTFLOAT) ;
    SOPTINT const      ncol = H->ncol ;
    int const       fortran = H->fortran ;
    int const           sym = H->sym ;
    int const      compress = (C->col == NULL) ? FALSE : TRUE ;
    int             Hformat = C->Hformat ;
    int const SparsityFixed = C->SparsityFixed ;
    int const       Aexists = C->Aexists ;
    int              status = 0 ;
    if ( Hformat == EMPTY )
    {
        if      ( H->p != NULL )       Hformat = 0 ;
        else if ( H->rows != NULL )    Hformat = 3 ;
        else if ( H->by_rows != NULL ) Hformat = 1 ;
        else                           Hformat = 2 ;
        C->Hformat = Hformat ;
    }

    if ( Hformat == 0 )
    {
        Hp = H->p ;
        Hi = H->i ;
        Hx = H->x ;
    }
    else if ( Hformat == 3 )
    {
        Hrows = H->rows ;
        Hcols = H->cols ;
        Hvals = H->vals ;
        /* If triples format is used & sym = TRUE, then only the elements on the
           main diagonal and on one side of the diagonal should be given. If
           fortran = TRUE, then indices should be in Fortran format where the
           first row and column are one. */
    }
    else if ( Hformat == 1 ) dense = H->by_rows ;
    else                     dense = H->by_cols ;

    /* malloc memory if necessary (determine max number of nnz for
       the triples*/
    if ( (C->rows == NULL) || !SparsityFixed )
    {
        if ( (Hformat == 3) && sym )
        {
            nnz = H->nnz ;
        }
        else
        {
            /* obtain nnz in sparse matrix format, then convert to triples */
            if ( Hformat == 0 )
            {
                nnz = Hp [ncol] ;
            }
            else if ( Hformat == 3 ) /* sym = FLASE */
            {
                nnz = H->nnz ; /* entire matrix is stored */
            }
            else /* dense matrix */
            {
                nnz = 0 ;
                le = (LONG) ncol * (LONG) ncol ;
                nnz = 0 ;
                for (li = 0; li < le; li++)
                {
                    if ( dense [li] ) nnz++ ;
                }
            }
            nnz = (ncol + nnz)/2 ; /* bound for nnz in triples format */
        }
        C->nnzsym = nnz ;
        /* malloc space for triples */
        if ( (C->rows != NULL) && !SparsityFixed )
        {
            sopt_free (C->rows) ;
            sopt_free (C->cols) ;
            sopt_free (C->vals) ;
            C->rows = NULL ;
        }
        if ( C->rows == NULL )
        {
            /* allocate for nnz of Hessian + constraint A + diagonal elements */
            li = nnz + Cdim + Annz ;
            C->rows = sopt_lmalloc (&status, li, sI) ;
            C->cols = sopt_lmalloc (&status, li, sI) ;
            C->vals = sopt_lmalloc (&status, li, sF) ;
            if ( status ) return (sopt_convert_error (location, status)) ;
            C->malloc_rows = C->rows ;
        }
    }
    SOPTINT   *rows = C->rows ;
    SOPTINT   *cols = C->cols ;
    SOPTFLOAT *vals = C->vals ;

    /* an easy special case: */
    if ( !Aexists && (Hformat == 3) && !compress )
    {
        nnz = H->nnz ;
        if ( initial )
        {
            if ( fortran )
            {
                sopt_copyi (rows, Hrows, nnz) ;
                sopt_copyi (cols, Hcols, nnz) ;
            }
            else /* add 1 to rows and columns to get Fortran indexing */
            {
                sopt_add2i (rows, Hrows, (SOPTINT) 1, nnz) ;
                sopt_add2i (cols, Hcols, (SOPTINT) 1, nnz) ;
            }
            C->nnz  = nnz ;
        }
        sopt_copyx (vals, Hvals, nnz) ;
        return (status) ;
    }

    /* malloc space for Smap if necessary
       Ccol [i] = column in C associated with column i in H (EMPTY if i does
                  not appear in C (set up at start of cg_descent)
       Smap [i] = index in Hx associated with Cx [i] (set up here) */
    if ( (initial&&(compress||(!compress && (C->Smap==NULL))))|| !SparsityFixed)
    {
        if ( C->Smap == NULL && SparsityFixed )
        {
            C->Smap = sopt_malloc (&status, C->nnzsym, sI) ;
            if ( status ) return (sopt_convert_error (location, status)) ;
        }
        Smap = C->Smap ;
        Ccol = C->col ;
        nnz = 0 ;
        /* depending on the format of H, the conversion to the triples of
           of C is different */
        if ( Hformat == 0 )
        {
            p = 0 ;
            for (j = 0; j < ncol; j++)
            {
                if ( compress )
                {
                    if ( (k = Ccol [j]) == EMPTY ) continue ;
                    SOPTINT const q = Hp [j+1] ;
                    for (p = Hp [j]; p < q; p++)
                    {
                        if ( ((i = Ccol [Hi [p]]) >= 0) && (i <= k) )
                        {
                            rows [nnz] = i + 1 ;
                            cols [nnz] = k + 1 ;
                            vals [nnz] = Hx [p] ;
                            if ( SparsityFixed ) Smap [nnz] = p ;
                            nnz++ ;
                        }
                    }
                }
                else /* no compression */
                {
                    SOPTINT const q = Hp [j+1] ;
                    for (; p < q; p++)
                    {
                        if ( (i = Hi [p]) <= j )
                        {
                            rows [nnz] = i + 1 ;/* adjust for Fortran indexing*/
                            cols [nnz] = j + 1 ;
                            vals [nnz] = Hx [p] ;
                            if ( SparsityFixed ) Smap [nnz] = p ;
                            nnz++ ;
                        }
                    }
                }
            }
        }
        else if ( Hformat == 3 )
        {
            l = H->nnz ;
            for (k = 0; k < l; k++)
            {
                if ( fortran )
                {
                    i = Hrows [k] - 1 ;
                    j = Hcols [k] - 1 ;
                }
                else
                {
                    i = Hrows [k] ;
                    j = Hcols [k] ;
                }
                if ( compress )
                {
                    if ( ((i = Ccol [i]) < 0)||((j = Ccol [j]) < 0) ) continue;
                }
                if ( sym || (j >= i) )
                {
                    rows [nnz] = i + 1 ;
                    cols [nnz] = j + 1 ;
                    vals [nnz] = Hvals [k] ;
                }
                else
                {
                    continue ;
                }
                if ( SparsityFixed )
                {
                    Smap [nnz] = k ;
                }
                nnz++ ;
            }
        }
        else /* dense H */
        {
            m = 0 ;
            for (j = 0; j < ncol; j++)
            {
                if ( compress )
                {
                    if ( (k = Ccol [j]) >= 0 )
                    {
                        for (i = 0; i < ncol; i++) 
                        {
                            if ( ((l = Ccol [i]) >= 0) && (l <= k)
                                                       && (t = dense[m+i])  )
                            {
                                    vals [nnz] = t ;
                                    if ( SparsityFixed ) Smap [nnz] = m + i ;
                                    rows [nnz] = l ;
                                    cols [nnz] = k ;
                                    nnz++ ;
                            }
                        }
                    }
                }
                else /* no bounds */
                {
                    for (i = 0; i < ncol; i++)
                    {
                        if ( t = dense [l = (i+m)] )
                        {
                            vals [nnz] = t ;
                            if ( SparsityFixed ) Smap [nnz] = l ;
                            rows [nnz] = i + 1 ;/* adjust for Fortran */
                            cols [nnz] = j + 1 ;
                            nnz++ ;
                        }
                    }
                }
                m += ncol ;
            }
        } /* end of dense H */
        C->nnz = nnz ;
        return (status) ;
    }

    /* C->rows and C->cols have been determined, evaluate vals */
    if      ( Hformat == 0 ) Hx = H->x ;
    else if ( Hformat == 3 ) Hx = H->vals ;
    else if ( Hformat == 1 ) Hx = H->by_rows ;
    else                     Hx = H->by_cols ;
    SOPTINT const q = C->nnz ;
    Smap = C->Smap ;
    for (p = 0; p < q; p++)
    {
        vals [p] = Hx [Smap [p]] ;
    }
    return (status) ;
}

/* ==========================================================================
   === sopt_check_matrix ====================================================
   ==========================================================================
    Check a sparse matrix for consistency. The following checks are
    implemented: that the row indices in each column are strictly increasing,
    that Ap [0] = 0 and Ap [j] <= Ap [j+1], and that the elements
    in Ax are all nonzero.
   ========================================================================== */
int sopt_check_matrix /* returned integer:
                         0 if no error was detected
                         SOPT_ERROR_IN_INPUT_MATRIX */
(
    SOPTINT   const  *Ap, /* column pointers */
    SOPTINT   const  *Ai, /* row indices */
    SOPTFLOAT const  *Ax, /* numerical entries */
    SOPTINT   const ncol  /* number of columns in matrix */
)
{
    SOPTINT j, p ;

    /* if there is no matrix, then return immediately */
    if ( (Ap == NULL) || (Ai == NULL) || (Ax == NULL) || (ncol == 0) )
    {
        return (0) ;
    }
    /* check that row indices in each column are in increasing order */
    if ( Ap [0] != 0 )
    {
        printf ("in check_matrix: Ap [0] != 0\n") ;
        return (SOPT_ERROR_IN_INPUT_MATRIX) ;
    }
    p = 0 ;
    for (j = 0; j < ncol; j++)
    {
        SOPTINT const q = Ap [j+1] ;
        if ( p > q )
        {
            printf ("\nIn check_matrix: Ap [%ld] > Ap [%ld]\n",
                     (LONG) j, (LONG) j+1) ;
        }
        if ( q > p )
        {
            SOPTINT oldAi = Ai [p] ;
            for (p++; p < q; p++)
            {
                if ( Ai [p] <= oldAi )
                {
                    printf ("\nIn check_matrix, the row indices for\n"
                    "column %ld are not strictly increasing. Either a\n"
                    "row index repeats or the row indices for the column\n"
                    "are not sorted in increasing order.\n\n", (LONG) j) ;
                    return (SOPT_ERROR_IN_INPUT_MATRIX) ;
                }
                oldAi = Ai [p] ;
            }
        }
    }
    SOPTINT const nnz = Ap [ncol] ;
    for (j = 0; j < nnz; j++)
    {
        if ( Ax [j] != SOPTZERO )
        {
            printf ("\nIn check_matrix, a matrix element in column %ld"
                    "was zero.\n", (LONG) j) ;
            return (SOPT_ERROR_IN_INPUT_MATRIX) ;
        }
    }
    return (0) ;
}

/* ====================================================================== */
/* === sopt_minsortx ==================================================== */
/* ======================================================================
       ________________________________________________________
      |                                                        |
      |       sort a SOPTFLOAT array in increasing order       |
      |                                                        |
      |    input:                                              |
      |                                                        |
      |         x     --SOPTFLOAT array of numbers of length n |
      |         w     --SOPTINT working array of length n      |
      |         n     --number of array elements to sort       |
      |                                                        |
      |    output:                                             |
      |                                                        |
      |         x     --original array (SOPTFLOAT *) length n  |
      |         y     --indices of x giving increasing order   |
      |                 (SOPTINT *) of length n                |
      |________________________________________________________| */

void sopt_minsortx
(
    SOPTINT         *y, /* n-by-1 (output) */
    SOPTFLOAT const *x, /* n-by-1 (input not modified) */
    SOPTINT         *w, /* n-by-1, (input, working array) */
    SOPTINT          n  /* number of elements to sort */
)
{
    SOPTINT *yi, *wi, i, j, k, l, m, p, q ;
    SOPTFLOAT s, t ;

    y [0] = 0 ;
    if ( n < 2 ) return ;
    if ( n < 3 )
    {
        if ( x [0] > x [1] )
        {
            y [0] = 1 ;
            y [1] = 0 ;
        }
        else y [1] = 1 ;
        return ;
    }

    j = k = 0 ;
    for (i = 1; i < n; i++)
    {
        if ( x [i] < x [j] )
        {
            w [k] = i ;
            k = i ;
        }
        y [i] = j = i ;
    }

    w [k] = n ;
    while ( k > 0 )
    {
        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = y [i] ;
            s = x [p] ;
            j = w [i] ;
            k = j ;
            if ( j == n )
            {
                y [i] = j ;
                l = j ;
                w [m] = p ;
                k += m - i ;
                yi = y+(i-m) ;
                for (m++; m < k; m++) w [m] = yi [m] ;
            }
            else
            {
                q = y [j] ;
                t = x [q] ;
                l = w [j] ;
                y [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        w [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            w [m] = p ;
                            k += m - i ;
                            yi = y+(i-m) ;
                            for (m++; m < k; m++) w [m] = yi [m] ;
                            break ;
                        }
                        q = y [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        w [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            w [m] = q ;
                            k = m + l - j ;
                            yi = y+(j-m) ;
                            for (m++; m < k; m++) w [m] = yi [m] ;
                            break ;
                        }
                        p = y [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n )
        {
            for (i = 0; i < n; i++) y [i] = w [i] ;
            return ;
        }

        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = w [i] ;
            s = x [p] ;
            j = y [i] ;
            k = j ;
            if ( j == n )
            {
                w [i] = j ;
                l = j ;
                y [m] = p ;
                k += m - i ;
                wi = w+(i-m) ;
                for (m++; m < k; m++) y [m] = wi [m] ;
            }
            else
            {
                q = w [j] ;
                t = x [q] ;
                l = y [j] ;
                w [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        y [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            y [m] = p ;
                            k += m - i ;
                            wi = w+(i-m) ;
                            for (m++; m < k; m++) y [m] = wi [m] ;
                            break ;
                        }
                        q = w [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        y [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            y [m] = q ;
                            k = m + l - j ;
                            wi = w+(j-m) ;
                            for (m++; m < k; m++) y [m] = wi [m] ;
                            break ;
                        }
                        p = w [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n ) return ;
    }
}

/* ====================================================================== */
/* === sopt_minsorti ==================================================== */
/* ======================================================================
       ________________________________________________________
      |                                                        |
      |       sort an SOPTINT array in increasing order        |
      |                                                        |
      |    input:                                              |
      |                                                        |
      |         x  --array of numbers (SOPTINT *) of length n  |
      |         w  --working array (SOPTINT *) of length n     |
      |         n  --number of array elements to sort          |
      |                                                        |
      |    output:                                             |
      |                                                        |
      |         x  --original array (SOPTINT *) of length n    |
      |         y  --indices of x giving increasing order      |
      |              (SOPTINT *) of length n                   |
      |________________________________________________________| */

void sopt_minsorti
(
    SOPTINT        *y, /* n-by-1 (output) */
    SOPTINT  const *x, /* n-by-1 (input not modified) */
    SOPTINT        *w, /* n-by-1, (input, working array) */
    SOPTINT         n  /* number of elements to sort */
)
{
    SOPTINT *yi, *wi, i, j, k, l, m, p, q ;
    SOPTINT s, t ;

    y [0] = 0 ;
    if ( n < 2 ) return ;

    j = k = 0 ;
    for (i = 1; i < n; i++)
    {
        if ( x [i] < x [j] )
        {
            w [k] = i ;
            k = i ;
        }
        y [i] = j = i ;
    }

    w [k] = n ;
    while ( k > 0 )
    {
        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = y [i] ;
            s = x [p] ;
            j = w [i] ;
            k = j ;
            if ( j == n )
            {
                y [i] = j ;
                l = j ;
                w [m] = p ;
                k += m - i ;
                yi = y+(i-m) ;
                for (m++; m < k; m++)
                {
                    w [m] = yi [m] ;
                }
            }
            else
            {
                q = y [j] ;
                t = x [q] ;
                l = w [j] ;
                y [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        w [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            w [m] = p ;
                            k += m - i ;
                            yi = y+(i-m) ;
                            for (m++; m < k; m++)
                            {
                                w [m] = yi [m] ;
                            }
                            break ;
                        }
                        q = y [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        w [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            w [m] = q ;
                            k = m + l - j ;
                            yi = y+(j-m) ;
                            for (m++; m < k; m++)
                            {
                                w [m] = yi [m] ;
                            }
                            break ;
                        }
                        p = y [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n )
        {
            for (i = 0; i < n; i++)
            {
                y [i] = w [i] ;
            }
            return ;
        }

        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = w [i] ;
            s = x [p] ;
            j = y [i] ;
            k = j ;
            if ( j == n )
            {
                w [i] = j ;
                l = j ;
                y [m] = p ;
                k += m - i ;
                wi = w+(i-m) ;
                for (m++; m < k; m++)
                {
                    y [m] = wi [m] ;
                }
            }
            else
            {
                q = w [j] ;
                t = x [q] ;
                l = y [j] ;
                w [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        y [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            y [m] = p ;
                            k += m - i ;
                            wi = w+(i-m) ;
                            for (m++; m < k; m++)
                            {
                                y [m] = wi [m] ;
                            }
                            break ;
                        }
                        q = w [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        y [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            y [m] = q ;
                            k = m + l - j ;
                            wi = w+(j-m) ;
                            for (m++; m < k; m++)
                            {
                                y [m] = wi [m] ;
                            }
                            break ;
                        }
                        p = w [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n ) return ;
    }
}
