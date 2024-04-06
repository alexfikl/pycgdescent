#ifndef _SuiteOPT_H_
#define _SuiteOPT_H_
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"
#endif

/* error macro */
#ifdef MATLAB_MEX_FILE
#define SOPTERROR(msg) \
{ \
    mexPrintf ("Error in %s, line %d\n", __FILE__, __LINE__) ; \
    mexErrMsgTxt (msg) ; \
}
#else
#define SOPTERROR(msg) \
{ \
    fprintf (stdout,"Error in %s, line %d (%s)\n", __FILE__, __LINE__,msg) ; \
    fflush (stdout) ; \
    ASSERT (0) ; \
    abort () ; \
}
#endif

/* status values returned by codes in SOPT */

#define SOPT_SPACE_MALLOCED                                      (-1)
#define SOPT_EXISTING_MATRIX_SPARSE                              (-2)
#define SOPT_NO_MATRIX_PROVIDED                                  (-3)
#define SOPT_SPARSE_MATRIX_CREATED                               (-4)
#define LSOPT                                                  (-100)
#define LPASA                                                  (-101)
#define LCG                                                    (-102)
#define LPPROJ                                                 (-103)
#define LNAPHEAP                                               (-104)
#define LCUTE                                                  (-105)
#define LSSM                                                   (-106)

#define SOPT_OUT_OF_MEMORY                                      (901)
#define SOPT_ERROR_IN_INPUT_MATRIX                              (902)
#define SOPT_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE                   (903)
#define SOPT_TRIPLES_FORMAT_ERROR                               (904)
#define SOPT_HESSIAN_NOT_COMPUTED                               (905)
#define SOPT_START_MESSAGES                                     (900)
#define SOPT_END_MESSAGES                                       (999)

/* status values returned by PASA */

#define PASA_ERROR_TOLERANCE_SATISFIED                            (0)
#define PASA_ERROR_DECAY_STAGNATES                                (1)
#define PASA_POLYHEDRON_INFEASIBLE                                (2)
#define PASA_INVALID_VARIABLE_BOUNDS                              (3)
#define PASA_INVALID_LINEAR_CONSTRAINT_BOUNDS                     (4)
#define PASA_INVALID_MATRIX_ELEMENT                               (5)
#define PASA_ITERATIONS_EXCEED_MAXITS_IN_GRAD_PROJ                (6)
#define PASA_ITERATIONS_EXCEED_MAXITS_IN_ACTIVE_GRAD_PROJ         (7)
#define PASA_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS_IN_GRAD_PROJ       (8)
#define PASA_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS_IN_ACTIVEGP        (9)
#define PASA_OUT_OF_MEMORY                                       (10)
#define PASA_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION_IN_GRAD_PROJ (11)
#define PASA_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION_IN_ACTIVEGP  (12)
#define PASA_MUST_USE_CHOLMOD                                    (13)
#define PASA_STARTING_FUNCTION_VALUE_INFINITE_OR_NAN             (14)
#define PASA_FUNCTION_NAN_OR_INF                                 (15)
#define PASA_LAMBDA_IS_NULL_BUT_USE_LAMBDA_IS_TRUE               (16)
#define PASA_MATRIX_INCOMPLETE                                   (17)
#define PASA_ERROR_IN_INPUT_MATRIX                               (18)
#define PASA_FUNCTION_VALUE_OR_GRADIENT_MISSING                  (19)
#define PASA_MATRIX_GIVEN_BUT_RHS_MISSING                        (20)
#define PASA_RHS_GIVEN_BUT_MATRIX_MISSING                        (21)
#define PASA_MISSING_OBJECTIVE                                   (22)
#define PASA_PROBLEM_OVERSPECIFIED                               (23)
#define PASA_BOTH_A_AND_A_EXIST                                  (24)
#define PASA_MISSING_HESSIAN_FOR_QP                              (25)
#define PASA_MISSING_HESSIAN_FOR_CG_QP                           (26)
#define PASA_PROBLEM_DIMENSION_NOT_GIVEN                         (27)
#define PASA_TOO_MANY_ASCENT_DIRECTIONS                          (28)
#define PASA_FACTORIZATION_FAILS_IN_CHOLMOD                      (29)
#define PASA_INTEGER_OVERFLOW_IN_CHOLMOD                         (30)
#define PASA_INVALID_DERIV_MODE_PARAMETER                        (31)
#define PASA_DERIV_MODE_USES_HESSIAN_BUT_NO_HESSIAN_PROVIDED     (32)
#define PASA_HESSIAN_NOT_COMPUTED                                (33)
#define PASA_HPROD_PLUS_HESSIAN                                  (34)
#define PASA_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE                    (35)
#define PASA_TRIPLES_FORMAT_ERROR                                (36)
#define PASA_START_MESSAGES                                       (0)
#define PASA_END_MESSAGES                                        (99)

/* status values returned by PPROJ */

#define PPROJ_SOLUTION_FOUND                                      (0)
#define PPROJ_ERROR_DECAY_STAGNATES                             (100)
#define PPROJ_OUT_OF_MEMORY                                     (101)
#define PPROJ_SSOR_NONASCENT                                    (102)
#define PPROJ_SSOR_MAX_ITS                                      (103)
#define PPROJ_MISSING_PRIOR_DATA                                (104)
#define PPROJ_START_GUESS_NEEDS_PRIOR_DATA                      (105)
#define PPROJ_INVALID_LINEAR_CONSTRAINT_BOUNDS                  (106)
#define PPROJ_DUAL_SOLVE_ERROR                                  (107)
#define PPROJ_START_GUESS_IS_1_BUT_LAMBDA_NULL                  (108)
#define PPROJ_START_GUESS_IS_2_BUT_CHOLMOD_FALSE                (109)
#define PPROJ_START_GUESS_IS_3_BUT_LAMBDA_NULL                  (110)
#define PPROJ_BOTH_NI_AND_NSING_POSITIVE                        (111)
#define PPROJ_OPTIMAL_COST_IS_MINUS_INFINITY                    (112)
#define PPROJ_NSING_START_GUESS_PROB                            (113)
#define PPROJ_FACTORIZATION_FAILS_IN_CHOLMOD                    (114)
#define PPROJ_INTEGER_OVERFLOW_IN_CHOLMOD                       (115)
#define PPROJ_OUT_OF_MEMORY_IN_CHOLMOD                          (116)
#define PPROJ_AN_ERROR_OCCURRED_IN_CHOLMOD                      (117)
#define PPROJ_ERROR_IN_INPUT_MATRIX                             (118)
#define PPROJ_PROBLEM_DIMENSION_NOT_GIVEN                       (119)
#define PPROJ_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE                  (120)
#define PPROJ_START_MESSAGES                                    (100)
#define PPROJ_END_MESSAGES                                      (199)

/* status values returned by CG_DESCENT */

#define CG_ERROR_TOLERANCE_SATISFIED                            (0)
#define CG_ITERATIONS_EXCEED_MAXITS                             (200)
#define CG_SLOPE_ALWAYS_NEGATIVE                                (201)
#define CG_LINE_SEARCH_STEPS_EXCEED_MAXSTEPS                    (202)
#define CG_SEARCH_DIRECTION_NOT_DESCENT_DIRECTION               (203)
#define CG_WOLFE_CONDITIONS_NOT_SATISFIED                       (204)
#define CG_DEBUGGER_IS_ON_AND_FUNCTION_VALUE_INCREASES          (205)
#define CG_NO_COST_OR_GRADIENT_IMPROVEMENT                      (206)
#define CG_OUT_OF_MEMORY                                        (207)
#define CG_QUADRATIC_OBJECTIVE_NO_LOWER_BOUND                   (208)
#define CG_STARTING_FUNCTION_VALUE_INFINITE_OR_NAN              (209)
#define CG_EXCESSIVE_UPDATING_OF_PERT_EPS                       (210)
#define CG_FUNCTION_NAN_OR_INF                                  (211)
#define CG_QP_LINEAR_TERM_GIVEN_BUT_HPROD_MISSING               (212)
#define CG_N_IS_EMPTY                                           (213)
#define CG_ERROR_IN_INPUT_MATRIX                                (214)
#define CG_MISSING_HESSIAN_FOR_QUADCOST                         (215)
#define CG_INVALID_DERIV_MODE_PARAMETER                         (216)
#define CG_DERIV_MODE_USES_HESSIAN_BUT_NO_HESSIAN_PROVIDED      (217)
#define CG_SYMMETRIC_SOLVER_FAILS                               (218)
#define CG_HESSIAN_NOT_COMPUTED                                 (219)
#define CG_HPROD_PLUS_HESSIAN                                   (220)
#define CG_VALUE_OR_GRAD_MISSING                                (221)
#define CG_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE                     (222)
#define CG_TRIPLES_FORMAT_ERROR                                 (223)
#define CG_MULTI_SOLVERS                                        (224)
#define CG_START_MESSAGES                                       (200)
#define CG_END_MESSAGES                                         (299)

/* status values returned by NAPHEAP */

#define NAPHEAP_STATUS_OK                                         (0)
#define NAPHEAP_STATUS_INFEASIBLE                               (300)
#define NAPHEAP_STATUS_UNBOUNDED                                (301)
#define NAPHEAP_STATUS_OUT_OF_MEMORY                            (302)
#define NAPHEAP_STATUS_INVALID_N                                (303)
#define NAPHEAP_STATUS_INVALID_BOUNDS                           (304)
#define NAPHEAP_STATUS_INVALID_B                                (305)
#define NAPHEAP_STATUS_INVALID_D                                (306)
#define NAPHEAP_STATUS_CONFLICTING_PARM_D                       (307)
#define NAPHEAP_STATUS_MISSING_PARM                             (308)
#define NAPHEAP_STATUS_MISSING_STAT                             (309)
#define NAPHEAP_START_MESSAGES                                  (300)
#define NAPHEAP_END_MESSAGES                                    (399)

/* status values returned by SSM */

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

/* status values returned by CUTE */

#define CUTE_OUT_OF_MEMORY                                      (801)
#define CUTE_ERROR_IN_INPUT_MATRIX                              (802)
#define CUTE_NROW_OR_NCOL_NOT_GIVEN_FOR_DENSE                   (803)
#define CUTE_START_MESSAGES                                     (800)
#define CUTE_END_MESSAGES                                       (899)

/* some constants used by SOPT */
#define SuiteOPTfalse 0
#define SuiteOPTtrue 1

#ifndef NULL
#define NULL 0
#endif

/* by default, sopt_timer (clock_gettime) is enabled */
#define SUITEOPT_TIMER

/* define the long version of integers */
#ifndef SOPTLONG

#ifdef _WIN64

#define SOPTLONG __int64
#define SOPT_LONG_MAX _I64_MAX
#define LONG SOPTLONG

#else

#define SOPTLONG long
#define SOPT_LONG_MAX LONG_MAX
#define LONG SOPTLONG

#endif

#endif

/* define the integer precision for the BLAS */
#define BLAS_INT SOPTLONG

/* If the compiler flag DLONG in Userconfig.mk is defined, then the integers
   in SuiteOPT are LONG integers; otherwise the integers are simply int's. */

#ifdef DLONG

#define SOPTFLOAT double
#define SOPTINT SOPTLONG
#define SuiteOPTinfint SOPT_LONG_MAX
/* When using long its, need to include "_l" when calling CHOLMOD routines */
#define CHOLMOD(name) cholmod_l_ ## name

#else

#define SOPTFLOAT double
#define SOPTINT int
#define SuiteOPTinfint INT_MAX
/* No need to include "_l" when calling CHOLMOD routines */
#define CHOLMOD(name) cholmod_ ## name

#endif

#define EMPTY (SOPTINT) -1
#define TRUE SuiteOPTtrue
#define FALSE SuiteOPTfalse

#define SOPTZERO ((SOPTFLOAT) 0)

/* ANSI C99 has a clean definition of IEEE infinity, defined in <math.h>.
   MATLAB has this as well.  With the ANSI (C90) version of C, there is no
   well-defined infinity, so DBL_MAX is used instead.

   You can override these defaults and define your own version of infinity
   with (for example):

   cc -ansi -DSuiteOPTinf=1e200 pproj.c ...
*/

/* infinite float */
#ifndef inf
#ifdef INFINITY
/* ANSI C99 (gcc -std=c99) */
#define inf INFINITY
#else
/* ANSI C90 (gcc -ansi) */
#define inf DBL_MAX
#endif
#endif

#define SOPTMAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SOPTMIN(a,b) ( ((a) < (b)) ? (a) : (b) )

/* debugging options */
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE
#define ASSERT(expression) (mxAssert ((expression), ""))
#else
#define ASSERT(expression) (assert (expression))
#endif
#else
#define ASSERT(expression)
#endif

    /* There are three different ways to input the matrix using the
       SOPT_matrix structure:

       1. A dense packed array containing the matrix entries, either by rows
          or by columns. If the matrix is stored by rows, then A_by_rows
          should point to the matrix start; if the matrix is stored by
          columns, then A_by cols should point to the matrix start. Both the
          nrow and ncol elements of the data structure must also be given.

       2. The nonzero elements in the matrix can be stored as a triple:
              rows  (row    indices of nonzero elements)
              cols  (column indices of nonzero elements)
              vals  (nonzero matrix entries)
              nnz (number of nonzeros)
          NOTE: When using triples, Fortran/MATLAB format should be used,
                where the first row and column are 1.

       3. The matrix can be stored in standard sparse matrix format:
              p  (column pointers, the number of entries in Ap is 1 plus the
                  number of columns in A. Ap [j] is the location of the first
                  nonzero in Ax associated with column j
              i  (the associated row indices of the nonzero matrix elements, in
                  each column, the row indices should be in increasing order)
              x  (the numerical nonzero matrix elements) */

typedef struct SOPT_matrix_struct
{
    SOPTINT         nrow ; /* number of rows in the matrix (optional) */
    SOPTINT         ncol ; /* number of columns in the matrix (optional) */

    /* NOTE: The inputs nrow and ncol may be optional when their
             values can be deduced from inputs to the problem data structure
             such as the number of linear constraints or the dimension of x  */

   /* dense matrix */
    SOPTFLOAT   *by_rows ; /* numerical entries in A by rows */
    SOPTFLOAT   *by_cols ; /* numerical entries in A by cols */

    /* triples format */
    SOPTINT        *rows ; /* size nnz (row indices) */
    SOPTINT        *cols ; /* size nnz (column indices) */
    SOPTFLOAT      *vals ; /* size nnz (nonzero matrix entries) */
    SOPTINT          nnz ; /* number of nonzero entries */
    int              sym ; /* TRUE  => matrix is symmetric, only elements on
                                       main diagonal & one side given (default)
                              FALSE => all nonzero matrix elements are given */
    int          fortran ; /* TRUE  => Fortran indexing is used where the first
                                       row and column are 1 not 0 (default)
                              FALSE => first row and column are zero */

    /* sparse matrix format */
    SOPTINT           *p ; /* size ncol+1, column pointers */
    SOPTINT           *i ; /* size p [ncol], row indices for sparse matrix */
    SOPTFLOAT         *x ; /* size p [ncol], numerical entries of matrix */

    SOPTINT       nnzmax ; /* not specified by user (SuiteOPT bookkeeping) */
} SOPT_matrix ;

/* A Cmatrix is an internal matrix structure, not intended for the user.
   A Cmatrix is a compression of a symmetric matrix.  The format of the
   compression can be either triples or sparse matrix. In the triples format,
   only elements on the main diagonal and above are retained. The compressed
   format is obtained by only retaining the rows and columns associated
   with the free set of variables (the keep set).  */
typedef struct SOPT_Cmatrix_struct
{
    SOPTINT             n ; /* number of columns in the matrix */
    SOPTINT       nnzfull ; /* max number nnz in C if sparsity fixed and
                               element above and below diagonal stored */
    SOPTINT        nnzsym ; /* max number nnz in C if sparsity fixed and
                               only elements on diagonal and above stored */
    /* triples format */
    SOPTINT         *rows ; /* size nnz (row indices) */
    SOPTINT         *cols ; /* size nnz (column indices) */
    SOPTFLOAT       *vals ; /* size nnz (nonzero matrix entries) */
    SOPTINT           nnz ; /* number of nonzero entries */

    /* sparse matrix format */
    SOPTINT            *p ; /* size ncol+1, column pointers */
    SOPTINT            *i ; /* size p [ncol], row indices for sparse matrix */
    SOPTFLOAT          *x ; /* size p [ncol], numerical entries of matrix */

    /* information concerning the compression: */
    SOPTINT          *col ; /* col [j] = column of C corresponding to column
                               column j of H, EMPTY <=> column deleted */
    SOPTINT        *order ; /* order [i] = column of H corresponding to
                               column i of C (basically the inverse of col) */
    SOPTINT          *map ; /* if the sparsity is fixed, then map is the map
                               from the elements of the compressed matrix
                               to the elements of the original H */
    SOPTINT         *Smap ; /* Smap is similar to map except that Smap is
                               used in SS while map is used in CG */
    int     SparsityFixed ; /* TRUE => same location in H of Hessian nonzeros
                               for each x */
    int           Aexists ; /* TRUE => linear constraints present */
    int           Hformat ; /* = 0, 1, 2, or 3 matrix format */
    SOPTINT     *malloc_p ; /* malloc'd pointer for C->p */
    SOPTINT  *malloc_rows ; /* malloc'd pointer for C->rows */
    SOPTFLOAT      *Hwork ; /* workspace for CG PRP+ */
    SOPTFLOAT        *rhs ; /* workspace right side of symmetric solve */
} SOPT_Cmatrix ;

/* By default, the BLAS are not used in SuiteOPT (they are used in
   SuiteSparse). To use the BLAS in SuiteOPT, comment out the next
   line and then below, specify whether your BLAS routines include an
   underscore, and set the threshhold dimension for using the BLAS.  */
#ifdef SUITEOPT_DISABLE_BLAS
#define NOBLAS
#endif

#ifndef NOBLAS

/* when BLAS are used, comment out the next statement if there is no
   underscore in the subroutine names */
#ifdef SUITEOPT_BLAS_UNDERSCORE
#define BLAS_UNDERSCORE
#endif

#ifdef BLAS_UNDERSCORE

#define SOPT_DGEMV dgemv_
#define SOPT_DAXPY daxpy_
#define SOPT_DDOT ddot_
#define SOPT_DSCAL dscal_
#define SOPT_DCOPY dcopy_
#define SOPT_IDAMAX idamax_

#else

#define SOPT_DGEMV dgemv
#define SOPT_DAXPY daxpy
#define SOPT_DDOT ddot
#define SOPT_DSCAL dscal
#define SOPT_DCOPY dcopy
#define SOPT_IDAMAX idamax

#endif

/* define the starting size of vectors when BLAS are used */
#define MATVEC_START 1000
#define DAXPY_START  1000
#define DDOT_START   1000
#define DSCAL_START  1000
#define DCOPY_START  1000
#define IDAMAX_START 1000

/* protypes for the BLAS routines that could be used */
void SOPT_DGEMV (char *trans, BLAS_INT *m, BLAS_INT *n, SOPTFLOAT *alpha,
        SOPTFLOAT *A, BLAS_INT *lda, SOPTFLOAT *X, BLAS_INT *incx,
        SOPTFLOAT *beta, SOPTFLOAT *Y, BLAS_INT *incy) ;

void SOPT_DAXPY (BLAS_INT *n, SOPTFLOAT *DA, SOPTFLOAT *DX, BLAS_INT *incx,
        SOPTFLOAT *DY, BLAS_INT *incy) ;

SOPTFLOAT SOPT_DDOT (BLAS_INT *n, SOPTFLOAT *DX, BLAS_INT *incx, SOPTFLOAT *DY,
        BLAS_INT *incy) ;

void SOPT_DSCAL (BLAS_INT *n, SOPTFLOAT *DA, SOPTFLOAT *DX, BLAS_INT *incx) ;

void SOPT_DCOPY (BLAS_INT *n, SOPTFLOAT *DX, BLAS_INT *incx, SOPTFLOAT *DY,
        BLAS_INT *incy) ;

BLAS_INT SOPT_IDAMAX (BLAS_INT *n, SOPTFLOAT *DX, BLAS_INT *incx) ;

#endif

/* prototypes for codes in SIOPT.c */
void sopt_free
(
    void * p
) ;

void * sopt_malloc
(
    int *status,
    SOPTINT   n,
    int    size
) ;

void * sopt_lmalloc
(
    int *status,
    LONG      n,
    int    size
) ;

void sopt_matrix_default
(
    SOPT_matrix *A
) ;

void sopt_Cmatrix_default
(
    SOPT_Cmatrix *C
) ;

void sopt_error
(
    int status,
    const char *file,
    int line,
    const char *message
) ;

double sopt_timer ( void ) ;

void sopt_print_TF
(
    int TF /* TRUE or FALSE */
) ;

void sopt_printA
(
    SOPTINT  ncol, /* number of cols in A */
    SOPTINT   *Ap, /* size ncol+1, column pointers */
    SOPTINT  *Anz, /* if NULL, A is packed; otherwise gives # nonzeros in cols*/
    SOPTINT   *Ai, /* size Ap [ncol], row indices for A */
    SOPTFLOAT *Ax, /* size Ap [ncol], numerical entries of A */
    char    *what  /* name of the matrix */
) ;

void sopt_printAMATLAB
(
    SOPTINT   const ncol, /* number of cols in A */
    SOPTINT   const  *Ap, /* size ncol+1, column pointers */
    SOPTINT   const  *Ai, /* size Ap [ncol], row indices for A */
    SOPTINT   const *Anz, /* if NULL, A packed; otherwise # nonzeros in cols */
    SOPTFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    char           *what  /* name of the matrix */
) ;

void sopt_printx
(
    SOPTFLOAT const *x, /* numerical entries in the vector */
    SOPTINT   const  n, /* dimension of the vector */
    char          *what  /* name of the vector */
) ;

void sopt_printxMATLAB
(
    SOPTFLOAT *x, /* numerical entries in the vector */
    SOPTINT    n, /* dimension of the vector */
    char   *what  /* name of the vector */
) ;

void sopt_printi
(
    SOPTINT *i, /* array of SOPTINT's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of array */
) ;

void sopt_printiMATLAB
(
    SOPTINT *i, /* array of SOPTINT's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of the array */
) ;

void sopt_print_int
(
    int     *i, /* array of int's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of array */
) ;

void sopt_print_intMATLAB
(
    int     *i, /* array of int's */
    SOPTINT  n, /* dimension of i */
    char *what  /* name of the array */
) ;

void sopt_add
(
    SOPTFLOAT       *x,  /* array to which s is added */
    SOPTFLOAT const  s,  /* scalar */
    SOPTINT   const  n   /* length of x */
) ;

void sopt_addi
(
    SOPTINT       *x,  /* array to which s is added */
    SOPTINT const  s,  /* scalar */
    SOPTINT const  n   /* length of x */
) ;

void sopt_add2i
(
    SOPTINT        *x,  /* array to which y is stored  */
    SOPTINT  const *y,  /* array to which s is added */
    SOPTINT  const  s,  /* scalar */
    SOPTINT  const  n   /* length of x */
) ;

void sopt_scale
(
    SOPTFLOAT       *x,  /* array to be scaled */
    SOPTFLOAT const *y,  /* array used for the scaling */
    SOPTFLOAT const  s,  /* scale */
    SOPTINT   const  n   /* length of x */
) ;

SOPTFLOAT sopt_scale_max
(
    SOPTFLOAT       *x,  /* scaled array */
    SOPTFLOAT const *y,  /* array used for the scaling */
    SOPTFLOAT const  s,  /* scale */
    SOPTINT   const  n   /* length of x */
) ;

void sopt_step
(
    SOPTFLOAT       *xnew, /* updated x vector */
    SOPTFLOAT const    *x, /* current x */
    SOPTFLOAT const    *d, /* search direction */
    SOPTFLOAT const alpha, /* stepsize */
    SOPTINT   const     n  /* dimension */
) ;

SOPTFLOAT sopt_step_max
(
    SOPTFLOAT       *xnew, /* updated x vector */
    SOPTFLOAT const    *x, /* current x */
    SOPTFLOAT const    *d, /* search direction */
    SOPTFLOAT const alpha, /* stepsize */
    SOPTINT   const     n  /* dimension */
) ;

void sopt_daxpy
(
    SOPTFLOAT       *x, /* input and output vector */
    SOPTFLOAT const *d, /* direction vector */
    SOPTFLOAT const  s, /* stepsize */
    SOPTINT   const  n  /* length of the vectors */
) ;

void sopt_copyx
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    SOPTINT   const  n  /* length of vectors */
) ;

void sopt_lcopyx
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    LONG      const  n  /* length of vectors */
) ;

void sopt_copyx_noblas
(
    SOPTFLOAT       *x, /* output of copy */
    SOPTFLOAT const *y, /* input of copy */
    SOPTINT   const  n  /* length of vectors */
) ;

void sopt_copyi_noblas
(
    SOPTINT       *x, /* output of copy */
    SOPTINT const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
) ;

void sopt_copyi
(
    SOPTINT       *x, /* output of copy */
    SOPTINT const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
) ;

void sopt_copy_int
(
    int           *x, /* output of copy */
    int     const *y, /* input of copy */
    SOPTINT const  n  /* length of vectors */
) ;

SOPTFLOAT sopt_dot
(
    SOPTFLOAT const *x, /* first vector */
    SOPTFLOAT const *y, /* second vector */
    SOPTINT   const  n  /* length of vectors */
) ;

void sopt_initx
(
    SOPTFLOAT      *x,  /* array to be initialized */
    SOPTFLOAT const s,  /* scalar */
    SOPTINT   const n   /* length of x */
) ;

void sopt_initi
(
    SOPTINT      *x,  /* array to be initialized */
    SOPTINT const s,  /* scalar */
    SOPTINT const n   /* length of x */
) ;

void sopt_init_int
(
    int          *x,  /* array to be initialized */
    int     const s,  /* scalar */
    SOPTINT const n   /* length of x */
) ;

SOPTFLOAT sopt_sup_normx
(
    SOPTFLOAT const *x, /* vector */
    SOPTINT   const  n  /* length of vector */
) ;

SOPTINT sopt_supi
(
    SOPTINT const *x, /* vector */
    SOPTINT const  n  /* length of vector */
) ;

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
) ;

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
) ;

int sopt_convert_to_sparse /* returned integer:
                                   0 if conversion successful
                                   OUT_OF_MEMORY
                                   ERROR_IN_INPUT_MATRIX */
(
    SOPT_matrix       *A, /* input matrix */
    SOPTINT         nrow, /* number of rows in A (can be empty) */
    SOPTINT         ncol, /* number of cols in A (can be empty) */
    int const order_cols, /* TRUE if row indices in each column of the sparse
                             matrix format should be put in increasing order */
    int const   no_zeros, /* if an element of Tx is zero, then delete it */
    int         location  /* location (PASA, CGDESCENT, PPROJ) where routine
                             was invoked */
) ;

int sopt_convert_error
(
    int  location,
    int    status
) ;

int sopt_convert_dense_to_sparse /* returned integer:
                                  0 if conversion successful
                                  SOPT_OUT_OF_MEMORY
                                  SOPT_ERROR_IN_INPUT_MATRIX */
(
    SOPT_matrix *A
) ;

int sopt_convert_triple_to_sparse /* returned integer:
                                   0 if conversion successful
                                   SOPT_OUT_OF_MEMORY
                                   SOPT_ERROR_IN_INPUT_MATRIX */
(
    SOPT_matrix *A, /* matrix structure where sparse matrix is stored */
    int order_cols, /* = TRUE => order columns in increasing order */
    int   no_zeros  /* = TRUE => delete all zeros, only store nonzero entries */
) ;

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
) ;

int sopt_convert_H_to_C_in_SS /* returned integer:
                              0 if conversion successful
                              SOPT_OUT_OF_MEMORY */
(
    SOPT_Cmatrix       *C, /* compressed matrix */
    SOPT_matrix        *H, /* original matrix */
    int const     initial, /* TRUE if the initial iteration in SS */
    int const        Cdim, /* max C dim = ncol + nrow */
    int const        Annz, /* extra nnz needed for A in SS */
    int const    location  /* = LPASA or LCG */

) ;

int sopt_check_matrix /* return 1 if an error was detected, otherwise return 0*/
(
    SOPTINT   const  *Ap, /* column pointers */
    SOPTINT   const  *Ai, /* row indices */
    SOPTFLOAT const  *Ax, /* numerical entries */
    SOPTINT   const ncol  /* number of columns in matrix */
) ;

void sopt_minsortx
(
    SOPTINT         *y, /* n-by-1 (output) */
    SOPTFLOAT const *x, /* n-by-1 (input not modified) */
    SOPTINT         *w, /* n-by-1, (input, working array) */
    SOPTINT          n  /* number of elements to sort */
) ;

void sopt_minsorti
(
    SOPTINT        *y, /* n-by-1 (output) */
    SOPTINT  const *x, /* n-by-1 (input not modified) */
    SOPTINT        *w, /* n-by-1, (input, working array) */
    SOPTINT         n  /* number of elements to sort */
) ;

#endif
