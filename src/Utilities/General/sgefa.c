/* src/Internal/Utilities_Nondistributable/General/sgefa.f -- translated by f2c (version 20200916).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sgefa_(real *a, integer *lda, integer *n, integer *ipvt, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, l;
    static real t;
    static integer kp1, nm1;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    saxpy_(integer *, real *, real *, integer *, real *, integer *);
    extern integer isamax_(integer *, real *, integer *);


/*     sgefa factors a real matrix by gaussian elimination. */

/*     sgefa is usually called by sgeco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for sgeco) = (1 + 9/n)*(time for sgefa) . */

/*     on entry */

/*        a       real(lda, n) */
/*                the matrix to be factored. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if  u(k,k) .eq. 0.0 .  this is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that sgesl or sgedi will divide by zero */
/*                     if called.  use  rcond  in sgeco for a reliable */
/*                     indication of singularity. */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas saxpy,sscal,isamax */

/*     internal variables */



/*     gaussian elimination with partial pivoting */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        find l = pivot index */

	i__2 = *n - k + 1;
	l = isamax_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        zero pivot implies this column already triangularized */

	if (a[l + k * a_dim1] == 0.f) {
	    goto L40;
	}

/*           interchange if necessary */

	if (l == k) {
	    goto L10;
	}
	t = a[l + k * a_dim1];
	a[l + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L10:

/*           compute multipliers */

	t = -1.f / a[k + k * a_dim1];
	i__2 = *n - k;
	sscal_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           row elimination with column indexing */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[l + j * a_dim1];
	    if (l == k) {
		goto L20;
	    }
	    a[l + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
L20:
	    i__3 = *n - k;
	    saxpy_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    if (a[*n + *n * a_dim1] == 0.f) {
	*info = *n;
    }
    return 0;
} /* sgefa_ */

