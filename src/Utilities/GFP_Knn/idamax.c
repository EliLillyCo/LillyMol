/* src/Internal/Utilities_Nondistributable/GFP_Tools/idamax.f -- translated by f2c (version 20200916).
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

integer idamax_(integer *n, doublereal *sx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, ii, ns;
    static doublereal xmag, smax;

/* ***BEGIN PROLOGUE  ISAMAX */
/*     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE */
/*     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS */
/*     From the book "Numerical Methods and Software" */
/*          by  D. Kahaner, C. Moler, S. Nash */
/*               Prentice Hall 1988 */
/* ***END PROLOGUE  ISAMAX */

/* ***FIRST EXECUTABLE STATEMENT  ISAMAX */
    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0;
    if (*n <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n <= 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENTS NOT EQUAL TO 1. */

    smax = abs(sx[1]);
    ns = *n * *incx;
    ii = 1;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	xmag = (d__1 = sx[i__], abs(d__1));
	if (xmag <= smax) {
	    goto L5;
	}
	ret_val = ii;
	smax = xmag;
L5:
	++ii;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENTS EQUAL TO 1. */

L20:
    smax = abs(sx[1]);
    i__2 = *n;
    for (i__ = 2; i__ <= i__2; ++i__) {
	xmag = (d__1 = sx[i__], abs(d__1));
	if (xmag <= smax) {
	    goto L30;
	}
	ret_val = i__;
	smax = xmag;
L30:
	;
    }
    return ret_val;
} /* idamax_ */

