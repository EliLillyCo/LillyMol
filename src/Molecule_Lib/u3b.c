/* u3b.f -- translated by f2c (version 20160102).
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

/* ccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc */
/*  w    - w(m) is weight for atom pair  c m           (given) */
/*  x    - x(i,m) are coordinates of atom c m in set x       (given) */
/*  y    - y(i,m) are coordinates of atom c m in set y       (given) */
/*  n    - n is number of atom pairs                         (given) */
/*  mode  - 0:calculate rms only                             (given) */
/*          1:calculate rms,u,t                              (takes longer) */
/*  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result) */
/*  u    - u(i,j) is   rotation  matrix for best superposition  (result) */
/*  t    - t(i)   is translation vector for best superposition  (result) */
/*  ier  - 0: a unique optimal superposition has been determined(result) */
/*       -1: superposition is not unique but optimal */
/*       -2: no result obtained because of negative weights w */
/*           or all weights equal to zero. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */ int u3b_(doublereal *w, doublereal *x, doublereal *y, 
	integer *n, integer *mode, doublereal *rms, doublereal *u, doublereal 
	*t, integer *ier)
{
    /* Initialized data */

    static doublereal sqrt3 = 1.73205080756888;
    static doublereal tol = .01;
    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal three = 3.;
    static integer ip[9] = { 1,2,4,2,3,5,4,5,6 };
    static integer ip2312[4] = { 2,3,1,2 };

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    static doublereal equiv_5[6], equiv_11[6], equiv_14[3];

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), cos(doublereal), 
	    sin(doublereal);

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */, b[9]	/* was [3][3] */, d__;
#define e (equiv_14)
    static doublereal g, h__;
    static integer i__, j, k, l, m;
    static doublereal p, r__[9]	/* was [3][3] */, e0;
#define e1 (equiv_14)
#define e2 (equiv_14 + 1)
#define e3 (equiv_14 + 2)
    static integer m1;
    static doublereal wc, xc[3], yc[3];
#define rr (equiv_5)
#define ss (equiv_11)
#define rr1 (equiv_5)
#define rr2 (equiv_5 + 1)
#define rr3 (equiv_5 + 2)
#define rr4 (equiv_5 + 3)
#define rr5 (equiv_5 + 4)
#define rr6 (equiv_5 + 5)
#define ss1 (equiv_11)
#define ss2 (equiv_11 + 1)
#define ss3 (equiv_11 + 2)
#define ss4 (equiv_11 + 3)
#define ss5 (equiv_11 + 4)
#define ss6 (equiv_11 + 5)
    static doublereal cof, det, cth, sth, spur, sigma, sqrth;

    /* Parameter adjustments */
    --t;
    u -= 4;
    y -= 4;
    x -= 4;
    --w;

    /* Function Body */
/* 156 "rms.for" */
    wc = zero;
    *rms = 0.f;
    e0 = zero;
    for (i__ = 1; i__ <= 3; ++i__) {
	xc[i__ - 1] = zero;
	yc[i__ - 1] = zero;
	t[i__] = 0.f;
	for (j = 1; j <= 3; ++j) {
	    d__ = zero;
	    if (i__ == j) {
		d__ = one;
	    }
	    u[i__ + j * 3] = d__;
	    a[i__ + j * 3 - 4] = d__;
/* L1: */
	    r__[i__ + j * 3 - 4] = zero;
	}
    }
    *ier = -1;
/* **** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y */
/* 170 "rms.for" */
    if (*n < 1) {
	return 0;
    }
/* 172 "rms.for" */
    *ier = -2;
    i__1 = *n;
    for (m = 1; m <= i__1; ++m) {
	if (w[m] < 0.f) {
	    return 0;
	}
	wc += w[m];
	for (i__ = 1; i__ <= 3; ++i__) {
	    xc[i__ - 1] += w[m] * x[i__ + m * 3];
/* L2: */
	    yc[i__ - 1] += w[m] * y[i__ + m * 3];
	}
    }
    if (wc <= zero) {
	return 0;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	xc[i__ - 1] /= wc;
/* **** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X */
/* 182 "rms.for" */
/* L3: */
	yc[i__ - 1] /= wc;
    }
/* 184 "rms.for" */
    i__1 = *n;
    for (m = 1; m <= i__1; ++m) {
	for (i__ = 1; i__ <= 3; ++i__) {
/* Computing 2nd power */
	    d__1 = x[i__ + m * 3] - xc[i__ - 1];
/* Computing 2nd power */
	    d__2 = y[i__ + m * 3] - yc[i__ - 1];
	    e0 += w[m] * (d__1 * d__1 + d__2 * d__2);
/* 187 "rms.for" */
	    d__ = w[m] * (y[i__ + m * 3] - yc[i__ - 1]);
	    for (j = 1; j <= 3; ++j) {
/* **** CALCULATE DETERMINANT OF R(I,J) */
/* 189 "rms.for" */
/* L4: */
		r__[i__ + j * 3 - 4] += d__ * (x[j + m * 3] - xc[j - 1]);
	    }
	}
    }
/* 191 "rms.for" */
    det = r__[0] * (r__[4] * r__[8] - r__[7] * r__[5]) - r__[3] * (r__[1] * 
	    r__[8] - r__[7] * r__[2]) + r__[6] * (r__[1] * r__[5] - r__[4] * 
	    r__[2]);
/* **** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R */
/* 194 "rms.for" */
    sigma = det;
/* 196 "rms.for" */
    m = 0;
    for (j = 1; j <= 3; ++j) {
	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++m;
/* ***************** EIGENVALUES ***************************************** */
/* **** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0 */
/* 200 "rms.for" */
/* L5: */
	    rr[m - 1] = r__[i__ * 3 - 3] * r__[j * 3 - 3] + r__[i__ * 3 - 2] *
		     r__[j * 3 - 2] + r__[i__ * 3 - 1] * r__[j * 3 - 1];
	}
    }
/* 203 "rms.for" */
    spur = (*rr1 + *rr3 + *rr6) / three;
    cof = (*rr3 * *rr6 - *rr5 * *rr5 + *rr1 * *rr6 - *rr4 * *rr4 + *rr1 * *
	    rr3 - *rr2 * *rr2) / three;
/* 205 "rms.for" */
    det *= det;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L6: */
	e[i__ - 1] = spur;
    }
/* **** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR */
/* 208 "rms.for" */
    if (spur <= zero) {
	goto L40;
    }
/* 210 "rms.for" */
    d__ = spur * spur;
    h__ = d__ - cof;
/* **** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER */
/* 212 "rms.for" */
    g = (spur * cof - det) / two - spur * h__;
/* 214 "rms.for" */
    if (h__ <= zero) {
	goto L8;
    }
    sqrth = sqrt(h__);
    d__ = h__ * h__ * h__ - g * g;
    if (d__ < zero) {
	d__ = zero;
    }
    d__ = atan2(sqrt(d__), -g) / three;
    cth = sqrth * cos(d__);
    sth = sqrth * sqrt3 * sin(d__);
    *e1 = spur + cth + cth;
    *e2 = spur - cth + sth;
    *e3 = spur - cth - sth;
/* .....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS */
/* 224 "rms.for" */
    if (*mode != 0) {
	goto L10;
    } else {
	goto L50;
    }
/* **************** EIGENVECTORS ***************************************** */
/* 226 "rms.for" */
L8:
    if (*mode != 0) {
	goto L30;
    } else {
	goto L50;
    }
/* 228 "rms.for" */
L10:
    for (l = 1; l <= 3; l += 2) {
	d__ = e[l - 1];
	*ss1 = (d__ - *rr3) * (d__ - *rr6) - *rr5 * *rr5;
	*ss2 = (d__ - *rr6) * *rr2 + *rr4 * *rr5;
	*ss3 = (d__ - *rr1) * (d__ - *rr6) - *rr4 * *rr4;
	*ss4 = (d__ - *rr3) * *rr4 + *rr2 * *rr5;
	*ss5 = (d__ - *rr1) * *rr5 + *rr2 * *rr4;
	*ss6 = (d__ - *rr1) * (d__ - *rr3) - *rr2 * *rr2;
	j = 1;
	if (abs(*ss1) >= abs(*ss3)) {
	    goto L12;
	}
	j = 2;
	if (abs(*ss3) >= abs(*ss6)) {
	    goto L13;
	}
L11:
	j = 3;
	goto L13;
L12:
	if (abs(*ss1) < abs(*ss6)) {
	    goto L11;
	}
L13:
	d__ = zero;
	j = (j - 1) * 3;
	for (i__ = 1; i__ <= 3; ++i__) {
	    k = ip[i__ + j - 1];
	    a[i__ + l * 3 - 4] = ss[k - 1];
/* L14: */
	    d__ += ss[k - 1] * ss[k - 1];
	}
	if (d__ > zero) {
	    d__ = one / sqrt(d__);
	}
	for (i__ = 1; i__ <= 3; ++i__) {
/* L15: */
	    a[i__ + l * 3 - 4] *= d__;
	}
    }
    d__ = a[0] * a[6] + a[1] * a[7] + a[2] * a[8];
    m1 = 3;
    m = 1;
    if (*e1 - *e2 > *e2 - *e3) {
	goto L16;
    }
    m1 = 1;
    m = 3;
L16:
    p = zero;
    for (i__ = 1; i__ <= 3; ++i__) {
	a[i__ + m1 * 3 - 4] -= d__ * a[i__ + m * 3 - 4];
/* L17: */
/* Computing 2nd power */
	d__1 = a[i__ + m1 * 3 - 4];
	p += d__1 * d__1;
    }
    if (p <= tol) {
	goto L19;
    }
    p = one / sqrt(p);
    for (i__ = 1; i__ <= 3; ++i__) {
/* L18: */
	a[i__ + m1 * 3 - 4] *= p;
    }
    goto L21;
L19:
    p = one;
    for (i__ = 1; i__ <= 3; ++i__) {
	if (p < (d__1 = a[i__ + m * 3 - 4], abs(d__1))) {
	    goto L20;
	}
	p = (d__1 = a[i__ + m * 3 - 4], abs(d__1));
	j = i__;
L20:
	;
    }
    k = ip2312[j - 1];
    l = ip2312[j];
/* Computing 2nd power */
    d__1 = a[k + m * 3 - 4];
/* Computing 2nd power */
    d__2 = a[l + m * 3 - 4];
    p = sqrt(d__1 * d__1 + d__2 * d__2);
    if (p <= tol) {
	goto L40;
    }
    a[j + m1 * 3 - 4] = zero;
    a[k + m1 * 3 - 4] = -(a[l + m * 3 - 4] / p);
    a[l + m1 * 3 - 4] = a[k + m * 3 - 4] / p;
L21:
    a[3] = a[7] * a[2] - a[1] * a[8];
    a[4] = a[8] * a[0] - a[2] * a[6];
/* ****************** ROTATION MATRIX ************************************ */
/* 282 "rms.for" */
    a[5] = a[6] * a[1] - a[0] * a[7];
/* 284 "rms.for" */
L30:
    for (l = 1; l <= 2; ++l) {
	d__ = zero;
	for (i__ = 1; i__ <= 3; ++i__) {
	    b[i__ + l * 3 - 4] = r__[i__ - 1] * a[l * 3 - 3] + r__[i__ + 2] * 
		    a[l * 3 - 2] + r__[i__ + 5] * a[l * 3 - 1];
/* 288 "rms.for" */
/* L31: */
/* Computing 2nd power */
	    d__1 = b[i__ + l * 3 - 4];
	    d__ += d__1 * d__1;
	}
	if (d__ > zero) {
	    d__ = one / sqrt(d__);
	}
	for (i__ = 1; i__ <= 3; ++i__) {
/* L32: */
	    b[i__ + l * 3 - 4] *= d__;
	}
    }
    d__ = b[0] * b[3] + b[1] * b[4] + b[2] * b[5];
    p = zero;
    for (i__ = 1; i__ <= 3; ++i__) {
	b[i__ + 2] -= d__ * b[i__ - 1];
/* L33: */
/* Computing 2nd power */
	d__1 = b[i__ + 2];
	p += d__1 * d__1;
    }
    if (p <= tol) {
	goto L35;
    }
    p = one / sqrt(p);
    for (i__ = 1; i__ <= 3; ++i__) {
/* L34: */
	b[i__ + 2] *= p;
    }
    goto L37;
L35:
    p = one;
    for (i__ = 1; i__ <= 3; ++i__) {
	if (p < (d__1 = b[i__ - 1], abs(d__1))) {
	    goto L36;
	}
	p = (d__1 = b[i__ - 1], abs(d__1));
	j = i__;
L36:
	;
    }
    k = ip2312[j - 1];
    l = ip2312[j];
/* Computing 2nd power */
    d__1 = b[k - 1];
/* Computing 2nd power */
    d__2 = b[l - 1];
    p = sqrt(d__1 * d__1 + d__2 * d__2);
    if (p <= tol) {
	goto L40;
    }
    b[j + 2] = zero;
    b[k + 2] = -(b[l - 1] / p);
    b[l + 2] = b[k - 1] / p;
L37:
    b[6] = b[1] * b[5] - b[4] * b[2];
    b[7] = b[2] * b[3] - b[5] * b[0];
    b[8] = b[0] * b[4] - b[3] * b[1];
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* ****************** TRANSLATION VECTOR ********************************* */
/* 320 "rms.for" */
/* L39: */
	    u[i__ + j * 3] = b[i__ - 1] * a[j - 1] + b[i__ + 2] * a[j + 2] + 
		    b[i__ + 5] * a[j + 5];
	}
    }
L40:
    for (i__ = 1; i__ <= 3; ++i__) {
/* ********************** RMS ERROR ************************************** */
/* 323 "rms.for" */
/* L41: */
	t[i__] = yc[i__ - 1] - u[i__ + 3] * xc[0] - u[i__ + 6] * xc[1] - u[
		i__ + 9] * xc[2];
    }
L50:
    for (i__ = 1; i__ <= 3; ++i__) {
	if (e[i__ - 1] < zero) {
	    e[i__ - 1] = zero;
	}
/* L51: */
	e[i__ - 1] = sqrt(e[i__ - 1]);
    }
    *ier = 0;
    if (*e2 <= *e1 * 1e-5) {
	*ier = -1;
    }
    d__ = *e3;
    if (sigma >= 0.f) {
	goto L52;
    }
    d__ = -d__;
    if (*e2 - *e3 <= *e1 * 1e-5) {
	*ier = -1;
    }
L52:
    d__ = d__ + *e2 + *e1;
    *rms = e0 - d__ - d__;
    if (*rms < 0.f) {
	*rms = 0.f;
    }
    return 0;
} /* u3b_ */

#undef ss6
#undef ss5
#undef ss4
#undef ss3
#undef ss2
#undef ss1
#undef rr6
#undef rr5
#undef rr4
#undef rr3
#undef rr2
#undef rr1
#undef ss
#undef rr
#undef e3
#undef e2
#undef e1
#undef e


