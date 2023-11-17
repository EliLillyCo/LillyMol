/* src/Internal/Utilities_Nondistributable/GFP_Tools/dq1da.f -- translated by f2c (version 20200916).
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

static integer c__4 = 4;
static integer c__1 = 1;
static doublereal c_b19 = 1.5;

/* Subroutine */ int dq1da_(D_fp f, doublereal *a, doublereal *b, doublereal *
	eps, doublereal *r__, doublereal *e, integer *kf, integer *iflag)
{
    static doublereal w[300]	/* was [50][6] */;
    static logical rst;
    static doublereal fmin, fmax;
    static integer nmax, nint;
    extern /* Subroutine */ int q1dax_(D_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, logical *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);

/* ***BEGIN PROLOGUE  Q1DA */
/* ***DATE WRITTEN   821018   (YYMMDD) */
/* ***REVISION DATE  870525   (YYMMDD) */
/* ***CATEGORY NO.  H2A1A1 */
/* ***KEYWORDS   ADAPTIVE  QUADRATURE, AUTOMATIC  QUADRATURE */
/* ***AUTHOR  KAHANER, DAVID K., SCIENTIFIC COMPUTING DIVISION, NBS. */
/* ***PURPOSE  Approximates one dimensional integrals of user defined */
/*            functions, easy to use. */

/* ***DESCRIPTION */
/*       Q1DA IS A SUBROUTINE FOR THE AUTOMATIC EVALUATION */
/*            OF THE DEFINITE INTEGRAL OF A USER DEFINED FUNCTION */
/*            OF ONE VARIABLE. */

/*       From the book "Numerical Methods and Software" */
/*            by  D. Kahaner, C. Moler, S. Nash */
/*               Prentice Hall 1988 */

/*       A R G U M E N T S   I N   T H E   C A L L   S E Q U E N C E */

/*       A */
/*       B     (INPUT) THE ENDPOINTS OF THE INTEGRATION INTERVAL */
/*       EPS   (INPUT) THE ACCURACY TO WHICH YOU WANT THE INTEGRAL */
/*                COMPUTED.  IF YOU WANT 2 DIGITS OF ACCURACY SET */
/*                EPS=.01, FOR 3 DIGITS SET EPS=.001, ETC. */
/*                EPS MUST BE POSITIVE. */
/*       R     (OUTPUT) Q1DA'S BEST ESTIMATE OF YOUR INTEGRAL */
/*       E     (OUTPUT) AN ESTIMATE OF ABS(INTEGRAL-R) */
/*       KF    (OUTPUT) THE COST OF THE INTEGRATION, MEASURED IN */
/*                   NUMBER OF EVALUATIONS OF YOUR INTEGRAND. */
/*                   KF WILL ALWAYS BE AT LEAST 30. */
/*       IFLAG (OUTPUT) TERMINATION FLAG...POSSIBLE VALUES ARE */
/*               0   NORMAL COMPLETION, E SATISFIES */
/*                        E<EPS  AND  E<EPS*ABS(R) */
/*               1   NORMAL COMPLETION, E SATISFIES */
/*                        E<EPS, BUT E>EPS*ABS(R) */
/*               2   NORMAL COMPLETION, E SATISFIES */
/*                        E<EPS*ABS(R), BUT E>EPS */
/*               3   NORMAL COMPLETION BUT EPS WAS TOO SMALL TO */
/*                     SATISFY ABSOLUTE OR RELATIVE ERROR REQUEST. */

/*               4   ABORTED CALCULATION BECAUSE OF SERIOUS ROUNDING */
/*                     ERROR.  PROBABLY E AND R ARE CONSISTENT. */
/*               5   ABORTED CALCULATION BECAUSE OF INSUFFICIENT STORAGE. */
/*                     R AND E ARE CONSISTENT. */
/*               6   ABORTED CALCULATION BECAUSE OF SERIOUS DIFFICULTIES */
/*                     MEETING YOUR ERROR REQUEST. */
/*               7   ABORTED CALCULATION BECAUSE EPS WAS SET <= 0.0 */

/*            NOTE...IF IFLAG=3, 4, 5 OR 6 CONSIDER USING Q1DAX INSTEAD. */

/*    W H E R E   I S   Y O U R   I N T E G R A N D ? */

/*        YOU MUST WRITE A FORTRAN FUNCTION, CALLED F, TO EVALUATE */
/*        THE INTEGRAND.  USUALLY THIS LOOKS LIKE... */
/*                 FUNCTION F(X) */
/*                    F=(EVALUATE THE INTEGRAND AT THE POINT X) */
/*                    RETURN */
/*                 END */


/*    T Y P I C A L   P R O B L E M   S E T U P */

/*          A=0.0 */
/*          B=1.0          (SET INTERVAL ENDPOINTS TO [0,1]) */
/*          EPS=0.001       (SET ACCURACY REQUEST FOR 3 DIGITS) */
/*          CALL Q1DA(A,B,EPS,R,E,KF,IFLAG) */
/*          END */
/*          FUNCTION F(X) */
/*            F=SIN(2.*X)-SQRT(X)     (FOR EXAMPLE) */
/*            RETURN */
/*          END */
/*      FOR THIS SAMPLE PROBLEM, THE OUTPUT IS */
/*  0.0    1.0     .001    .041406750    .69077E-07    30    0 */

/*    R E M A R K   I. */

/*           A SMALL AMOUT OF RANDOMIZATION IS BUILT INTO THIS PROGRAM. */
/*           CALLING Q1DA A FEW TIMES IN SUCCESSION WILL GIVE DIFFERENT */
/*           BUT HOPEFULLY CONSISTENT RESULTS. */

/*   R E M A R K   II. */

/*           THIS ROUTINE IS DESIGNED FOR INTEGRATION OVER A FINITE */
/*           INTERVAL.  THUS THE INPUT ARGUMENTS A AND B MUST BE */
/*           VALID REAL NUMBERS ON YOUR COMPUTER.  IF YOU WANT TO DO */
/*           AN INTEGRAL OVER AN INFINITE INTERVAL SET A OR B OR BOTH */
/*           LARGE ENOUGH SO THAT THE INTERVAL [A,B] CONTAINS MOST OF */
/*           THE INTEGRAND.  CARE IS NECESSARY, HOWEVER.  FOR EXAMPLE, */
/*           TO INTEGRATE EXP(-X*X) ON THE ENTIRE REAL LINE ONE COULD */
/*           TAKE A=-20., B=20. OR SIMILAR VALUES TO GET GOOD RESULTS. */
/*           IF YOU TOOK A=-1.E10 AND B=+1.E10 TWO BAD THINGS WOULD */
/*           OCCUR. FIRST, YOU WILL CERTAINLY GET AN ERROR MESSAGE FROM */
/*           THE EXP ROUTINE, AS ITS ARGUMENT IS TOO SMALL.  OTHER */
/*           THINGS COULD HAPPEN TOO, FOR EXAMPLE AN UNDERFLOW. */
/*           SECOND, EVEN IF THE ARITHMETIC WORKED PROPERLY Q1DA WILL */
/*           SURELY GIVE AN INCORRECT ANSWER, BECAUSE ITS FIRST TRY */
/*           AT SAMPLING THE INTEGRAND IS BASED ON YOUR SCALING AND */
/*           IT IS VERY UNLIKELY TO SELECT EVALUATION POINTS IN THE */
/*           INFINITESMALLY SMALL INTERVAL [-20,20] WHERE ALL THE */
/*           INTEGRAND IS CONCENTRATED, WHEN A, B ARE SO LARGE. */

/*    M O R E   F L E X I B I L I T Y */

/*           Q1DA IS AN EASY TO USE DRIVER FOR ANOTHER PROGRAM, Q1DAX. */
/*           Q1DAX PROVIDES SEVERAL OPTIONS WHICH ARE NOT AVAILABLE */
/*                WITH Q1DA. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  Q1DAX */
/* ***END PROLOGUE  Q1DA */


/* ***FIRST EXECUTUTABLE STATEMENT  Q1DA */
    nint = 1;
    rst = FALSE_;
    nmax = 50;
    q1dax_((D_fp)f, a, b, eps, r__, e, &nint, &rst, w, &nmax, &fmin, &fmax, 
	    kf, iflag);
    return 0;
} /* dq1da_ */

/* Subroutine */ int q1dax_(D_fp f, doublereal *a, doublereal *b, doublereal *
	eps, doublereal *r__, doublereal *e, integer *nint, logical *rst, 
	doublereal *w, integer *nmax, doublereal *fmin, doublereal *fmax, 
	integer *kf, integer *iflag)
{
    /* System generated locals */
    integer w_dim1, w_offset, i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer c__, i__;
    static doublereal t, eb, te, xm, tr, te1, te2, tr1, tr2, rab;
    static integer loc;
    static doublereal fmn, rav, fmx;
    extern doublereal uni_(void);
    static doublereal rabs;
    extern /* Subroutine */ int gl15t_(D_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal fminl;
    static integer iroff;
    static doublereal fmaxl, fminr, fmaxr, uflow;
    extern doublereal d1mach_(integer *);
    static integer mxtry;
    static doublereal epmach;
    extern integer idamax_(integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  Q1DAX */
/*     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE */
/*     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS */
/* ***END PROLOGUE  Q1DAX */


/* ***FIRST EXECUTABLE STATEMENT  Q1DAX */
    /* Parameter adjustments */
    w_dim1 = *nmax;
    w_offset = 1 + w_dim1;
    w -= w_offset;

    /* Function Body */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);
    mxtry = *nmax / 2;
/*          In case there is no more room, we can toss out easy intervals, */
/*             at most MXTRY times. */
    if (*a == *b) {
	*r__ = 0.;
	*e = 0.;
	*nint = 0;
	*iflag = 0;
	*kf = 1;
	*fmin = (*f)(a);
	*fmax = *fmin;
	goto L20;
    }
    if (*rst) {
	if (*iflag < 3) {
/* Computing MAX */
/* Computing MAX */
	    d__3 = *eps, d__4 = epmach * 50.;
	    d__1 = uflow * 100., d__2 = max(d__3,d__4) * abs(*r__);
	    eb = max(d__1,d__2);
	    i__1 = *nint;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if ((d__1 = w[i__ + w_dim1 * 3], abs(d__1)) > eb * (w[i__ + (
			w_dim1 << 1)] - w[i__ + w_dim1]) / (*b - *a)) {
		    w[i__ + w_dim1 * 3] = (d__1 = w[i__ + w_dim1 * 3], abs(
			    d__1));
		} else {
		    w[i__ + w_dim1 * 3] = -(d__1 = w[i__ + w_dim1 * 3], abs(
			    d__1));
		}
/* L19: */
	    }
	    goto L15;
	} else {
	    goto L20;
	}
    }
    *kf = 0;
    if (*eps <= 0.f || *nint <= 0 || *nint >= *nmax) {
	*iflag = 7;
	goto L20;
    }
    if (*nint == 1) {
	w[w_dim1 + 1] = *a;
	w[(w_dim1 << 1) + 2] = *b;
	w[w_dim1 * 5 + 1] = *a;
	w[w_dim1 * 6 + 1] = *b;
	w[w_dim1 * 5 + 2] = *a;
	w[w_dim1 * 6 + 2] = *b;
/*          SELECT FIRST SUBDIVISION RANDOMLY */
	w[(w_dim1 << 1) + 1] = *a + (*b - *a) / 2. * (uni_() * 2 + 7.) / 8.;
	w[w_dim1 + 2] = w[(w_dim1 << 1) + 1];
	*nint = 2;
    } else {
	if (w[w_dim1 + 1] != *a || w[*nint + 1 + w_dim1] != *b) {
	    *iflag = 8;
	    goto L20;
	}
	w[w_dim1 * 5 + 1] = *a;
	i__1 = *nint;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    w[i__ + (w_dim1 << 1)] = w[i__ + 1 + w_dim1];
	    w[i__ + w_dim1 * 5] = w[i__ + w_dim1];
	    w[i__ + w_dim1 * 6] = w[i__ + (w_dim1 << 1)];
/* L89: */
	}
    }

    *iflag = 0;
    iroff = 0;
    rabs = 0.f;
    i__1 = *nint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = w[i__ + w_dim1 * 5];
	d__2 = w[i__ + w_dim1 * 6];
	gl15t_((D_fp)f, &w[i__ + w_dim1], &w[i__ + (w_dim1 << 1)], &d__1, &
		d__2, &w[i__ + (w_dim1 << 2)], &w[i__ + w_dim1 * 3], &rab, &
		rav, &fmn, &fmx);
	*kf += 15;
	if (i__ == 1) {
	    *r__ = w[i__ + (w_dim1 << 2)];
	    *e = w[i__ + w_dim1 * 3];
	    rabs += rab;
	    *fmin = fmn;
	    *fmax = fmx;
	} else {
	    *r__ += w[i__ + (w_dim1 << 2)];
	    *e += w[i__ + w_dim1 * 3];
	    rabs += rab;
	    *fmax = max(*fmax,fmx);
	    *fmin = min(*fmin,fmn);
	}
/* L3: */
    }
    i__1 = *nmax;
    for (i__ = *nint + 1; i__ <= i__1; ++i__) {
	w[i__ + w_dim1 * 3] = 0.;
/* L10: */
    }
L15:

/*   MAIN SUBPROGRAM LOOP */

    if (epmach * 100. * rabs >= abs(*r__) && *e < *eps) {
	goto L20;
    }
/* Computing MAX */
/* Computing MAX */
    d__3 = *eps, d__4 = epmach * 50.;
    d__1 = uflow * 100., d__2 = max(d__3,d__4) * abs(*r__);
    eb = max(d__1,d__2);
    if (*e <= eb) {
	goto L20;
    }
    if (*nint < *nmax) {
	++(*nint);
	c__ = *nint;
    } else {
	c__ = 0;
L16:
	if (c__ == *nmax || mxtry <= 0) {
	    *iflag = 5;
	    goto L20;
	}
	++c__;
	if (w[c__ + w_dim1 * 3] > 0.) {
	    goto L16;
	}
/*            Found an interval to throw out */
	--mxtry;
    }
    loc = idamax_(nint, &w[w_dim1 * 3 + 1], &c__1);
    xm = w[loc + w_dim1] + (w[loc + (w_dim1 << 1)] - w[loc + w_dim1]) / 2.;
/* Computing MAX */
    d__3 = (d__1 = w[loc + w_dim1], abs(d__1)), d__4 = (d__2 = w[loc + (
	    w_dim1 << 1)], abs(d__2));
    if (max(d__3,d__4) > (epmach * 100. + 1.) * (abs(xm) + uflow * 1e3)) {
	d__1 = w[loc + w_dim1 * 5];
	d__2 = w[loc + w_dim1 * 6];
	gl15t_((D_fp)f, &w[loc + w_dim1], &xm, &d__1, &d__2, &tr1, &te1, &rab,
		 &rav, &fminl, &fmaxl);
	*kf += 15;
	if (te1 < eb * (xm - w[loc + w_dim1]) / (*b - *a)) {
	    te1 = -te1;
	}
	d__1 = w[loc + w_dim1 * 5];
	d__2 = w[loc + w_dim1 * 6];
	gl15t_((D_fp)f, &xm, &w[loc + (w_dim1 << 1)], &d__1, &d__2, &tr2, &
		te2, &rab, &rav, &fminr, &fmaxr);
	*kf += 15;
/* Computing MIN */
	d__1 = min(*fmin,fminl);
	*fmin = min(d__1,fminr);
/* Computing MAX */
	d__1 = max(*fmax,fmaxl);
	*fmax = max(d__1,fmaxr);
	if (te2 < eb * (w[loc + (w_dim1 << 1)] - xm) / (*b - *a)) {
	    te2 = -te2;
	}
	te = (d__1 = w[loc + w_dim1 * 3], abs(d__1));
	tr = w[loc + (w_dim1 << 2)];
	w[c__ + w_dim1 * 3] = te2;
	w[c__ + (w_dim1 << 2)] = tr2;
	w[c__ + w_dim1] = xm;
	w[c__ + (w_dim1 << 1)] = w[loc + (w_dim1 << 1)];
	w[c__ + w_dim1 * 5] = w[loc + w_dim1 * 5];
	w[c__ + w_dim1 * 6] = w[loc + w_dim1 * 6];
	w[loc + w_dim1 * 3] = te1;
	w[loc + (w_dim1 << 2)] = tr1;
	w[loc + (w_dim1 << 1)] = xm;
	*e = *e - te + (abs(te1) + abs(te2));
	*r__ = *r__ - tr + (tr1 + tr2);
	if ((d__1 = abs(te1) + abs(te2) - te, abs(d__1)) < te * .001f) {
	    ++iroff;
	    if (iroff >= 10) {
		*iflag = 4;
		goto L20;
	    }
	}
    } else {
	if (eb > w[loc + w_dim1 * 3]) {
	    w[loc + w_dim1 * 3] = 0.;
	} else {
	    *iflag = 6;
	    goto L20;
	}
    }
    goto L15;

/*        ALL EXITS FROM HERE */

L20:
    if (*iflag >= 4) {
	return 0;
    }
    *iflag = 3;
    t = *eps * abs(*r__);
    if (*e > *eps && *e > t) {
	return 0;
    }
    *iflag = 2;
    if (*e > *eps && *e < t) {
	return 0;
    }
    *iflag = 1;
    if (*e < *eps && *e > t) {
	return 0;
    }
    *iflag = 0;
    return 0;
} /* q1dax_ */

/* Subroutine */ int gl15t_(D_fp f, doublereal *a, doublereal *b, doublereal *
	xl, doublereal *xr, doublereal *r__, doublereal *ae, doublereal *ra, 
	doublereal *rasc, doublereal *fmin, doublereal *fmax)
{
    /* Initialized data */

    static doublereal wg[4] = { .129484966168869693270611432679082,
	    .27970539148927666790146777142378,
	    .381830050505118944950369775488975,
	    .417959183673469387755102040816327 };
    static doublereal xgk[8] = { .991455371120812639206854697526329,
	    .949107912342758524526189684047851,
	    .864864423359769072789712788640926,
	    .741531185599394439863864773280788,
	    .58608723546769113029414483825873,
	    .405845151377397166906606412076961,
	    .207784955007898467600689403773245,0. };
    static doublereal wgk[8] = { .02293532201052922496373200805897,
	    .063092092629978553290700663189204,
	    .104790010322250183839876322541518,
	    .140653259715525918745189590510238,
	    .16900472663926790282658342659855,
	    .190350578064785409913256402421014,
	    .204432940075298892414161999234649,
	    .209482141084727828012999174891714 };
    static doublereal epmach = 0.;
    static doublereal uflow = 0.;

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer j;
    static doublereal u, fc, sl, sr, fv1[7], fv2[7];
    static integer jtw;
    static doublereal absc, resg, resk, phiu, fsum, fval1, fval2;
    static integer jtwm1;
    static doublereal hlgth, centr, reskh;
    extern doublereal d1mach_(integer *);
    static doublereal dhlgth;

/* ***BEGIN PROLOGUE  GL15T */
/* ***REFER TO  Q1DA,Q1DAX */
/* ***END PROLOGUE  GL15T */


/*           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1) */
/*           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR */
/*           CORRESPONDING WEIGHTS ARE GIVEN. */

/*           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE */
/*                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT */
/*                    GAUSS RULE */
/*                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY */
/*                    ADDED TO THE 7-POINT GAUSS RULE */

/*           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE */

/*           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE */





/*           LIST OF MAJOR VARIABLES */
/*           ----------------------- */

/*           CENTR  - MID POINT OF THE INTERVAL */
/*           HLGTH  - HALF-LENGTH OF THE INTERVAL */
/*           ABSC   - ABSCISSA */
/*           FVAL*  - FUNCTION VALUE */
/*           RESG   - R OF THE 7-POINT GAUSS FORMULA */
/*           RESK   - R OF THE 15-POINT KRONROD FORMULA */
/*           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B), */
/*                    I.E. TO I/(B-A) */

/*           MACHINE DEPENDENT CONSTANTS */
/*           --------------------------- */

/*           EPMACH IS THE LARGEST RELATIVE SPACING. */
/*           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE. */

/* ***FIRST EXECUTABLE STATEMENT  GL15T */
    if (epmach == 0.) {
	epmach = d1mach_(&c__4);
	uflow = d1mach_(&c__1);
    }

    if (*xl < *xr) {
	sl = (real) (*xl);
	sr = (real) (*xr);
    } else {
	sl = (real) (*xr);
	sr = (real) (*xl);
    }
    hlgth = (*b - *a) * .5;
    centr = *a + hlgth;
    dhlgth = abs(hlgth);

/*           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO */
/*           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR. */

    u = (centr - *xr) / (*xr - *xl);
    phiu = *xr - (*xr - *xl) * u * u * (2. * u + 3.);
    if (phiu <= sl || phiu >= sr) {
	phiu = centr;
    }
    *fmin = (*f)(&phiu);
    *fmax = *fmin;
    fc = *fmin * (-6. * u * (u + 1.));
    resg = fc * wg[3];
    resk = fc * wgk[7];
    *ra = abs(resk);
    for (j = 1; j <= 3; ++j) {
	jtw = j << 1;
	absc = hlgth * xgk[jtw - 1];
	u = (centr - absc - *xr) / (*xr - *xl);
	phiu = *xr - (*xr - *xl) * u * u * (2. * u + 3.);
	if (phiu <= sl || phiu >= sr) {
	    phiu = centr;
	}
	fval1 = (*f)(&phiu);
	*fmax = max(*fmax,fval1);
	*fmin = min(*fmin,fval1);
	fval1 *= -6. * u * (u + 1.);
	u = (centr + absc - *xr) / (*xr - *xl);
	phiu = *xr - (*xr - *xl) * u * u * (2. * u + 3.);
	if (phiu <= sl || phiu >= sr) {
	    phiu = centr;
	}
	fval2 = (*f)(&phiu);
	*fmax = max(*fmax,fval2);
	*fmin = min(*fmin,fval2);
	fval2 *= -6. * u * (u + 1.);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*ra += wgk[jtw - 1] * (abs(fval1) + abs(fval2));
/* L10: */
    }
    for (j = 1; j <= 4; ++j) {
	jtwm1 = (j << 1) - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	u = (centr - absc - *xr) / (*xr - *xl);
	phiu = *xr - (*xr - *xl) * u * u * (2. * u + 3.);
	if (phiu <= sl || phiu >= sr) {
	    phiu = centr;
	}
	fval1 = (*f)(&phiu);
	*fmax = max(*fmax,fval1);
	*fmin = min(*fmin,fval1);
	fval1 *= -6. * u * (u + 1.);
	u = (centr + absc - *xr) / (*xr - *xl);
	phiu = *xr - (*xr - *xl) * u * u * (2. * u + 3.);
	if (phiu <= sl || phiu >= sr) {
	    phiu = centr;
	}
	fval2 = (*f)(&phiu);
	*fmax = max(*fmax,fval2);
	*fmin = min(*fmin,fval2);
	fval2 *= -6. * u * (u + 1.);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*ra += wgk[jtwm1 - 1] * (abs(fval1) + abs(fval2));
/* L15: */
    }
    reskh = resk * .5;
    *rasc = wgk[7] * (d__1 = fc - reskh, abs(d__1));
    for (j = 1; j <= 7; ++j) {
	*rasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, abs(d__1)) + (d__2 
		= fv2[j - 1] - reskh, abs(d__2)));
/* L20: */
    }
    *r__ = resk * hlgth;
    *ra *= dhlgth;
    *rasc *= dhlgth;
    *ae = (d__1 = (resk - resg) * hlgth, abs(d__1));
    if (*rasc != 0. && *ae != 0.) {
/* Computing MIN */
	d__3 = *ae * 200. / *rasc;
	d__1 = 1., d__2 = pow_dd(&d__3, &c_b19);
	*ae = *rasc * min(d__1,d__2);
    }
    if (*ra > uflow / (epmach * 50.)) {
/* Computing MAX */
	d__1 = epmach * 50. * *ra;
	*ae = max(d__1,*ae);
    }
    return 0;
} /* gl15t_ */

