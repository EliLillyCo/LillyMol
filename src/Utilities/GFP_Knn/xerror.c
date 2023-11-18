/* src/Internal/Utilities_Nondistributable/GFP_Tools/xerror.f -- translated by f2c (version 20200916).
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

static integer c__0 = 0;
static real c_b6 = 0.f;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__4 = 4;

/* Subroutine */ int xerror_(char *messg, integer *nmessg, integer *nerr, 
	integer *level, ftnlen messg_len)
{
    extern /* Subroutine */ int xerrwv_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, real *, real *, 
	    ftnlen);

/* ***BEGIN PROLOGUE  XERROR */
/* ***DATE WRITTEN   790801   (YYMMDD) */
/* ***REVISION DATE  870930   (YYMMDD) */
/* ***CATEGORY NO.  R3C */
/* ***KEYWORDS  ERROR,XERROR PACKAGE */
/* ***AUTHOR  JONES, R. E., (SNLA) */
/* ***PURPOSE  Processes an error (diagnostic) message. */
/* ***DESCRIPTION */
/*    From the book "Numerical Methods and Software" */
/*       by  D. Kahaner, C. Moler, S. Nash */
/*           Prentice Hall 1988 */
/*     Abstract */
/*        XERROR processes a diagnostic message. It is a stub routine */
/*        written for the book above. Actually, XERROR is a sophisticated */
/*        error handling package with many options, and is described */
/*        in the reference below. Our version has the same calling sequence */
/*        but only prints an error message and either returns (if the */
/*        input value of ABS(LEVEL) is less than 2) or stops (if the */
/*        input value of ABS(LEVEL) equals 2). */

/*     Description of Parameters */
/*      --Input-- */
/*        MESSG - the Hollerith message to be processed. */
/*        NMESSG- the actual number of characters in MESSG. */
/*                (this is ignored in this stub routine) */
/*        NERR  - the error number associated with this message. */
/*                NERR must not be zero. */
/*                (this is ignored in this stub routine) */
/*        LEVEL - error category. */
/*                =2 means this is an unconditionally fatal error. */
/*                =1 means this is a recoverable error.  (I.e., it is */
/*                   non-fatal if XSETF has been appropriately called.) */
/*                =0 means this is a warning message only. */
/*                =-1 means this is a warning message which is to be */
/*                   printed at most once, regardless of how many */
/*                   times this call is executed. */
/*                 (in this stub routine */
/*                       LEVEL=2 causes a message to be printed and then a */
/*                                         stop. */
/*                       LEVEL<2 causes a message to be printed and then a */
/*                                         return. */

/*     Examples */
/*        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2) */
/*        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.', */
/*                    43,2,1) */
/*        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F */
/*    1ULLY COLLAPSED.',65,3,0) */
/*        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1) */

/* ***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR- */
/*                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES, */
/*                 1982. */
/* ***ROUTINES CALLED  XERRWV */
/* ***END PROLOGUE  XERROR */
/* ***FIRST EXECUTABLE STATEMENT  XERROR */
    xerrwv_(messg, nmessg, nerr, level, &c__0, &c__0, &c__0, &c__0, &c_b6, &
	    c_b6, messg_len);
    return 0;
} /* xerror_ */

/* Subroutine */ int xerrwv_(char *messg, integer *nmessg, integer *nerr, 
	integer *level, integer *ni, integer *i1, integer *i2, integer *nr, 
	real *r1, real *r2, ftnlen messg_len)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 6, 0, 0, 0 };
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };


/* ***BEGIN PROLOGUE  XERRWV */
/* ***DATE WRITTEN   800319   (YYMMDD) */
/* ***REVISION DATE  870930   (YYMMDD) */
/* ***CATEGORY NO.  R3C */
/* ***KEYWORDS  ERROR,XERROR PACKAGE */
/* ***AUTHOR  JONES, R. E., (SNLA) */
/* ***PURPOSE  Processes error message allowing 2 integer and two real */
/*            values to be included in the message. */
/* ***DESCRIPTION */
/*    From the book "Numerical Methods and Software" */
/*       by  D. Kahaner, C. Moler, S. Nash */
/*           Prentice Hall 1988 */
/*     Abstract */
/*        XERRWV prints a diagnostic error message. */
/*        In addition, up to two integer values and two real */
/*        values may be printed along with the message. */
/*        A stub routine for the book above. The actual XERRWV is described */
/*        in the reference below and contains many other options. */

/*     Description of Parameters */
/*      --Input-- */
/*        MESSG - the Hollerith message to be processed. */
/*        NMESSG- the actual number of characters in MESSG. */
/*                (ignored in this stub) */
/*        NERR  - the error number associated with this message. */
/*                NERR must not be zero. */
/*                (ignored in this stub) */
/*        LEVEL - error category. */
/*                =2 means this is an unconditionally fatal error. */
/*                =1 means this is a recoverable error.  (I.e., it is */
/*                   non-fatal if XSETF has been appropriately called.) */
/*                =0 means this is a warning message only. */
/*                =-1 means this is a warning message which is to be */
/*                   printed at most once, regardless of how many */
/*                   times this call is executed. */
/*                  (in this stub LEVEL=2 causes an error message to be */
/*                                          printed followed by a stop, */
/*                                LEVEL<2 causes an error message to be */
/*                                          printed followed by a return.) */
/*        NI    - number of integer values to be printed. (0 to 2) */
/*        I1    - first integer value. */
/*        I2    - second integer value. */
/*        NR    - number of real values to be printed. (0 to 2) */
/*        R1    - first real value. */
/*        R2    - second real value. */

/*     Examples */
/*        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2, */
/*    1   1,NUM,0,0,0.,0.) */
/*        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM ( */
/*    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN) */

/* ***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR- */
/*                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES, */
/*                 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  XERRWV */
/* ***FIRST EXECUTABLE STATEMENT  XERRWV */
    s_wsle(&io___1);
    do_lio(&c__9, &c__1, messg, messg_len);
    e_wsle();
    if (*ni == 2) {
	s_wsle(&io___2);
	do_lio(&c__3, &c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*i2), (ftnlen)sizeof(integer));
	e_wsle();
    } else if (*ni == 1) {
	s_wsle(&io___3);
	do_lio(&c__3, &c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (*nr == 2) {
	s_wsle(&io___4);
	do_lio(&c__4, &c__1, (char *)&(*r1), (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&(*r2), (ftnlen)sizeof(real));
	e_wsle();
    } else if (*nr == 1) {
	s_wsle(&io___5);
	do_lio(&c__4, &c__1, (char *)&(*r1), (ftnlen)sizeof(real));
	e_wsle();
    }
    if (abs(*level) < 2) {
	return 0;
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* xerrwv_ */

