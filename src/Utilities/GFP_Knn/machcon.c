/* src/Internal/Utilities_Nondistributable/GFP_Tools/machcon.f -- translated by f2c (version 20200916).
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

static integer c__25 = 25;
static integer c__1 = 1;
static integer c__2 = 2;

doublereal r1mach_(integer *i__)
{
    /* Initialized data */

    static struct {
	real e_1[5];
	integer fill_2[1];
	} equiv_4 = { 1.18e-38f, 3.4e38f, 5.95e-8f, 1.19e-7f, .30102999566f };


    /* System generated locals */
    real ret_val;

    /* Local variables */
#define log10 ((integer *)&equiv_4 + 4)
#define large ((integer *)&equiv_4 + 1)
#define rmach ((real *)&equiv_4)
#define small ((integer *)&equiv_4)
#define diver ((integer *)&equiv_4 + 3)
#define right ((integer *)&equiv_4 + 2)
    extern /* Subroutine */ int xerror_(char *, integer *, integer *, integer 
	    *, ftnlen);

/* ***BEGIN PROLOGUE  R1MACH */
/* ***DATE WRITTEN   790101   (YYMMDD) */
/* ***REVISION DATE  831014   (YYMMDD) */
/* ***CATEGORY NO.  R1 */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  FOX, P. A., (BELL LABS) */
/*           HALL, A. D., (BELL LABS) */
/*           SCHRYER, N. L., (BELL LABS) */
/* ***PURPOSE  Returns single precision machine dependent constants */
/* ***DESCRIPTION */
/*     From the book, "Numerical Methods and Software" by */
/*                D. Kahaner, C. Moler, S. Nash */
/*                Prentice Hall, 1988 */


/*     R1MACH can be used to obtain machine-dependent parameters */
/*     for the local machine environment.  It is a function */
/*     subroutine with one (input) argument, and can be called */
/*     as follows, for example */

/*          A = R1MACH(I) */

/*     where I=1,...,5.  The (output) value of A above is */
/*     determined by the (input) value of I.  The results for */
/*     various values of I are discussed below. */

/*  Single-Precision Machine Constants */
/*  R1MACH(1) = B**(EMIN-1), the smallest positive magnitude. */
/*  R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*  R1MACH(3) = B**(-T), the smallest relative spacing. */
/*  R1MACH(4) = B**(1-T), the largest relative spacing. */
/*  R1MACH(5) = LOG10(B) */
/* ***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR */
/*                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE- */
/*                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978, */
/*                 PP. 177-188. */
/* ***ROUTINES CALLED  XERROR */
/* ***END PROLOGUE  R1MACH */





/*     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5). */

/*      DATA RMACH(1) / O"00014000000000000000" / */
/*      DATA RMACH(2) / O"37767777777777777777" / */
/*      DATA RMACH(3) / O"16404000000000000000" / */
/*      DATA RMACH(4) / O"16414000000000000000" / */
/*      DATA RMACH(5) / O"17164642023241175720" / */

/*     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES */

/*     DATA RMACH(1) / X'9000400000000000' / */
/*     DATA RMACH(2) / X'6FFF7FFFFFFFFFFF' / */
/*     DATA RMACH(3) / X'FFA3400000000000' / */
/*     DATA RMACH(4) / X'FFA4400000000000' / */
/*     DATA RMACH(5) / X'FFD04D104D427DE8' / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES. */

/*     DATA RMACH(1) / 00564000000000000000B / */
/*     DATA RMACH(2) / 37767777777777777776B / */
/*     DATA RMACH(3) / 16414000000000000000B / */
/*     DATA RMACH(4) / 16424000000000000000B / */
/*     DATA RMACH(5) / 17164642023241175720B / */

/*     MACHINE CONSTANTS FOR THE CRAY 1 */

/*     DATA RMACH(1) / 200034000000000000000B / */
/*     DATA RMACH(2) / 577767777777777777776B / */
/*     DATA RMACH(3) / 377224000000000000000B / */
/*     DATA RMACH(4) / 377234000000000000000B / */
/*     DATA RMACH(5) / 377774642023241175720B / */


/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA RMACH(1) / Z00100000 / */
/*     DATA RMACH(2) / Z7FFFFFFF / */
/*     DATA RMACH(3) / Z3B100000 / */
/*     DATA RMACH(4) / Z3C100000 / */
/*     DATA RMACH(5) / Z41134413 / */

/*     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS) */


/*     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR). */

/*     DATA RMACH(1) / "000400000000 / */
/*     DATA RMACH(2) / "377777777777 / */
/*     DATA RMACH(3) / "146400000000 / */
/*     DATA RMACH(4) / "147400000000 / */
/*     DATA RMACH(5) / "177464202324 / */


/*     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP, */
/*       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO */
/*       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.) */

/*     DATA SMALL(1) / X'00800000' / */
/*     DATA LARGE(1) / X'7F7FFFFF' / */
/*     DATA RIGHT(1) / X'33800000' / */
/*     DATA DIVER(1) / X'34000000' / */
/*     DATA LOG10(1) / X'3E9A209B' / */


/*     MACHINE CONSTANTS FOR THE VAX 11/780 */
/*    (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*  *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS*** */

/*     DATA SMALL(1) /       128 / */
/*     DATA LARGE(1) /    -32769 / */
/*     DATA RIGHT(1) /     13440 / */
/*     DATA DIVER(1) /     13568 / */
/*     DATA LOG10(1) / 547045274 / */

/*  ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS*** */
/*     DATA SMALL(1) / Z00000080 / */
/*     DATA LARGE(1) / ZFFFF7FFF / */
/*     DATA RIGHT(1) / Z00003480 / */
/*     DATA DIVER(1) / Z00003500 / */
/*     DATA LOG10(1) / Z209B3F9A / */


/* ***FIRST EXECUTABLE STATEMENT  R1MACH */
    if (*i__ < 1 || *i__ > 5) {
	xerror_("R1MACH -- I OUT OF BOUNDS", &c__25, &c__1, &c__2, (ftnlen)25)
		;
    }

    ret_val = rmach[*i__ - 1];
    return ret_val;

} /* r1mach_ */

#undef right
#undef diver
#undef small
#undef rmach
#undef large
#undef log10


doublereal d1mach_(integer *i__)
{
    /* Initialized data */

    static struct {
	doublereal e_1[5];
	doublereal fill_2[1];
	} equiv_4 = { 2.23e-308, 1.79e308, 1.11e-16, 2.22e-16, 
		.301029995663981195 };


    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
#define log10 ((integer *)&equiv_4 + 8)
#define dmach ((doublereal *)&equiv_4)
#define large ((integer *)&equiv_4 + 2)
#define small ((integer *)&equiv_4)
#define diver ((integer *)&equiv_4 + 6)
#define right ((integer *)&equiv_4 + 4)
    extern /* Subroutine */ int xerror_(char *, integer *, integer *, integer 
	    *, ftnlen);

/* ***BEGIN PROLOGUE  D1MACH */
/* ***DATE WRITTEN   750101   (YYMMDD) */
/* ***REVISION DATE  831014   (YYMMDD) */
/* ***CATEGORY NO.  R1 */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  FOX, P. A., (BELL LABS) */
/*           HALL, A. D., (BELL LABS) */
/*           SCHRYER, N. L., (BELL LABS) */
/* ***PURPOSE  Returns double precision machine dependent constants */
/* ***DESCRIPTION */
/*     From the book, "Numerical Methods and Software" by */
/*                D. Kahaner, C. Moler, S. Nash */
/*                Prentice Hall, 1988 */


/*     D1MACH can be used to obtain machine-dependent parameters */
/*     for the local machine environment.  It is a function */
/*     subprogram with one (input) argument, and can be called */
/*     as follows, for example */

/*          D = D1MACH(I) */

/*     where I=1,...,5.  The (output) value of D above is */
/*     determined by the (input) value of I.  The results for */
/*     various values of I are discussed below. */

/*  Double-precision machine constants */
/*  D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude. */
/*  D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*  D1MACH( 3) = B**(-T), the smallest relative spacing. */
/*  D1MACH( 4) = B**(1-T), the largest relative spacing. */
/*  D1MACH( 5) = LOG10(B) */
/* ***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A */
/*                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL */
/*                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188. */
/* ***ROUTINES CALLED  XERROR */
/* ***END PROLOGUE  D1MACH */





/*     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5). */

/*      DATA SMALL(1) / O"00604000000000000000" / */
/*      DATA SMALL(2) / O"00000000000000000000" / */

/*      DATA LARGE(1) / O"37767777777777777777" / */
/*      DATA LARGE(2) / O"37167777777777777777" / */

/*      DATA RIGHT(1) / O"15604000000000000000" / */
/*      DATA RIGHT(2) / O"15000000000000000000" / */

/*      DATA DIVER(1) / O"15614000000000000000" / */
/*      DATA DIVER(2) / O"15010000000000000000" / */

/*      DATA LOG10(1) / O"17164642023241175717" / */
/*      DATA LOG10(2) / O"16367571421742254654" / */

/*     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES */

/*     DATA SMALL(1) / X'9000400000000000' / */
/*     DATA SMALL(2) / X'8FD1000000000000' / */

/*     DATA LARGE(1) / X'6FFF7FFFFFFFFFFF' / */
/*     DATA LARGE(2) / X'6FD07FFFFFFFFFFF' / */

/*     DATA RIGHT(1) / X'FF74400000000000' / */
/*     DATA RIGHT(2) / X'FF45000000000000' / */

/*     DATA DIVER(1) / X'FF75400000000000' / */
/*     DATA DIVER(2) / X'FF46000000000000' / */

/*     DATA LOG10(1) / X'FFD04D104D427DE7' / */
/*     DATA LOG10(2) / X'FFA17DE623E2566A' / */


/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES. */

/*     DATA SMALL(1) / 00564000000000000000B / */
/*     DATA SMALL(2) / 00000000000000000000B / */

/*     DATA LARGE(1) / 37757777777777777777B / */
/*     DATA LARGE(2) / 37157777777777777777B / */

/*     DATA RIGHT(1) / 15624000000000000000B / */
/*     DATA RIGHT(2) / 00000000000000000000B / */

/*     DATA DIVER(1) / 15634000000000000000B / */
/*     DATA DIVER(2) / 00000000000000000000B / */

/*     DATA LOG10(1) / 17164642023241175717B / */
/*     DATA LOG10(2) / 16367571421742254654B / */

/*     MACHINE CONSTANTS FOR THE CRAY 1 */

/*     DATA SMALL(1) / 201354000000000000000B / */
/*     DATA SMALL(2) / 000000000000000000000B / */

/*     DATA LARGE(1) / 577767777777777777777B / */
/*     DATA LARGE(2) / 000007777777777777774B / */

/*     DATA RIGHT(1) / 376434000000000000000B / */
/*     DATA RIGHT(2) / 000000000000000000000B / */

/*     DATA DIVER(1) / 376444000000000000000B / */
/*     DATA DIVER(2) / 000000000000000000000B / */

/*     DATA LOG10(1) / 377774642023241175717B / */
/*     DATA LOG10(2) / 000007571421742254654B / */


/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 / */
/*     DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF / */
/*     DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 / */
/*     DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 / */
/*     DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF / */

/*     MACHINE CONSTATNS FOR THE IBM PC FAMILY (D. KAHANER NBS) */


/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR). */

/*     DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 / */
/*     DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 / */
/*     DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 / */
/*     DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 / */
/*     DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR). */

/*     DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 / */
/*     DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 / */
/*     DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 / */
/*     DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 / */
/*     DATA LOG10(1),LOG10(2) / "177464202324, "476747767461 / */


/*     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP, */
/*       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO */
/*       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.) */

/*      DATA SMALL(1),SMALL(2) / X'00100000', X'00000000' / */
/*      DATA LARGE(1),LARGE(2) / X'7FEFFFFF', X'FFFFFFFF' / */
/*      DATA RIGHT(1),RIGHT(2) / X'3CA00000', X'00000000' / */
/*      DATA DIVER(1),DIVER(2) / X'3CB00000', X'00000000' / */
/*      DATA LOG10(1),LOG10(2) / X'3FD34413', X'509F79FF' / */


/*     MACHINE CONSTANTS FOR VAX 11/780 */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS*** */

/*     DATA SMALL(1), SMALL(2) /        128,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /       9344,           0 / */
/*     DATA DIVER(1), DIVER(2) /       9472,           0 / */
/*     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 / */

/*    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS*** */
/*     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB / */

/*   MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING) */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS*** */

/*     DATA SMALL(1), SMALL(2) /         16,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /      15552,           0 / */
/*     DATA DIVER(1), DIVER(2) /      15568,           0 / */
/*     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 / */

/*    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS*** */
/*     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F / */


/* ***FIRST EXECUTABLE STATEMENT  D1MACH */
    if (*i__ < 1 || *i__ > 5) {
	xerror_("D1MACH -- I OUT OF BOUNDS", &c__25, &c__1, &c__2, (ftnlen)25)
		;
    }

    ret_val = dmach[*i__ - 1];
    return ret_val;

} /* d1mach_ */

#undef right
#undef diver
#undef small
#undef large
#undef dmach
#undef log10


integer i1mach_(integer *i__)
{
    /* Initialized data */

    static struct {
	integer e_1[16];
	} equiv_0 = { 5, 6, 0, 6, 32, 4, 2, 31, 2147483647, 2, 24, -125, 127, 
		53, -1021, 1023 };


    /* System generated locals */
    integer ret_val;

    /* Local variables */
#define imach ((integer *)&equiv_0)
    extern /* Subroutine */ int xerror_(char *, integer *, integer *, integer 
	    *, ftnlen);
#define output ((integer *)&equiv_0 + 3)

/* ***BEGIN PROLOGUE  I1MACH */
/* ***DATE WRITTEN   750101   (YYMMDD) */
/* ***REVISION DATE  840405   (YYMMDD) */
/* ***CATEGORY NO.  R1 */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  FOX, P. A., (BELL LABS) */
/*           HALL, A. D., (BELL LABS) */
/*           SCHRYER, N. L., (BELL LABS) */
/* ***PURPOSE  Returns integer machine dependent constants */
/* ***DESCRIPTION */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   These machine constant routines must be activated for */
/*   a particular environment. */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*     I1MACH can be used to obtain machine-dependent parameters */
/*     for the local machine environment.  It is a function */
/*     subroutine with one (input) argument, and can be called */
/*     as follows, for example */

/*          K = I1MACH(I) */

/*     where I=1,...,16.  The (output) value of K above is */
/*     determined by the (input) value of I.  The results for */
/*     various values of I are discussed below. */

/*  I/O unit numbers. */
/*    I1MACH( 1) = the standard input unit. */
/*    I1MACH( 2) = the standard output unit. */
/*    I1MACH( 3) = the standard punch unit. */
/*    I1MACH( 4) = the standard error message unit. */

/*  Words. */
/*    I1MACH( 5) = the number of bits per integer storage unit. */
/*    I1MACH( 6) = the number of characters per integer storage unit. */

/*  Integers. */
/*    assume integers are represented in the S-digit, base-A form */

/*               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) */

/*               where 0 .LE. X(I) .LT. A for I=0,...,S-1. */
/*    I1MACH( 7) = A, the base. */
/*    I1MACH( 8) = S, the number of base-A digits. */
/*    I1MACH( 9) = A**S - 1, the largest magnitude. */

/*  Floating-Point Numbers. */
/*    Assume floating-point numbers are represented in the T-digit, */
/*    base-B form */
/*               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*               where 0 .LE. X(I) .LT. B for I=1,...,T, */
/*               0 .LT. X(1), and EMIN .LE. E .LE. EMAX. */
/*    I1MACH(10) = B, the base. */

/*  Single-Precision */
/*    I1MACH(11) = T, the number of base-B digits. */
/*    I1MACH(12) = EMIN, the smallest exponent E. */
/*    I1MACH(13) = EMAX, the largest exponent E. */

/*  Double-Precision */
/*    I1MACH(14) = T, the number of base-B digits. */
/*    I1MACH(15) = EMIN, the smallest exponent E. */
/*    I1MACH(16) = EMAX, the largest exponent E. */

/*  To alter this function for a particular environment, */
/*  the desired set of DATA statements should be activated by */
/*  removing the C from column 1.  Also, the values of */
/*  I1MACH(1) - I1MACH(4) should be checked for consistency */
/*  with the local operating system. */
/* ***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A */
/*                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL */
/*                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188. */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  I1MACH */



/*     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5). */

/*      DATA IMACH( 1) /    5 / */
/*      DATA IMACH( 2) /    6 / */
/*      DATA IMACH( 3) /    7 / */
/*      DATA IMACH( 4) /    6 / */
/*      DATA IMACH( 5) /   60 / */
/*      DATA IMACH( 6) /   10 / */
/*      DATA IMACH( 7) /    2 / */
/*      DATA IMACH( 8) /   48 / */
/*      DATA IMACH( 9) / O"00007777777777777777" / */
/*      DATA IMACH(10) /    2 / */
/*      DATA IMACH(11) /   48 / */
/*      DATA IMACH(12) / -974 / */
/*      DATA IMACH(13) / 1070 / */
/*      DATA IMACH(14) /   96 / */
/*      DATA IMACH(15) / -927 / */
/*      DATA IMACH(16) / 1070 / */

/*     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES */

/*     DATA IMACH( 1) /      5 / */
/*     DATA IMACH( 2) /      6 / */
/*     DATA IMACH( 3) /      7 / */
/*     DATA IMACH( 4) /      6 / */
/*     DATA IMACH( 5) /     64 / */
/*     DATA IMACH( 6) /      8 / */
/*     DATA IMACH( 7) /      2 / */
/*     DATA IMACH( 8) /     47 / */
/*     DATA IMACH( 9) / X'00007FFFFFFFFFFF' / */
/*     DATA IMACH(10) /      2 / */
/*     DATA IMACH(11) /     47 / */
/*     DATA IMACH(12) / -28625 / */
/*     DATA IMACH(13) /  28718 / */
/*     DATA IMACH(14) /     94 / */
/*     DATA IMACH(15) / -28625 / */
/*     DATA IMACH(16) /  28718 / */


/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES. */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    7 / */
/*     DATA IMACH( 4) /6LOUTPUT/ */
/*     DATA IMACH( 5) /   60 / */
/*     DATA IMACH( 6) /   10 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   48 / */
/*     DATA IMACH( 9) / 00007777777777777777B / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   47 / */
/*     DATA IMACH(12) / -929 / */
/*     DATA IMACH(13) / 1070 / */
/*     DATA IMACH(14) /   94 / */
/*     DATA IMACH(15) / -929 / */
/*     DATA IMACH(16) / 1069 / */

/*     MACHINE CONSTANTS FOR THE CRAY 1 */

/*     DATA IMACH( 1) /   100 / */
/*     DATA IMACH( 2) /   101 / */
/*     DATA IMACH( 3) /   102 / */
/*     DATA IMACH( 4) /   101 / */
/*     DATA IMACH( 5) /    64 / */
/*     DATA IMACH( 6) /     8 / */
/*     DATA IMACH( 7) /     2 / */
/*     DATA IMACH( 8) /    63 / */
/*     DATA IMACH( 9) /  777777777777777777777B / */
/*     DATA IMACH(10) /     2 / */
/*     DATA IMACH(11) /    47 / */
/*     DATA IMACH(12) / -8189 / */
/*     DATA IMACH(13) /  8190 / */
/*     DATA IMACH(14) /    94 / */
/*     DATA IMACH(15) / -8099 / */
/*     DATA IMACH(16) /  8190 / */


/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA IMACH( 1) /   5 / */
/*     DATA IMACH( 2) /   6 / */
/*     DATA IMACH( 3) /   7 / */
/*     DATA IMACH( 4) /   6 / */
/*     DATA IMACH( 5) /  32 / */
/*     DATA IMACH( 6) /   4 / */
/*     DATA IMACH( 7) /  16 / */
/*     DATA IMACH( 8) /  31 / */
/*     DATA IMACH( 9) / Z7FFFFFFF / */
/*     DATA IMACH(10) /  16 / */
/*     DATA IMACH(11) /   6 / */
/*     DATA IMACH(12) / -64 / */
/*     DATA IMACH(13) /  63 / */
/*     DATA IMACH(14) /  14 / */
/*     DATA IMACH(15) / -64 / */
/*     DATA IMACH(16) /  63 / */

/*     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS) */

/*               NOTE! I1MACH(3) IS NOT WELL DEFINED AND IS SET TO ZERO. */


/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR). */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    5 / */
/*     DATA IMACH( 4) /    6 / */
/*     DATA IMACH( 5) /   36 / */
/*     DATA IMACH( 6) /    5 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   27 / */
/*     DATA IMACH(12) / -128 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   54 / */
/*     DATA IMACH(15) / -101 / */
/*     DATA IMACH(16) /  127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR). */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    5 / */
/*     DATA IMACH( 4) /    6 / */
/*     DATA IMACH( 5) /   36 / */
/*     DATA IMACH( 6) /    5 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   27 / */
/*     DATA IMACH(12) / -128 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   62 / */
/*     DATA IMACH(15) / -128 / */
/*     DATA IMACH(16) /  127 / */


/*     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP, */
/*       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO */
/*       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.) */

/*      DATA IMACH( 1) /    5 / */
/*      DATA IMACH( 2) /    6 / */
/*      DATA IMACH( 3) /    6 / */
/*      DATA IMACH( 4) /    0 / */
/*      DATA IMACH( 5) /   32 / */
/*      DATA IMACH( 6) /    4 / */
/*      DATA IMACH( 7) /    2 / */
/*      DATA IMACH( 8) /   31 / */
/*      DATA IMACH( 9) / 2147483647 / */
/*      DATA IMACH(10) /    2 / */
/*      DATA IMACH(11) /   24 / */
/*      DATA IMACH(12) / -125 / */
/*      DATA IMACH(13) /  128 / */
/*      DATA IMACH(14) /   53 / */
/*      DATA IMACH(15) / -1021 / */
/*      DATA IMACH(16) /  1024 / */


/*     MACHINE CONSTANTS FOR THE VAX 11/780 */

/*     DATA IMACH(1) /    5 / */
/*     DATA IMACH(2) /    6 / */
/*     DATA IMACH(3) /    5 / */
/*     DATA IMACH(4) /    6 / */
/*     DATA IMACH(5) /   32 / */
/*     DATA IMACH(6) /    4 / */
/*     DATA IMACH(7) /    2 / */
/*     DATA IMACH(8) /   31 / */
/*     DATA IMACH(9) /2147483647 / */
/*     DATA IMACH(10)/    2 / */
/*     DATA IMACH(11)/   24 / */
/*     DATA IMACH(12)/ -127 / */
/*     DATA IMACH(13)/  127 / */
/*     DATA IMACH(14)/   56 / */
/*     DATA IMACH(15)/ -127 / */
/*     DATA IMACH(16)/  127 / */

/* ***FIRST EXECUTABLE STATEMENT  I1MACH */
    if (*i__ < 1 || *i__ > 16) {
	xerror_("I1MACH -- I OUT OF BOUNDS", &c__25, &c__1, &c__2, (ftnlen)25)
		;
    }

    ret_val = imach[*i__ - 1];
    return ret_val;

} /* i1mach_ */

#undef output
#undef imach


