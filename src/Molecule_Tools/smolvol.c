/* Molecule_Tools/smolvol.f -- translated by f2c (version 20200916).
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

/* Common Block Declarations */

struct {
    doublereal x_chain__[301], y_chain__[301], z_chain__[301], radii[301], 
	    refradius;
    integer n_spheres__;
} rawdata_;

#define rawdata_1 rawdata_

struct {
    doublereal nx[11250], ny[11250], nz[11250], magnitude[22500]	/* 
	    was [11250][2] */, rho_sq__[11250];
    integer n_planes__[301], plane_list__[271803]	/* was [301][301][3] 
	    */, point[22500]	/* was [11250][2] */, planef[11250], planecnt;
} planeinfo_;

#define planeinfo_1 planeinfo_

struct {
    doublereal areag[22500]	/* was [11250][2] */;
    integer numarcsg[22500]	/* was [11250][2] */, trueedgesg[22500]	/* 
	    was [11250][2] */;
} global_;

#define global_1 global_

struct {
    doublereal radius, rsq, o_radius__, o_rsq__, o_rcu__;
    integer sphere, numplanes, tripcnt;
} sphereinfo_;

#define sphereinfo_1 sphereinfo_

struct {
    doublereal atomic_volume__[301], atomic_area__[301];
} volumeoutput_;

#define volumeoutput_1 volumeoutput_

struct {
    doublereal xa[3], xb[3];
    integer deadplane;
    logical wipeout, interior, remove;
} conclusion_;

#define conclusion_1 conclusion_

struct {
    doublereal nxs[6000], nys[6000], nzs[6000], v_mags__[6000], rho_sqs__[
	    6000];
    integer indexs[6000], planes[6000], points[6000], ith, jth;
} spheresplanes_;

#define spheresplanes_1 spheresplanes_

struct {
    doublereal apexx[6000], apexy[6000], apexz[6000];
    integer apexf[12000]	/* was [6000][2] */, n_edges__[300], 
	    edge_list__[5400000]	/* was [300][6000][3] */, apexcnt;
} apexinfo_;

#define apexinfo_1 apexinfo_

struct {
    integer twin[6000];
} twinlist_;

#define twinlist_1 twinlist_

struct {
    integer plane_order__[18000]	/* was [6000][3] */;
} debug_;

#define debug_1 debug_

struct {
    doublereal ni[3], vi_mag__, rhoi_sq__, curvei, o_rhoi_sq__, r_rhoi__, 
	    vi_rhoi__;
    integer edgesi[12000]	/* was [6000][2] */, indexi, planei, pointi, 
	    planenum, numedges;
} thisplane_;

#define thisplane_1 thisplane_

struct {
    doublereal vertexx[300], vertexy[300], vertexz[300], ngx[300], ngy[300], 
	    ngz[300], tangeo[300];
    integer arc_list__[600]	/* was [300][2] */, arc_map__[300], arccnt, 
	    vertexcnt;
} arcinfo_;

#define arcinfo_1 arcinfo_

struct {
    doublereal etad[2];
} diangleinfo_;

#define diangleinfo_1 diangleinfo_

struct {
    doublereal eab[3], ecd[3], abxcd[3], ra[3], rc[3], rmagab, rmagcd;
    integer apexa, apexb, apexc, apexd, edgei, edgej, ith, jth, kth;
} edgepair_;

#define edgepair_1 edgepair_

struct {
    doublereal lambda[2], lambdaf, lambdas, nu1;
    integer first, second, movement;
} formove_;

#define formove_1 formove_

struct {
    integer numstart[5], numend[5], locstart[30000]	/* was [6000][5] */, 
	    locend[30000]	/* was [6000][5] */, edge_num__[18003]	/* 
	    was [6001][3] */;
} sortedout_;

#define sortedout_1 sortedout_

struct {
    integer local_list__[12000]	/* was [6000][2] */, numarcs;
} localarcinfo_;

#define localarcinfo_1 localarcinfo_

struct {
    doublereal ra1[3], axb[3];
} foraccum_;

#define foraccum_1 foraccum_

struct {
    integer festoon[9000]	/* was [300][30] */, festarcs[30], festooncnt;
} festooninfo_;

#define festooninfo_1 festooninfo_

/* Table of constant values */

static doublereal c_b4 = 1.;
static integer c__1 = 1;
static doublereal c_b31 = .5;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__5 = 5;

/*     volume.f - volume determination code */

/*     Author: Lawrence R. Dodd <dodd@roebling.poly.edu> */
/*             Doros N. Theodorou <doros@pylos.cchem.berkeley.edu> */
/*     Maintainer: Lawrence R. Dodd <dodd@roebling.poly.edu> */
/*     Created: March 21, 1990 */
/*     Version: 2.1 */
/*     Date: 1994/09/29 22:16:57 */
/*     Keywrds: volume and area determination */
/*     Time-stamp: <94/09/29 18:15:01 dodd> */
/*     Copyright (c) 1990, 1991, 1992, 1993, 1994 */
/*     by Lawrence R. Dodd and Doros N. Theodorou. */
/*     This program is free software; you can redistribute it and/or */
/*     modify it under the terms of the GNU General Public License as */
/*     published by the Free Software Foundation; either version 2 of the */
/*     License, or (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/*     General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */
/* ---------------------------------------------------------------------C */
/*                     Plane Sphere Intersections                      C */
/* ---------------------------------------------------------------------C */
/*     This program will find the total and individual volume and      C */
/*     exposed surface area of an arbitrary collection of spheres of   C */
/*     arbitrary radii cut by an arbitrary collection of planes        C */
/*     analytically by analyzing the plane/sphere intersections.       C */
/* ---------------------------------------------------------------------C */
/*     Algorithm by: Doros N. Theodorou and Lawrence R. Dodd           C */
/*     Coded by: L.R. Dodd                                             C */
/* ---------------------------------------------------------------------C */
/*     Created on: March 21, 1990                                      C */
/*       Phase 1 Completed on: March 23, 1990                          C */
/*       Phase 2 Completed on: April 16, 1990                          C */
/*       Phase 3 Completed on: May   17, 1990                          C */
/*       Phase 4 Completed on: June   5, 1990                          C */
/*       Phase 5 Completed on: July  26, 1990                          C */
/* ---------------------------------------------------------------------C */
/*     Reference:                                                      C */
/*                                                                     C */
/*       "Analytical treatment of the volume and surface area of       C */
/*       molecules formed by an arbitrary collection of unequal        C */
/*       spheres intersected by planes"                                C */
/*                                                                     C */
/*     L.R. Dodd and D.N. Theodorou                                    C */
/*     MOLECULAR PHYSICS, Volume 72, Number 6, 1313-1345, April 1991   C */
/* ---------------------------------------------------------------------C */
/*     Acknowlegement:                                                 C */
/*                                                                     C */
/*     LRD wishes to thank his mentor DNT for a stimulating and        C */
/*     enjoyable post-doctoral experience.                             C */
/* ---------------------------------------------------------------------C */
/*     General Notes On Program:                                       C */
/*                                                                     C */
/*     This program has been written with an eye towards both          C */
/*     efficiency and clarity. On a philosophical note, many believe   C */
/*     that these ideals are mutually exclusive but in general they    C */
/*     are not. There are, however, a few instances where one ideal    C */
/*     has been given more prominence over the other. The comments in  C */
/*     the program, together with the associated journal article,      C */
/*     should help to explain any apparent logical leaps in the        C */
/*     algorithm.                                                      C */
/*                                                                     C */
/*     The program was intended to be used as a subroutine called      C */
/*     repeatly by some main program. In this case the subroutine      C */
/*     "VOLUME" is called by some main routine which has placed the    C */
/*     necessary information in common block /Raw Data/. The answers   C */
/*     are returned in common block /Volume Output/. I must apologize  C */
/*     for the poor input/output for the program. For example, the     C */
/*     area/volume of each sphere is not placed in /Volume Output/.    C */
/*                                                                     C */
/*     This program was developed on a Sun SPARCstation 330 using Sun  C */
/*     FORTRAN 1.3.1 (all trademarks of Sun Microsystems, Inc.). We    C */
/*     have used some of extensions to the ANSI standard including:    C */
/*                                                                     C */
/*         o  long variable names (i.e., more than six characters)     C */
/*         o  variable names containing the characters '_dollar_' and '_'     C */
/*         o  END DO used in place of the CONTINUE statement           C */
/*         o  DO-WHILE used in place of IF-GOTO constructs             C */
/*         o  excessive number of continuation lines in some FORMATs   C */
/*         o  generic intrinsic function calls (e.g., SIN for DSIN)    C */
/*         o  IMPLICIT NONE statement (needed in development)          C */
/*                                                                     C */
/*     The advantage of using non-standard FORTRAN is that it makes it C */
/*     considerably easier to follow the flow of a program. There are  C */
/*     no extraneous statement labels in this program that may have    C */
/*     obscured the logic (not a single GOTO was used). The previews   C */
/*     of the new F90 standard appear to adopt many of the features    C */
/*     already implemented in VMS, Sun, Cray, and IBM FORTRAN.         C */
/*                                                                     C */
/*     Note that this algorithm is completely parallelizable.          C */
/*                                                                     C */
/*                           Larry Dodd                                C */
/*                           dodd@mycenae.cchem.berkeley.edu           C */
/*                                                                     C */
/*                           Department of Chemical Engineering        C */
/*                           College of Chemistry                      C */
/*                           University of California at Berkeley      C */
/*                           Berkeley, California 94720-9989           C */
/*                           (415) 643-7691 (LRD)                      C */
/*                           (415) 643-8523 (DNT)                      C */
/*                           (415) 642-5927 (Lab)                      C */
/*                                                                     C */
/*                            dodd@mycenae.cchem.berkeley.edu          C */
/*                           doros@mycenae.cchem.berkeley.edu          C */
/*                                                                     C */
/* ---------------------------------------------------------------------C */
/*     Note:                                                           C */
/*       Plane_Ordering of common block /Debug/ is, as the name        C */
/*       implies, for debugging purposes only as is routine ORDERING.  C */
/*       The information contain therein is not necessary for solving  C */
/*       the sphere plane problem but proved incredibly useful during  C */
/*       program development.                                          C */
/* ---------------------------------------------------------------------C */
/* --------------------------------------------------------------------C */
/*     Main Subroutine                                                C */
/* --------------------------------------------------------------------C */
/* Subroutine */ int volume_(logical *error_encountered__)
{
    extern /* Subroutine */ int find_planes__(void), all_spheres__(logical *);

    find_planes__();
    all_spheres__(error_encountered__);
    return 0;
} /* volume_ */

/* ---------------------------------------------------------------------C */
/*     Find_Planes                                                     C */
/* ---------------------------------------------------------------------C */
/*     Purpose:                                                        C */
/*        This routine finds the planes formed by the intersection of  C */
/*     all possible sphere pairs.  The plane information is stored for C */
/*     both spheres forming a plane because eventually the problem     C */
/*     will be decoupled completely. For each sphere, the information  C */
/*     stored for the plane includes the unit normal vector pointing   C */
/*     away from the cutout of the sphere, the magnitude of the vector C */
/*     connecting the sphere center normal to the plane, a variable    C */
/*     called Point, which is +1 if the sphere center is contained in  C */
/*     the cutout for this plane and -1 if the sphere center is not    C */
/*     contained the plane's cutout (if the plane goes through the     C */
/*     center of the sphere then Point is zero), the radius of the     C */
/*     circle inscribed by the plane on the sphere, and a flag         C */
/*     indicating whether the plane has been removed.                  C */
/* ---------------------------------------------------------------------C */
/*     Called by main program.                                         C */
/*     Output stored in /Plane Info/.                                  C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int find_planes__(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal rjsq_risq__, ri, rj, xi, yi, zi, dij;
    static integer n_spheresm1__;
    static doublereal delr, rdij;
    static integer sphi, sphj;
    static doublereal sumr, del_x__, del_y__, del_z__;
    static integer plane;
    static doublereal ri_sq__, termi, termj, dij_sq__, atermi, atermj;

/*     < LOCAL */
/*     < Initialize. */
    plane = 0;
/*     write statements */
/*      Write(0,*) ' Number of Spheres= ', n_Spheres */
/*      Do I= 0, n_Spheres-1 */
/*        write (0, *) 'x=', x_Chain(I),'  y=', y_Chain(I) ,'  z=' */
/*     >, z_Chain(I) */
/*      End Do */
/*     end of write statement */
    if (rawdata_1.n_spheres__ == 0) {
	s_stop("no spheres", (ftnlen)10);
    }
/* JW      Do I= 0, n_Spheres */
    i__1 = rawdata_1.n_spheres__ - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	planeinfo_1.n_planes__[i__] = 0;
    }
/*     < Sweep through all pairs of spheres to find planes. */
    n_spheresm1__ = rawdata_1.n_spheres__ - 1;
    i__1 = n_spheresm1__;
    for (sphi = 0; sphi <= i__1; ++sphi) {
/*       < Has sphere SphI been removed? */
	if (planeinfo_1.n_planes__[sphi] != -1) {
/*         < Pull out coordinates and radius of sphere SphI, ri and Ri */
	    xi = rawdata_1.x_chain__[sphi];
	    yi = rawdata_1.y_chain__[sphi];
	    zi = rawdata_1.z_chain__[sphi];
	    ri = rawdata_1.radii[sphi];
	    ri_sq__ = ri * ri;
/*         < Initialize */
	    sphj = sphi;
/*         < Loop through all SphJ > SphI until SphJ= n_Spheres */
/*         < or until sphere SphI has been removed */
/*         < Skip any spheres SphJ that have been removed */
/* JW          Do While( SphJ .Lt. n_Spheres ) */
	    while(sphj < rawdata_1.n_spheres__ - 1) {
		++sphj;
/*           < Has sphere SphJ been removed? */
		if (planeinfo_1.n_planes__[sphj] != -1) {
/*             < Find the center to center vector, rij= ri - rj */
		    del_x__ = rawdata_1.x_chain__[sphj] - xi;
		    del_y__ = rawdata_1.y_chain__[sphj] - yi;
		    del_z__ = rawdata_1.z_chain__[sphj] - zi;
/*             < Find the squared magnitude of rij, |rij|^2 */
		    dij_sq__ = del_x__ * del_x__ + del_y__ * del_y__ + 
			    del_z__ * del_z__;
/*             < Find the sum of the radii, Rj + Ri */
		    rj = rawdata_1.radii[sphj];
		    sumr = rj + ri;
/*             < Do we have an intersection? Note that we make the */
/*             < criterion for an intersection a little tighter by */
/*             < putting a negative on the epsilon.  That is, dij^2 has */
/*             < to be less than (Ri + Rj)^2 MINUS epsilon. */
		    if (dij_sq__ - sumr * sumr < -1e-8) {
/*               < The spheres intersect. Is one sphere completely */
/*               < embedded within the other sphere?  Find the */
/*               < difference of the radii, Rj - Ri */
			delr = rj - ri;
/*               < Is one sphere inside the other? Note that we make */
/*               < the criterion for embedment a little looser by */
/*               < making the epsilon positive.  That is, dij^2 has to */
/*               < be less than (Ri - Rj)^2 PLUS epsilon. */
			if (dij_sq__ - delr * delr < 1e-8) {
/*                 < One sphere is completely embedded within the other */
			    if (delr > 0.) {
/*                   < Sphere SphI is in sphere SphJ. Go on to next */
/*                   < SphI. This sphere SphI will not show up again */
/*                   < Turn-off all of sphere SphI's planes and mark */
/*                   < them as having died in trying to form the */
/*                   < planes. */
				i__2 = planeinfo_1.n_planes__[sphi];
				for (k = 1; k <= i__2; ++k) {
				    i__ = planeinfo_1.plane_list__[sphi + (k 
					    + 301) * 301 - 90601];
				    j = planeinfo_1.plane_list__[sphi + (k + 
					    602) * 301 - 90601];
				    planeinfo_1.planef[i__ - 1] = -1;
				    global_1.numarcsg[i__ + j * 11250 - 11251]
					     = -2;
				    global_1.trueedgesg[i__ + j * 11250 - 
					    11251] = -2;
				}
/*                   < Turn-off sphere SphI itself */
				planeinfo_1.n_planes__[sphi] = -1;
/*                   < Get out of inner loop */
/* JW                    SphJ= n_Spheres */
				sphj = rawdata_1.n_spheres__ - 1;
			    } else {
/*                   < Sphere SphJ is in sphere SphI. Go on to next */
/*                   < SphJ.  Next time SphJ shows up we will skip it */
/*                   < Turn-off all of sphere SphJ's planes and mark */
/*                   < them as having died in trying to form the */
/*                   < planes. */
				i__2 = planeinfo_1.n_planes__[sphj];
				for (k = 1; k <= i__2; ++k) {
				    i__ = planeinfo_1.plane_list__[sphj + (k 
					    + 301) * 301 - 90601];
				    j = planeinfo_1.plane_list__[sphj + (k + 
					    602) * 301 - 90601];
				    planeinfo_1.planef[i__ - 1] = -1;
				    global_1.numarcsg[i__ + j * 11250 - 11251]
					     = -2;
				    global_1.trueedgesg[i__ + j * 11250 - 
					    11251] = -2;
				}
/*                   < Turn-off sphere SphJ itself */
				planeinfo_1.n_planes__[sphj] = -1;
			    }
			} else {
/*                 < We have a non-trivial intersection of spheres SphI */
/*                 < and SphJ. Increment the plane count for both */
/*                 < spheres. */
			    ++planeinfo_1.n_planes__[sphi];
			    ++planeinfo_1.n_planes__[sphj];
/*                 < Increment the plane number */
			    ++plane;
/*                 < Flag this plane as existing */
			    planeinfo_1.planef[plane - 1] = 1;
/*                 < Store the plane associated with spheres SphI and SphJ */
/*                 < Also store location in normal vector array */
			    planeinfo_1.plane_list__[sphi + (
				    planeinfo_1.n_planes__[sphi] + 301) * 301 
				    - 90601] = plane;
			    planeinfo_1.plane_list__[sphi + (
				    planeinfo_1.n_planes__[sphi] + 602) * 301 
				    - 90601] = 1;
			    planeinfo_1.plane_list__[sphi + (
				    planeinfo_1.n_planes__[sphi] + 903) * 301 
				    - 90601] = sphj;
			    planeinfo_1.plane_list__[sphj + (
				    planeinfo_1.n_planes__[sphj] + 301) * 301 
				    - 90601] = plane;
			    planeinfo_1.plane_list__[sphj + (
				    planeinfo_1.n_planes__[sphj] + 602) * 301 
				    - 90601] = 2;
			    planeinfo_1.plane_list__[sphj + (
				    planeinfo_1.n_planes__[sphj] + 903) * 301 
				    - 90601] = sphi;
/*                 < Find the unit normal plane vectors, n.  These are */
/*                 < the unit normal vector pointing outside cutout for */
/*                 < each sphere.  It is defined as */
/*                 < */
/*                 <      ni= rij/|rij| */
/*                 <      nj= rji/|rji|= -rij/|rij|= -ni */
/*                 < */
/*                 < where rij is defined as rj - ri and is the vector */
/*                 < connecting the center of sphere SphI to the center */
/*                 < of sphere SphJ and where |rij| is simply dij the */
/*                 < center-to-center distance.  These vectors are */
/*                 < uneffected by location of the plane of intersection */
/*                 < of the two spheres.  The vectors connecting the */
/*                 < respective sphere centers to the normal point on */
/*                 < the plane of intesection are called vi and vj.  We */
/*                 < do not form the v's explicitly since they may or */
/*                 < may not point out of the cutouts but we do */
/*                 < calculate the magnitudes of the v's and record */
/*                 < where the sphere centers lie relative to their */
/*                 < respective cutouts.  We will define a variable */
/*                 < called Point to do this the values of which will be */
/*                 < */
/*                 <   If Point(Pi) = -1      ===> */
/*                 <       SphI's center is NOT in its own cutout */
/*                 <       ri is on the 'negative' side of plane (P) */
/*                 <       ni= -vi/|vi| */
/*                 <   Else If Point(Pi) +1   ===> */
/*                 <       SphI's center is in its own cutout */
/*                 <       ri is on the 'positive' side of plane (P) */
/*                 <       ni= vi/|vi| */
/*                 <   Else If Point(Pi) = 0 ===> */
/*                 <       Plane goes through SphI's center */
/*                 <       rk and ri are the same point */
/*                 <       |vi|= 0 and ni := -vj */
/*                 < */
/*                 < Both are positive or one (and only one) is negative */
			    dij = sqrt(dij_sq__);
			    rdij = 1. / dij;
/*                 < Find the unit normal vector for the leading sphere */
/*                 < (i.e., the sphere with the smaller number).  The */
/*                 < unit normal vector for the other sphere is simply */
/*                 < the negative of the stored unit normal vector.  We */
/*                 < will save on memory by just storing the unit */
/*                 < normal for one side of the plane and then finding */
/*                 < the other unit normal vector when needed by */
/*                 < multiplying this one by minus one (see Find_Apices). */
			    planeinfo_1.nx[plane - 1] = del_x__ * rdij;
			    planeinfo_1.ny[plane - 1] = del_y__ * rdij;
			    planeinfo_1.nz[plane - 1] = del_z__ * rdij;
/*                 < Initialize the global lists for consistency checks. */
			    global_1.numarcsg[plane - 1] = -999;
			    global_1.numarcsg[plane + 11249] = -999;
			    global_1.trueedgesg[plane - 1] = -999;
			    global_1.trueedgesg[plane + 11249] = -999;
/*                 < Determine the proper vector magnitudes, the */
/*                 < location of the plane of intersection and the */
/*                 < radius of the inscribed circle. */
			    rjsq_risq__ = delr * sumr;
			    termi = (dij_sq__ - rjsq_risq__) * .5;
			    termj = (dij_sq__ + rjsq_risq__) * .5;
			    atermi = abs(termi);
			    atermj = abs(termj);
			    if (atermi < 1e-8) {
/*                   < If |Termi| is less than Eps then we know that */
/*                   < dij^2 = (Rj^2 - Ri^2) implying that the triangle */
/*                   < formed by dij, Rj and Ri is a right triangle with */
/*                   < Rj as the hypotenuse.  This means that the plane */
/*                   < of intersection (P) goes through SphI's center */
/*                   < ri and that rk (the normal point on the plane) */
/*                   < coincide. */
/*                   < Topology Note */
/*                   < */
/*                   < Plane (P) goes through SphI's Center, ri, the */
/*                   < normal vector vi is zero and so the magnitude of */
/*                   < vi is zero and the magnitude of vj is dij.  The */
/*                   < radius of the circle inscribed by plane on each */
/*                   < sphere is Ri.  Finally, Point(i) is defined as */
/*                   < zero and Point(j) is +1. */
/*                   < */
/*                   <    Plane (P) */
/*                   <           | */
/*                   <           | */
/*                   <           | */
/*                   <           |   ni, vi= 0 */
/*                   <     ---   o------------> */
/*                   <      |    |      nj, vj */
/*                   < rho= Ri   |   <----------o */
/*                   <      |    | */
/*                   <     --- --x--------------x----- */
/*                   <           ri              rj */
/*                   <           rk */
/*                   < */
/*                   <           |               | */
/*                   <           |               | */
/*                   <           |----- dij -----| */
/*                   <           |               | */
				planeinfo_1.point[plane - 1] = 0;
				planeinfo_1.point[plane + 11249] = 1;
				planeinfo_1.magnitude[plane - 1] = 0.;
				planeinfo_1.magnitude[plane + 11249] = dij;
				planeinfo_1.rho_sq__[plane - 1] = ri_sq__;
			    } else if (atermj < 1e-8) {
/*                   < If |Termj| is less than Eps then we know that */
/*                   < dij^2 = (Ri^2 - Rj^2) implying that the triangle */
/*                   < formed by dij, Rj and Ri is a right triangle with */
/*                   < Ri as the hypotenuse.  This means that the plane */
/*                   < of intersection (P) goes through SphJ's center */
/*                   < rj and that rk (the normal point on the plane) */
/*                   < coincide. */
/*                   < Topology Note */
/*                   < */
/*                   < Plane (P) goes through SphJ's Center, rj, the */
/*                   < normal vector vj is zero and so the magnitude of */
/*                   < vj is zero and the magnitude of vi is dij.  The */
/*                   < radius of the circle inscribed by plane on each */
/*                   < sphere is Rj. Finally, Point(j) is defined as */
/*                   < zero and Point(i) is +1. */
/*                   < */
/*                   <                     Plane (P) */
/*                   <                          | */
/*                   <                          | */
/*                   <                          | */
/*                   <                 nj, vj=0 | */
/*                   <               <----------o   --- */
/*                   <              ni, vi      |    | */
/*                   <           o--------->    |   Rj = rho */
/*                   <                          |    | */
/*                   <         --x--------------x-- --- */
/*                   <           ri             rj */
/*                   <                          rk */
/*                   < */
/*                   <           |               | */
/*                   <           |               | */
/*                   <           |----- dij -----| */
/*                   <           |               | */
				planeinfo_1.point[plane - 1] = 1;
				planeinfo_1.point[plane + 11249] = 0;
				planeinfo_1.magnitude[plane - 1] = dij;
				planeinfo_1.magnitude[plane + 11249] = 0.;
				planeinfo_1.rho_sq__[plane - 1] = rj * rj;
			    } else {
/*                   < The plane of intersection does not pass through */
/*                   < the center of either sphere.  We have a three */
/*                   < possible cases: the plane is between the sphere */
/*                   < centers meaning that the center of each sphere is */
/*                   < in its own cutout, the center of sphere SphI is */
/*                   < not in its own cutout or the center of sphere */
/*                   < SphJ is not in its own cutout. */
/*                   < */
/*                   < Topology Note */
/*                   < */
/*                   < Three possible locations of plane (P) */
/*                   < */
/*                   <     P'             P''             P''' */
/*                   <    -|-----x--------|---------x-----|- */
/*                   <     rk'   ri       rk''      rj    rk''' */
/*                   < */
/*                   <           |                  | */
/*                   <           |--- dij= |rij| ---| */
/*                   <           |                  | */
/*                   < */
/*                   <   1. within  center to center line ('') */
/*                   <   2. outside center to center line */
/*                   <     a) such that ri is NOT in SphI's cutout (') */
/*                   <     b) such that rj is NOT in SphJ's cutout (''') */
/*                   < */
/*                   <    If Termi > 0 then the triangle formed by dij, */
/*                   < Ri and Rj is such that the angle between dij and */
/*                   < Ri, angle Thetai, is less than Pi/2 radians. */
/*                   < Therefore, rk is between ri and rj and the */
/*                   < magnitude of the vector vi is simply */
/*                   < */
/*                   <   |vi|= Ri.cos(Thetai) */
/*                   <       = Termi/dij */
/*                   < */
/*                   < This is case (1) in the above diagram. */
/*                   < */
/*                   < */
/*                   <    On the other hand, if Termi < 0 then the */
/*                   < triangle formed by dij, Ri and Rj is such that */
/*                   < the angle Thetai is greater than Pi/2 radians. */
/*                   < Therefore, rk is not between ri and rj and sphere */
/*                   < SphI's center is not in its own cutout.  In this */
/*                   < case the |vi| is Ri.cos(Thetai') where Thetai' is */
/*                   < the complementary angle of Theta */
/*                   < */
/*                   <   Thetai'= Pi - Thetai */
/*                   < */
/*                   <   cos(Thetai')= cos[Pi/2 + (Pi/2 - Thetai)] */
/*                   <               = -sin[Pi/2 - Thetai] */
/*                   <               = -cos(Thetai) */
/*                   < */
/*                   <   |vi|= Ri.cos(Thetai') */
/*                   <       = -Ri.cos(Thetai) */
/*                   <       = -Termi/dij */
/*                   < */
/*                   < This is case (2a) in the above diagram. */
/*                   < */
/*                   < So therefore, in general we may write */
/*                   < */
/*                   <   |vi|= |Termi|/dij */
/*                   < */
/*                   < Compiler Note */
/*                   << Sign(1,x)= +1 for x > 0 */
/*                   << Sign(1,x)= +1 for x = 0 */
/*                   << Sign(1,x)= -1 for x < 0 */
				d__1 = d_sign(&c_b4, &termi);
				planeinfo_1.point[plane - 1] = i_dnnt(&d__1);
				d__1 = d_sign(&c_b4, &termj);
				planeinfo_1.point[plane + 11249] = i_dnnt(&
					d__1);
				planeinfo_1.magnitude[plane - 1] = atermi * 
					rdij;
				planeinfo_1.magnitude[plane + 11249] = atermj 
					* rdij;
				planeinfo_1.rho_sq__[plane - 1] = ri_sq__ - 
					planeinfo_1.magnitude[plane - 1] * 
					planeinfo_1.magnitude[plane - 1];
			    }
			}
		    }
		}
	    }
	}
    }
/*     < Store the number of planes. */
    planeinfo_1.planecnt = plane;
/* iaw  Write(0,*) ' Number of planes= ', PlaneCnt */
/*     < Initialize the global areas for the consistency checks. */
    i__1 = planeinfo_1.planecnt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	global_1.areag[i__ - 1] = 0.;
	global_1.areag[i__ + 11249] = 0.;
    }
    return 0;
/* L100: */
} /* find_planes__ */

/* ---------------------------------------------------------------------C */
/*     All_Spheres                                                     C */
/* ---------------------------------------------------------------------C */
/*     Purpose:                                                        C */
/*        This routine calls the routine Each_Sphere for every sphere  C */
/*     in this system.                                                 C */
/* ---------------------------------------------------------------------C */
/*     Called by main program.                                         C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int all_spheres__(logical *error_encountered__)
{
    /* Format strings */
    static char fmt_200[] = "(1x,i4,3(\002 |\002,f7.3,2i4),\002 NOT EQUAL"
	    " \002)";
    static char fmt_999[] = "(/1x,\002 ********** WARNING ********** \002/"
	    "1x,\002    \002,i5,\002 Planes NOT EQUAL \002/1x,\002 **********"
	    " WARNING ********** \002/)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer notequal;
    static doublereal rcu;
    extern /* Subroutine */ int each_sphere__(doublereal *, doublereal *, 
	    logical *);
    static doublereal area[301];
    static integer trueplanecnt;
    static doublereal diffa;
    static integer diffn, difft, planei;
    static doublereal refvol, volume[301], refarea, asphere, vsphere;

    /* Fortran I/O blocks */
    static cilist io___39 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_999, 0 };


/*     passed */
/*     < LOCAL */
    *error_encountered__ = FALSE_;
/*     < Find the volume and surface area of of the reference sphere in */
/*     < the units of the radii (e.g., Angstroms). */
    refarea = rawdata_1.refradius * 12.566370614359172 * rawdata_1.refradius;
    refvol = refarea * rawdata_1.refradius / 3.;
/*     < Loop through all the spheres. */
/* JW      Do Sphere= 0, n_Spheres */
/*      write (6, *) 'HERE' */
    i__1 = rawdata_1.n_spheres__ - 1;
    for (sphereinfo_1.sphere = 0; sphereinfo_1.sphere <= i__1; 
	    ++sphereinfo_1.sphere) {
/*       < Find the number of planes intersecting sphere */
	sphereinfo_1.numplanes = planeinfo_1.n_planes__[sphereinfo_1.sphere];
/*       < How many planes do we have for this sphere? */
	if (sphereinfo_1.numplanes > 0) {
/*         < We have at least one plane. */
	    each_sphere__(&vsphere, &asphere, error_encountered__);
	    if (*error_encountered__) {
		return 0;
	    }
	    volume[sphereinfo_1.sphere] = vsphere;
	    area[sphereinfo_1.sphere] = asphere;
	} else if (sphereinfo_1.numplanes == 0) {
/*         < We have no planes intersecting this sphere.  Its volume */
/*         < and area are that of an uncut sphere. */
	    volume[sphereinfo_1.sphere] = 1.;
	    area[sphereinfo_1.sphere] = 1.;
	    sphereinfo_1.rsq = rawdata_1.radii[sphereinfo_1.sphere] * 
		    rawdata_1.radii[sphereinfo_1.sphere];
	    rcu = sphereinfo_1.rsq * rawdata_1.radii[sphereinfo_1.sphere];
	} else {
/*         < The value of n_Planes for this sphere is negative */
/*         < indicating that this sphere has been removed because it is */
/*         < completely embedded within another sphere. */
	    volume[sphereinfo_1.sphere] = 0.;
	    area[sphereinfo_1.sphere] = 0.;
	}
    }
/*      write (6, *) 'THERE' */
/*     < Find the total volume and surface area of this system of */
/*     < spheres in units of reference sphere volume and area. */
/* JW      Do Sphere= 0, n_Spheres */
    i__1 = rawdata_1.n_spheres__ - 1;
    for (sphereinfo_1.sphere = 0; sphereinfo_1.sphere <= i__1; 
	    ++sphereinfo_1.sphere) {
	sphereinfo_1.rsq = rawdata_1.radii[sphereinfo_1.sphere] * 
		rawdata_1.radii[sphereinfo_1.sphere];
	volumeoutput_1.atomic_area__[sphereinfo_1.sphere] = area[
		sphereinfo_1.sphere] * sphereinfo_1.rsq * refarea;
	rcu = sphereinfo_1.rsq * rawdata_1.radii[sphereinfo_1.sphere];
	volumeoutput_1.atomic_volume__[sphereinfo_1.sphere] = volume[
		sphereinfo_1.sphere] * rcu * refvol;
    }
/*      write (6, *) 'THERE1' */
/*     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/*     < NOTE: This section of code is only for debugging purposes and */
/*     < is probably not needed. In fact if external (or artificial) */
/*     < planes are defined in file volume.art spurious warning */
/*     < messages will be issued. */
/*     < Find out if the head and tail of each plane match. */
    notequal = 0;
    i__1 = planeinfo_1.planecnt;
    for (planei = 1; planei <= i__1; ++planei) {
	diffa = global_1.areag[planei - 1] - global_1.areag[planei + 11249];
	diffn = global_1.numarcsg[planei - 1] - global_1.numarcsg[planei + 
		11249];
	difft = global_1.trueedgesg[planei - 1] - global_1.trueedgesg[planei 
		+ 11249];
	if (abs(diffa) > 1e-6 || diffn != 0 || difft != 0) {
	    ++notequal;
	    s_wsfe(&io___39);
	    do_fio(&c__1, (char *)&planei, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&global_1.areag[planei - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&global_1.numarcsg[planei - 1], (ftnlen)
		    sizeof(integer));
	    do_fio(&c__1, (char *)&global_1.trueedgesg[planei - 1], (ftnlen)
		    sizeof(integer));
	    do_fio(&c__1, (char *)&global_1.areag[planei + 11249], (ftnlen)
		    sizeof(doublereal));
	    do_fio(&c__1, (char *)&global_1.numarcsg[planei + 11249], (ftnlen)
		    sizeof(integer));
	    do_fio(&c__1, (char *)&global_1.trueedgesg[planei + 11249], (
		    ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&diffa, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&diffn, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&difft, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
/*      write (6, *) 'THERE2' */
    trueplanecnt = 0;
    i__1 = planeinfo_1.planecnt;
    for (planei = 1; planei <= i__1; ++planei) {
	if (global_1.numarcsg[planei - 1] >= 0 || global_1.trueedgesg[planei 
		- 1] >= 0 || global_1.numarcsg[planei - 1] == -3) {
	    ++trueplanecnt;
	}
    }
/*      write (6, *) 'THERE3' */
/*     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    if (notequal != 0) {
	s_wsfe(&io___41);
	do_fio(&c__1, (char *)&notequal, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
/* L100: */
/* L150: */
/* L151: */
} /* all_spheres__ */

/* ---------------------------------------------------------------------C */
/*     Each Sphere                                                     C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int each_sphere__(doublereal *volume, doublereal *area, 
	logical *error_encountered__)
{
    /* Builtin functions */
    double d_int(doublereal *), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal volume_cp__;
    extern /* Subroutine */ int arc_polys__(doublereal *, doublereal *, 
	    doublereal *, logical *), festoonery_(doublereal *, logical *), 
	    find_apices__(void);
    static doublereal area_sphseg__, areat, realr, volume_sphsec__, area_d__, 
	    area_l__, ascale, vscale, area_sp__, realrcu, realrsq, volumet;

/*     < PASSED */
/*     < LOCAL */
/*     < Initialize some stuff. */
    *error_encountered__ = FALSE_;
    *volume = 0.;
    volumet = 0.;
    volume_cp__ = 0.;
    volume_sphsec__ = 0.;
    *area = 0.;
    areat = 0.;
    area_l__ = 0.;
    area_sp__ = 0.;
    area_d__ = 0.;
    area_sphseg__ = 0.;
/*     < Find some sphere information. */
    sphereinfo_1.radius = rawdata_1.radii[sphereinfo_1.sphere];
    sphereinfo_1.rsq = sphereinfo_1.radius * sphereinfo_1.radius;
    sphereinfo_1.o_radius__ = 1. / sphereinfo_1.radius;
    sphereinfo_1.o_rsq__ = sphereinfo_1.o_radius__ * sphereinfo_1.o_radius__;
    sphereinfo_1.o_rcu__ = sphereinfo_1.o_rsq__ * sphereinfo_1.o_radius__;
/*     < For all plane pairs find the plane-plane intersections and */
/*     < hence the plane-plane lines that will become the edges of the */
/*     < arc-polygons */
    find_apices__();
/*     < Need to check to see if this sphere is a WIPEOUT */
    if (sphereinfo_1.numplanes > 0) {
/*       < Find the arc-polygons. */
	arc_polys__(&volume_cp__, &area_d__, &area_sphseg__, 
		error_encountered__);
	if (*error_encountered__) {
	    return 0;
	}
/*       < Find the festoons. */
	festoonery_(&area_sp__, error_encountered__);
	if (*error_encountered__) {
	    return 0;
	}
/*       < Find the lateral surface area of the spherical sector */
/*       < reduced by the sphere area 4 Pi R^3. */
	area_l__ = area_sp__ + area_sphseg__ - area_d__;
/*       < The lateral surface area may be off by an integer multiple */
/*       < of a sphere area.  Therefore, we need to adjust Area_L by */
/*       < an integer (positive or negative) such that Area_L is */
/*       < between zero and unity. */
	volume_sphsec__ = area_l__ - d_int(&area_l__) + (.5 - d_sign(&c_b31, &
		area_l__));
/*       < The volume of the spherical sector is equal to 1/3 of the */
/*       < product of the lateral area and the sphere radius. */
/*       < Therefore, the reduced spherical sector volume is exactly */
/*       < the reduced lateral surface area.  Find the total volume of */
/*       < this sphere reduced by 4/3 Pi R^3 and total exposed area of */
/*       < this sphere reduced by 4 Pi R^2. */
	*volume = volume_cp__ + volume_sphsec__;
	*area = volume_sphsec__;
/*       < Find the scaling factors for this sphere. */
	realr = rawdata_1.refradius * rawdata_1.radii[sphereinfo_1.sphere];
	realrsq = realr * realr;
	realrcu = realrsq * realr;
	ascale = realrsq * 12.566370614359172;
	vscale = realrcu * 4.1887902047863905;
/*       < Find the true volume and exposed area of this sphere in */
/*       < cubic angstroms. */
	volumet = *volume * vscale;
	areat = *area * ascale;
    }
    return 0;
} /* each_sphere__ */

/* ---------------------------------------------------------------------C */
/*     Find_Apices                                                     C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int find_apices__(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int plane_pair__(void);
    static doublereal sgn;
    static integer apexa, apexb, planei, indexi;

/*     < LOCAL */
/*     < Initialize some quantities. */
    apexb = 0;
    i__1 = sphereinfo_1.numplanes;
    for (i__ = 1; i__ <= i__1; ++i__) {
	apexinfo_1.n_edges__[i__ - 1] = 0;
    }
/*     < Store the information about this sphere's planes. */
    i__1 = sphereinfo_1.numplanes;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       < Find plane I, (P) */
	planei = planeinfo_1.plane_list__[sphereinfo_1.sphere + (i__ + 301) * 
		301 - 90601];
/*       < Index of plane (P) (i.e., 1st or 2nd sphere?) */
	indexi = planeinfo_1.plane_list__[sphereinfo_1.sphere + (i__ + 602) * 
		301 - 90601];
/*       < Has this plane been removed for this sphere? */
	if (planeinfo_1.planef[planei - 1] == -1) {
/*         < Yes, it has.  Turn off all its edges. */
	    apexinfo_1.n_edges__[i__ - 1] = -1;
/*         << Store order of plane */
	    spheresplanes_1.indexs[i__ - 1] = indexi;
/*         << True plane number */
	    spheresplanes_1.planes[i__ - 1] = planei;
	} else {
/*         << Store order of plane */
	    spheresplanes_1.indexs[i__ - 1] = indexi;
/*         << True plane number */
	    spheresplanes_1.planes[i__ - 1] = planei;
/*         << Sign of plane (P) */
	    spheresplanes_1.points[i__ - 1] = planeinfo_1.point[planei + 
		    indexi * 11250 - 11251];
/*         << Normal vector for plane (P) pointing out of cutout. */
/*         << Note that if Indexi is unity then the sphere we are */
/*         << currently on was the leading sphere in forming the */
/*         << plane and the unit normal vector stored is the proper */
/*         << sign.  If, however, Indexi is two then the current */
/*         << sphere was not the leading sphere in the creation of */
/*         << this plane (it was the second sphere in the double */
/*         << loop forming planes, see Find_Planes).  Therefore, the */
/*         << unit normal vector is anti-parallel to the direction */
/*         << it should be point and so we need to multiply its */
/*         << coordinates by minus one.  The variable Sgn is +1 if */
/*         << Indexi is unity and -1 if Indexi is two. */
	    sgn = (doublereal) (3 - (indexi << 1));
	    spheresplanes_1.nxs[i__ - 1] = planeinfo_1.nx[planei - 1] * sgn;
	    spheresplanes_1.nys[i__ - 1] = planeinfo_1.ny[planei - 1] * sgn;
	    spheresplanes_1.nzs[i__ - 1] = planeinfo_1.nz[planei - 1] * sgn;
/*         << Magnitude of vector connecting r0 to rk, |vk| */
	    spheresplanes_1.v_mags__[i__ - 1] = planeinfo_1.magnitude[planei 
		    + indexi * 11250 - 11251];
/*         << Squared radius of circle enscribed by plane I on sphere */
	    spheresplanes_1.rho_sqs__[i__ - 1] = planeinfo_1.rho_sq__[planei 
		    - 1];
	}
    }
/*     ### Loop Over All Possible Pairs Of Planes. */
/*     ------------------------------------------- */
/*     Note- center of sphere is denoted by r0 */
/*           n_Edges(plane): */
/*                -1 ==> plane removed; */
/*                 0 ==> no apices or edges for plane; */
/*                >0 ==> number of pairs of apices or edges on plane */
/*                       including some apex pairs or edges that */
/*                       have been turned off. */
/*     < Initialize variables */
    conclusion_1.remove = FALSE_;
    conclusion_1.interior = FALSE_;
    conclusion_1.wipeout = FALSE_;
    i__1 = sphereinfo_1.numplanes - 1;
    for (spheresplanes_1.ith = 1; spheresplanes_1.ith <= i__1; 
	    ++spheresplanes_1.ith) {
/*       < Has this plane been turned off? */
	if (apexinfo_1.n_edges__[spheresplanes_1.ith - 1] != -1) {
/*         < Initialize counter */
	    spheresplanes_1.jth = spheresplanes_1.ith;
	    while(spheresplanes_1.jth < sphereinfo_1.numplanes) {
/*           < Increment counter */
		++spheresplanes_1.jth;
/*           < Has this plane been turned off? */
		if (apexinfo_1.n_edges__[spheresplanes_1.jth - 1] != -1) {
/*             < Determine information about this pair of the planes */
/*             << Data passed through common block /Spheres Planes/ */
/*             << Results passed back through common block /Conclusion/ */
		    plane_pair__();
/*             < Do we have a wipeout of the sphere? */
/*             < How did the planes intersect (if they did)? */
/*             < Do we have a plane removal? */
		    if (conclusion_1.wipeout) {
/*               < We have a wipeout of this sphere. Turn off all of */
/*               < its planes and mark them as having died in trying to */
/*               < form the plane-plane edges. */
			i__2 = sphereinfo_1.numplanes;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    planeinfo_1.planef[spheresplanes_1.planes[i__ - 1]
				     - 1] = -1;
			    global_1.numarcsg[spheresplanes_1.planes[i__ - 1] 
				    - 1] = -2;
			    global_1.numarcsg[spheresplanes_1.planes[i__ - 1] 
				    + 11249] = -2;
			    global_1.trueedgesg[spheresplanes_1.planes[i__ - 
				    1] - 1] = -2;
			    global_1.trueedgesg[spheresplanes_1.planes[i__ - 
				    1] + 11249] = -2;
			}
			planeinfo_1.n_planes__[sphereinfo_1.sphere] = -1;
			sphereinfo_1.numplanes = -1;
			return 0;
		    } else if (conclusion_1.interior) {
			conclusion_1.interior = FALSE_;
/*               < Planes intersect inside sphere */
/*               < Increment the apex pair or edge count for each plane */
			++apexinfo_1.n_edges__[spheresplanes_1.ith - 1];
			++apexinfo_1.n_edges__[spheresplanes_1.jth - 1];
/*               < Increment the number of the apices */
/*               < Apex= ApexB= ApexA + 1 */
			apexa = apexb + 1;
			apexb += 2;
/*               < Store coordinates of first apex, apex A */
			apexinfo_1.apexx[apexa - 1] = conclusion_1.xa[0];
			apexinfo_1.apexy[apexa - 1] = conclusion_1.xa[1];
			apexinfo_1.apexz[apexa - 1] = conclusion_1.xa[2];
/*               < Store coordinates of second apex, apex B */
			apexinfo_1.apexx[apexb - 1] = conclusion_1.xb[0];
			apexinfo_1.apexy[apexb - 1] = conclusion_1.xb[1];
			apexinfo_1.apexz[apexb - 1] = conclusion_1.xb[2];
/*               < For each apex, flag that it exists as part of */
/*               < an apex pair (i.e., edge) */
			apexinfo_1.apexf[apexa - 1] = 2;
			apexinfo_1.apexf[apexa + 5999] = -1;
			apexinfo_1.apexf[apexb - 1] = 2;
			apexinfo_1.apexf[apexb + 5999] = -1;
/*               < Initialize the Twin list. */
			twinlist_1.twin[apexa - 1] = 0;
			twinlist_1.twin[apexb - 1] = 0;
/*               < Store (in order) apices associated with plane i */
/*               < (Pi):(A,B) stored anti-clockwise looking from */
/*               < outside. This is the edge formed by the intersection */
/*               < of plane i with plane j within the sphere. */
			apexinfo_1.edge_list__[spheresplanes_1.ith + (
				apexinfo_1.n_edges__[spheresplanes_1.ith - 1] 
				+ 6000) * 300 - 1800301] = apexa;
			apexinfo_1.edge_list__[spheresplanes_1.ith + (
				apexinfo_1.n_edges__[spheresplanes_1.ith - 1] 
				+ 12000) * 300 - 1800301] = apexb;
			apexinfo_1.edge_list__[spheresplanes_1.ith + (
				apexinfo_1.n_edges__[spheresplanes_1.ith - 1] 
				+ 18000) * 300 - 1800301] = 
				spheresplanes_1.jth;
/*               < Store (in order) apices associated with plane j */
/*               < (Pj):(B,A) stored anti-clockwise looking from */
/*               < outside. This is the edge formed by the intersection */
/*               < of plane i with plane j within the sphere. */
			apexinfo_1.edge_list__[spheresplanes_1.jth + (
				apexinfo_1.n_edges__[spheresplanes_1.jth - 1] 
				+ 6000) * 300 - 1800301] = apexb;
			apexinfo_1.edge_list__[spheresplanes_1.jth + (
				apexinfo_1.n_edges__[spheresplanes_1.jth - 1] 
				+ 12000) * 300 - 1800301] = apexa;
			apexinfo_1.edge_list__[spheresplanes_1.jth + (
				apexinfo_1.n_edges__[spheresplanes_1.jth - 1] 
				+ 18000) * 300 - 1800301] = 
				spheresplanes_1.ith;
/*               < Store (in order) the planes around apex A */
/*               < (A):(Pj,Pi) stored anti-clockwise from outside */
			debug_1.plane_order__[apexa - 1] = 
				spheresplanes_1.jth;
			debug_1.plane_order__[apexa + 5999] = 
				spheresplanes_1.ith;
			debug_1.plane_order__[apexa + 11999] = 0;
/*               < Store (in order) the planes around apex B */
/*               < (B):(Pi,Pj) stored anti-clockwise from outside */
			debug_1.plane_order__[apexb - 1] = 
				spheresplanes_1.ith;
			debug_1.plane_order__[apexb + 5999] = 
				spheresplanes_1.jth;
			debug_1.plane_order__[apexb + 11999] = 0;
		    } else if (conclusion_1.remove) {
			conclusion_1.remove = FALSE_;
/*               < Planes intersect outside sphere such */
/*               < that we have a removal of a plane */
/*               < REMOVE all vestiges of dead plane */
/*               < and leave inner loop */
/*               < Turn-off all of dead plane's apices (i.e., edges) */
			i__2 = apexinfo_1.n_edges__[conclusion_1.deadplane - 
				1];
			for (i__ = 1; i__ <= i__2; ++i__) {
			    apexinfo_1.apexf[apexinfo_1.edge_list__[
				    conclusion_1.deadplane + (i__ + 6000) * 
				    300 - 1800301] - 1] = -2;
			    apexinfo_1.apexf[apexinfo_1.edge_list__[
				    conclusion_1.deadplane + (i__ + 12000) * 
				    300 - 1800301] - 1] = -2;
			}
/*               < Turn-off dead plane itself */
			apexinfo_1.n_edges__[conclusion_1.deadplane - 1] = -1;
			global_1.numarcsg[spheresplanes_1.planes[
				conclusion_1.deadplane - 1] + 
				spheresplanes_1.indexs[conclusion_1.deadplane 
				- 1] * 11250 - 11251] = -2;
			global_1.trueedgesg[spheresplanes_1.planes[
				conclusion_1.deadplane - 1] + 
				spheresplanes_1.indexs[conclusion_1.deadplane 
				- 1] * 11250 - 11251] = -2;
/*               < Get out of inner loop if necessary.  That is, if the */
/*               < plane that has been killed is the Jth plane then do */
/*               < nothing to the inner loop counter.  However, if the */
/*               < dead plane is Ith then alter the inner loop counter */
/*               < to be greater than NumPlanes */
			spheresplanes_1.jth += sphereinfo_1.numplanes * (
				spheresplanes_1.jth - conclusion_1.deadplane);
		    }
		}
	    }
	}
    }
/*     < Store the number of apices. */
    apexinfo_1.apexcnt = apexb;
    return 0;
} /* find_apices__ */

/* ---------------------------------------------------------------------C */
/*     Plane_Pair                                                      C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int plane_pair__(void)
{
    extern /* Subroutine */ int parallel_(void), intersect_(doublereal *);
    static doublereal magsq, nixnj[3];

/*     < LOCAL */
/*     < Find cross-product of n and n', ni x nj. */
    nixnj[0] = spheresplanes_1.nys[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nzs[spheresplanes_1.jth - 1] - 
	    spheresplanes_1.nzs[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nys[spheresplanes_1.jth - 1];
    nixnj[1] = spheresplanes_1.nzs[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nxs[spheresplanes_1.jth - 1] - 
	    spheresplanes_1.nxs[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nzs[spheresplanes_1.jth - 1];
    nixnj[2] = spheresplanes_1.nxs[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nys[spheresplanes_1.jth - 1] - 
	    spheresplanes_1.nys[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nxs[spheresplanes_1.jth - 1];
/*     < Find the squared-magnitude of cross-product. */
    magsq = nixnj[0] * nixnj[0] + nixnj[1] * nixnj[1] + nixnj[2] * nixnj[2];
/*     < Are the two planes parallel?  Note- If the cross product of */
/*     < the normal vectors is the zero vector then the planes are */
/*     < parallel. */
    if (magsq < 1e-8) {
/*       < Planes Are Parallel */
	parallel_();
    } else {
/*       < Planes Intersect */
	intersect_(nixnj);
    }
    return 0;
} /* plane_pair__ */

/* ---------------------------------------------------------------------C */
/*     Parallel Planes                                                 C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int parallel_(void)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);

    /* Local variables */
    extern /* Subroutine */ int diff_side__(integer *), same_side__(integer *)
	    ;
    static integer pipj, ni_dollar_nj__;
    extern /* Subroutine */ int center_(integer *);

/*     < LOCAL */
/*     < Need to find dot-product of ni and nj, ni.nj := ni_dollar_nj. */
/*     < Note that we represent the dot product as an integer where */
/*     <   ni.nj= +1 means the planes normal vectors are parallel */
/*     <   ni.nj= -1 means the planes normal vectors are anti-parallel. */
    d__1 = spheresplanes_1.nxs[spheresplanes_1.ith - 1] * spheresplanes_1.nxs[
	    spheresplanes_1.jth - 1] + spheresplanes_1.nys[
	    spheresplanes_1.ith - 1] * spheresplanes_1.nys[
	    spheresplanes_1.jth - 1] + spheresplanes_1.nzs[
	    spheresplanes_1.ith - 1] * spheresplanes_1.nzs[
	    spheresplanes_1.jth - 1];
    ni_dollar_nj__ = i_dnnt(&d__1);
/*     ### Find The Product Of The Points. */
/*     ----------------------------------- */
    pipj = spheresplanes_1.points[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.points[spheresplanes_1.jth - 1];
/*     ### Do Either Of The Planes Go Through The Sphere Center? */
/*     --------------------------------------------------------- */
    if (pipj == 0) {
/*       < One or both planes go through the center */
	center_(&ni_dollar_nj__);
    } else {
/*       < Neither of the planes go through the sphere center */
	if (ni_dollar_nj__ * pipj == 1) {
/*         < Planes are on the same side of sphere center */
	    same_side__(&ni_dollar_nj__);
	} else {
/*         < Planes are on different sides of sphere center */
	    diff_side__(&ni_dollar_nj__);
	}
    }
    return 0;
} /* parallel_ */

/* ---------------------------------------------------------------------C */
/*     Parallel Planes Through Center                                  C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int center_(integer *ni_dollar_nj__)
{
    /* System generated locals */
    integer i__1, i__2;

/*     < PASSED */
    if (*ni_dollar_nj__ == -1) {
/*       < (ni_dollar_nj = -1) ==> */
/*       < Possible WIPEOUT */
/* Computing MAX */
	i__1 = spheresplanes_1.points[spheresplanes_1.ith - 1], i__2 = 
		spheresplanes_1.points[spheresplanes_1.jth - 1];
	if (max(i__1,i__2) == 0) {
/*         < (Pointi,Pointj)= (-1,0), (0,-1) or (0,0) */
/*         < WIPEOUT of sphere */
	    conclusion_1.wipeout = TRUE_;
	}
    } else {
/*       < (ni_dollar_nj = +1) ==> */
/*       < REMOVAL of a plane */
	if (spheresplanes_1.points[spheresplanes_1.ith - 1] == 0) {
/*         < REMOVE one plane */
	    if (spheresplanes_1.points[spheresplanes_1.jth - 1] == 1) {
/*           < (Pointi,Pointj)= (0,+1) */
/*           < REMOVE plane not going through center, planej */
		conclusion_1.remove = TRUE_;
		conclusion_1.deadplane = spheresplanes_1.jth;
	    } else {
/*           < (Pointi,Pointj)= (0,-1) */
/*           < REMOVE plane through center, planei */
/*           <   or */
/*           < (Pointi,Pointj)= (0,0) */
/*           < REMOVE either plane (e.g., planei) */
		conclusion_1.remove = TRUE_;
		conclusion_1.deadplane = spheresplanes_1.ith;
	    }
	} else {
/*         < Pointj= 0 */
/*         < REMOVE one plane */
	    if (spheresplanes_1.points[spheresplanes_1.ith - 1] == 1) {
/*           < (Pointi,Pointj)= (+1,0) */
/*           < REMOVE plane through center, planei */
		conclusion_1.remove = TRUE_;
		conclusion_1.deadplane = spheresplanes_1.ith;
	    } else {
/*           < (Pointi,Pointj)= (-1,0) */
/*           < REMOVE plane not going through center, planej */
		conclusion_1.remove = TRUE_;
		conclusion_1.deadplane = spheresplanes_1.jth;
	    }
	}
    }
    return 0;
} /* center_ */

/* ---------------------------------------------------------------------C */
/*     Parallel Planes On The Same Side Of Center                      C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int same_side__(integer *ni_dollar_nj__)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

/*     < PASSED */
/*     < We know (ni.nj)*(Pi.Pj)= +1 */
    if (*ni_dollar_nj__ == -1) {
/*       < (ni_dollar_nj = -1) ==> Points are of opposite sign */
/*       < Possible WIPEOUT */
/*       < Does the plane with larger |v| have Point= -1? */
/* Computing MAX */
	d__1 = spheresplanes_1.v_mags__[spheresplanes_1.ith - 1], d__2 = 
		spheresplanes_1.v_mags__[spheresplanes_1.jth - 1];
/* Computing MIN */
	d__3 = spheresplanes_1.points[spheresplanes_1.ith - 1] * 
		spheresplanes_1.v_mags__[spheresplanes_1.ith - 1], d__4 = 
		spheresplanes_1.points[spheresplanes_1.jth - 1] * 
		spheresplanes_1.v_mags__[spheresplanes_1.jth - 1];
	if (max(d__1,d__2) == -min(d__3,d__4)) {
/*         < Point= -1 for plane with Max(|v|) or |vi|=|vj| */
/*         < WIPEOUT of sphere */
	    conclusion_1.wipeout = TRUE_;
	}
    } else {
/*       < (ni_dollar_nj = +1) ==> Points are of same sign */
/*       < REMOVE one plane */
	if ((spheresplanes_1.v_mags__[spheresplanes_1.ith - 1] - 
		spheresplanes_1.v_mags__[spheresplanes_1.jth - 1]) * 
		spheresplanes_1.points[spheresplanes_1.ith - 1] > 0.) {
/*         < REMOVE plane i */
	    conclusion_1.remove = TRUE_;
	    conclusion_1.deadplane = spheresplanes_1.ith;
	} else {
/*         < REMOVE plane j */
	    conclusion_1.remove = TRUE_;
	    conclusion_1.deadplane = spheresplanes_1.jth;
	}
    }
    return 0;
} /* same_side__ */

/* ---------------------------------------------------------------------C */
/*     Parallel Planes Different Side Of Center                        C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int diff_side__(integer *ni_dollar_nj__)
{
/*     < PASSED */
/*     < We know (ni.nj)*(Pi.Pj)= -1 */
    if (*ni_dollar_nj__ == -1) {
/*       < (ni_dollar_nj= -1) ==> Points are of same sign */
/*       < Possible WIPEOUT */
	if (spheresplanes_1.points[spheresplanes_1.ith - 1] == -1) {
/*         < WIPEOUT of sphere */
	    conclusion_1.wipeout = TRUE_;
	}
    } else {
/*       < (ni_dollar_nj= +1) ==> Points are of opposite sign */
/*       < REMOVE one plane */
	if (spheresplanes_1.points[spheresplanes_1.ith - 1] == 1) {
/*         < REMOVE plane i */
	    conclusion_1.remove = TRUE_;
	    conclusion_1.deadplane = spheresplanes_1.ith;
	} else {
/*         < REMOVE plane j */
	    conclusion_1.remove = TRUE_;
	    conclusion_1.deadplane = spheresplanes_1.jth;
	}
    }
    return 0;
} /* diff_side__ */

/* ---------------------------------------------------------------------C */
/*     Intersecting Planes                                             C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int intersect_(doublereal *nixnj)
{
    static doublereal b[2], f[3], beta, ni_dollar_nj__, alpha;
    extern /* Subroutine */ int inside_(doublereal *, doublereal *, 
	    doublereal *);
    static doublereal mag_fsq__;
    extern /* Subroutine */ int outside_(doublereal *, doublereal *, 
	    doublereal *);

/*     < PASSED */
/*     < LOCAL */
/*     ### Find Point On Line Of Intersection Of Two Planes. */
/*     ----------------------------------------------------- */

/*         F is a vector from the sphere center (r0) to point rH on line */
/*     l(Pi,Pj) formed by intersection of planes ("F" stands for */
/*     "Foot").  This vector is perpendicular the line of intersection */
/*     and hence the shortest possible vector connecting r0 to l(Pi,Pj). */
/*     Therefore we have */

/*               ni.(F - vi)= 0 */
/*               nj.(F - vj)= 0 */
/*                 (nixnj).F= 0 */

/*     with ni.vi= Pointi*|vi|, nj.vj= Pointj*|vj| we have */

/*               ni.F= ni.vi = Pointi*|vi| */
/*               nj.F= nj.vj = Pointj*|vj| */
/*          (nixnj).F= 0 */

/*     we can solve this system for vector F using Cramer's Rule */

/*               |  ni   |       | Pointi*|vi| |    | bi | */
/*               |  nj   | . F = | Pointj*|vj| | := | bj | */
/*               | nixnj |       |      0      |    | 0  | */

/*     We propose */

/*               F= Alpha*ni + Beta*nj + Gamma*(nixnj) */

/*     therefore */

/*            ni.F= Alpha + Beta*(ni.nj) + Gamma*0 := bi */
/*            nj.F= Alpha*(ni.nj) + Beta + Gamma*0 := bj */
/*       (nixnj).F= Alpha*0 + Beta*0 + Gamma*(nixnj).(nixnj) := 0 */

/*     with a solution of */

/*           Alpha= [bi - Beta*(ni.nj)] */
/*            Beta= [bj - bi*(ni.nj)]/[1 - (ni.nj)^2] */
/*           Gamma= 0 */
    /* Parameter adjustments */
    --nixnj;

    /* Function Body */
    b[0] = spheresplanes_1.points[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.v_mags__[spheresplanes_1.ith - 1];
    b[1] = spheresplanes_1.points[spheresplanes_1.jth - 1] * 
	    spheresplanes_1.v_mags__[spheresplanes_1.jth - 1];
/*     ### Need To Find Dot-Product of ni and nj, ni.nj := ni_dollar_nj */
/*     --------------------------------------------------------- */
    ni_dollar_nj__ = spheresplanes_1.nxs[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nxs[spheresplanes_1.jth - 1] + 
	    spheresplanes_1.nys[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nys[spheresplanes_1.jth - 1] + 
	    spheresplanes_1.nzs[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.nzs[spheresplanes_1.jth - 1];
/*     ### Find The Constants. */
/*     ----------------------- */
    beta = (b[1] - b[0] * ni_dollar_nj__) / (1. - ni_dollar_nj__ * 
	    ni_dollar_nj__);
    alpha = b[0] - beta * ni_dollar_nj__;
/*     ### Find The Vector F. */
/*     ---------------------- */
    f[0] = alpha * spheresplanes_1.nxs[spheresplanes_1.ith - 1] + beta * 
	    spheresplanes_1.nxs[spheresplanes_1.jth - 1];
    f[1] = alpha * spheresplanes_1.nys[spheresplanes_1.ith - 1] + beta * 
	    spheresplanes_1.nys[spheresplanes_1.jth - 1];
    f[2] = alpha * spheresplanes_1.nzs[spheresplanes_1.ith - 1] + beta * 
	    spheresplanes_1.nzs[spheresplanes_1.jth - 1];
/*     ### Find The Square-Magnitude Of F, |F|^2 */
/*     ----------------------------------------- */
    mag_fsq__ = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
/*     ### Do the planes intersect outside the sphere? */
/*     ### This will be true if |F|^2 is greater then */
/*     ### the square of the radius of the sphere. */
/*     ----------------------------------------------- */
    if (mag_fsq__ > sphereinfo_1.rsq) {
/*       < Planes Intersect Outside Sphere */
	outside_(b, &ni_dollar_nj__, &mag_fsq__);
    } else {
/*       < Planes Intersect Inside Sphere */
	inside_(&nixnj[1], b, &ni_dollar_nj__);
    }
    return 0;
} /* intersect_ */

/* ---------------------------------------------------------------------C */
/*     Planes Intersecting Outside Sphere                              C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int outside_(doublereal *b, doublereal *ni_dollar_nj__, 
	doublereal *mag_fsq__)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ int centerout_(doublereal *);
    static integer pipj;
    static doublereal scalar;

/*     < PASSED */
/*     < LOCAL */
/*     ### Find The Product Of The Points. */
/*     ----------------------------------- */
    /* Parameter adjustments */
    --b;

    /* Function Body */
    pipj = spheresplanes_1.points[spheresplanes_1.ith - 1] * 
	    spheresplanes_1.points[spheresplanes_1.jth - 1];
/*     ### Do either of the planes go through the sphere center? */
/*     --------------------------------------------------------- */
    if (pipj == 0) {
/*       < One of the planes goes through the center */
	centerout_(ni_dollar_nj__);
    } else {
/*       < Neither plane goes through center */
/*       < Form  [vj x F].[vi x F]/(|vi|*|vj|) */
/*       Note: */
/*       There are three vectors eminenting from the sphere center: */
/*       F= (rH - r0), vi= (rki - r0) and vj= (rkj - r0). */
/*       The vectors (vi x F) and (vj x F) will be parallel or */
/*       anti parallel (i.e., their dot product will be positive */
/*       or negative) depending on the order of the vectors. If F is */
/*       found between vi and vj then Scalar (below) will be negative */
/*       and this means r0 is between the planes (Pi) and (Pj). */

/*       Topology Note - */

/*               rki  [vi]                       rki  [vi] */
/*              /                               / */
/*             /                               / */
/*            /                               / */
/*           /                               / */
/*       r0 <---------->rH  [F]          r0 <---------->rkj [vj] */
/*           \                               \ */
/*            \           Scalar < 0          \           Scalar > 0 */
/*             \                               \ */
/*              \                               \ */
/*               rkj  [vj]                       rH  [F] */
	scalar = *mag_fsq__ * *ni_dollar_nj__ * pipj - 
		spheresplanes_1.v_mags__[spheresplanes_1.ith - 1] * 
		spheresplanes_1.v_mags__[spheresplanes_1.jth - 1];
/*       < Is Sphere Center Between The Planes? */
	if (scalar < 0.) {
/*         < Sphere center is between planes */
	    if (pipj == 1) {
/*           < Possible WIPEOUT */
		if (spheresplanes_1.points[spheresplanes_1.ith - 1] == -1) {
/*             < WIPEOUT of sphere */
		    conclusion_1.wipeout = TRUE_;
		}
	    } else {
/*           < REMOVE one plane */
		if (spheresplanes_1.points[spheresplanes_1.jth - 1] == 1) {
/*             < REMOVE plane j */
		    conclusion_1.remove = TRUE_;
		    conclusion_1.deadplane = spheresplanes_1.jth;
		} else {
/*             < REMOVE plane i */
		    conclusion_1.remove = TRUE_;
		    conclusion_1.deadplane = spheresplanes_1.ith;
		}
	    }
	} else {
/*         < Sphere center is not between planes */
	    if (pipj == 1) {
/*           < REMOVE one plane */
/*           < The plane  with the larger |v| if its Point is +1 */
/*           < or the plane with the smaller |v| if its Point is -1. */
		if ((spheresplanes_1.v_mags__[spheresplanes_1.ith - 1] - 
			spheresplanes_1.v_mags__[spheresplanes_1.jth - 1]) * 
			spheresplanes_1.points[spheresplanes_1.ith - 1] > 0.) 
			{
/*             < REMOVE Plane i */
		    conclusion_1.remove = TRUE_;
		    conclusion_1.deadplane = spheresplanes_1.ith;
		} else {
/*             < REMOVE Plane j */
		    conclusion_1.remove = TRUE_;
		    conclusion_1.deadplane = spheresplanes_1.jth;
		}
	    } else {
/*           < Possible WIPEOUT */
/*           < Does the plane with larger |v| have Point= -1? */
/* Computing MAX */
		d__1 = spheresplanes_1.v_mags__[spheresplanes_1.ith - 1], 
			d__2 = spheresplanes_1.v_mags__[spheresplanes_1.jth - 
			1];
		if (max(d__1,d__2) == -min(b[1],b[2])) {
/*             < Point= -1 for plane with Max(|v|) or |vi|=|vj| */
/*             < WIPEOUT of sphere */
		    conclusion_1.wipeout = TRUE_;
		}
	    }
	}
    }
    return 0;
} /* outside_ */

/* ---------------------------------------------------------------------C */
/*     Planes Intersecting Outside Sphere One Going Through Center     C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int centerout_(doublereal *ni_dollar_nj__)
{
    /* System generated locals */
    integer i__1, i__2;

/*     < PASSED */
/*     ### One Of The Planes Goes Through The Center Of The Sphere. */
/*     ------------------------------------------------------------ */
    if (*ni_dollar_nj__ > 0.) {
	if ((spheresplanes_1.v_mags__[spheresplanes_1.ith - 1] - 
		spheresplanes_1.v_mags__[spheresplanes_1.jth - 1]) * 
		spheresplanes_1.points[spheresplanes_1.ith - 1] > 0.) {
/*         < REMOVE plane i */
	    conclusion_1.remove = TRUE_;
	    conclusion_1.deadplane = spheresplanes_1.ith;
	} else {
/*         < REMOVE plane j */
	    conclusion_1.remove = TRUE_;
	    conclusion_1.deadplane = spheresplanes_1.jth;
	}
    } else /* if(complicated condition) */ {
/* Computing MIN */
	i__1 = spheresplanes_1.points[spheresplanes_1.ith - 1], i__2 = 
		spheresplanes_1.points[spheresplanes_1.jth - 1];
	if (min(i__1,i__2) == -1) {
/*       < WIPEOUT of sphere */
	    conclusion_1.wipeout = TRUE_;
	}
    }
    return 0;
} /* centerout_ */

/* ---------------------------------------------------------------------C */
/*     Planes Intersecting Inside Sphere                               C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int inside_(doublereal *nixnj, doublereal *b, doublereal *
	ni_dollar_nj__)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal beta, alpha, group, gammaa;

/*     < PASSED */
/*     < LOCAL */
/*     ### Find Where The Line Of Intersection Cuts The Sphere Surface. */
/*     ---------------------------------------------------------------- */
/*     Note - */

/*     Solve the system */

/*           ni.(x - vi)= 0 */
/*           ni.(x - vj)= 0 */

/*     for x= (Point on Line - r0) such that |x|= Sphere Radius or */

/*           ni.x= Pointi*|vi| */
/*           nj.x= Pointj*|vj| */
/*            x.x= R^2 */

/*     System is quadratic; there are two solutions, xA and xB */

/*             xA= Alpha*ni + Beta*nj + GammaA*(ni x nj) */
/*             xB= Alpha*ni + Beta*nj + GammaB*(ni x nj) */
    /* Parameter adjustments */
    --b;
    --nixnj;

    /* Function Body */
    conclusion_1.interior = TRUE_;
    group = 1. / (1. - *ni_dollar_nj__ * *ni_dollar_nj__);
    beta = (b[2] - b[1] * *ni_dollar_nj__) * group;
    alpha = b[1] - beta * *ni_dollar_nj__;
    gammaa = -sqrt(group * (sphereinfo_1.rsq - (alpha * b[1] + beta * b[2])));
    group = gammaa * nixnj[1];
    conclusion_1.xa[0] = alpha * spheresplanes_1.nxs[spheresplanes_1.ith - 1] 
	    + beta * spheresplanes_1.nxs[spheresplanes_1.jth - 1];
    conclusion_1.xb[0] = conclusion_1.xa[0] - group;
    conclusion_1.xa[0] += group;
    group = gammaa * nixnj[2];
    conclusion_1.xa[1] = alpha * spheresplanes_1.nys[spheresplanes_1.ith - 1] 
	    + beta * spheresplanes_1.nys[spheresplanes_1.jth - 1];
    conclusion_1.xb[1] = conclusion_1.xa[1] - group;
    conclusion_1.xa[1] += group;
    group = gammaa * nixnj[3];
    conclusion_1.xa[2] = alpha * spheresplanes_1.nzs[spheresplanes_1.ith - 1] 
	    + beta * spheresplanes_1.nzs[spheresplanes_1.jth - 1];
    conclusion_1.xb[2] = conclusion_1.xa[2] - group;
    conclusion_1.xa[2] += group;
/*     < Order the points A and B for each plane such that when looking */
/*     < at each plane from the outside of the cutout (i.e., the part of */
/*     < the sphere that remains) the points are read in an */
/*     < anti-clockwise direction.  Thus, the vector connecting A and B */
/*     < points in an anti-clockwise direction relative to each plane. */
/*     < This is accomplished by looking at (nixnj).(xA-xB).  If this is */
/*     < negative then A and B are stored as (A,B) for plane (Pi) and */
/*     < (B,A) for plane (Pj).  Otherwise if it is positive then A and B */
/*     < are stored as (B,A) for plane (Pi) and (A,B) for plane (Pj) */
/*     < */
/*     < The line of intersection is */
/*     <      (xA - xB)= (GammaA - GammaB)*(ni x nj) */
/*     < Therefore the dot product is */
/*     <       (ni x nj).(xA - xB)= (GammaA - GammaB)*|nixnj|^2 */
/*     < Sign[(ni x nj).(xA - xB)]= Sign(GammaA - GammaB) */
/*     < */
/*     < But we know GammaB= - GammaA therefore */
/*     < */
/*     < Sign[(ni x nj).(xA - xB)]= Sign(GammaA*2) */
/*     <                          = Sign(GammaA) */
/*     < */
/*     < But GammaA is defined as (-Sqrt) so (nixnj).(xA-xB) is by */
/*     < definition less than zero.  Therefore we store */
/*     < */
/*     <     (Pi):(A,B) and (Pj):(B,A) */
/*     <     (A):(Pj,Pi) and (B):(Pi,Pj) */
    return 0;
} /* inside_ */

/* ---------------------------------------------------------------------C */
/*     Arc-Polygons                                                    C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int arc_polys__(doublereal *volume_cp__, doublereal *
	area_d__, doublereal *area_sphseg__, logical *error_encountered__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int build_ap__(doublereal *, logical *);
    static integer n_sphseg__, i__, j, edgei;
    static doublereal eta_ap__, o_rhoi__;
    extern /* Subroutine */ int find_ap__(void);
    static doublereal pvi_mag__;

/*     < PASSED */
/*     < LOCAL */
    *error_encountered__ = FALSE_;
/*     ### Initialize Some Things For This Sphere. */
/*     ------------------------------------------- */
/*     < Initialize the triple plane intersection counter, */
/*     < the arc counter and the vertex counter. */
    sphereinfo_1.tripcnt = 0;
    arcinfo_1.arccnt = 0;
    arcinfo_1.vertexcnt = 0;
/*     < Initialize the Cone-Pyramid volume, the diangle area, and the */
/*     < area of the spherical segments. */
    *volume_cp__ = 0.;
    *area_d__ = 0.;
    *area_sphseg__ = 0.;
    n_sphseg__ = 0;
/*     ### Single Loop Over All Planes For This Sphere. */
/*     ------------------------------------------------ */
    i__1 = sphereinfo_1.numplanes;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       < Find the number of apex pairs (i.e., edges) for plane I */
/*       <    -1 ==> plane does not exist */
/*       <     0 ==> plane exists as a spherical cap only */
/*       <    >0 ==> plane contains one or more edges */
	thisplane_1.numedges = apexinfo_1.n_edges__[i__ - 1];
	if (thisplane_1.numedges == 0) {
/*         < This plane does not intersect with other planes.  We have */
/*         < a spherical cap or SPHERICAL SEGMENT for this plane. */
/*         < Eta_AP (arc-polygon area divided by rhoi_sq/2) is 2Pi and */
/*         < we add this to the reduced volume of the cone-pyramid. */
	    pvi_mag__ = spheresplanes_1.points[i__ - 1] * 
		    spheresplanes_1.v_mags__[i__ - 1];
	    *volume_cp__ += pvi_mag__ * 6.2831853071795862 * 
		    spheresplanes_1.rho_sqs__[i__ - 1];
/*         < Find the reduced area bordered by this spherical segment. */
/*         < This will be added to the area of the final spherical */
/*         < polygon later. */
	    ++n_sphseg__;
	    *area_sphseg__ += pvi_mag__;
	    thisplane_1.indexi = spheresplanes_1.indexs[i__ - 1];
	    thisplane_1.planei = spheresplanes_1.planes[i__ - 1];
	    global_1.numarcsg[thisplane_1.planei + thisplane_1.indexi * 11250 
		    - 11251] = -3;
	    global_1.trueedgesg[thisplane_1.planei + thisplane_1.indexi * 
		    11250 - 11251] = -3;
	} else if (thisplane_1.numedges > 0) {
/*         < We have an ARC-POLYGON */
/*         < Pull out this plane's information */
	    thisplane_1.indexi = spheresplanes_1.indexs[i__ - 1];
	    thisplane_1.planei = spheresplanes_1.planes[i__ - 1];
	    thisplane_1.pointi = spheresplanes_1.points[i__ - 1];
/*         < Pull out its unit normal vector. */
	    thisplane_1.ni[0] = spheresplanes_1.nxs[i__ - 1];
	    thisplane_1.ni[1] = spheresplanes_1.nys[i__ - 1];
	    thisplane_1.ni[2] = spheresplanes_1.nzs[i__ - 1];
/*         < Pull out its distance from the sphere center and inscribed */
/*         < radius. */
	    thisplane_1.vi_mag__ = spheresplanes_1.v_mags__[i__ - 1];
	    thisplane_1.rhoi_sq__ = spheresplanes_1.rho_sqs__[i__ - 1];
/*         < Find the sphere radius and |vi| reduced by the inscribed */
/*         < radius of the plane. */
	    thisplane_1.o_rhoi_sq__ = 1. / thisplane_1.rhoi_sq__;
	    o_rhoi__ = sqrt(thisplane_1.o_rhoi_sq__);
	    thisplane_1.r_rhoi__ = sphereinfo_1.radius * o_rhoi__;
	    thisplane_1.vi_rhoi__ = thisplane_1.vi_mag__ * o_rhoi__;
/*         < Precalculate one-half the square curvature of the inscribed */
/*         < circle. */
	    thisplane_1.curvei = thisplane_1.o_rhoi_sq__ * .5;
/*         < Store all the edges for this plane */
	    j = 0;
	    i__2 = thisplane_1.numedges;
	    for (edgei = 1; edgei <= i__2; ++edgei) {
		if (apexinfo_1.apexf[apexinfo_1.edge_list__[i__ + (edgei + 
			6000) * 300 - 1800301] - 1] != -2) {
		    ++j;
		    thisplane_1.edgesi[j - 1] = apexinfo_1.edge_list__[i__ + (
			    edgei + 6000) * 300 - 1800301];
		    thisplane_1.edgesi[j + 5999] = apexinfo_1.edge_list__[i__ 
			    + (edgei + 12000) * 300 - 1800301];
		}
	    }
	    thisplane_1.numedges = j;
/*         < Store the local plane number.  This is the number of this */
/*         < plane for this sphere. */
	    thisplane_1.planenum = i__;
/*         < Find this Arc-Polygon.  Intersect all the edges on this */
/*         < plane face clipping back or discarding the unneeded parts. */
	    find_ap__();
/*         < Has this plane been killed? */
	    if (thisplane_1.numedges <= 0) {
		global_1.numarcsg[thisplane_1.planei + thisplane_1.indexi * 
			11250 - 11251] = -1;
		global_1.trueedgesg[thisplane_1.planei + thisplane_1.indexi * 
			11250 - 11251] = -1;
	    } else {
/*           < Form the Arc-Polygon. Connect all the edges properly and */
/*           < find its area and the area of angle of its diangles. */
/*           < Initialize the diangle accumulators for this plane. */
		diangleinfo_1.etad[0] = 0.;
		diangleinfo_1.etad[1] = 0.;
		build_ap__(&eta_ap__, error_encountered__);
		if (*error_encountered__) {
		    return 0;
		}
/*           < Accumulate the reduced volume of the cone-pyramid.  The */
/*           < sign convention is to add the cone-pyramid volume for a */
/*           < plane if the sphere center is contained in that plane's */
/*           < cutout and subtract it if not.  Therefore, we */
/*           < premultiply by Pointi for this plane.  Note that if the */
/*           < plane goes through the sphere center then Pointi and */
/*           < vi_Mag is zero the plane makes no contribution to the */
/*           < cone-pyramid volume. */
		*volume_cp__ += thisplane_1.pointi * eta_ap__ * 
			thisplane_1.vi_mag__ * thisplane_1.rhoi_sq__;
/*           < Find the total diangle area for this plane and add it to */
/*           < the total for this sphere. The sign convention for the */
/*           < diangles is to add them positively if the sphere center */
/*           < is contained in the plane's cutout and subtract them if */
/*           < not.  Therefore, we premultiply by Pointi for this */
/*           < plane.  Later this area will be subtracted to find the */
/*           < lateral area. */
		*area_d__ += thisplane_1.pointi * (diangleinfo_1.etad[0] - 
			thisplane_1.vi_mag__ * sphereinfo_1.o_radius__ * 
			diangleinfo_1.etad[1]);
	    }
	}
    }
/*     < Find the cone-pyramid volume reduced by 4/3 Pi R^3. */
    *volume_cp__ = *volume_cp__ * sphereinfo_1.o_rcu__ * .039788735772973836;
/*     < Find the diangle area for this sphere reduced by 4 Pi R^2. */
    *area_d__ *= .15915494309189535;
/*     < Find the area of the spherical segments reduced by 4 Pi R^2. */
    *area_sphseg__ = ((doublereal) n_sphseg__ + sphereinfo_1.o_radius__ * *
	    area_sphseg__) * .5;
    return 0;
} /* arc_polys__ */

/* ---------------------------------------------------------------------C */
/*     Find The Arc-Polygon                                            C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int find_ap__(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal maxabxcd;
    extern /* Subroutine */ int three_plane__(void);
    static doublereal magsq;
    extern /* Subroutine */ int cross_(doublereal *, doublereal *, doublereal 
	    *), parallel_edges__(void), aparallel_edges__(void);

/*     < LOCAL */
/*     ### Go Over All Pairs Of Edges (Line Segments). */
/*     ----------------------------------------------- */
    i__1 = thisplane_1.numedges - 1;
    for (edgepair_1.edgei = 1; edgepair_1.edgei <= i__1; ++edgepair_1.edgei) {
/*       < Pull out apices A and B of EdgeI */
	edgepair_1.apexa = thisplane_1.edgesi[edgepair_1.edgei - 1];
	edgepair_1.apexb = thisplane_1.edgesi[edgepair_1.edgei + 5999];
	i__2 = thisplane_1.numedges;
	for (edgepair_1.edgej = edgepair_1.edgei + 1; edgepair_1.edgej <= 
		i__2; ++edgepair_1.edgej) {
/*         < Pull out and store apex A's coordinates.  Note that */
/*         < this needs to remain in the inner loop since the */
/*         < coordinates of ApexA or ApexB may be changed within */
/*         < the inner loop. */
	    edgepair_1.ra[0] = apexinfo_1.apexx[edgepair_1.apexa - 1];
	    edgepair_1.ra[1] = apexinfo_1.apexy[edgepair_1.apexa - 1];
	    edgepair_1.ra[2] = apexinfo_1.apexz[edgepair_1.apexa - 1];
/*         < Find vector connecting A to B (need to pull out rB) */
/*         < This is the line segment formed by another plane on */
/*         < the plane face of plane PlaneNum. */
	    edgepair_1.eab[0] = apexinfo_1.apexx[edgepair_1.apexb - 1] - 
		    edgepair_1.ra[0];
	    edgepair_1.eab[1] = apexinfo_1.apexy[edgepair_1.apexb - 1] - 
		    edgepair_1.ra[1];
	    edgepair_1.eab[2] = apexinfo_1.apexz[edgepair_1.apexb - 1] - 
		    edgepair_1.ra[2];
/*         < Find the magnitude of this vector. */
	    edgepair_1.rmagab = 1. / sqrt(edgepair_1.eab[0] * edgepair_1.eab[
		    0] + edgepair_1.eab[1] * edgepair_1.eab[1] + 
		    edgepair_1.eab[2] * edgepair_1.eab[2]);
/*         < Find the unit vector AB. */
	    edgepair_1.eab[0] *= edgepair_1.rmagab;
	    edgepair_1.eab[1] *= edgepair_1.rmagab;
	    edgepair_1.eab[2] *= edgepair_1.rmagab;
/*         < Pull out apices C and D of EdgeJ */
	    edgepair_1.apexc = thisplane_1.edgesi[edgepair_1.edgej - 1];
	    edgepair_1.apexd = thisplane_1.edgesi[edgepair_1.edgej + 5999];
/*         < Pull out and store apex C's coordinates */
	    edgepair_1.rc[0] = apexinfo_1.apexx[edgepair_1.apexc - 1];
	    edgepair_1.rc[1] = apexinfo_1.apexy[edgepair_1.apexc - 1];
	    edgepair_1.rc[2] = apexinfo_1.apexz[edgepair_1.apexc - 1];
/*         < Find vector connecting C to D (need to pull out rD) */
/*         < This is the line segment formed by another plane on */
/*         < the plane face of plane PlaneNum. */
	    edgepair_1.ecd[0] = apexinfo_1.apexx[edgepair_1.apexd - 1] - 
		    edgepair_1.rc[0];
	    edgepair_1.ecd[1] = apexinfo_1.apexy[edgepair_1.apexd - 1] - 
		    edgepair_1.rc[1];
	    edgepair_1.ecd[2] = apexinfo_1.apexz[edgepair_1.apexd - 1] - 
		    edgepair_1.rc[2];
/*         < Find the reciprocal magnitude of this vector. */
	    edgepair_1.rmagcd = 1. / sqrt(edgepair_1.ecd[0] * edgepair_1.ecd[
		    0] + edgepair_1.ecd[1] * edgepair_1.ecd[1] + 
		    edgepair_1.ecd[2] * edgepair_1.ecd[2]);
/*         < Find the unit vector CD. */
	    edgepair_1.ecd[0] *= edgepair_1.rmagcd;
	    edgepair_1.ecd[1] *= edgepair_1.rmagcd;
	    edgepair_1.ecd[2] *= edgepair_1.rmagcd;
/*         < Find the cross product of the two vectors. */
	    cross_(edgepair_1.eab, edgepair_1.ecd, edgepair_1.abxcd);
/*         < Find the squared-magnitude of the cross-product */
	    magsq = edgepair_1.abxcd[0] * edgepair_1.abxcd[0] + 
		    edgepair_1.abxcd[1] * edgepair_1.abxcd[1] + 
		    edgepair_1.abxcd[2] * edgepair_1.abxcd[2];
/*         < We must determine how these two edges lie relative to one */
/*         < another and relative to the cutout.  First, we must */
/*         < determine if the edges are parallel (or anti-parallel). */
/*         < We do this by examining the squared-magnitude of the the */
/*         < cross-product of the two edge vectors.  If this is less */
/*         < than some Epsilon then we find the dot-product of the two */
/*         < edge vectors to determine if they are parallel or */
/*         < anti-parallel. If the squared-magnitude is greater than */
/*         < some Epsilon then the vectors intersect.  We must then */
/*         < largest, and by definition, non-zero element of the */
/*         < cross-product for this will be the determinate of the */
/*         < linear system we will solve to find the point of */
/*         < intersection of the two edge vectors. */
	    if (magsq < 1e-9) {
		if (edgepair_1.eab[0] * edgepair_1.ecd[0] + edgepair_1.eab[1] 
			* edgepair_1.ecd[1] + edgepair_1.eab[2] * 
			edgepair_1.ecd[2] < 0.) {
/*             < The two edges are anti-parallel. We should check to */
/*             < see if they coincide indicating that the arc-polygon */
/*             < is a line and we should kill this plane face.  If they */
/*             < do not coincide we then need to see if they point in */
/*             < an anti-clockwise fashion.  If not then the plane */
/*             < needs to be killed. */
		    aparallel_edges__();
/*             < Has this plane been killed? */
		    if (thisplane_1.numedges == -1) {
			return 0;
		    }
		} else {
/*             < The two edges are parallel. The edge further from the */
/*             < convex body should be removed. The plane can not be */
/*             < killed in this case so we do not test for it. */
		    parallel_edges__();
		}
	    } else {
/*           < The edges are not parallel (or anti-parallel for that */
/*           < matter) and therefore they intersect.  We must find a */
/*           < non-zero element of the cross-product.  This will be the */
/*           < determinate of the linear system we need to solve to */
/*           < find the point of intersection.  Find the largest */
/*           < element of the cross-product. */
/* Computing MAX */
		d__1 = abs(edgepair_1.abxcd[0]), d__2 = abs(edgepair_1.abxcd[
			1]), d__1 = max(d__1,d__2), d__2 = abs(
			edgepair_1.abxcd[2]);
		maxabxcd = max(d__1,d__2);
		if (abs(edgepair_1.abxcd[0]) == maxabxcd) {
		    edgepair_1.kth = 1;
		    edgepair_1.ith = 2;
		    edgepair_1.jth = 3;
		} else if (abs(edgepair_1.abxcd[1]) == maxabxcd) {
		    edgepair_1.kth = 2;
		    edgepair_1.ith = 3;
		    edgepair_1.jth = 1;
		} else {
		    edgepair_1.kth = 3;
		    edgepair_1.ith = 1;
		    edgepair_1.jth = 2;
		}
		three_plane__();
/*           < Has this plane been killed? */
		if (thisplane_1.numedges == -1) {
		    return 0;
		}
	    }
	}
    }
    return 0;
} /* find_ap__ */

/* ---------------------------------------------------------------------C */
/*     Parallel Edges                                                  C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int parallel_edges__(void)
{
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal rac[3];
    static integer elim;
    extern /* Subroutine */ int coincidental_(void);
    static doublereal acxab[3], magsq;
    extern /* Subroutine */ int cross_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal acxab_dollar_ni__;

/*     < LOCAL */
/*     < Form the vector connecting the start of the AB vector */
/*     < with the start of the CD vector. */
    rac[0] = edgepair_1.rc[0] - edgepair_1.ra[0];
    rac[1] = edgepair_1.rc[1] - edgepair_1.ra[1];
    rac[2] = edgepair_1.rc[2] - edgepair_1.ra[2];
/*     < Find the cross product of this vector with the AB vector */
    cross_(rac, edgepair_1.eab, acxab);
    magsq = acxab[0] * acxab[0] + acxab[1] * acxab[1] + acxab[2] * acxab[2];
/*     < Find out if these two edges coincide. */
    if (magsq < 1e-8) {
/*       < They coincide. */
	coincidental_();
    } else {
/*       < Find the dot product of the cross product ACxAB with the */
/*       < normal vector of plane PlaneNum, ni */
	acxab_dollar_ni__ = acxab[0] * thisplane_1.ni[0] + acxab[1] * 
		thisplane_1.ni[1] + acxab[2] * thisplane_1.ni[2];
/*       < If (ACxAB).ni is positive then edge AB is on the "inside" */
/*       < (closer to the cutout) and edge CD should be removed. */
/*       < Otherwise, (ACxAB).ni negative, edge AB is on the "outside" */
/*       < and should be removed. */
/*       <   AB Inside --> ACxAB.ni > 0 --> ELIM= EdgeJ */
/*       <   CD Inside --> ACxAB.ni < 0 --> ELIM= EdgeI */
	elim = edgepair_1.edgei + (edgepair_1.edgej - edgepair_1.edgei) * (
		integer) (d_sign(&c_b31, &acxab_dollar_ni__) + 1.);
/*       < Remove the edge ELIM. */
	apexinfo_1.apexf[thisplane_1.edgesi[elim - 1] - 1] = -1;
	apexinfo_1.apexf[thisplane_1.edgesi[elim + 5999] - 1] = -1;
    }
    return 0;
} /* parallel_edges__ */

/* ---------------------------------------------------------------------C */
/*     Parallel Coincidental Edges                                     C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int coincidental_(void)
{
/*     < We have two edges on a plane face that are parallel and */
/*     < coincidental.  This means that three planes are coming */
/*     < together to form a three plane line.  On two plane faces we */
/*     < will have coincidental parallel edges and on a third plane */
/*     < faces we will have coincidental anti-parallel edges.  One of */
/*     < the parallel coincidental edges will need to be perserved on */
/*     < the plane faces.  For the present calculations we should */
/*     < arbitrarily kill one of the edges if both are still alive */
/*     < since they are numerically identical. However, for purposes of */
/*     < bookkeeping they are different. This will need to be rectified */
/*     < for later work dealing with festoon formation using twin */
/*     < information. */
    if (apexinfo_1.apexf[edgepair_1.apexa - 1] != -1) {
	if (apexinfo_1.apexf[edgepair_1.apexc - 1] != -1) {
/*         < Both the AB edge and the CD edge are alive.  Have either */
/*         < of these edges already been involved in a parallel */
/*         < coincidental plane face? */
	    if (twinlist_1.twin[edgepair_1.apexa - 1] != 0) {
/*           < The AB edge was previously involved in a parallel */
/*           < coincidence on this or another plane face. This means */
/*           < that the AB edge should continue to live. Kill the CD */
/*           < edge letting it adopt AB as its twin. */
		apexinfo_1.apexf[edgepair_1.apexc - 1] = -1;
		apexinfo_1.apexf[edgepair_1.apexd - 1] = -1;
		twinlist_1.twin[edgepair_1.apexc - 1] = edgepair_1.apexa;
		twinlist_1.twin[edgepair_1.apexd - 1] = edgepair_1.apexb;
	    } else if (twinlist_1.twin[edgepair_1.apexc - 1] != 0) {
/*           < The CD edge was previously involved in a parallel */
/*           < coincidence on this or another plane face. This means */
/*           < that the CD edge should continue to live. Kill the AB */
/*           < edge letting it adopt CD as its twin. */
		apexinfo_1.apexf[edgepair_1.apexa - 1] = -1;
		apexinfo_1.apexf[edgepair_1.apexb - 1] = -1;
		twinlist_1.twin[edgepair_1.apexa - 1] = edgepair_1.apexc;
		twinlist_1.twin[edgepair_1.apexb - 1] = edgepair_1.apexd;
	    } else {
/*           < Neither of the edges was previously involved in a */
/*           < parallel coincidence. We must mark each apex with its */
/*           < twin and then arbitrarily kill one of the edges (i.e., */
/*           < edge AB) . The array Twin will contain these twin */
/*           < mappings.  For instance, ApexA is identical to ApexC. We */
/*           < give Edge CD the negative apex numbers of AB to indicate */
/*           < that AB was killed in this confrontation and CD was kept */
/*           < alive. */
		twinlist_1.twin[edgepair_1.apexa - 1] = edgepair_1.apexc;
		twinlist_1.twin[edgepair_1.apexb - 1] = edgepair_1.apexd;
		twinlist_1.twin[edgepair_1.apexc - 1] = -edgepair_1.apexa;
		twinlist_1.twin[edgepair_1.apexd - 1] = -edgepair_1.apexb;
		apexinfo_1.apexf[edgepair_1.apexa - 1] = -1;
		apexinfo_1.apexf[edgepair_1.apexb - 1] = -1;
	    }
	} else if (twinlist_1.twin[edgepair_1.apexc - 1] > 0) {
/*         < Edge AB is still alive while edge CD is dead.  We also */
/*         < know that Edge CD was previously involved in a parallel */
/*         < coincidence but was improperly and arbitrarily killed. If */
/*         < either the start or end of Edge AB is still a two-plane */
/*         < intersection and if the corresponding start or end of CD's */
/*         < twin is still a two-plane intersection then the vertices */
/*         < should be transferred.  These are assigned in RenumArc for */
/*         < use in forming the vertices. */
	    if (apexinfo_1.apexf[edgepair_1.apexa - 1] == 2 && 
		    apexinfo_1.apexf[twinlist_1.twin[edgepair_1.apexc - 1] - 
		    1] == 2) {
		apexinfo_1.apexf[edgepair_1.apexa + 5999] = apexinfo_1.apexf[
			twinlist_1.twin[edgepair_1.apexc - 1] + 5999];
	    }
	    if (apexinfo_1.apexf[edgepair_1.apexb - 1] == 2 && 
		    apexinfo_1.apexf[twinlist_1.twin[edgepair_1.apexd - 1] - 
		    1] == 2) {
		apexinfo_1.apexf[edgepair_1.apexb + 5999] = apexinfo_1.apexf[
			twinlist_1.twin[edgepair_1.apexd - 1] + 5999];
	    }
	}
    } else if (apexinfo_1.apexf[edgepair_1.apexc - 1] != -1 && 
	    twinlist_1.twin[edgepair_1.apexa - 1] > 0) {
/*       < Edge CD is still alive while edge AB is dead.  We also know */
/*       < that Edge AB was previously involved in a parallel */
/*       < coincidence but was improperly and arbitrarily killed. If */
/*       < either the start or end of Edge CD is still a two-plane */
/*       < intersection and if the corresponding start or end of AB's */
/*       < twin is still a two-plane intersection then the vertices */
/*       < should be transferred.  These are assigned in RenumArc for */
/*       < use in forming the vertices. */
	if (apexinfo_1.apexf[edgepair_1.apexc - 1] == 2 && apexinfo_1.apexf[
		twinlist_1.twin[edgepair_1.apexa - 1] - 1] == 2) {
	    apexinfo_1.apexf[edgepair_1.apexc + 5999] = apexinfo_1.apexf[
		    twinlist_1.twin[edgepair_1.apexa - 1] + 5999];
	}
	if (apexinfo_1.apexf[edgepair_1.apexd - 1] == 2 && apexinfo_1.apexf[
		twinlist_1.twin[edgepair_1.apexb - 1] - 1] == 2) {
	    apexinfo_1.apexf[edgepair_1.apexd + 5999] = apexinfo_1.apexf[
		    twinlist_1.twin[edgepair_1.apexb - 1] + 5999];
	}
    }
    return 0;
} /* coincidental_ */

/* ---------------------------------------------------------------------C */
/*     AParallel Edges                                                 C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int aparallel_edges__(void)
{
    static doublereal rac[3], acxab[3], magsq;
    extern /* Subroutine */ int cross_(doublereal *, doublereal *, doublereal 
	    *);

/*     < LOCAL */
/*     < Form the vector connecting the start of the AB vector */
/*     < with the start of the CD vector. */
    rac[0] = edgepair_1.rc[0] - edgepair_1.ra[0];
    rac[1] = edgepair_1.rc[1] - edgepair_1.ra[1];
    rac[2] = edgepair_1.rc[2] - edgepair_1.ra[2];
/*     < Find the cross product of this vector with the AB vector */
    cross_(rac, edgepair_1.eab, acxab);
    magsq = acxab[0] * acxab[0] + acxab[1] * acxab[1] + acxab[2] * acxab[2];
/*     < If the magnitude of the cross-product is zero then we know */
/*     < that AB and CD are anti-parallel and coincide.  The */
/*     < arc-polygon is contained between these two vectors and since */
/*     < AB and CD coincide it means that the arc-polygon does not */
/*     < exist at all.  Turn off this plane face. */
    if (magsq < 1e-8) {
/*       < Turn-off plane */
	thisplane_1.numedges = -1;
	apexinfo_1.n_edges__[thisplane_1.planenum - 1] = -1;
    } else if (acxab[0] * thisplane_1.ni[0] + acxab[1] * thisplane_1.ni[1] + 
	    acxab[2] * thisplane_1.ni[2] > 0.) {
/*       < The anti-parallel edges do not coincide but they do violate */
/*       < the laws of convexity. The triple product with the normal */
/*       < vector of this plane is positive. The anti-parallel edges */
/*       < are positioned clockwise on this plane face and so there is */
/*       < no arc-polygon and this plane must be killed. */
	thisplane_1.numedges = -1;
	apexinfo_1.n_edges__[thisplane_1.planenum - 1] = -1;
    }
    return 0;
} /* aparallel_edges__ */

/* ---------------------------------------------------------------------C */
/*     Three Plane Intersection                                        C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int three_plane__(void)
{
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer k, l;
    static logical edgealive[2], elimination;
    extern /* Subroutine */ int find_lambdas__(doublereal *, doublereal *), 
	    move_(void);
    static integer leading, lagging;
    static doublereal abxcd_dollar_ni__;

/*     < LOCAL */
/*     < Find the Lambdas for the point of intersection. */
    find_lambdas__(formove_1.lambda, &formove_1.nu1);
/*     < Now that we have the Lambda's we can proceed onto the moving of */
/*     < apices and the elimination of edges.  This whole algorithm is */
/*     < based on two planes intersecting on another plane face and the */
/*     < ordering of the edges formed in a right-hand sense. */
/*     < The first step is to find out which of the two edges is the */
/*     < leading edge vector in a right-handed sense.  That is, which of */
/*     < the two edge vector when crossed with the other and dotted with */
/*     < the normal vector of the plane yields a positive number.  Form */
/*     < this triple product. */
    abxcd_dollar_ni__ = edgepair_1.abxcd[0] * thisplane_1.ni[0] + 
	    edgepair_1.abxcd[1] * thisplane_1.ni[1] + edgepair_1.abxcd[2] * 
	    thisplane_1.ni[2];
/*     < If ABxCD.ni is positive then K will be zero while if ABxCD.ni */
/*     < is negative then K will be unity: */
/*     <   AB Leading --> ABxCD.ni > 0 --> K= 0 --> L= 0 */
/*     <   CD Leading --> ABxCD.ni < 0 --> K= 1 --> L= EdgeJ - EdgeI */
/*     < Thus, below First will equal EdgeI in the former case and */
/*     < EdgeJ in the later case. */
    k = (integer) (1. - d_sign(&c_b31, &abxcd_dollar_ni__));
    l = k * (edgepair_1.edgej - edgepair_1.edgei);
    formove_1.first = edgepair_1.edgei + l;
    formove_1.second = edgepair_1.edgej - l;
    leading = k + 1;
    lagging = 2 - k;
/*     < Retrieve the proper Lambda's based on the value of K. */
    formove_1.lambdaf = formove_1.lambda[leading - 1];
    formove_1.lambdas = formove_1.lambda[lagging - 1];
/*     < Find out who is alive and who is dead. */
    edgealive[0] = apexinfo_1.apexf[edgepair_1.apexa - 1] != -1;
    edgealive[1] = apexinfo_1.apexf[edgepair_1.apexc - 1] != -1;
/*     < Now that we know the values of Lambda corresponding to the */
/*     < leading and lagging edges we can proceed. */
/*     < */
/*     < When two planes intersect they form a line or an edge that */
/*     < would appear on both plane faces.  When three planes intersect */
/*     < they form a point and so the face of each plane involved would */
/*     < have two edges and a point appearing on it.  So our interest is */
/*     < in the intersection of two lines or edges.  More generally our */
/*     < interest is in the intersection of two rays eminating from the */
/*     < edge vectors.  The picture for the intersection of the two edge */
/*     < rays is the following (this is the "template") */
/*     < */
/*     <                        o           o */
/*     <                         \         / */
/*     <             +--------->  \       / */
/*     <             |             \     / */
/*     <        These sections      \   / */
/*     <        need to be           \ /       Region of */
/*     <        discarded             X rT     desired */
/*     <             |               / \       cutout */
/*     <             +--------->    /   \ */
/*     <                           /     \ */
/*     <                          /       \ */
/*     <                     mu1 /         \ mu2 */
/*     <                        V_         _V */
/*     <                   Leading      Lagging */
/*     <                    Edge         Edge */
/*     <                    Ray          Ray */
/*     < */
/*     < The point rT is the triple plane intersection.  The parts of */
/*     < the real edge vectors that lie on the sections to the left will */
/*     < need to be discarded and if the entire edge vector lies in */
/*     < these sections then the entire edge should be eliminated.  So */
/*     < the matter of interest is where the edge vectors lie on this */
/*     < coordinate system relative to the point rT. */
/*     < */
/*     < Each apex of the edge vectors (i.e., the start and terminus of */
/*     < each edge vector) will be transformed from r-space or Cartesian */
/*     < space to mu-space where the two rays eminating from the edge */
/*     < vectors will be the mu's with rT as the origin.  We do this by */
/*     < defining */
/*     < */
/*     < (mu1,mu2)= (   -LambdaF,0) for the start of the leading edge */
/*     < (mu1,mu2)= (1 - LambdaF,0) for the terminus of the leading edge */
/*     < (mu1,mu2)= (0,   -LambdaS) for the start of the lagging edge */
/*     < (mu1,mu2)= (0,1 - LambdaS) for the terminus of the lagging edge */
/*     < */
/*     < The rules for dealing with these edges are then */
/*     < */
/*     <    o  any apex with mu1 > +Epsilon should be moved to mu1= 0 */
/*     <       this is a clipping of the edge that apex is on and is a */
/*     <       result of a three-plane intersection. */
/*     <    o  any apex with mu2 < -Epsilon should be moved to mu2= 0 */
/*     <       this is a clipping of the edge that apex is on and is a */
/*     <       result of a three-plane intersection. */
/*     <    o  any edge with (mu1,mu2)= (0,0) for both apices should be */
/*     <       eliminated this is a requirement for convexity of the */
/*     <       final arc-polygon. */
/*     < */
/*     < Epsilon (a very small positive number) is included in the two */
/*     < inequalities because when two edges intersect to form a */
/*     < three-plane intersection, rT, they will be both be clipped back */
/*     < to form two edges of the final arc-polygon for the present */
/*     < plane face.  Each of these clipped edges (call them 1 and 2) */
/*     < will appear separately on subsequent plane faces and will */
/*     < terminate or originate from the triple-plane intersection rT. */
/*     < One of the edges on each of these subsequent plane faces will */
/*     < be the previously unknown third edge that forms the */
/*     < triple-plane intersection.  Thus, it will be found that clipped */
/*     < edge vector 1 will terminate or originate at the point rT */
/*     < somewhere on the third edge (thus its edge vector 1's Lambda */
/*     < will be ~1 or ~0 while the third edge vector's Lambda will be */
/*     < in (0,1)) and we will need to clip this third edge at the */
/*     < triple point intersection.  Finally, when we arrive at the */
/*     < third plane face of the three planes involved in the three */
/*     < plane intersection, the clipped edge vector 2 will terminate or */
/*     < originate at rT on the start or terminus of the now clipped */
/*     < third edge vector.  In this case the one Lambda will be ~0 and */
/*     < another will be ~1.  From above then the three planes i, j and */
/*     < k will form a triad with edges 1, 2 and 3 meeting at point rT */
/*     < */
/*     <                                (1) */
/*     <                      Plane     / */
/*     <                      Face j   / */
/*     <                              /    Plane */
/*     <                  (3)--------< rT  Face i */
/*     <                              \ */
/*     <                      Plane    \ */
/*     <                      Face k    \ */
/*     <                                (2) */
/*     < */
/*     < Try to imagine point rT above the plane of the screen (or paper */
/*     < or whatever media you are using to read this program) and the */
/*     < three plane faces i, j and k all tilting up converging at it. */
/*     < */
/*     < Thus edge vectors 1 and 2 clipped each other on plane face i */
/*     < forming the point rT but edge vector 3 was unknown on that */
/*     < plane face since it was formed by the intersection of planes j */
/*     < and k and when executing plane face j edge vector 1 will be */
/*     < found to originate on edge vector 3 (Lambda1~0 and Lambda3 in */
/*     < (0,1)).  We will clip edge vector 3 and then when examining */
/*     < plane face k we will find that edge vector 2 terminates at the */
/*     < start of edge vector 3 (Lambda2~1 and Lambda3~0).  To make a */
/*     < long story come to an end, the Epsilon is needed to handle the */
/*     < cases of ~0 and ~1. */
/*     < */
/*     < The other possible values of Lambda cover the cases where one */
/*     < edge is farther from the arc-polygon than the other vector and */
/*     < is either eliminated entirely or partially clipped. It also */
/*     < covers the cases where nothing should be done.  That is, the */
/*     < two edges are entirely contained in the rays closest to the */
/*     < cutout.  Another interesting case that is covered is when we */
/*     < get a four plane intersection (yes, they can occur) in which */
/*     < case the two of the edge vectors either have their starts or */
/*     < their termini touching.  One of these two edge vectors will */
/*     < reside on one of the rays closest to the cutout and will be */
/*     < kept.  In this case, the other edge vector will be eliminated. */
/*     < It is also possible for neither of the edge vectors to appear */
/*     < on the rays closest to the cutout and thus both will be */
/*     < eliminated. */
/*     < We are ready to test the values of the Lambdas.  Initialize the */
/*     < flags that will let us determine whether or not we have a */
/*     < double clipping, that is, a true intersection. */
    elimination = FALSE_;
    formove_1.movement = 0;
/*     < ---------------------------- */
/*     < Start with the LEADING EDGE. */
/*     < ---------------------------- */
    if (-formove_1.lambdaf > -5e-7) {
/*       < Since (-LambdaF > -Eps) we can add one to both sides and */
/*       < write the following */
/*       < */
/*       <        1 - LambdaF > 1 - Eps  > Eps */
/*       < */
/*       < This means that the leading edge begins and ends in the */
/*       < undesired region and thus should be eliminated. */
/*       < */
/*       <   Lagging */
/*       <     o */
/*       <      \ */
/*       <       \ */
/*       <        \ o */
/*       <         %   <---  rT */
/*       <        / \ */
/*       <       /   \ */
/*       <      /     \ */
/*       <     V_     _V */
/*       <   Leading */
/*       < */
/*       < Kill both of the apices of the First edge. */
/*       < ELIMINATION */
	apexinfo_1.apexf[thisplane_1.edgesi[formove_1.first - 1] - 1] = -1;
	apexinfo_1.apexf[thisplane_1.edgesi[formove_1.first + 5999] - 1] = -1;
	elimination = TRUE_;
    } else if (1. - formove_1.lambdaf > 5e-7 && edgealive[leading - 1]) {
/*       < We know that (-LambdaF < -Eps) and now we also know that */
/*       < (1 - LambdaF > Eps) so we have */
/*       < */
/*       <   Lagging */
/*       <     o */
/*       <      \ */
/*       <       \   o */
/*       <        \ / */
/*       <         %   <---  rT */
/*       <        / \ */
/*       <       /   \ */
/*       <      /     \ */
/*       <     V_     _V */
/*       <   Leading */
/*       < */
/*       < This means that the leading edge starts at least Eps into */
/*       < the desired region but its terminus is in the undesired */
/*       < region for the leading edge.  We have to move (clip) the end */
/*       < of the first edge back to rT. */
	formove_1.movement = 1;
    }
/*     < ------------------------ */
/*     < Now do the LAGGING EDGE. */
/*     < ------------------------ */
    if (1. - formove_1.lambdas < 5e-7) {
/*       < Since (1 - LambdaS < Eps) we can subtract one from both */
/*       < sides and write the following */
/*       < */
/*       <        - LambdaS < Eps - 1 < -Eps */
/*       < */
/*       < This means that the lagging edge should be eliminated. */
/*       < */
/*       <   Lagging */
/*       <     o       o */
/*       <      \     / */
/*       <       \   / */
/*       <        \ / */
/*       <         %   <---  rT */
/*       <        /_V */
/*       <       / */
/*       <      / */
/*       <     V_ */
/*       <   Leading */
/*       < */
/*       < Kill both apices of the Second edge. */
/*       < ELIMINATION */
	if (elimination) {
/*         < Both edges have been entirely eliminated from this plane */
/*         < face indicating that they were clockwise to each other and */
/*         < so the cutout was not between them.  This means that this */
/*         < plane should be killed. */
	    thisplane_1.numedges = -1;
	    apexinfo_1.n_edges__[thisplane_1.planenum - 1] = -1;
	} else {
/*         < The leading edge was not killed entirely.  Just kill the */
/*         < lagging edge. */
	    apexinfo_1.apexf[thisplane_1.edgesi[formove_1.second - 1] - 1] = 
		    -1;
	    apexinfo_1.apexf[thisplane_1.edgesi[formove_1.second + 5999] - 1] 
		    = -1;
	}
    } else if (-formove_1.lambdas < -5e-7 && edgealive[lagging - 1]) {
/*       < We know that (1 - LambdaS > Eps) and now we also know that */
/*       < (-LambdaS < -Eps) so we have */
/*       < */
/*       <   Lagging */
/*       <     o       o */
/*       <      \     / */
/*       <       \   / */
/*       <        \ / */
/*       <         %   <---  rT */
/*       <        / \ */
/*       <       /  _V */
/*       <      / */
/*       <     V_ */
/*       <   Leading */
/*       < */
/*       < This means that the terminus of the lagging edge is at least */
/*       < Eps into the desired region but the start of the lagging */
/*       < edge is in the undesired region.  We have to move (clip) the */
/*       < start of the lagging edge up to rT. */
	formove_1.movement += 2;
    }
/*     < Do we have a movement of an apex? */
    if (formove_1.movement != 0) {
	move_();
    }
    return 0;
} /* three_plane__ */

/* ---------------------------------------------------------------------C */
/*     Move                                                            C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int move_(void)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int ordering_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer p2, p3, apex1, apex2;

/*     < LOCAL */
/*     < Pull out the location of second apex of leading edge and the */
/*     < location of first apex of lagging edge. */
    apex1 = thisplane_1.edgesi[formove_1.first + 5999];
    apex2 = thisplane_1.edgesi[formove_1.second - 1];
    if (formove_1.movement == 3) {
/*       < Both apices have to be moved.  We have an intersection of */
/*       < the two edges on the plane face. */
/*       < */
/*       <   Apex2 */
/*       <     o       o */
/*       <      \     / */
/*       <       \   / */
/*       <        \ / */
/*       <         X   <---  rT */
/*       <        / \ */
/*       <       /   \ */
/*       <      /     \ */
/*       <     V_     _V */
/*       <   Apex1 */
/*       < */
/*       < Replace coordinates of Apex1 with those of rT. */
	apexinfo_1.apexx[apex1 - 1] = edgepair_1.ra[0] + formove_1.nu1 * 
		edgepair_1.eab[0];
	apexinfo_1.apexy[apex1 - 1] = edgepair_1.ra[1] + formove_1.nu1 * 
		edgepair_1.eab[1];
	apexinfo_1.apexz[apex1 - 1] = edgepair_1.ra[2] + formove_1.nu1 * 
		edgepair_1.eab[2];
/*       < Move Apex2 to the new Apex1. */
	apexinfo_1.apexx[apex2 - 1] = apexinfo_1.apexx[apex1 - 1];
	apexinfo_1.apexy[apex2 - 1] = apexinfo_1.apexy[apex1 - 1];
	apexinfo_1.apexz[apex2 - 1] = apexinfo_1.apexz[apex1 - 1];
/*       < Increment the triple plane intersection counter. */
	++sphereinfo_1.tripcnt;
/*       < Flag this apex as resulting from 3-plane intersection then */
/*       < store the number of the 3-plane intersection for future */
/*       < reference. */
	apexinfo_1.apexf[apex1 - 1] = 3;
	apexinfo_1.apexf[apex1 + 5999] = sphereinfo_1.tripcnt;
	apexinfo_1.apexf[apex2 - 1] = 3;
	apexinfo_1.apexf[apex2 + 5999] = sphereinfo_1.tripcnt;
	ordering_(&apex1, &apex2, &thisplane_1.planenum, &p2, &p3);
	debug_1.plane_order__[apex1 - 1] = thisplane_1.planenum;
	debug_1.plane_order__[apex1 + 5999] = p2;
	debug_1.plane_order__[apex1 + 11999] = p3;
	debug_1.plane_order__[apex2 - 1] = thisplane_1.planenum;
	debug_1.plane_order__[apex2 + 5999] = p2;
	debug_1.plane_order__[apex2 + 11999] = p3;
    } else if (formove_1.movement == 1) {
	if (abs(formove_1.lambdas) < 1e-8) {
/*         < Only the terminus of the leading edge needs to be moved. */
/*         < Since Movement is +1 this means that the start of the */
/*         < leading edge is at least Eps into the desired region and */
/*         < we have this picture. */
/*         < */
/*         <             o */
/*         <            / */
/*         <           / */
/*         <          / */
/*         <         %   <---  rT AND Apex2 */
/*         <        / \ */
/*         <       /   \ */
/*         <      /     \ */
/*         <     V_     _V */
/*         <   Apex1 */
/*         < */
/*         < That is, we are sure that the beginning of the leading */
/*         < edge is in the desired section and therefore it needs to */
/*         < be clipped back and that the entire leading edge need not */
/*         < be eliminated.  So only Apex1 has to be moved and the */
/*         < lagging edge starts on the leading edge.  This is a */
/*         < degenerate intersection.  Move Apex1 to start of lagging */
/*         < edge and adopt all of the characteristics of Apex2. */
	    apexinfo_1.apexx[apex1 - 1] = apexinfo_1.apexx[apex2 - 1];
	    apexinfo_1.apexy[apex1 - 1] = apexinfo_1.apexy[apex2 - 1];
	    apexinfo_1.apexz[apex1 - 1] = apexinfo_1.apexz[apex2 - 1];
	    apexinfo_1.apexf[apex1 - 1] = 3;
	    apexinfo_1.apexf[apex1 + 5999] = apexinfo_1.apexf[apex2 + 5999];
	    debug_1.plane_order__[apex1 - 1] = debug_1.plane_order__[apex2 - 
		    1];
	    debug_1.plane_order__[apex1 + 5999] = debug_1.plane_order__[apex2 
		    + 5999];
	    debug_1.plane_order__[apex1 + 11999] = debug_1.plane_order__[
		    apex2 + 11999];
	} else {
/*         < Only Apex1 has to be moved but the start of the lagging */
/*         < does not make actual contact with the leading edge. */
/*         < */
/*         <             o */
/*         <            / */
/*         <           / */
/*         <          / */
/*         <  rT --> % */
/*         <        / */
/*         <       /    o Apex2 */
/*         <      /      \ */
/*         <     V_       \ */
/*         <   Apex1       \ */
/*         <               _V */
/*         < */
/*         < Replace coordinates of Apex1 with those of rT. */
	    apexinfo_1.apexx[apex1 - 1] = edgepair_1.ra[0] + formove_1.nu1 * 
		    edgepair_1.eab[0];
	    apexinfo_1.apexy[apex1 - 1] = edgepair_1.ra[1] + formove_1.nu1 * 
		    edgepair_1.eab[1];
	    apexinfo_1.apexz[apex1 - 1] = edgepair_1.ra[2] + formove_1.nu1 * 
		    edgepair_1.eab[2];
/*         < Increment the triple plane intersection counter. */
	    ++sphereinfo_1.tripcnt;
/*         < Let Apex1 adopt all of the characteristics of Apex2 except */
/*         < for its TripCnt number. */
	    apexinfo_1.apexf[apex1 - 1] = 3;
	    apexinfo_1.apexf[apex1 + 5999] = sphereinfo_1.tripcnt;
/*         < Assign plane ordering to the new end of the leading edge as */
/*         < if the lagging edge truly did pass through the leading */
/*         < edge and we had real contact.  We will not modify the plane */
/*         < ordering of Apex2 however. */
	    ordering_(&apex1, &apex2, &thisplane_1.planenum, &p2, &p3);
	    debug_1.plane_order__[apex1 - 1] = thisplane_1.planenum;
	    debug_1.plane_order__[apex1 + 5999] = p2;
	    debug_1.plane_order__[apex1 + 11999] = p3;
	}
	if (apexinfo_1.apexf[apex2 - 1] == 2) {
/*         < RARE CASE: The start of the lagging edge is still a */
/*         < two-plane intersection even though it lies on the leading */
/*         < edge (more or less). This means that three planes are */
/*         < intersecting on the sphere surface.  Act as if these two */
/*         < edges intersected normally forming a three-plane point. */
	    apexinfo_1.apexx[apex2 - 1] = apexinfo_1.apexx[apex1 - 1];
	    apexinfo_1.apexy[apex2 - 1] = apexinfo_1.apexy[apex1 - 1];
	    apexinfo_1.apexz[apex2 - 1] = apexinfo_1.apexz[apex1 - 1];
	    apexinfo_1.apexf[apex2 - 1] = 3;
	    apexinfo_1.apexf[apex2 + 5999] = apexinfo_1.apexf[apex1 + 5999];
	    debug_1.plane_order__[apex2 - 1] = debug_1.plane_order__[apex1 - 
		    1];
	    debug_1.plane_order__[apex2 + 5999] = debug_1.plane_order__[apex1 
		    + 5999];
	    debug_1.plane_order__[apex2 + 11999] = debug_1.plane_order__[
		    apex1 + 11999];
	}
    } else {
/*       < Movement must equal 2. Only Apex2 has to be moved. */
/*       < Does the end of the leading edge sit on the lagging edge? */
	if ((d__1 = 1. - formove_1.lambdaf, abs(d__1)) < 1e-8) {
/*         < Only the start of the lagging needs to be moved.  Since */
/*         < Movement is +2 this means that the terminus of the lagging */
/*         < edge is at least Eps into the desired region and we have */
/*         < this picture. */
/*         < */
/*         <   Apex2 */
/*         <     o       o */
/*         <      \     / */
/*         <       \   / */
/*         <        \ / */
/*         <         V_  <---  rT AND Apex1 */
/*         <          \ */
/*         <           \ */
/*         <            \ */
/*         <            _V */
/*         < */
/*         < That is, we are sure that the end of the lagging edge is */
/*         < in the desired section and therefore it needs to be */
/*         < clipped back and that the entire lagging edge need not be */
/*         < eliminated.  So only Apex2 has to be moved and the leading */
/*         < edge ends on the lagging edge.  This is a degenerate */
/*         < intersection.  Move Apex2 to end of the leading edge and */
/*         < have it adopt all of the characteristics of Apex1. */
	    apexinfo_1.apexx[apex2 - 1] = apexinfo_1.apexx[apex1 - 1];
	    apexinfo_1.apexy[apex2 - 1] = apexinfo_1.apexy[apex1 - 1];
	    apexinfo_1.apexz[apex2 - 1] = apexinfo_1.apexz[apex1 - 1];
	    apexinfo_1.apexf[apex2 - 1] = 3;
	    apexinfo_1.apexf[apex2 + 5999] = apexinfo_1.apexf[apex1 + 5999];
	    debug_1.plane_order__[apex2 - 1] = debug_1.plane_order__[apex1 - 
		    1];
	    debug_1.plane_order__[apex2 + 5999] = debug_1.plane_order__[apex1 
		    + 5999];
	    debug_1.plane_order__[apex2 + 11999] = debug_1.plane_order__[
		    apex1 + 11999];
	} else {
/*         < The end of the leading edge cuts the lagging edge but does */
/*         < not make actual contact with it. */
/*         < */
/*         <               o */
/*         <   Apex2      / */
/*         <     o       / */
/*         <      \     / */
/*         <       \   V_ Apex1 */
/*         <        \ */
/*         <         %  <---  rT */
/*         <          \ */
/*         <           \ */
/*         <            \ */
/*         <            _V */
/*         < */
/*         < Replace coordinates of Apex2 with those of rT. */
	    apexinfo_1.apexx[apex2 - 1] = edgepair_1.ra[0] + formove_1.nu1 * 
		    edgepair_1.eab[0];
	    apexinfo_1.apexy[apex2 - 1] = edgepair_1.ra[1] + formove_1.nu1 * 
		    edgepair_1.eab[1];
	    apexinfo_1.apexz[apex2 - 1] = edgepair_1.ra[2] + formove_1.nu1 * 
		    edgepair_1.eab[2];
/*         < Increment the triple plane intersection counter. */
	    ++sphereinfo_1.tripcnt;
/*         < Let Apex2 adopt all of the characteristics of Apex1 except */
/*         < for the TripCnt of Apex1. */
	    apexinfo_1.apexf[apex2 - 1] = 3;
	    apexinfo_1.apexf[apex2 + 5999] = sphereinfo_1.tripcnt;
/*         < Assign plane ordering to the new start of the lagging edge */
/*         < as if the leading edge truly did pass through the lagging */
/*         < edge and we had real contact.  We will not modify the plane */
/*         < ordering of Apex1 however. */
	    ordering_(&apex1, &apex2, &thisplane_1.planenum, &p2, &p3);
	    debug_1.plane_order__[apex2 - 1] = thisplane_1.planenum;
	    debug_1.plane_order__[apex2 + 5999] = p2;
	    debug_1.plane_order__[apex2 + 11999] = p3;
	}
	if (apexinfo_1.apexf[apex1 - 1] == 2) {
/*         < RARE CASE: The end of the leading edge is still a */
/*         < two-plane intersection even though it lies on the lagging */
/*         < edge (more or less). This means that three planes are */
/*         < intersecting on the sphere surface. Act as if these two */
/*         < edges intersected normally forming a three-plane point. */
	    apexinfo_1.apexx[apex1 - 1] = apexinfo_1.apexx[apex2 - 1];
	    apexinfo_1.apexy[apex1 - 1] = apexinfo_1.apexy[apex2 - 1];
	    apexinfo_1.apexz[apex1 - 1] = apexinfo_1.apexz[apex2 - 1];
	    apexinfo_1.apexf[apex1 - 1] = 3;
	    apexinfo_1.apexf[apex1 + 5999] = apexinfo_1.apexf[apex2 + 5999];
	    debug_1.plane_order__[apex1 - 1] = debug_1.plane_order__[apex2 - 
		    1];
	    debug_1.plane_order__[apex1 + 5999] = debug_1.plane_order__[apex2 
		    + 5999];
	    debug_1.plane_order__[apex1 + 11999] = debug_1.plane_order__[
		    apex2 + 11999];
	}
    }
    return 0;
} /* move_ */

/* ---------------------------------------------------------------------C */
/*     Find The Lambdas For The Point Of Intersection                  C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int find_lambdas__(doublereal *lambda, doublereal *nu1)
{
    static doublereal raci, racj, rdeter;

/*     < PASSED */
/*     < LOCAL */
/*     < The point of intersection between the two vectors */
/*     < rAB= rB-rA  and rCD= rD-rC is given by point rT */
/*     < */
/*     <   rT= rA + Lambda1*rAB */
/*     <   rT= rC + Lambda2*rCD */
/*     < */
/*     < where the Lambda's are constants. Thus we write */
/*     < */
/*     <   Lambda1*rAB - Lambda2*rCD= (rC - rA) */
/*     <                            = rAC */
/*     < */
/*     < therefore */
/*     < */
/*     < [1]  Lambda1*rAB(1) - Lambda2*rCD(1)= (rC(1) - rA(1))= rACx */
/*     < [2]  Lambda1*rAB(2) - Lambda2*rCD(2)= (rC(2) - rA(2))= rACy */
/*     < [3]  Lambda1*rAB(3) - Lambda2*rCD(3)= (rC(3) - rA(3))= rACz */
/*     < */
/*     < Three equations in two unknowns.  We express this as the */
/*     < general two by two system Ax= b where */
/*     < */
/*     <         [ rAB(i) -rCD(i) ]     [Lambda1]    [rAC(i)] */
/*     <      A= [ rAB(j) -rCD(j) ], x= [Lambda2], b=[rAC(j)] */
/*     < */
/*     < with (i,j) an element of {(1,2), (3,1), (2,3)}.  This system */
/*     < is non-singular iff the determinate of A is non-zero, that is */
/*     < A is non-singular iff */
/*     < */
/*     <   -[rABxrCD](k)= -[rAB(i)*rCD(j) - rAB(j)*rCD(i)] < > 0 */
/*     < */
/*     < where k is an element of {3,2,1}.  Solution of this system */
/*     < using Cramer's rule is given by */
/*     < */
/*     <    Lambda1= Deter(1)(k)/Deter(k) */
/*     <    Lambda2= Deter(2)(k)/Deter(k) */
/*     < */
/*     < where Deter(k) is the determinate of A and Deter(m)(k) is the */
/*     < determinate of A with the mth column replaced with b, all for */
/*     < a given value of k. */
/*     < */
/*     <    We will use the unit vectors so the expressions for */
/*     < Lambda1, Lambda2 are */
/*     < */
/*     <               1   rAC(i)*eCD(j) - rAC(j)*eCD(i) */
/*     <    Lambda1= ----- ----------------------------- */
/*     <             |rAB|       [eAB x eCD](k) */
/*     < */
/*     <               1   rAC(i)*eAB(j) - rAC(j)*eAB(i) */
/*     <    Lambda2= ----- ----------------------------- */
/*     <             |rCD|       [eAB x eCD](k) */
/*     < */
/*     < Note that the permissible (Kth,Ith,Jth) triplets are cyclic */
/*     < permutations of (1,2,3), that is, (1,2,3), (2,3,1), and */
/*     < (3,1,2). First find the reciprocal of Kth determinate */
    /* Parameter adjustments */
    --lambda;

    /* Function Body */
    rdeter = 1. / edgepair_1.abxcd[edgepair_1.kth - 1];
/*     < Find the Ith and Jth right-hand side. */
    raci = edgepair_1.rc[edgepair_1.ith - 1] - edgepair_1.ra[edgepair_1.ith - 
	    1];
    racj = edgepair_1.rc[edgepair_1.jth - 1] - edgepair_1.ra[edgepair_1.jth - 
	    1];
/*     < Form the Lambda's.  We will also store Lambda1*|rAB| for use */
/*     < in the expression */
/*     < */
/*     <   rT= rA + Lambda1*rAB */
/*     <     = rA + Lambda1*eAB*|rAB| */
/*     <     = rA + (Lambda1*|rAB|)*eAB */
/*     <     = rA + Nu1*eAB */
    *nu1 = (raci * edgepair_1.ecd[edgepair_1.jth - 1] - racj * edgepair_1.ecd[
	    edgepair_1.ith - 1]) * rdeter;
    lambda[1] = *nu1 * edgepair_1.rmagab;
    lambda[2] = (raci * edgepair_1.eab[edgepair_1.jth - 1] - racj * 
	    edgepair_1.eab[edgepair_1.ith - 1]) * rdeter * edgepair_1.rmagcd;
    return 0;
} /* find_lambdas__ */

/* ---------------------------------------------------------------------C */
/*     Find The Proper Plane Ordering For The Point Of Intersection    C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int ordering_(integer *apex1, integer *apex2, integer *
	planenum, integer *p2, integer *p3)
{
/*     < PASSED */
/*     < The ordering of the planes involved in the three-plane */
/*     < intersection is defined via this algorithm */
/*     < */
/*     <   o Locate the present plane in the ordered plane list of */
/*     <     the end of the leading edge. */
/*     <   o Pull out the plane immediately following this plane in */
/*     <     the cyclic permutation. */
/*     <   o Locate the present plane in the ordered plane list of */
/*     <     the start of the lagging edge. */
/*     <   o Pull out the plane immediately preceding this plane */
/*     <     in the cyclic permutation. */
/*     <   o Form the triple plane order by combining the present */
/*     <     plane with the two other planes located. */
    if (debug_1.plane_order__[*apex1 - 1] == *planenum) {
	*p2 = debug_1.plane_order__[*apex1 + 5999];
    } else if (debug_1.plane_order__[*apex1 + 5999] == *planenum) {
	*p2 = debug_1.plane_order__[*apex1 + 11999];
	if (*p2 == 0) {
	    *p2 = debug_1.plane_order__[*apex1 - 1];
	}
    } else {
	*p2 = debug_1.plane_order__[*apex1 - 1];
    }
    if (debug_1.plane_order__[*apex2 - 1] == *planenum) {
	*p3 = debug_1.plane_order__[*apex2 + 11999];
	if (*p3 == 0) {
	    *p3 = debug_1.plane_order__[*apex2 + 5999];
	}
    } else if (debug_1.plane_order__[*apex2 + 5999] == *planenum) {
	*p3 = debug_1.plane_order__[*apex2 - 1];
    } else {
	*p3 = debug_1.plane_order__[*apex2 + 5999];
    }
    return 0;
} /* ordering_ */

/* ---------------------------------------------------------------------C */
/*     Build The Arc-Polygon                                           C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int build_ap__(doublereal *eta_ap__, logical *
	error_encountered__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_sign(integer *, integer *), s_wsle(cilist *), do_lio(integer *, 
	    integer *, char *, ftnlen), e_wsle(void);

    /* Local variables */
    static doublereal arcareas, polyarea;
    static logical twoplane;
    static integer i__;
    extern /* Subroutine */ int form_arcs__(doublereal *);
    static integer trueedges;
    extern /* Subroutine */ int form_polygon__(doublereal *, integer *, 
	    logical *);
    static integer loc_i__, flag_i__;

    /* Fortran I/O blocks */
    static cilist io___113 = { 0, 0, 0, 0, 0 };
    static cilist io___114 = { 0, 0, 0, 0, 0 };
    static cilist io___115 = { 0, 0, 0, 0, 0 };
    static cilist io___116 = { 0, 0, 0, 0, 0 };
    static cilist io___117 = { 0, 0, 0, 0, 0 };
    static cilist io___118 = { 0, 0, 0, 0, 0 };
    static cilist io___119 = { 0, 0, 0, 0, 0 };
    static cilist io___120 = { 0, 0, 0, 0, 0 };
    static cilist io___121 = { 0, 0, 0, 0, 0 };
    static cilist io___122 = { 0, 0, 0, 0, 0 };


/*     < PASSED */
/*     < LOCAL */
    *error_encountered__ = FALSE_;
/*     < Initialize the arc-polygon areas. */
    arcareas = 0.;
    polyarea = 0.;
/*     < Sort the two-plane and the three-plane apices and the ends */
/*     < and starts. At the conclusion of this section we will have */
/*     < */
/*     <   NumStart(-1)= number of dead apex starts */
/*     <   NumStart(+2)= number of 2-plane apex starts */
/*     <   NumStart(+3)= number of 3-plane apex starts */
/*     < */
/*     <   LocStart(1:NumStart(-1),-1)= location of dead apex starts */
/*     <   LocStart(1:NumStart(+2),+2)= location of 2-plane apex starts */
/*     <   LocStart(1:NumStart(+3),+3)= location of 3-plane apex starts */
/*     < */
/*     < and likewise for the apex ends. Note that number of edge starts */
/*     < must equal the number of edge ends. We also will have an */
/*     < updated list of edges containing only living (non-dead) edges. */
/*     < */
/*     <   Edge_Num(*,-1)= edge numbers of all dead edges */
/*     <   Edge_Num(*,+1)= edge numbers of all alive (non-dead) edges */
/*     < */
/*     < where '*' represents '1:(NumStart(+2) + NumStart(+3))' */
/*     < Initialize the counters */
    sortedout_1.numstart[0] = 0;
    sortedout_1.numstart[3] = 0;
    sortedout_1.numstart[4] = 0;
    sortedout_1.numend[0] = 0;
    sortedout_1.numend[3] = 0;
    sortedout_1.numend[4] = 0;
/*     < Loop through the edges */
    i__1 = thisplane_1.numedges;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       < START */
/*       << Pull out the Edge I's start and its flag */
/*       << Flag_I will be equal to -1, +2 or +3 depending on */
/*       << whether the edge is dead, a 2-plane or a 3-plane */
	loc_i__ = thisplane_1.edgesi[i__ - 1];
	flag_i__ = apexinfo_1.apexf[loc_i__ - 1];
/*       << Increment the proper flag counter */
	++sortedout_1.numstart[flag_i__ + 1];
/*       << Store the location of this type of apex */
	sortedout_1.locstart[sortedout_1.numstart[flag_i__ + 1] + flag_i__ * 
		6000 + 5999] = loc_i__;
/*       < END */
/*       << Pull out the Edge I's end and its flag */
/*       << Flag_I will be equal to -1, +2 or +3 depending on */
/*       << whether the edge is dead, a 2-plane or a 3-plane */
	loc_i__ = thisplane_1.edgesi[i__ + 5999];
	flag_i__ = apexinfo_1.apexf[loc_i__ - 1];
/*       << Increment the proper flag counter */
	++sortedout_1.numend[flag_i__ + 1];
/*       << Store the location of this type of apex */
	sortedout_1.locend[sortedout_1.numend[flag_i__ + 1] + flag_i__ * 6000 
		+ 5999] = loc_i__;
/*       < EDGE */
/*       << Store the number of this type of edge */
/*       << (NB, the types here are -1 and +1) */
	sortedout_1.edge_num__[sortedout_1.numstart[3] + sortedout_1.numstart[
		4] + i_sign(&c__1, &flag_i__) * 6001 + 6001] = i__;
    }
/*     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/*     < NOTE: These lines test whether or not the number of two-plane */
/*     < START apices match the number of two-plane END apices and if */
/*     < the number of three-plane START apices match the number of */
/*     < three-plane END apices.  This was used in debugging the code. */
/*     < However, it is possible that because of numerical precision, */
/*     < two or more planes may just graze a sphere and cause one or */
/*     < more two-plane apices to be mistaken for three-plane apices or */
/*     < vise-versa.  This happens rarely but when it does the program */
/*     < should be stopped.  The solution is to perturb the coordinates */
/*     < or radii of the spheres involved by a tiny amount (< 0.01 */
/*     < Angstroms).  This is a patch but it works and the change in */
/*     < radii has almost no effect on the final answers.  One could */
/*     < replace the STOP statement with something like DEAD= .True. */
/*     < RETURN where the variable DEAD is tested in the calling */
/*     < routine and all subsequent calling routines like: If( DEAD ) */
/*     < RETURN In this way one could "bubble" up to the top routine, */
/*     < perturb the radii, and restart the calculation.  I never */
/*     < automated this but it was effectively what I did. */
/*     < */
/*     < Finally, for the same reason of numerical imprecision, a */
/*     < similar error can occur in the routine Bin_Arcs.  Perturbing */
/*     < coordinates slightly solves that problem too. */
/*     < LRD, Thu Dec 12 08:07:11 PST 1991. */
    if (sortedout_1.numstart[3] != sortedout_1.numend[3]) {
	s_wsle(&io___113);
	do_lio(&c__9, &c__1, " (+2) <> ", (ftnlen)9);
	do_lio(&c__3, &c__1, (char *)&sortedout_1.numstart[3], (ftnlen)sizeof(
		integer));
	do_lio(&c__3, &c__1, (char *)&sortedout_1.numend[3], (ftnlen)sizeof(
		integer));
	e_wsle();
	s_wsle(&io___114);
	do_lio(&c__3, &c__1, (char *)&thisplane_1.numedges, (ftnlen)sizeof(
		integer));
	e_wsle();
	s_wsle(&io___115);
	do_lio(&c__3, &c__5, (char *)&sortedout_1.numstart[0], (ftnlen)sizeof(
		integer));
	e_wsle();
	s_wsle(&io___116);
	do_lio(&c__3, &c__5, (char *)&sortedout_1.numend[0], (ftnlen)sizeof(
		integer));
	e_wsle();
	s_wsle(&io___117);
	do_lio(&c__3, &c__1, (char *)&sphereinfo_1.sphere, (ftnlen)sizeof(
		integer));
	e_wsle();
	*error_encountered__ = TRUE_;
	return 0;
    }
    if (sortedout_1.numstart[4] != sortedout_1.numend[4]) {
	s_wsle(&io___118);
	do_lio(&c__9, &c__1, " (+3) <> ", (ftnlen)9);
	do_lio(&c__3, &c__1, (char *)&sortedout_1.numstart[4], (ftnlen)sizeof(
		integer));
	do_lio(&c__3, &c__1, (char *)&sortedout_1.numend[4], (ftnlen)sizeof(
		integer));
	e_wsle();
	s_wsle(&io___119);
	do_lio(&c__3, &c__1, (char *)&thisplane_1.numedges, (ftnlen)sizeof(
		integer));
	e_wsle();
	s_wsle(&io___120);
	do_lio(&c__3, &c__5, (char *)&sortedout_1.numstart[0], (ftnlen)sizeof(
		integer));
	e_wsle();
	s_wsle(&io___121);
	do_lio(&c__3, &c__5, (char *)&sortedout_1.numend[0], (ftnlen)sizeof(
		integer));
	e_wsle();
	s_wsle(&io___122);
	do_lio(&c__3, &c__1, (char *)&sphereinfo_1.sphere, (ftnlen)sizeof(
		integer));
	e_wsle();
	*error_encountered__ = TRUE_;
	return 0;
    }
/*     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/*     < If we have any 2-plane starts at all we should form the arcs */
    localarcinfo_1.numarcs = sortedout_1.numstart[3];
    trueedges = 0;
    twoplane = localarcinfo_1.numarcs > 0;
    if (twoplane) {
	form_arcs__(&arcareas);
    }
/*     < If we have any 3-plane starts at all we should form the */
/*     < inscribed polygon on the plane face */
    if (sortedout_1.numstart[3] + sortedout_1.numstart[4] > 1) {
/*       < The number of living (non-dead) edges is equal to sum of the */
/*       < number of 2-plane edge starts and 3-plane edge starts.  We */
/*       < will use only this subset of edges as ennumerated in the */
/*       < Edge_Num(*,1) for all future work on this plane face. */
	trueedges = sortedout_1.numstart[3] + sortedout_1.numstart[4];
/*       < Form the polygon */
	form_polygon__(&polyarea, &trueedges, &twoplane);
    }
/*     < Find the reduced area for this arc-polygon. Eta_AP is defined */
/*     < as the area of the arc-polygon divided by rhoi_sq/2, where */
/*     < rhoi_sq is the radius of the inscribed circle due to this */
/*     < plane cutting the sphere. */
    *eta_ap__ = arcareas + polyarea;
/*     < Store these globally for purposes of debugging. */
    global_1.areag[thisplane_1.planei + thisplane_1.indexi * 11250 - 11251] = 
	    *eta_ap__ * .15915494309189535;
    global_1.numarcsg[thisplane_1.planei + thisplane_1.indexi * 11250 - 11251]
	     = localarcinfo_1.numarcs;
    global_1.trueedgesg[thisplane_1.planei + thisplane_1.indexi * 11250 - 
	    11251] = trueedges;
    return 0;
} /* build_ap__ */

/* ---------------------------------------------------------------------C */
/*     Form The Arcs Of The Arc-Polygon                                C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int form_arcs__(doublereal *arcareas)
{
    /* Initialized data */

    static doublereal thetastart[6000] = { 6.2831853071795862 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    static doublereal thetaend[6000];
    extern /* Subroutine */ int renumarc_(integer *, integer *);
    static integer i__, j, k;
    static doublereal ra[3], rb[3];
    static integer iend[6000];
    extern /* Subroutine */ int sort_(doublereal *, integer *, integer *);
    static integer loc_a__, loc_b__;
    extern /* Subroutine */ int angle_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal theta;
    static integer istart[6000];

/*     < PASSED */
/*     < LOCAL */
/*     < Loop through all the 2-plane apices (starts and ends) and */
/*     < find the anti-clockwise angle between it and the first */
/*     < 2-plane apex start. This angle will be between 0 and 2 Pi */
/*     < radians. */
/*     < Pull out the reference 2-plane apex start, call it A */
/*     < Apex A's angle is defined as 2Pi */
    loc_a__ = sortedout_1.locstart[18000];
/*     < Pull out A's coordinates. Note that these are vectors */
/*     < connecting the sphere center to the point A on the edge. */
    ra[0] = apexinfo_1.apexx[loc_a__ - 1];
    ra[1] = apexinfo_1.apexy[loc_a__ - 1];
    ra[2] = apexinfo_1.apexz[loc_a__ - 1];
/*     < Pull out the first 2-plane apex end, call it B */
    loc_b__ = sortedout_1.locend[18000];
/*     < Pull out B's coordinates. Note that these are vectors */
/*     < connecting the sphere center to the point B on the edge. */
    rb[0] = apexinfo_1.apexx[loc_b__ - 1];
    rb[1] = apexinfo_1.apexy[loc_b__ - 1];
    rb[2] = apexinfo_1.apexz[loc_b__ - 1];
/*     < Find the angle between A and B */
/*     < HERE IS THE PROBLEM */
    angle_(ra, rb, thetaend);
/*     < Loop through all the remaining apices and find all the angles */
/*     < relative to the first 2-plane apex start, apex A.  We will do */
/*     < the starts and the ends at the same time.  Note that we may */
/*     < have only one pair of 2-plane apices.  There may in fact be */
/*     < just one edge on this entire plane face so that these two */
/*     < apices are connected via a chord.  Whichever the case may be, */
/*     < the following code still applies (the next do-loop would not */
/*     < be executed and the sorting routine would return immediately). */
    i__1 = localarcinfo_1.numarcs;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*       < Pull out the next 2-plane apex start */
	loc_b__ = sortedout_1.locstart[i__ + 17999];
/*       < Pull out B's coordinates. Note that these are vectors */
/*       < connecting the sphere center to the point B on the edge. */
	rb[0] = apexinfo_1.apexx[loc_b__ - 1];
	rb[1] = apexinfo_1.apexy[loc_b__ - 1];
	rb[2] = apexinfo_1.apexz[loc_b__ - 1];
/*       < Find the angle between A and B */
	angle_(ra, rb, &thetastart[i__ - 1]);
/*       < Pull out the next 2-plane apex end */
	loc_b__ = sortedout_1.locend[i__ + 17999];
/*       < Pull out B's coordinates. Note that these are vectors */
/*       < connecting the sphere center to the point B on the edge. */
	rb[0] = apexinfo_1.apexx[loc_b__ - 1];
	rb[1] = apexinfo_1.apexy[loc_b__ - 1];
	rb[2] = apexinfo_1.apexz[loc_b__ - 1];
/*       < Find the angle between A and B */
	angle_(ra, rb, &thetaend[i__ - 1]);
    }
/*     < We must form and order the arcs of this plane face. The method */
/*     < we will use to form these arcs is as follows */
/*     < */
/*     <      o sort the theta-starts into ascending order (the angles */
/*     <        between each 2-plane apex start and the reference */
/*     <        2-plane apex start, A). Note that the reference */
/*     <        2-plane apex start will be the first element of this */
/*     <        ordered list since its theta-start is defined as zero. */
/*     < */
/*     <      o sort the theta-ends into ascending order (the angles */
/*     <        between each 2-plane apex end and the reference */
/*     <        2-plane apex start, A) */
/*     < */
/*     <      o match the two ordered lists; the largest theta-end */
/*     <        should form an arc with the largest theta-start */
/*     <        (i.e., the reference apex defined as having an */
/*     <        angle of 2Pi) and in general the Kth largest */
/*     <        theta-end should pair with the Kth largest */
/*     <        theta-start. */
/*     < Sort the angles for the 2-plane starts into ascending order. */
    sort_(thetastart, istart, &localarcinfo_1.numarcs);
/*     < Sort the angles for the 2-plane ends into ascending order. */
    sort_(thetaend, iend, &localarcinfo_1.numarcs);
/*     < Initialize the accumulator for the arc areas. */
    *arcareas = 0.;
/*     < The arrays IStart and IEnd contain the ordered indices for the */
/*     < 2-plane apex starts and ends, respectively.  We may now sweep */
/*     < through the indices, matching the ascending starts and ends. */
    i__1 = localarcinfo_1.numarcs;
    for (k = 1; k <= i__1; ++k) {
/*       < Retrieve the Index of the arc pair. The Kth 2-plane */
/*       < apex end forms an arc with the Kth 2-plane apex start. */
	i__ = istart[k - 1];
	j = iend[k - 1];
/*       < Store this arc as a 2-plane end and 2-plane start */
	localarcinfo_1.local_list__[k - 1] = sortedout_1.locend[j + 17999];
	localarcinfo_1.local_list__[k + 5999] = sortedout_1.locstart[i__ + 
		17999];
/*       < Find the true anti-clockwise angle for the arc in radians */
	theta = thetastart[i__ - 1] - thetaend[j - 1];
/*       < Find the area of the circle segment formed by the arc and the */
/*       < chord connecting the ends of the arc.  Note that since we are */
/*       < dealing with the true arc angles the circular segment area */
/*       < formula gives the true circular segment area.  We will sum */
/*       < these arc areas and multiply by the squared radius of the */
/*       < inscribed circle's radius at the end. */
/*       < There might be a way to get Sin of Theta from (rAxrB).n - LRD */
	*arcareas += theta - sin(theta);
/*       < Find all the arc information for use later with the festoons. */
	renumarc_(&sortedout_1.locend[j + 17999], &sortedout_1.locstart[i__ + 
		17999]);
    }
    return 0;
} /* form_arcs__ */

/* ---------------------------------------------------------------------C */
/*     Angle Between Two 2-Plane Apices                                C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int angle_(doublereal *ra, doublereal *rb, doublereal *theta)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), acos(doublereal);

    /* Local variables */
    static doublereal w, rab[3], axb[3], cos_w__;
    extern /* Subroutine */ int cross_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal rab_dollar_rab__;

/*     < PASSED */
/*     < LOCAL */
/*     < Find the difference point B and point A */
    /* Parameter adjustments */
    --rb;
    --ra;

    /* Function Body */
    rab[0] = rb[1] - ra[1];
    rab[1] = rb[2] - ra[2];
    rab[2] = rb[3] - ra[3];
/*     < Find the squared magnitude of this vector */
    rab_dollar_rab__ = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];
/*     < Find the cosine of the angle between rA and rB */
    cos_w__ = 1. - thisplane_1.curvei * rab_dollar_rab__;
/*     < Find the angle. Make sure |Cos_w| < 1 */
    if (abs(cos_w__) > 1.) {
	w = (.5 - d_sign(&c_b31, &cos_w__)) * 3.14159265358979323846;
    } else {
	w = acos(cos_w__);
    }
/*     < Determine whether w is the real anti-clockwise angle or not. */
/*     < Find the cross product of rA and rB, rA x rB */
    cross_(&ra[1], &rb[1], axb);
/*     < Find the dot product of this vector with the */
/*     < normal vector of the plane ni, (rA x rB).ni and then */
/*     < find the true anti-clockwise angle between A and B */
/*     < */
/*     <   (rA x rB).ni > 0 ==> Theta= w */
/*     <   (rA x rB).ni = 0 ==> Theta= w */
/*     <   (rA x rB).ni < 0 ==> Theta= 2Pi - w */
    d__1 = axb[0] * thisplane_1.ni[0] + axb[1] * thisplane_1.ni[1] + axb[2] * 
	    thisplane_1.ni[2];
    *theta = (w - 3.14159265358979323846) * d_sign(&c_b4, &d__1) + 
	    3.14159265358979323846;
    return 0;
} /* angle_ */

/* ---------------------------------------------------------------------C */
/*     Renumber Arcs                                                   C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int renumarc_(integer *start, integer *terminus)
{
    extern /* Subroutine */ int diangle_(integer *, integer *);

/*     < PASSED */
/*     < Increment the arc counter */
    ++arcinfo_1.arccnt;
/*     < Every arc for a sphere is made up of two two-plane apices both */
/*     < are not necessarily created by the same pair of planes.  Each */
/*     < of these two-plane apices will appear on exactly two plane */
/*     < faces. We will renumber all of this apices defining them as */
/*     < vertices of the ultimate spherical polygon.  We must perserve */
/*     < the non-uniqueness of this apices, that is, a two-plane apex */
/*     < that is an arc start on one plane face is the end of another */
/*     < arc on another plane face.  We wish to perserve this */
/*     < connection.  Therefore, as we assign new vertex numbers to */
/*     < these arc starts and termini we must check to see that they do */
/*     < not already have a new vertex number. */
/*     < START */
/*     < Have we done this vertex yet? */
    if (apexinfo_1.apexf[*start + 5999] == -1) {
/*       < We have NOT done this vertex before. */
/*       < Assign the apex a new vertex number. */
	++arcinfo_1.vertexcnt;
/*       < Copy the coordinates. */
	arcinfo_1.vertexx[arcinfo_1.vertexcnt - 1] = apexinfo_1.apexx[*start 
		- 1];
	arcinfo_1.vertexy[arcinfo_1.vertexcnt - 1] = apexinfo_1.apexy[*start 
		- 1];
	arcinfo_1.vertexz[arcinfo_1.vertexcnt - 1] = apexinfo_1.apexz[*start 
		- 1];
/*       < Mark this apex as done for next time. */
	apexinfo_1.apexf[*start + 5999] = arcinfo_1.vertexcnt;
/*       < Define this arc as starting with this vertex. */
	arcinfo_1.arc_list__[arcinfo_1.arccnt - 1] = arcinfo_1.vertexcnt;
    } else {
/*       < We have come across this apex before.  Pull out its new */
/*       < vertex number and define this arc as starting with it. */
	arcinfo_1.arc_list__[arcinfo_1.arccnt - 1] = apexinfo_1.apexf[*start 
		+ 5999];
    }
/*     < Later we will need to be able to establish the connectiveness */
/*     < of the arcs for this sphere so it will be incredibly useful to */
/*     < have a map from the vertex numbers of the arc starts to the */
/*     < arc numbers. */
    arcinfo_1.arc_map__[arcinfo_1.arc_list__[arcinfo_1.arccnt - 1] - 1] = 
	    arcinfo_1.arccnt;
/*     < TERMINUS */
/*     < Now repeat for the arc terminus. */
    if (apexinfo_1.apexf[*terminus + 5999] == -1) {
/*       < We have NOT done this vertex before. */
/*       < Assign the apex a new vertex number. */
	++arcinfo_1.vertexcnt;
/*       < Copy the coordinates. */
	arcinfo_1.vertexx[arcinfo_1.vertexcnt - 1] = apexinfo_1.apexx[*
		terminus - 1];
	arcinfo_1.vertexy[arcinfo_1.vertexcnt - 1] = apexinfo_1.apexy[*
		terminus - 1];
	arcinfo_1.vertexz[arcinfo_1.vertexcnt - 1] = apexinfo_1.apexz[*
		terminus - 1];
/*       < Mark this apex as done for next time. */
	apexinfo_1.apexf[*terminus + 5999] = arcinfo_1.vertexcnt;
/*       < Define this arc as terminating with this vertex. */
	arcinfo_1.arc_list__[arcinfo_1.arccnt + 299] = arcinfo_1.vertexcnt;
    } else {
/*       < We have come across this apex before.  Pull out its new */
/*       < vertex number and define this arc as terminating with it. */
	arcinfo_1.arc_list__[arcinfo_1.arccnt + 299] = apexinfo_1.apexf[*
		terminus + 5999];
    }
/*     < Now that we have the arc start and terminus we need to find */
/*     < the geodesic quarter-angle, the diangle area, and the geodesic */
/*     < unit normal for this arc. */
    diangle_(&arcinfo_1.arc_list__[arcinfo_1.arccnt - 1], &
	    arcinfo_1.arc_list__[arcinfo_1.arccnt + 299]);
    return 0;
} /* renumarc_ */

/* ---------------------------------------------------------------------C */
/*     Diangle Information                                             C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int diangle_(integer *start, integer *terminus)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *), acos(
	    doublereal);

    /* Local variables */
    static doublereal cos_theta__, xs[3], xt[3], sxt[3], rmag, cosx, denom;
    extern /* Subroutine */ int cross_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal cos_chi__, cot_chi__;

/*     < PASSED */
/*     < LOCAL */
/*     < Pull out the coordinates of the vector connecting the sphere */
/*     < center to the start of this arc, ArcCnt. */
    xs[0] = arcinfo_1.vertexx[*start - 1];
    xs[1] = arcinfo_1.vertexy[*start - 1];
    xs[2] = arcinfo_1.vertexz[*start - 1];
/*     < Pull out the coordinates of the vector connecting the sphere */
/*     < center to the terminus of this arc, ArcCnt. */
    xt[0] = arcinfo_1.vertexx[*terminus - 1];
    xt[1] = arcinfo_1.vertexy[*terminus - 1];
    xt[2] = arcinfo_1.vertexz[*terminus - 1];
/*     < Find the cosine of the angle between these vectors.  Both */
/*     < vectors are of magnitude R, the radius of the sphere.  The dot */
/*     < product is defined as, */
/*     < */
/*     <        xA.xB= |xA||xB| Cos(Theta) */
/*     <             = (R*R) Cos(Theta) */
/*     < */
/*     < Therefore, Cos(Theta)= xA.xB/(R*R) */
    cos_theta__ = (xs[0] * xt[0] + xs[1] * xt[1] + xs[2] * xt[2]) * 
	    sphereinfo_1.o_rsq__;
/*     < Find out if Cos(Theta) is close is -1. */
    if ((d__1 = cos_theta__ + 1., abs(d__1)) < 1e-8) {
/*       < The geodesic angle between the arc start and end is Pi */
/*       < radians.  This means that the tangent of the geodesic */
/*       < quarter-angle is unity, the normal vector of the plane is */
/*       < the geodesic unit normal vector, and there is no diangle. */
	arcinfo_1.tangeo[arcinfo_1.arccnt - 1] = 1.;
	arcinfo_1.ngx[arcinfo_1.arccnt - 1] = thisplane_1.ni[0];
	arcinfo_1.ngy[arcinfo_1.arccnt - 1] = thisplane_1.ni[1];
	arcinfo_1.ngz[arcinfo_1.arccnt - 1] = thisplane_1.ni[2];
/*       < Return now. */
	return 0;
    }
/*     < Find the square of the tangent of the geodesic quarter-angle. */
/*     < The square tangent of the quarter angle is related to the */
/*     < cosine and sine half-angle via, */
/*     < */
/*     <   Tan^2(Theta/4)= Sin^2(Theta/2)/(1 + Cos(Theta/2))^2 */
/*     < */
/*     < and the half-angles are related the whole angles via, */
/*     < */
/*     <   Sin^2(Theta/2)= (1 - Cos(Theta))/2 */
/*     <     Cos(Theta/2)= Sqrt(1 + Cos(Theta))/Sqrt(2) */
/*     < */
/*     < thus, after multiplying numerator and denominator by 2, we have */
/*     < */
/*     <                           1 - Cos(Theta) */
/*     <   Tan^2(Theta/4)= --------------------------------- */
/*     <                   [Sqrt(2) + Sqrt(1 + Cos(Theta)]^2 */
/*     < */
    denom = sqrt(cos_theta__ + 1.) + 1.4142135623730950488;
    arcinfo_1.tangeo[arcinfo_1.arccnt - 1] = (1. - cos_theta__) / (denom * 
	    denom);
/*     < Find the cross-product of these two vectors. */
    cross_(xs, xt, sxt);
/*     < Find the reciprocal magnitude of the cross product. */
    rmag = 1. / sqrt(sxt[0] * sxt[0] + sxt[1] * sxt[1] + sxt[2] * sxt[2]);
/*     < Find the geodesic unit normal vector, nG, defined as the unit */
/*     < vector normal the the plane that xA and xB lie on. */
    arcinfo_1.ngx[arcinfo_1.arccnt - 1] = sxt[0] * rmag;
    arcinfo_1.ngy[arcinfo_1.arccnt - 1] = sxt[1] * rmag;
    arcinfo_1.ngz[arcinfo_1.arccnt - 1] = sxt[2] * rmag;
/*     < Find the cosine of the diangle angle.  This is defined by the */
/*     < dot product of the geodesic unit normal, nG, with the normal */
/*     < of the plane, ni, */
/*     < */
/*     < nG.ni= |nG|*|ni| Cos(Chi)= Cos(Chi) */
    cos_chi__ = arcinfo_1.ngx[arcinfo_1.arccnt - 1] * thisplane_1.ni[0] + 
	    arcinfo_1.ngy[arcinfo_1.arccnt - 1] * thisplane_1.ni[1] + 
	    arcinfo_1.ngz[arcinfo_1.arccnt - 1] * thisplane_1.ni[2];
/*     < Find out if Cos(Chi) is close to +/-1. */
    if ((d__1 = 1. - abs(cos_chi__), abs(d__1)) < 1e-8) {
/*       < The great plane passing through points A and B is the plane */
/*       < itself. Note that the formula for the diangle would have */
/*       < |vi| times the cotangent of Chi. Chi in this case is zero */
/*       < (or 2Pi) and so Cot(Chi) is infinity.  Therefore, the */
/*       < diangle is zero in this case.  Return now. */
	return 0;
    }
/*     < Find the cotangent of this angle.  Note that Cot(Chi)'s sign */
/*     < over 0 to Pi is do strictly to Cos(Chi) thus the sign is */
/*     < perserved in the formula below. */
    cot_chi__ = cos_chi__ / sqrt(1. - cos_chi__ * cos_chi__);
/*     < Find the area of the diangle via Gibson's formulae [eqn. A10 */
/*     < of appendix of Gibson and Scheraga, Molecular Physics, 62 */
/*     < (1987) 1247-1265].  We have to consider only type one and type */
/*     < two diangles (z is always positive with Chi being both acute */
/*     < or obtuse). AD_2Rsq is the area of diangle over twice the */
/*     < sphere radius squared.  The quantities are defined as follows */
/*     < */
/*     <    R_rhoi := R/rhoi= R/Sqrt(R^2 - z^2) */
/*     <      vi_R := |vi|/R */
/*     <   vi_rhoi := |vi|/rhoi= |vi|/Sqrt(R^2 - z^2) */
/*     < */
/*     < where z is Gibson's z defined as the distance between the */
/*     < plane and the geodesic plane and in our terminology is |vi| */
/*     < and hence Sqrt(R^2 - z^2) is really rhoi the radius of the */
/*     < inscribed circle.  We break up the formula for the diangle */
/*     < area into two contributions.  We sum these for all the arcs on */
/*     < this plane and then we will multiply the second part by */
/*     < -|vk|/R and add it to the first part.  Eta is equal to the */
/*     < area of the diangle over twice the sphere radius squared. */
    cosx = thisplane_1.r_rhoi__ * cos_chi__;
/* Computing MIN */
    d__2 = 1., d__3 = abs(cosx);
    d__1 = min(d__2,d__3);
    cosx = d_sign(&d__1, &cosx);
    diangleinfo_1.etad[0] += acos(cosx);
    cosx = thisplane_1.vi_rhoi__ * cot_chi__;
/* Computing MIN */
    d__2 = 1., d__3 = abs(cosx);
    d__1 = min(d__2,d__3);
    cosx = d_sign(&d__1, &cosx);
    diangleinfo_1.etad[1] += acos(cosx);
    return 0;
} /* diangle_ */

/* ---------------------------------------------------------------------C */
/*     Form Polygon                                                    C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int form_polygon__(doublereal *area, integer *trueedges, 
	logical *twoplane)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, a1tripnum, a1;
    extern /* Subroutine */ int accumulate_(integer *, integer *);
    static integer eend, edgei, estart;

/*     < PASSED */
/*     < LOCAL */
/*     < The basic algorithm used to find the area of the inscribed */
/*     < polygon will be the following */
/*     < */
/*     < If there are no two-plane apices then */
/*     < */
/*     <     o pick a reference 3-plane apex START, call it A1 */
/*     <     o loop through all the REMAINING edges on this plane face */
/*     <     o find the area of the triangle formed by the reference */
/*     <       apex A1 and the start and end of the edge */
/*     <       (A1 Estart Eend) */
/*     <     o sum all of these areas */
/*     < */
/*     <   Since there are no arcs and we simply have an inscribed */
/*     <   polygon on this plane face. We loop through only the */
/*     <   REMAINING edges of the plane face comparing the triple plane */
/*     <   number of the end of each edge with A1's triple plane number. */
/*     <   If they match then we skip this edge. */
/*     < */
/*     < Else if there are two-plane apices then */
/*     < */
/*     <     o pick a reference 2-plane apex END, call it A1 */
/*     <     o loop through ALL the edges on this plane face */
/*     <     o find the area of the triangle formed by the reference */
/*     <       apex A1 and the start and end of the edge */
/*     <       (A1 Estart Eend) */
/*     <     o sum all of these areas */
/*     < */
/*     <   We loop through ALL edges, however, we check the end of each */
/*     <   edge to see whether the edge end is a 2-plane apex and if so */
/*     <   whether it is A1 or not.  If it is A1 we skip this edge. */
/*     < */
/*     <   At the end of the sweep through all edges the total area will */
/*     <   not contain the areas of the triangles formed by the pseudo */
/*     <   edges connecting the starts and ends of the arcs.  It will, */
/*     <   however, contain the area for the arc that A1 starts.  We */
/*     <   will need to do the following */
/*     < */
/*     <     o loop through all the remaining arcs on this plane face */
/*     <     o find the area of the triangle formed by the reference */
/*     <       apex A1 and the start and end of each arc */
/*     <       (A1 Astart Aend) */
/*     <     o sum all of these areas */
/*     < Initialize the a x b vector */
    foraccum_1.axb[0] = 0.;
    foraccum_1.axb[1] = 0.;
    foraccum_1.axb[2] = 0.;
/*     < Find out if our reference apex will be 2-plane or 3-plane */
    if (*twoplane) {
/*       < Pull out the 2-plane apex end that forms the start of the */
/*       < first arc we found above. */
	a1 = localarcinfo_1.local_list__[0];
/*       < Pull out the coordinates of our reference apex */
	foraccum_1.ra1[0] = apexinfo_1.apexx[a1 - 1];
	foraccum_1.ra1[1] = apexinfo_1.apexy[a1 - 1];
	foraccum_1.ra1[2] = apexinfo_1.apexz[a1 - 1];
/*       < Loop through ALL the edges */
	i__1 = *trueedges;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*         < Pull out living (not dead) edge number */
	    edgei = sortedout_1.edge_num__[i__ + 12002];
/*         < Pull out the end of edge I */
	    eend = thisplane_1.edgesi[edgei + 5999];
/*         < Is it our reference apex? */
/*         < Note by definition this edge is alive. */
	    if (eend != a1) {
/*           < The end of edge I, Eend is not a 2-plane and is not our */
/*           < reference apex A1. We must form the area of the triangle */
/*           < (A1 Estart Eend) */
/*           < Pull out the start of edge I */
		estart = thisplane_1.edgesi[edgei - 1];
/*           < Accumulate the cross-product of the vector from A1 to the */
/*           < edge start and the vector from A1 to the edge end. */
		accumulate_(&estart, &eend);
	    }
	}
/*       < Loop through the remaining arcs of this plane face. */
/*       < These will be our pseudo-edges. */
	i__1 = localarcinfo_1.numarcs;
	for (i__ = 2; i__ <= i__1; ++i__) {
/*         < Accumulate the cross-product of the vector from A1 to the */
/*         < arc start, Local_List(I,1), and the vector from A1 */
/*         < to the arc end, Local_List(I,2). */
	    accumulate_(&localarcinfo_1.local_list__[i__ - 1], &
		    localarcinfo_1.local_list__[i__ + 5999]);
	}
    } else {
/*       < Pull out the 3-plane apex start that forms the start of the */
/*       < first living (non-dead) edge in the edge list.  Note that we */
/*       < know that its not dead because we are selecting it from the */
/*       < Edge_Num list. */
	a1 = thisplane_1.edgesi[sortedout_1.edge_num__[12003] - 1];
/*       < Pull out the coordinates of our reference apex */
	foraccum_1.ra1[0] = apexinfo_1.apexx[a1 - 1];
	foraccum_1.ra1[1] = apexinfo_1.apexy[a1 - 1];
	foraccum_1.ra1[2] = apexinfo_1.apexz[a1 - 1];
/*       < Pull out the triple plane number for our reference apex */
/*       < We know that ApexF(A1,1)= +3 because we have no 2-plane */
/*       < apices and this edge is the first living (non-dead) edge */
/*       < (i.e., we know that ApexF(A1,1) is not -1). */
	a1tripnum = apexinfo_1.apexf[a1 + 5999];
/*       < Loop through the REMAINING living (non-dead) edges */
	i__1 = *trueedges;
	for (i__ = 2; i__ <= i__1; ++i__) {
/*         < Pull out living (not dead) edge number */
	    edgei = sortedout_1.edge_num__[i__ + 12002];
/*         < Pull out the end of edge I */
	    eend = thisplane_1.edgesi[edgei + 5999];
/*         < Is it our reference apex? */
/*         < Note we know that this edge is alive (non-dead) */
	    if (apexinfo_1.apexf[eend + 5999] != a1tripnum) {
/*           < The end of edge I, Eend does not have the same triple */
/*           < plane number as the triple plane number of our reference */
/*           < apex A1.  We must form the area of the triangle */
/*           < (A1 Estart Eend) */
/*           < Pull out the start of edge I */
		estart = thisplane_1.edgesi[edgei - 1];
/*           < Accumulate the cross-product of the vector from A1 to the */
/*           < edge start and the vector from A1 to the edge end. */
		accumulate_(&estart, &eend);
	    }
	}
    }
/*     < Dot the accumulator vector, (a x b), with this plane's unit */
/*     < normal vector ni.  This will give the sum of the areas of all */
/*     < the parallelograms on this plane face.  The total area of the */
/*     < inscribed polygon is the sum of the areas of all the triangles */
/*     < formed.  Thus the total area is one-half the area of all the */
/*     < parallelograms formed. */
    *area = thisplane_1.o_rhoi_sq__ * (foraccum_1.axb[0] * thisplane_1.ni[0] 
	    + foraccum_1.axb[1] * thisplane_1.ni[1] + foraccum_1.axb[2] * 
	    thisplane_1.ni[2]);
    return 0;
} /* form_polygon__ */

/* ---------------------------------------------------------------------C */
/*     Accumulate Cross Products                                       C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int accumulate_(integer *estart, integer *eend)
{
    static doublereal ra[3], rb[3], sxe[3];
    extern /* Subroutine */ int cross_(doublereal *, doublereal *, doublereal 
	    *);

/*     < PASSED */
/*     < LOCAL */
/*     < Find the vector connecting A1 to Estart, call it rA */
    ra[0] = apexinfo_1.apexx[*estart - 1] - foraccum_1.ra1[0];
    ra[1] = apexinfo_1.apexy[*estart - 1] - foraccum_1.ra1[1];
    ra[2] = apexinfo_1.apexz[*estart - 1] - foraccum_1.ra1[2];
/*     < Find the vector connecting A1 to Eend, call it rB */
    rb[0] = apexinfo_1.apexx[*eend - 1] - foraccum_1.ra1[0];
    rb[1] = apexinfo_1.apexy[*eend - 1] - foraccum_1.ra1[1];
    rb[2] = apexinfo_1.apexz[*eend - 1] - foraccum_1.ra1[2];
/*     < Find the cross-product of these two vectors, call it sxe */
/*     < (for 'start-vector cross-product end-vector') */
    cross_(ra, rb, sxe);
/*     < The area of the parallelogram formed by these two vectors, rA */
/*     < and rB, is given by the dot of the cross-product with a unit */
/*     < normal vector pointing in the same direction as the */
/*     < cross-product.  We have formed the cross-product in a */
/*     < right-handed sense and crossed the start vector into the end */
/*     < vector and so the cross-product is normal to the plane face */
/*     < pointing away from the cutout.  Therefore the plane's unit */
/*     < normal vector ni points in the same direction as the cross */
/*     < product, sxe.  Furthermore, this same unit normal vector can be */
/*     < used for all the cross-products of this plane face and we will */
/*     < save time by performing the dot product at the end of the loop */
/*     < on the sum of the cross-products.  Mathematically this means */
/*     < */
/*     <    Area of Parallelograms= Sum_{j=1,n} [(sj x ej).ni] */
/*     <                          = [Sum_{j=1,n} (sj x ej)].ni */
/*     <                          = (a x b).ni */
/*     < */
/*     < where (a x b)= Sum_{j=1,n} (sj x ej).  Note that axb was */
/*     < initialized earlier. */
    foraccum_1.axb[0] += sxe[0];
    foraccum_1.axb[1] += sxe[1];
    foraccum_1.axb[2] += sxe[2];
    return 0;
} /* accumulate_ */

/* ---------------------------------------------------------------------C */
/*     Festoonery                                                      C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int festoonery_(doublereal *area_sp__, logical *
	error_encountered__)
{
    extern /* Subroutine */ int bin_arcs__(logical *), spherepoly_(doublereal 
	    *);

/*     < PASSED */
    *error_encountered__ = FALSE_;
/*     < How many arcs do we have for this sphere?  Note that none of */
/*     < them will be spherical caps so that an arc is necessariy */
/*     < connected to another arc. */
    if (arcinfo_1.arccnt <= 2) {
/*       < Initialize the area of the spherical polygon. */
	*area_sp__ = 0.;
    } else {
/*       < There are three or more arcs and so we must bin them into */
/*       < festoons.  The spherical polygon may be just a simple */
/*       < spherical triangle.  Note that four arcs could form two */
/*       < linear spherical polygons and not necessarily a spherical */
/*       < rectangle. */
	bin_arcs__(error_encountered__);
	if (*error_encountered__) {
	    return 0;
	}
/*       < Find the area bordered by each spherical polygon. */
	spherepoly_(area_sp__);
    }
    return 0;
} /* festoonery_ */

/* ---------------------------------------------------------------------C */
/*     Bin Arcs Into Festoons                                          C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int bin_arcs__(logical *error_encountered__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer usedarcs, arc, cnt, arcref;
    static logical notused[300];

    /* Fortran I/O blocks */
    static cilist io___165 = { 0, 0, 0, 0, 0 };


/*     < LOCAL */
    *error_encountered__ = FALSE_;
/*     < Initialize some variables. */
    usedarcs = 0;
    arcref = 0;
    festooninfo_1.festooncnt = 0;
    i__1 = arcinfo_1.arccnt;
    for (arc = 1; arc <= i__1; ++arc) {
	notused[arc - 1] = TRUE_;
    }
/*     < We have to put the arcs into festoon bins.  It takes three or */
/*     < more arcs to make a festoon.  We have to analyze the */
/*     < connectivity of the arcs for this sphere to determine the */
/*     < festoons. */
    while(usedarcs < arcinfo_1.arccnt) {
/*       < Define the reference arc starting this festoon. */
	++arcref;
/*       < Has this been used? */
	if (notused[arcref - 1]) {
/*         < Increment the festoon count. */
	    ++festooninfo_1.festooncnt;
/*         < Initialize the arc count for this festoon. */
	    cnt = 1;
/*         < Store the first arc. */
	    festooninfo_1.festoon[cnt + festooninfo_1.festooncnt * 300 - 301] 
		    = arcref;
/*         < Find the arc that reference arc's terminus starts. */
	    arc = arcinfo_1.arc_map__[arcinfo_1.arc_list__[arcref + 299] - 1];
/*         < Tag this arc as used. */
	    notused[arc - 1] = FALSE_;
/*         < Connect the arcs until the reference arc appears again. */
	    while(arc != arcref) {
/*           < Increment the arc count for this festoon. */
		++cnt;
/*           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/*           < Test to see if we are endlessly looping. */
/*           < See note in routine Build_AP concerning two-plane versus */
/*           < three-plane apices. */
/*           < LRD - Saturday, February 08, 1992, 10:14:54 PST */
		if (cnt > 300) {
		    s_wsle(&io___165);
		    do_lio(&c__9, &c__1, " Perturb some coordinates ", (
			    ftnlen)26);
		    e_wsle();
		    *error_encountered__ = TRUE_;
		    return 0;
		}
/*           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/*           < Store the arc. */
		festooninfo_1.festoon[cnt + festooninfo_1.festooncnt * 300 - 
			301] = arc;
/*           < Find the arc that the last arc's terminus starts. */
		arc = arcinfo_1.arc_map__[arcinfo_1.arc_list__[arc + 299] - 1]
			;
/*           < Tag this arc as used. */
		notused[arc - 1] = FALSE_;
	    }
/*         < Increment the number of arcs used. */
	    usedarcs += cnt;
/*         < Store the number of arcs belonging to this festoon. */
	    festooninfo_1.festarcs[festooninfo_1.festooncnt - 1] = cnt;
	}
    }
    return 0;
} /* bin_arcs__ */

/* ---------------------------------------------------------------------C */
/*     Find The Area Of All The Spherical Polygons                     C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int spherepoly_(doublereal *area_sp__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal), acos(
	    doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static doublereal x[3], cos_omega__, cos_theta__, ng[3], xt[3];
    static integer arc, cnt;
    static doublereal ng_dollar_x__, rmag;
    static integer fest;
    extern /* Subroutine */ int lune_(doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *);
    static integer vert;
    static doublereal denom, alune;
    extern /* Subroutine */ int euler_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), cross_(doublereal *, doublereal *, doublereal *)
	    ;
    static doublereal tandia, areaet;
    static integer arcref;
    static doublereal tanarc, tanold;
    static integer arclast;
    static logical notdone;
    static doublereal tanlast;
    static integer numarcs;

    /* Fortran I/O blocks */
    static cilist io___186 = { 0, 0, 0, 0, 0 };
    static cilist io___189 = { 0, 0, 0, 0, 0 };


/*     < PASSED */
/*     < LOCAL */
/*     < The method for finding the area of the spherical polygon is */
/*     < analogous to finding the area of a polygon on a plane. There */
/*     < we pick a vertex of the polygon and find the area of triangles */
/*     < formed by each adjacent pair of vertices of the polygon and */
/*     < the reference vertex.  So for a polygon with n sides the area */
/*     < is defined as a sum of n-2 triangles.  For a spherical polygon */
/*     < we do the same thing except now the triangles are spherical */
/*     < triangles and in particular Euler triangles. */
/*     < Initialize the area of the spherical polygon. */
    *area_sp__ = 0.;
/*     < Single Loop over all festoons (i.e., a single loop over all */
/*     < the spherical polygons on this sphere). */
    i__1 = festooninfo_1.festooncnt;
    for (fest = 1; fest <= i__1; ++fest) {
/*       < Pull out the number of arcs for this festoon. */
	numarcs = festooninfo_1.festarcs[fest - 1];
/*       < The festoon has to contain at least three arcs in order to */
/*       < be non-trivial. If it contains two arcs then the spherical */
/*       < polygon has degenerated down to a line. It can not contain */
/*       < only one arc because this would be a spherical segment and */
/*       < would have be considered earlier. */
	if (numarcs > 2) {
/*         < Pull out the arc number of the first arc for this festoon. */
/*         < This first arc will be our reference arc. */
	    arcref = festooninfo_1.festoon[fest * 300 - 300];
/*         < Pull out the square tangent of the geodesic quarter-angle */
/*         < for this reference arc. We will disguise it as diagonal's. */
	    tanold = arcinfo_1.tangeo[arcref - 1];
/*         < Pull out the vertex number of the vertex starting this */
/*         < reference arc. */
	    vert = arcinfo_1.arc_list__[arcref - 1];
/*         < Pull out the vector connecting the sphere center to the */
/*         < vertex starting this reference arc, call this vector x. */
	    x[0] = arcinfo_1.vertexx[vert - 1];
	    x[1] = arcinfo_1.vertexy[vert - 1];
	    x[2] = arcinfo_1.vertexz[vert - 1];
/*         < Do a single loop over all the remaining arcs of this */
/*         < festoon but two. Note that if the total number of arcs for */
/*         < this festoon is three then this loop will not be executed. */
	    cnt = 1;
	    notdone = TRUE_;
	    while(cnt < numarcs - 2) {
/*           < Increment counter. */
		++cnt;
/*           < Pull out the arc number of arc Cnt. */
		arc = festooninfo_1.festoon[cnt + fest * 300 - 301];
/*           < Pull out the square tangent of the geodesic quarter-angle */
/*           < for the present arc Cnt. */
		tanarc = arcinfo_1.tangeo[arc - 1];
/*           < Pull out the number of the terminating vertex of arc Cnt. */
		vert = arcinfo_1.arc_list__[arc + 299];
/*           < Pull out the vector connecting the sphere center to the */
/*           < end of arc Arc.  Call it xT; T is for terminus. */
		xt[0] = arcinfo_1.vertexx[vert - 1];
		xt[1] = arcinfo_1.vertexy[vert - 1];
		xt[2] = arcinfo_1.vertexz[vert - 1];
/*           < Find the square tangent of the geodesic quarter-angle */
/*           < between the vector connecting the sphere center to the */
/*           < terminus of this arc (vertex Vert) and the vector for */
/*           < the reference arc, x.  This is the diagonal of the Euler */
/*           < triangle to be considered.  First, find the cosine of */
/*           < the angle between these vectors.  Both vectors are of */
/*           < magnitude R, the radius of the sphere. Therefore, */
/*           < Cos(Theta)= x.xT/(R*R) */
		cos_theta__ = (x[0] * xt[0] + x[1] * xt[1] + x[2] * xt[2]) * 
			sphereinfo_1.o_rsq__;
/*           < Is Cos(Theta) close to -1? */
		if ((d__1 = cos_theta__ + 1., abs(d__1)) < 1e-8) {
/*             < The diagonal of this triangle is exactly Pi radians */
/*             < long.  This means that the lengths of the other two */
/*             < sides (one or both being arcs of the spherical */
/*             < polygon) sum to Pi radians also.  Furthermore, the */
/*             < next triangle will also contain this diagonal so the */
/*             < these two triangles have reduced to a spherical lune. */
/*             < Pull out the unit normal geodesic vector for this arc */
/*             < and dot it with the unit normal geodesic vector for */
/*             < the next arc that will share this diagonal. */
		    lune_(xt, &ng_dollar_x__, &alune, &cnt, &fest, &arc);
/*             < Sum the area adding or subtracting the spherical */
/*             < lune's area when appropriate. */
		    *area_sp__ -= d_sign(&alune, &ng_dollar_x__);
/*             < The variable Cnt has been incremented in Lune because */
/*             < we have used two arcs and the reference vertex to form */
/*             < the lune.  If Cnt is still less or equal to than */
/*             < (NumArcs - 2) then we need to get the geodesic angle */
/*             < for the diagonal.  Otherwise we are done with this */
/*             < spherical polygon. */
		    if (cnt <= numarcs - 2) {
/*               < In order to find the area of the next spherical */
/*               < triangle, which may be the last, we need the square */
/*               < of the tangent of the geodesic quarter-angle for the */
/*               < diagonal formed by the reference vertex and the */
/*               < terminus of the second arc.  Actually, since the */
/*               < diagonal is of length Pi, the second arc and this */
/*               < diagonal have to sum to a length of Pi.  This could */
/*               < be used as a short cut. Find the Cos(Theta). */
			cos_theta__ = (x[0] * xt[0] + x[1] * xt[1] + x[2] * 
				xt[2]) * sphereinfo_1.o_rsq__;
/*               < The square tangent of the quarter angle is related to */
/*               < the cosine angle, assign it directly to TanOld. */
			denom = sqrt(cos_theta__ + 1.) + 
				1.4142135623730950488;
			tanold = (1. - cos_theta__) / (denom * denom);
		    } else {
/*               < This spherical polygon is complete.  We will need to */
/*               < skip over the last triangle normally done outside */
/*               < this loop because we have already incorporated it */
/*               < into the lune. */
			notdone = FALSE_;
		    }
		} else {
/*             < The square tangent of the quarter angle is related to */
/*             < the cosine angle via, */
/*             < */
/*             <                            1 - Cos(Theta) */
/*             <   Tan^2(Theta/4)= ---------------------------------- */
/*             <                   [Sqrt(2) + Sqrt(1 + Cos(Theta))]^2 */
/*             < */
		    denom = sqrt(cos_theta__ + 1.) + 1.4142135623730950488;
		    tandia = (1. - cos_theta__) / (denom * denom);
/*             < Find out if the square of the tangent of the geodesic */
/*             < quarter-angle of the arc is unity. */
		    if ((d__1 = 1. - tanarc, abs(d__1)) < 1e-8) {
/*               < If so then the arc is of length Pi and the other two */
/*               < sides sum to Pi also. This means that the triangle */
/*               < reduces to a lune with the arc forming one side. The */
/*               < other side is formed by the diagonals going from the */
/*               < arc end to the arc start through the reference */
/*               < vertex. This is also a geodesic plane. We must find */
/*               < the angle between these two planes. First find the */
/*               < geodesic unit normal vector for the geodesic going */
/*               < from the end of the arc Arc through the reference */
/*               < vertex to the start of the arc Arc. */
/*               < Find the cross product of xT into x the reference */
/*               < vertex. Call it nG= (xT) x (x). */
			cross_(xt, x, ng);
/*               < Find the magnitude of nG */
			rmag = 1. / sqrt(ng[0] * ng[0] + ng[1] * ng[1] + ng[2]
				 * ng[2]);
/*               < Find the dot product between this geodesic normal */
/*               < and the geodesic unit normal for arc Arc.  Multiply */
/*               < by the magnitude RMag. This then is the cosine of */
/*               < the complement of the angle between the two geodesic */
/*               < planes.  The angle we wish, call it alpha, is Pi/2 */
/*               < minus the angle between the two normal vectors, */
/*               < cos(w)= -cos(alpha) */
			cos_omega__ = (arcinfo_1.ngx[arc - 1] * ng[0] + 
				arcinfo_1.ngy[arc - 1] * ng[1] + 
				arcinfo_1.ngz[arc - 1] * ng[2]) * rmag;
/*               < The area of a lune is given as A= 2R^2 Alpha. Thus, */
/*               < the area of the lune divided by 4R^2 is Alpha/2. We */
/*               < will assign this to the variable AreaET since a lune */
/*               < is simply a special case of an Euler triangle. */
			areaet = acos(-cos_omega__) * .5;
			s_wsle(&io___186);
			do_lio(&c__9, &c__1, " Lune: Side ", (ftnlen)12);
			e_wsle();
		    } else {
/*               < Find the area of this Euler triangle. The area is */
/*               < reduced by 4 Pi R^2. */
			euler_(&tanold, &tanarc, &tandia, &areaet);
		    }
/*             < Find the directionality of this triangle.  The */
/*             < directionality of the triangle represents its location */
/*             < relative to the cutout.  If it is anti-clockwise then */
/*             < the triangle is out of the cutout and so the area of */
/*             < this area should be subtracted from that of a sphere. */
/*             < Otherwise the area should be added to find the area of */
/*             < the cutout.  We will establish directionality by */
/*             < considering the three points on the surface of the */
/*             < sphere making up the triangle and the sphere center as */
/*             < a tetrahedron.  The directionality of the base of this */
/*             < tetrahedron is found by finding the triple product */
/*             < (rAxrB).rC where rA is the vector connecting the */
/*             < sphere center to the point A on the sphere surface. */
/*             < If B and C are both on an arc then the unit geodesic */
/*             < normal, nG, is equal to rBxrC/|rBxrC| and we need */
/*             < merely dot this with the third point, rA.  In this */
/*             < case rA is the reference vertex of the reference arc. */
/*             < So find nG(Arc).x and if this is positive then it */
/*             < means that the triangle is directed in an */
/*             < anti-clockwise manner and we must subtract this area. */
		    ng_dollar_x__ = arcinfo_1.ngx[arc - 1] * x[0] + 
			    arcinfo_1.ngy[arc - 1] * x[1] + arcinfo_1.ngz[arc 
			    - 1] * x[2];
/*             < Accumulate the area of the triangle to the area of the */
/*             < spherical polygon adding or subtracting the when */
/*             < appropriate. */
		    *area_sp__ -= d_sign(&areaet, &ng_dollar_x__);
/*             < Define the square tangent of the diagonal as the old */
/*             < square tangent. */
		    tanold = tandia;
		}
	    }
/*         < Check if a lune has already used the last two arcs of the */
/*         < spherical polygon. */
	    if (notdone) {
/*           < At this point only the last two arcs of the spherical */
/*           < polygon have not been used yet.  Pull out the number of */
/*           < the penultimate arc. */
		arc = festooninfo_1.festoon[numarcs - 1 + fest * 300 - 301];
/*           < Pull out the square tangent of the geodesic quarter-angle */
/*           < for the penultimate arc. */
		tanarc = arcinfo_1.tangeo[arc - 1];
/*           < Pull out the number of the final arc. */
		arclast = festooninfo_1.festoon[numarcs + fest * 300 - 301];
/*           < Pull out the square tangent of the geodesic quarter-angle */
/*           < for the final arc. */
		tanlast = arcinfo_1.tangeo[arclast - 1];
/*           < Find out if the square of the tangent of the geodesic */
/*           < quarter-angle of the penultimate arc is unity. */
		if ((d__1 = 1. - tanarc, abs(d__1)) < 1e-8 || (d__2 = 1. - 
			tanlast, abs(d__2)) < 1e-8) {
/*             < Either the last or the penultimate arc is of length */
/*             < Pi.  This means that the reference vertex is on the */
/*             < geodesic plane of the other arc.  Thus the Euler */
/*             < triangle has reduced to a lune. The angle of the lune */
/*             < is found by taking the dot product of the unit normal */
/*             < geodesic vectors of each plane and then adding or */
/*             < subtracting Pi/2. */
		    cos_omega__ = arcinfo_1.ngx[arc - 1] * arcinfo_1.ngx[
			    arclast - 1] + arcinfo_1.ngy[arc - 1] * 
			    arcinfo_1.ngy[arclast - 1] + arcinfo_1.ngz[arc - 
			    1] * arcinfo_1.ngz[arclast - 1];
/*             < The area of a lune is given as A= 2R^2 Alpha. Thus, */
/*             < the area of the lune divided by 4R^2 is Alpha/2. We */
/*             < will assign this to the variable AreaET since a lune */
/*             < is simply a special case of an Euler triangle. */
		    areaet = acos(-cos_omega__) * .5;
		    s_wsle(&io___189);
		    do_lio(&c__9, &c__1, " Lune: Side ", (ftnlen)12);
		    e_wsle();
		} else {
/*             < Find the area of this last Euler triangle for this */
/*             < spherical polygon. */
		    euler_(&tanold, &tanarc, &tanlast, &areaet);
		}
/*           < Find the directionality of this triangle. */
		ng_dollar_x__ = arcinfo_1.ngx[arc - 1] * x[0] + arcinfo_1.ngy[
			arc - 1] * x[1] + arcinfo_1.ngz[arc - 1] * x[2];
/*           < Sum the area adding or subtracting the Euler triangle's */
/*           < area when appropriate. */
		*area_sp__ -= d_sign(&areaet, &ng_dollar_x__);
	    }
	}
    }
/*     < Divide the reduced spherical polygon area by Pi making it */
/*     < reduced by 4 Pi R^2 the surface area of a sphere. */
    *area_sp__ *= .31830988618379069;
    return 0;
} /* spherepoly_ */

/* ---------------------------------------------------------------------C */
/*     Find The Area Of An Euler Triangle                              C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int euler_(doublereal *asq, doublereal *bsq, doublereal *csq,
	 doublereal *area)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal);

    /* Local variables */
    static doublereal sum, x1sq, x2sq, x3sq, absq, group, sqrtarg;

/*     < PASSED */
/*     < LOCAL */
/*     < Define some useful quantities. */
    sum = *asq + *bsq;
    absq = *asq * *bsq;
    group = absq + 1.;
/*     < Find the three squares. */
    x1sq = sum - *csq * group;
    x1sq *= x1sq;
    x2sq = group - *csq * sum;
    x2sq *= x2sq;
    x3sq = *csq + 1.;
    x3sq = absq * 4. * x3sq * x3sq;
/*     < The square of the tangent of the quarter-angle for the */
/*     < spherical excess is given by Tan^2(Epsilon/4). Epsilon is the */
/*     < spherical excess defined as the sum of the three angles of the */
/*     < Euler triangle minus Pi radians.  The sum of the angles of an */
/*     < Euler triangle must be greater than Pi but less than 3Pi thus */
/*     < the spherical excess must be between zero and 2Pi. The reduced */
/*     < area is given as A/4R^2= Epsilon/4, which is less than Pi/2. */
/* Computing MAX */
    d__1 = 0., d__2 = (x3sq - x1sq) / (x2sq - x3sq);
    sqrtarg = max(d__1,d__2);
    *area = atan(sqrt(sqrtarg));
    return 0;
} /* euler_ */

/* ---------------------------------------------------------------------C */
/*     Lune                                                            C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int lune_(doublereal *xt, doublereal *ng_dollar_xt__, 
	doublereal *area, integer *cnt, integer *fest, integer *arc)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double acos(doublereal);

    /* Local variables */
    static doublereal cos_omega__, ng[3];
    static integer vert;

    /* Fortran I/O blocks */
    static cilist io___197 = { 0, 0, 0, 0, 0 };


/*     < PASSED */
/*     < LOCAL */
    /* Parameter adjustments */
    --xt;

    /* Function Body */
    s_wsle(&io___197);
    do_lio(&c__9, &c__1, " Lune: Diagonal ", (ftnlen)16);
    e_wsle();
/*     < Pull out the unit normal geodesic vector */
    ng[0] = arcinfo_1.ngx[*arc - 1];
    ng[1] = arcinfo_1.ngy[*arc - 1];
    ng[2] = arcinfo_1.ngz[*arc - 1];
/*     < Increment the arc number. */
    ++(*cnt);
/*     < Pull out the arc number of arc Cnt. */
    *arc = festooninfo_1.festoon[*cnt + *fest * 300 - 301];
/*     < Find the cosine of the angle between nG and the unit */
/*     < geodesic norm for the second arc, cos(w)= nG.nG(Arc). */
    cos_omega__ = ng[0] * arcinfo_1.ngx[*arc - 1] + ng[1] * arcinfo_1.ngy[*
	    arc - 1] + ng[2] * arcinfo_1.ngz[*arc - 1];
/*     < This angle is the complement of the angle between the two */
/*     < geodesic planes: one forming the first arc and the other */
/*     < forming the second arc. Find the angle between the planes. */
/*     < The area of a lune is defined as A= 2R^2 Alpha so the area */
/*     < reduced by 4R^2 is Alpha/2. */
    *area = acos(-cos_omega__) * .5;
/*     < Pull out the number of the terminating vertex of this */
/*     < second arc Cnt. */
    vert = arcinfo_1.arc_list__[*arc + 299];
/*     < Pull out the vector connecting the sphere center to the */
/*     < end of arc Arc.  Call it xT; T is for terminus. */
    xt[1] = arcinfo_1.vertexx[vert - 1];
    xt[2] = arcinfo_1.vertexy[vert - 1];
    xt[3] = arcinfo_1.vertexz[vert - 1];
/*     < The area of this lune should be subtracted from total */
/*     < spherical polygon area if it directed in an anti-clockwise */
/*     < manner on the sphere surface and added if directed in a */
/*     < clockwise manner.  The directionality is determined by dot */
/*     < product of geodesic unit normal of the first arc with the */
/*     < vector connecting the sphere center with the terminus of the */
/*     < second arc, Vert. Note that we can not use the vector x */
/*     < connecting the sphere center to the reference vertex as we */
/*     < normally do since this lies on the geodesic plane of the first */
/*     < arc and the second arc. */
    *ng_dollar_xt__ = ng[0] * xt[1] + ng[1] * xt[2] + ng[2] * xt[3];
    return 0;
} /* lune_ */

/* =====================================================================C */
/*                          UTILITY ROUTINES                           C */
/* =====================================================================C */
/* ---------------------------------------------------------------------C */
/*     Cross Product Of Two 3x1 Vectors                                C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int cross_(doublereal *a, doublereal *b, doublereal *axb)
{
/*     < PASSED */
    /* Parameter adjustments */
    --axb;
    --b;
    --a;

    /* Function Body */
    axb[1] = a[2] * b[3] - a[3] * b[2];
    axb[2] = a[3] * b[1] - a[1] * b[3];
    axb[3] = a[1] * b[2] - a[2] * b[1];
    return 0;
} /* cross_ */

/* ---------------------------------------------------------------------C */
/*     Sort                                                            C */
/* ---------------------------------------------------------------------C */
/* Subroutine */ int sort_(doublereal *a, integer *index, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static doublereal hold;
    static integer holdi;

/*     < PASSED */
/*     < LOCAL */
/*     < Initialize the index list */
    /* Parameter adjustments */
    --index;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	index[i__] = i__;
    }
/*     < Sort A(Index) into ascending order */
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
/*       < Pull out Jth value */
	hold = a[index[j]];
	holdi = index[j];
	i__ = j - 1;
	while(i__ > 0 && a[index[i__]] > hold) {
	    index[i__ + 1] = index[i__];
	    --i__;
	}
	index[i__ + 1] = holdi;
    }
    return 0;
} /* sort_ */

