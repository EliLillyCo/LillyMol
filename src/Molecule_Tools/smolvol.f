c     volume.f - volume determination code
c     
c     Author: Lawrence R. Dodd <dodd@roebling.poly.edu>
c             Doros N. Theodorou <doros@pylos.cchem.berkeley.edu>
c     Maintainer: Lawrence R. Dodd <dodd@roebling.poly.edu>
c     Created: March 21, 1990
c     Version: 2.1
c     Date: 1994/09/29 22:16:57
c     Keywrds: volume and area determination
c     Time-stamp: <94/09/29 18:15:01 dodd>

c     Copyright (c) 1990, 1991, 1992, 1993, 1994
c     by Lawrence R. Dodd and Doros N. Theodorou.

c     This program is free software; you can redistribute it and/or
c     modify it under the terms of the GNU General Public License as
c     published by the Free Software Foundation; either version 2 of the
c     License, or (at your option) any later version.
c
c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     General Public License for more details.
c    
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

C---------------------------------------------------------------------C
C                     Plane Sphere Intersections                      C
C---------------------------------------------------------------------C
C     This program will find the total and individual volume and      C
C     exposed surface area of an arbitrary collection of spheres of   C
C     arbitrary radii cut by an arbitrary collection of planes        C
C     analytically by analyzing the plane/sphere intersections.       C
C---------------------------------------------------------------------C
C     Algorithm by: Doros N. Theodorou and Lawrence R. Dodd           C
C     Coded by: L.R. Dodd                                             C
C---------------------------------------------------------------------C
C     Created on: March 21, 1990                                      C
C       Phase 1 Completed on: March 23, 1990                          C
C       Phase 2 Completed on: April 16, 1990                          C
C       Phase 3 Completed on: May   17, 1990                          C
C       Phase 4 Completed on: June   5, 1990                          C
C       Phase 5 Completed on: July  26, 1990                          C
C---------------------------------------------------------------------C
C     Reference:                                                      C
C                                                                     C
C       "Analytical treatment of the volume and surface area of       C
C       molecules formed by an arbitrary collection of unequal        C
C       spheres intersected by planes"                                C
C                                                                     C
C     L.R. Dodd and D.N. Theodorou                                    C
C     MOLECULAR PHYSICS, Volume 72, Number 6, 1313-1345, April 1991   C
C---------------------------------------------------------------------C
C     Acknowlegement:                                                 C
C                                                                     C
C     LRD wishes to thank his mentor DNT for a stimulating and        C
C     enjoyable post-doctoral experience.                             C
C---------------------------------------------------------------------C
C     General Notes On Program:                                       C
C                                                                     C
C     This program has been written with an eye towards both          C
C     efficiency and clarity. On a philosophical note, many believe   C
C     that these ideals are mutually exclusive but in general they    C
C     are not. There are, however, a few instances where one ideal    C
C     has been given more prominence over the other. The comments in  C
C     the program, together with the associated journal article,      C
C     should help to explain any apparent logical leaps in the        C
C     algorithm.                                                      C
C                                                                     C
C     The program was intended to be used as a subroutine called      C
C     repeatly by some main program. In this case the subroutine      C
C     "VOLUME" is called by some main routine which has placed the    C
C     necessary information in common block /Raw Data/. The answers   C
C     are returned in common block /Volume Output/. I must apologize  C
C     for the poor input/output for the program. For example, the     C
C     area/volume of each sphere is not placed in /Volume Output/.    C
C                                                                     C
C     This program was developed on a Sun SPARCstation 330 using Sun  C
C     FORTRAN 1.3.1 (all trademarks of Sun Microsystems, Inc.). We    C
C     have used some of extensions to the ANSI standard including:    C
C                                                                     C
C         o  long variable names (i.e., more than six characters)     C
C         o  variable names containing the characters '_dollar_' and '_'     C
C         o  END DO used in place of the CONTINUE statement           C
C         o  DO-WHILE used in place of IF-GOTO constructs             C
C         o  excessive number of continuation lines in some FORMATs   C
C         o  generic intrinsic function calls (e.g., SIN for DSIN)    C
C         o  IMPLICIT NONE statement (needed in development)          C
C                                                                     C
C     The advantage of using non-standard FORTRAN is that it makes it C
C     considerably easier to follow the flow of a program. There are  C
C     no extraneous statement labels in this program that may have    C
C     obscured the logic (not a single GOTO was used). The previews   C
C     of the new F90 standard appear to adopt many of the features    C
C     already implemented in VMS, Sun, Cray, and IBM FORTRAN.         C
C                                                                     C
C     Note that this algorithm is completely parallelizable.          C
C                                                                     C
C                           Larry Dodd                                C
C                           dodd@mycenae.cchem.berkeley.edu           C
C                                                                     C
C                           Department of Chemical Engineering        C
C                           College of Chemistry                      C
C                           University of California at Berkeley      C
C                           Berkeley, California 94720-9989           C
C                           (415) 643-7691 (LRD)                      C
C                           (415) 643-8523 (DNT)                      C
C                           (415) 642-5927 (Lab)                      C
C                                                                     C
C                            dodd@mycenae.cchem.berkeley.edu          C
C                           doros@mycenae.cchem.berkeley.edu          C
C                                                                     C
C---------------------------------------------------------------------C
C     Note:                                                           C
C       Plane_Ordering of common block /Debug/ is, as the name        C
C       implies, for debugging purposes only as is routine ORDERING.  C
C       The information contain therein is not necessary for solving  C
C       the sphere plane problem but proved incredibly useful during  C
C       program development.                                          C
C---------------------------------------------------------------------C

C--------------------------------------------------------------------C
C     Main Subroutine                                                C
C--------------------------------------------------------------------C

      Subroutine Volume (error_encountered)

      Implicit None

      logical error_encountered

      Call Find_Planes

      Call All_Spheres (error_encountered)

      Return
      End

C---------------------------------------------------------------------C
C     Find_Planes                                                     C
C---------------------------------------------------------------------C
C     Purpose:                                                        C
C        This routine finds the planes formed by the intersection of  C
C     all possible sphere pairs.  The plane information is stored for C
C     both spheres forming a plane because eventually the problem     C
C     will be decoupled completely. For each sphere, the information  C
C     stored for the plane includes the unit normal vector pointing   C
C     away from the cutout of the sphere, the magnitude of the vector C
C     connecting the sphere center normal to the plane, a variable    C
C     called Point, which is +1 if the sphere center is contained in  C
C     the cutout for this plane and -1 if the sphere center is not    C
C     contained the plane's cutout (if the plane goes through the     C
C     center of the sphere then Point is zero), the radius of the     C
C     circle inscribed by the plane on the sphere, and a flag         C
C     indicating whether the plane has been removed.                  C
C---------------------------------------------------------------------C
C     Called by main program.                                         C
C     Output stored in /Plane Info/.                                  C
C---------------------------------------------------------------------C

      Subroutine Find_Planes

      Implicit None

      Integer*4 n_Max, Max_Plane

        Parameter ( n_Max= 300, Max_Plane= n_Max*n_Max/8 )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Common /Raw Data/ x_Chain, y_Chain, z_Chain, Radii, RefRadius
     >     , n_Spheres

        Integer*4 n_Spheres

        Real*8 RefRadius

        Real*8 x_Chain(0:n_Max), y_Chain(0:n_Max), z_Chain(0:n_Max)
     >       , Radii(0:n_Max)

      Common /Plane Info/ nx, ny, nz, Magnitude, rho_sq
     >       , n_Planes, Plane_List, Point, PlaneF, PlaneCnt

        Integer*4 PlaneCnt

        Integer*4 n_Planes(0:n_Max), Plane_List(0:n_Max,0:n_Max,3)
     >     , Point(Max_Plane,2), PlaneF(Max_Plane)

        Real*8 nx(Max_Plane), ny(Max_Plane), nz(Max_Plane)
     >       , Magnitude(Max_Plane,2), rho_sq(Max_Plane)

      Common /Global/ AreaG, NumArcsG, TrueEdgesG

        Real*8 AreaG(Max_Plane,2)

        Integer*4 NumArcsG(Max_Plane,2), TrueEdgesG(Max_Plane,2)

C     < LOCAL

        Integer*4 I, J, SphI, SphJ, K, Plane
        integer*4 n_Spheresm1

        Real*8 xi, yi, zi
     >       , del_x, del_y, del_z
     >       , dij_sq, dij, Rdij
     >       , Ri, Rj, Ri_sq, SumR, DelR, Rjsq_Risq
     >       , Termi, Termj, aTermi, aTermj


C     < Initialize.

      Plane= 0

C     write statements
C      Write(0,*) ' Number of Spheres= ', n_Spheres
C      Do I= 0, n_Spheres-1
C        write (0, *) 'x=', x_Chain(I),'  y=', y_Chain(I) ,'  z='
C     >, z_Chain(I)
C      End Do
C     end of write statement      

      if (n_Spheres .eq. 0) stop 'no spheres'

CJW      Do I= 0, n_Spheres
      Do I= 0, n_Spheres-1
        n_Planes(I)= 0
      End Do

C     < Sweep through all pairs of spheres to find planes.

      n_Spheresm1 = n_Spheres - 1
      Do SphI= 0, n_Spheresm1

C       < Has sphere SphI been removed?

        If( n_Planes(SphI) .Ne. -1 ) Then

C         < Pull out coordinates and radius of sphere SphI, ri and Ri

          xi= x_Chain(SphI)
          yi= y_Chain(SphI)
          zi= z_Chain(SphI)

          Ri= Radii(SphI)
          Ri_sq= Ri*Ri

C         < Initialize

          SphJ= SphI

C         < Loop through all SphJ > SphI until SphJ= n_Spheres
C         < or until sphere SphI has been removed
C         < Skip any spheres SphJ that have been removed

CJW          Do While( SphJ .Lt. n_Spheres )

          Do While( SphJ .Lt. n_Spheres-1)

            SphJ= SphJ + 1

C           < Has sphere SphJ been removed?

            If( n_Planes(SphJ) .Ne. -1 ) Then

C             < Find the center to center vector, rij= ri - rj

              del_x= x_Chain(SphJ) - xi
              del_y= y_Chain(SphJ) - yi
              del_z= z_Chain(SphJ) - zi

C             < Find the squared magnitude of rij, |rij|^2

              dij_sq= del_x*del_x + del_y*del_y + del_z*del_z

C             < Find the sum of the radii, Rj + Ri

              Rj= Radii(SphJ)
              SumR= Rj + Ri

C             < Do we have an intersection? Note that we make the
C             < criterion for an intersection a little tighter by
C             < putting a negative on the epsilon.  That is, dij^2 has
C             < to be less than (Ri + Rj)^2 MINUS epsilon.

              If( dij_sq - SumR*SumR .Lt. -Diminutive ) Then

C               < The spheres intersect. Is one sphere completely
C               < embedded within the other sphere?  Find the
C               < difference of the radii, Rj - Ri

                DelR= Rj - Ri

C               < Is one sphere inside the other? Note that we make
C               < the criterion for embedment a little looser by
C               < making the epsilon positive.  That is, dij^2 has to
C               < be less than (Ri - Rj)^2 PLUS epsilon.

                If( dij_sq - DelR*DelR .Lt. Diminutive ) Then

C                 < One sphere is completely embedded within the other

                  If( DelR .Gt. 0.0D0 ) Then

C                   < Sphere SphI is in sphere SphJ. Go on to next
C                   < SphI. This sphere SphI will not show up again

C                   < Turn-off all of sphere SphI's planes and mark
C                   < them as having died in trying to form the
C                   < planes.

                    Do K= 1, n_Planes(SphI)
                      I= Plane_List(SphI,K,1)
                      J= Plane_List(SphI,K,2)
                      PlaneF(I)= -1
                      NumArcsG(I,J)= -2
                      TrueEdgesG(I,J)= -2
                    End Do

C                   < Turn-off sphere SphI itself

                    n_Planes(SphI)= -1

C                   < Get out of inner loop

CJW                    SphJ= n_Spheres
                    SphJ= n_Spheres-1

                  Else

C                   < Sphere SphJ is in sphere SphI. Go on to next
C                   < SphJ.  Next time SphJ shows up we will skip it

C                   < Turn-off all of sphere SphJ's planes and mark
C                   < them as having died in trying to form the
C                   < planes.

                    Do K= 1, n_Planes(SphJ)
                      I= Plane_List(SphJ,K,1)
                      J= Plane_List(SphJ,K,2)
                      PlaneF(I)= -1
                      NumArcsG(I,J)= -2
                      TrueEdgesG(I,J)= -2
                    End Do

C                   < Turn-off sphere SphJ itself

                    n_Planes(SphJ)= -1

                  End If

                Else

C                 < We have a non-trivial intersection of spheres SphI
C                 < and SphJ. Increment the plane count for both
C                 < spheres.

                  n_Planes(SphI)= n_Planes(SphI) + 1
                  n_Planes(SphJ)= n_Planes(SphJ) + 1

C                 < Increment the plane number

                  Plane= Plane + 1

C                 < Flag this plane as existing

                  PlaneF(Plane)= +1

C                 < Store the plane associated with spheres SphI and SphJ
C                 < Also store location in normal vector array

                  Plane_List(SphI,n_Planes(SphI),1)= Plane
                  Plane_List(SphI,n_Planes(SphI),2)= 1
                  Plane_List(SphI,n_Planes(SphI),3)= SphJ

                  Plane_List(SphJ,n_Planes(SphJ),1)= Plane
                  Plane_List(SphJ,n_Planes(SphJ),2)= 2
                  Plane_List(SphJ,n_Planes(SphJ),3)= SphI

C                 < Find the unit normal plane vectors, n.  These are
C                 < the unit normal vector pointing outside cutout for
C                 < each sphere.  It is defined as
C                 <
C                 <      ni= rij/|rij|
C                 <      nj= rji/|rji|= -rij/|rij|= -ni
C                 <
C                 < where rij is defined as rj - ri and is the vector
C                 < connecting the center of sphere SphI to the center
C                 < of sphere SphJ and where |rij| is simply dij the
C                 < center-to-center distance.  These vectors are
C                 < uneffected by location of the plane of intersection
C                 < of the two spheres.  The vectors connecting the
C                 < respective sphere centers to the normal point on
C                 < the plane of intesection are called vi and vj.  We
C                 < do not form the v's explicitly since they may or
C                 < may not point out of the cutouts but we do
C                 < calculate the magnitudes of the v's and record
C                 < where the sphere centers lie relative to their
C                 < respective cutouts.  We will define a variable
C                 < called Point to do this the values of which will be
C                 <
C                 <   If Point(Pi) = -1      ===>
C                 <       SphI's center is NOT in its own cutout
C                 <       ri is on the 'negative' side of plane (P)
C                 <       ni= -vi/|vi|
C                 <   Else If Point(Pi) +1   ===>
C                 <       SphI's center is in its own cutout
C                 <       ri is on the 'positive' side of plane (P)
C                 <       ni= vi/|vi|
C                 <   Else If Point(Pi) = 0 ===>
C                 <       Plane goes through SphI's center
C                 <       rk and ri are the same point
C                 <       |vi|= 0 and ni := -vj
C                 <
C                 < Both are positive or one (and only one) is negative

                  dij= Sqrt( dij_sq )
                  Rdij= 1.0D0/dij

C                 < Find the unit normal vector for the leading sphere
C                 < (i.e., the sphere with the smaller number).  The
C                 < unit normal vector for the other sphere is simply
C                 < the negative of the stored unit normal vector.  We
C                 < will save on memory by just storing the unit
C                 < normal for one side of the plane and then finding
C                 < the other unit normal vector when needed by
C                 < multiplying this one by minus one (see Find_Apices).

                  nx(Plane)= del_x*Rdij
                  ny(Plane)= del_y*Rdij
                  nz(Plane)= del_z*Rdij

C                 < Initialize the global lists for consistency checks.

                  NumArcsG(Plane,1)= -999
                  NumArcsG(Plane,2)= -999
                  TrueEdgesG(Plane,1)= -999
                  TrueEdgesG(Plane,2)= -999

C                 < Determine the proper vector magnitudes, the
C                 < location of the plane of intersection and the
C                 < radius of the inscribed circle.

                  Rjsq_Risq= DelR*SumR

                  Termi= 0.50D0*(dij_sq - Rjsq_Risq)
                  Termj= 0.50D0*(dij_sq + Rjsq_Risq)

                  aTermi= Abs( Termi )
                  aTermj= Abs( Termj )

                  If( aTermi .Lt. Diminutive ) Then

C                   < If |Termi| is less than Eps then we know that
C                   < dij^2 = (Rj^2 - Ri^2) implying that the triangle
C                   < formed by dij, Rj and Ri is a right triangle with
C                   < Rj as the hypotenuse.  This means that the plane
C                   < of intersection (P) goes through SphI's center
C                   < ri and that rk (the normal point on the plane)
C                   < coincide.

C                   < Topology Note
C                   <
C                   < Plane (P) goes through SphI's Center, ri, the
C                   < normal vector vi is zero and so the magnitude of
C                   < vi is zero and the magnitude of vj is dij.  The
C                   < radius of the circle inscribed by plane on each
C                   < sphere is Ri.  Finally, Point(i) is defined as
C                   < zero and Point(j) is +1.
C                   <
C                   <    Plane (P)
C                   <           |
C                   <           |
C                   <           |
C                   <           |   ni, vi= 0
C                   <     ---   o------------>
C                   <      |    |      nj, vj
C                   < rho= Ri   |   <----------o
C                   <      |    |
C                   <     --- --x--------------x-----
C                   <           ri              rj
C                   <           rk
C                   <
C                   <           |               |
C                   <           |               |
C                   <           |----- dij -----|
C                   <           |               |

                    Point(Plane,1)=  0
                    Point(Plane,2)= +1

                    Magnitude(Plane,1)= 0.0D0
                    Magnitude(Plane,2)= dij

                    rho_sq(Plane)= Ri_sq

                  Else If( aTermj .Lt. Diminutive ) Then

C                   < If |Termj| is less than Eps then we know that
C                   < dij^2 = (Ri^2 - Rj^2) implying that the triangle
C                   < formed by dij, Rj and Ri is a right triangle with
C                   < Ri as the hypotenuse.  This means that the plane
C                   < of intersection (P) goes through SphJ's center
C                   < rj and that rk (the normal point on the plane)
C                   < coincide.

C                   < Topology Note
C                   <
C                   < Plane (P) goes through SphJ's Center, rj, the
C                   < normal vector vj is zero and so the magnitude of
C                   < vj is zero and the magnitude of vi is dij.  The
C                   < radius of the circle inscribed by plane on each
C                   < sphere is Rj. Finally, Point(j) is defined as
C                   < zero and Point(i) is +1.
C                   <
C                   <                     Plane (P)
C                   <                          |
C                   <                          |
C                   <                          |
C                   <                 nj, vj=0 |
C                   <               <----------o   ---
C                   <              ni, vi      |    |
C                   <           o--------->    |   Rj = rho
C                   <                          |    |
C                   <         --x--------------x-- ---
C                   <           ri             rj
C                   <                          rk
C                   <
C                   <           |               |
C                   <           |               |
C                   <           |----- dij -----|
C                   <           |               |

                    Point(Plane,1)= +1
                    Point(Plane,2)=  0

                    Magnitude(Plane,1)= dij
                    Magnitude(Plane,2)= 0.0D0

                    rho_sq(Plane)= Rj*Rj

                  Else

C                   < The plane of intersection does not pass through
C                   < the center of either sphere.  We have a three
C                   < possible cases: the plane is between the sphere
C                   < centers meaning that the center of each sphere is
C                   < in its own cutout, the center of sphere SphI is
C                   < not in its own cutout or the center of sphere
C                   < SphJ is not in its own cutout.
C                   <
C                   < Topology Note
C                   <
C                   < Three possible locations of plane (P)
C                   <
C                   <     P'             P''             P'''
C                   <    -|-----x--------|---------x-----|-
C                   <     rk'   ri       rk''      rj    rk'''
C                   <
C                   <           |                  |
C                   <           |--- dij= |rij| ---|
C                   <           |                  |
C                   <
C                   <   1. within  center to center line ('')
C                   <   2. outside center to center line
C                   <     a) such that ri is NOT in SphI's cutout (')
C                   <     b) such that rj is NOT in SphJ's cutout (''')
C                   <
C                   <    If Termi > 0 then the triangle formed by dij,
C                   < Ri and Rj is such that the angle between dij and
C                   < Ri, angle Thetai, is less than Pi/2 radians.
C                   < Therefore, rk is between ri and rj and the
C                   < magnitude of the vector vi is simply
C                   <
C                   <   |vi|= Ri.cos(Thetai)
C                   <       = Termi/dij
C                   <
C                   < This is case (1) in the above diagram.
C                   <
C                   <
C                   <    On the other hand, if Termi < 0 then the
C                   < triangle formed by dij, Ri and Rj is such that
C                   < the angle Thetai is greater than Pi/2 radians.
C                   < Therefore, rk is not between ri and rj and sphere
C                   < SphI's center is not in its own cutout.  In this
C                   < case the |vi| is Ri.cos(Thetai') where Thetai' is
C                   < the complementary angle of Theta
C                   <
C                   <   Thetai'= Pi - Thetai
C                   <
C                   <   cos(Thetai')= cos[Pi/2 + (Pi/2 - Thetai)]
C                   <               = -sin[Pi/2 - Thetai]
C                   <               = -cos(Thetai)
C                   <
C                   <   |vi|= Ri.cos(Thetai')
C                   <       = -Ri.cos(Thetai)
C                   <       = -Termi/dij
C                   <
C                   < This is case (2a) in the above diagram.
C                   <
C                   < So therefore, in general we may write
C                   <
C                   <   |vi|= |Termi|/dij
C                   <
C                   < Compiler Note
C                   << Sign(1,x)= +1 for x > 0
C                   << Sign(1,x)= +1 for x = 0
C                   << Sign(1,x)= -1 for x < 0

                    Point(Plane,1)= NInt( Sign(1.0D0,Termi) )
                    Point(Plane,2)= NInt( Sign(1.0D0,Termj) )

                    Magnitude(Plane,1)= aTermi*Rdij
                    Magnitude(Plane,2)= aTermj*Rdij

                    rho_sq(Plane)= Ri_sq
     >                    - Magnitude(Plane,1)*Magnitude(Plane,1)

                  End If

                End If

              End If

            End If
          End Do

        End If
      End Do

C     < Store the number of planes.

      PlaneCnt= Plane

ciaw  Write(0,*) ' Number of planes= ', PlaneCnt

C     < Initialize the global areas for the consistency checks.

      Do I= 1, PlaneCnt
        AreaG(I,1)= 0.0D0
        AreaG(I,2)= 0.0D0
      End Do

      Return
  100 Format(1X,' Warning: Sphere number ', I3
     >     ,' is EMBEDDED in sphere number ', I3)
      End



C---------------------------------------------------------------------C
C     All_Spheres                                                     C
C---------------------------------------------------------------------C
C     Purpose:                                                        C
C        This routine calls the routine Each_Sphere for every sphere  C
C     in this system.                                                 C
C---------------------------------------------------------------------C
C     Called by main program.                                         C
C---------------------------------------------------------------------C

      Subroutine All_Spheres (Error_Encountered)

      Implicit None

      Integer*4 n_Max, Max_Plane

        Parameter ( n_Max= 300, Max_Plane= n_Max*n_Max/8 )

C     passed

      logical Error_Encountered

      Real*8 Pi, Pi4

        Parameter ( Pi= 3.14159265358979323846D0, Pi4= 4.0D0*Pi )

      Common /Raw Data/ x_Chain, y_Chain, z_Chain, Radii, RefRadius
     >     , n_Spheres

        Integer*4 n_Spheres

        Real*8 RefRadius

        Real*8 x_Chain(0:n_Max), y_Chain(0:n_Max), z_Chain(0:n_Max)
     >       , Radii(0:n_Max)

      Common /Plane Info/ nx, ny, nz, Magnitude, rho_sq
     >       , n_Planes, Plane_List, Point, PlaneF, PlaneCnt

        Integer*4 PlaneCnt

        Integer*4 n_Planes(0:n_Max), Plane_List(0:n_Max,0:n_Max,3)
     >     , Point(Max_Plane,2), PlaneF(Max_Plane)

        Real*8 nx(Max_Plane), ny(Max_Plane), nz(Max_Plane)
     >       , Magnitude(Max_Plane,2), rho_sq(Max_Plane)

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Global/ AreaG, NumArcsG, TrueEdgesG

        Real*8 AreaG(Max_Plane,2)

        Integer*4 NumArcsG(Max_Plane,2), TrueEdgesG(Max_Plane,2)

      Common /Volume Output/ atomic_volume, atomic_area

      Real*8 atomic_volume(0:n_Max), atomic_area(0:n_Max)

C     < LOCAL

        Integer*4 Planei, DiffN, DiffT, NOTEQUAL, TruePlaneCnt

        Real*8 DiffA, VSphere, ASphere, Rcu
     >       , RefArea, RefVol

        Real*8 Volume(0:n_Max), Area(0:n_Max)

        Error_Encountered = .false.

C     < Find the volume and surface area of of the reference sphere in
C     < the units of the radii (e.g., Angstroms).

      RefArea= Pi4*RefRadius*RefRadius
      RefVol= RefArea*RefRadius/3.0D0

C     < Loop through all the spheres.

CJW      Do Sphere= 0, n_Spheres
C      write (6, *) 'HERE'
      Do Sphere= 0, n_Spheres-1
C       < Find the number of planes intersecting sphere

        NumPlanes= n_Planes(Sphere)

C       < How many planes do we have for this sphere?

        If( NumPlanes .Gt. 0 ) Then

C         < We have at least one plane.

          Call Each_Sphere( VSphere, ASphere , Error_Encountered)
          if (Error_Encountered) return

          Volume(Sphere)= VSphere
          Area(Sphere)= ASphere

        Else If( NumPlanes .Eq. 0 ) Then

C         < We have no planes intersecting this sphere.  Its volume
C         < and area are that of an uncut sphere.

          Volume(Sphere)= 1.0D0
          Area(Sphere)= 1.0D0

          Rsq= Radii(Sphere)*Radii(Sphere)
          Rcu= Rsq*Radii(Sphere)

        Else

C         < The value of n_Planes for this sphere is negative
C         < indicating that this sphere has been removed because it is
C         < completely embedded within another sphere.

          Volume(Sphere)= 0.0D0
          Area(Sphere)= 0.0D0

        End If

      End Do
C      write (6, *) 'THERE'

C     < Find the total volume and surface area of this system of
C     < spheres in units of reference sphere volume and area. 

CJW      Do Sphere= 0, n_Spheres
      Do Sphere= 0, n_Spheres-1
        Rsq= Radii(Sphere)*Radii(Sphere)
        atomic_area(Sphere) = Area(Sphere) * Rsq * RefArea
        Rcu= Rsq*Radii(Sphere)
        atomic_volume(Sphere) = Volume(Sphere) * Rcu * RefVol
      End Do

C      write (6, *) 'THERE1'

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     < NOTE: This section of code is only for debugging purposes and
C     < is probably not needed. In fact if external (or artificial)
C     < planes are defined in file volume.art spurious warning
C     < messages will be issued.

C     < Find out if the head and tail of each plane match.

      NOTEQUAL= 0

      Do Planei= 1, PlaneCnt

        DiffA= AreaG(Planei,1) - AreaG(Planei,2)
        DiffN= NumArcsG(Planei,1) - NumArcsG(Planei,2)
        DiffT= TrueEdgesG(Planei,1) - TrueEdgesG(Planei,2)

        If( Abs(DiffA) .Gt. 1.0D-6 .Or. DiffN .Ne. 0
     >                             .Or. DiffT .Ne. 0 ) Then

          NOTEQUAL= NOTEQUAL + 1
          Write(0,200) Planei
     >         , AreaG(Planei,1), NumArcsG(Planei,1)
     >         , TrueEdgesG(Planei,1)
     >         , AreaG(Planei,2), NumArcsG(Planei,2)
     >         , TrueEdgesG(Planei,2)
     >         , DiffA, DiffN, DiffT

        End If

      End Do
C      write (6, *) 'THERE2'

      TruePlaneCnt= 0
      Do Planei= 1, PlaneCnt
        If( NumArcsG(Planei,1) .Ge. 0 .Or.
     >       TrueEdgesG(Planei,1) .Ge. 0 .Or.
     >       NumArcsG(Planei,1) .Eq. -3  )
     >             TruePlaneCnt= TruePlaneCnt + 1
      End Do
C      write (6, *) 'THERE3'

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      If( NOTEQUAL .Ne. 0 ) Write(0,999) NOTEQUAL

      Return
  100 Format(1X, 50('@')//
     >    '            Heads             Tails        Heads - Tails'/
     >    '         -------------    -------------    -------------'/
     >    '   P#    Area    A   E    Area    A   E    Area    A   E'/
     >    '  ========================================================')
  150 Format(1X, I4, 11X, '1', F12.6
     >     , 11X, '0', 11X, '1', F12.6, 4(11X,'0'))
  151 Format(1X, I4, 9(11X, '0'))
  200 Format(1X,I4, 3(' |', F7.3,2I4), ' NOT EQUAL ' )
  999 Format(/1X,' ********** WARNING ********** '
     >       /1X,'    ', I5,' Planes NOT EQUAL '
     >       /1X,' ********** WARNING ********** '/)
      End


C---------------------------------------------------------------------C
C     Each Sphere                                                     C
C---------------------------------------------------------------------C

      Subroutine Each_Sphere( Volume, Area , Error_Encountered)

      Implicit None

      Integer*4 n_Max

        Parameter ( n_Max= 300 )

      Real*8 Pi, Pi4, Pi4_3

        Parameter ( Pi= 3.14159265358979323846D0, Pi4= 4.0D0*Pi
     >                                          , Pi4_3= Pi4/3.0D0 )

C     < PASSED

        Real*8 Volume, Area
        logical Error_Encountered

      Common /Raw Data/ x_Chain, y_Chain, z_Chain, Radii, RefRadius
     >     , n_Spheres

        Integer*4 n_Spheres

        Real*8 RefRadius

        Real*8 x_Chain(0:n_Max), y_Chain(0:n_Max), z_Chain(0:n_Max)
     >       , Radii(0:n_Max)

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

C     < LOCAL

        Real*8 Volume_CP, Area_D, Area_SphSeg, Area_SP
     >       , Area_L, Volume_SphSec
     >       , RealR, RealRsq, RealRcu, AScale, VScale, VolumeT, AreaT

C     < Initialize some stuff.

      Error_Encountered = .false.

      Volume= 0.0D0
      VolumeT= 0.0D0
      Volume_CP= 0.0D0
      Volume_SphSec= 0.0D0
      Area= 0.0D0
      AreaT= 0.0D0
      Area_L= 0.0D0
      Area_SP= 0.0D0
      Area_D= 0.0D0
      Area_SphSeg= 0.0D0

C     < Find some sphere information.

      Radius= Radii(Sphere)
      Rsq= Radius*Radius
      O_Radius= 1.0D0/Radius
      O_Rsq= O_Radius*O_Radius
      O_Rcu= O_Rsq*O_Radius

C     < For all plane pairs find the plane-plane intersections and
C     < hence the plane-plane lines that will become the edges of the
C     < arc-polygons

      Call Find_Apices

C     < Need to check to see if this sphere is a WIPEOUT

      If( NumPlanes .Gt. 0 ) Then

C       < Find the arc-polygons.

        Call Arc_Polys( Volume_CP, Area_D, Area_SphSeg ,
     .                  Error_Encountered)
        if (Error_Encountered) return

C       < Find the festoons.

        Call Festoonery( Area_SP , Error_Encountered)
        if (Error_Encountered) return

C       < Find the lateral surface area of the spherical sector
C       < reduced by the sphere area 4 Pi R^3.

        Area_L= Area_SP + Area_SphSeg - Area_D

C       < The lateral surface area may be off by an integer multiple
C       < of a sphere area.  Therefore, we need to adjust Area_L by
C       < an integer (positive or negative) such that Area_L is
C       < between zero and unity.

        Volume_SphSec= (Area_L - DInt(Area_L)) +
     >          (0.50D0 - Sign(0.50D0,Area_L))

C       < The volume of the spherical sector is equal to 1/3 of the
C       < product of the lateral area and the sphere radius.
C       < Therefore, the reduced spherical sector volume is exactly
C       < the reduced lateral surface area.  Find the total volume of
C       < this sphere reduced by 4/3 Pi R^3 and total exposed area of
C       < this sphere reduced by 4 Pi R^2.

        Volume= Volume_CP + Volume_SphSec
        Area= Volume_SphSec

C       < Find the scaling factors for this sphere.

        RealR= RefRadius*Radii(Sphere)
        RealRsq= RealR*RealR
        RealRcu= RealRsq*RealR

        AScale= Pi4*RealRsq
        VScale= Pi4_3*RealRcu

C       < Find the true volume and exposed area of this sphere in
C       < cubic angstroms.

        VolumeT= Volume*VScale
        AreaT= Area*AScale

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Find_Apices                                                     C
C---------------------------------------------------------------------C

      Subroutine Find_Apices

      Implicit None

      Integer*4 n_Max, Max_Plane, Max_Apex

        Parameter ( n_Max= 300, Max_Plane= n_Max*n_Max/8
     >            , Max_Apex= 20*n_Max )

      Common /Plane Info/ nx, ny, nz, Magnitude, rho_sq
     >       , n_Planes, Plane_List, Point, PlaneF, PlaneCnt

        Integer*4 PlaneCnt

        Integer*4 n_Planes(0:n_Max), Plane_List(0:n_Max,0:n_Max,3)
     >     , Point(Max_Plane,2), PlaneF(Max_Plane)

        Real*8 nx(Max_Plane), ny(Max_Plane), nz(Max_Plane)
     >       , Magnitude(Max_Plane,2), rho_sq(Max_Plane)

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Conclusion/ xA, xB, DeadPlane, WIPEOUT, INTERIOR, REMOVE

        Integer*4 DeadPlane

        Logical WIPEOUT, INTERIOR, REMOVE

        Real*8 xA(3), xB(3)

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Twin List/ Twin

        Integer*4 Twin(Max_Apex)

      Common /Debug/ Plane_Order

        Integer*4 Plane_Order(Max_Apex,3)

      Common /Global/ AreaG, NumArcsG, TrueEdgesG

        Real*8 AreaG(Max_Plane,2)

        Integer*4 NumArcsG(Max_Plane,2), TrueEdgesG(Max_Plane,2)

C     < LOCAL

        Integer*4 I, Planei, Indexi, ApexA, ApexB

        Real*8 Sgn

C     < Initialize some quantities.

      ApexB= 0
      Do I= 1, NumPlanes
        n_Edges(I)= 0
      End Do

C     < Store the information about this sphere's planes.

      Do I= 1, NumPlanes

C       < Find plane I, (P)

        Planei= Plane_List(Sphere,I,1)

C       < Index of plane (P) (i.e., 1st or 2nd sphere?)

        Indexi= Plane_List(Sphere,I,2)

C       < Has this plane been removed for this sphere?

        If( PlaneF(Planei) .Eq. -1 ) Then

C         < Yes, it has.  Turn off all its edges.

          n_Edges(I)= -1

C         << Store order of plane

          IndexS(I)= Indexi

C         << True plane number

          PlaneS(I)= Planei

        Else

C         << Store order of plane

          IndexS(I)= Indexi

C         << True plane number

          PlaneS(I)= Planei

C         << Sign of plane (P)

          PointS(I)= Point(Planei,Indexi)

C         << Normal vector for plane (P) pointing out of cutout.
C         << Note that if Indexi is unity then the sphere we are
C         << currently on was the leading sphere in forming the
C         << plane and the unit normal vector stored is the proper
C         << sign.  If, however, Indexi is two then the current
C         << sphere was not the leading sphere in the creation of
C         << this plane (it was the second sphere in the double
C         << loop forming planes, see Find_Planes).  Therefore, the
C         << unit normal vector is anti-parallel to the direction
C         << it should be point and so we need to multiply its
C         << coordinates by minus one.  The variable Sgn is +1 if
C         << Indexi is unity and -1 if Indexi is two.

          Sgn= Dble(3 - 2*Indexi)
          nxS(I)= nx(Planei)*Sgn
          nyS(I)= ny(Planei)*Sgn
          nzS(I)= nz(Planei)*Sgn

C         << Magnitude of vector connecting r0 to rk, |vk|

          v_MagS(I)= Magnitude(Planei,Indexi)

C         << Squared radius of circle enscribed by plane I on sphere

          rho_sqS(I)= rho_sq(Planei)

        End If

      End Do


C     ### Loop Over All Possible Pairs Of Planes.
C     -------------------------------------------
C     Note- center of sphere is denoted by r0
C           n_Edges(plane):
C                -1 ==> plane removed;
C                 0 ==> no apices or edges for plane;
C                >0 ==> number of pairs of apices or edges on plane
C                       including some apex pairs or edges that
C                       have been turned off.

C     < Initialize variables

      REMOVE= .False.
      INTERIOR= .False.
      WIPEOUT= .False.

      Do Ith= 1, NumPlanes - 1

C       < Has this plane been turned off?

        If( n_Edges(Ith) .Ne. -1 ) Then

C         < Initialize counter

          Jth= Ith

          Do While( Jth .Lt. NumPlanes )

C           < Increment counter

            Jth= Jth + 1

C           < Has this plane been turned off?

            If( n_Edges(Jth) .Ne. -1 ) Then

C             < Determine information about this pair of the planes
C             << Data passed through common block /Spheres Planes/
C             << Results passed back through common block /Conclusion/

              Call Plane_Pair

C             < Do we have a wipeout of the sphere?
C             < How did the planes intersect (if they did)?
C             < Do we have a plane removal?

              If( WIPEOUT ) Then

C               < We have a wipeout of this sphere. Turn off all of
C               < its planes and mark them as having died in trying to
C               < form the plane-plane edges.

                Do I= 1, NumPlanes
                  PlaneF(PlaneS(I))= -1
                  NumArcsG(PlaneS(I),1)= -2
                  NumArcsG(PlaneS(I),2)= -2
                  TrueEdgesG(PlaneS(I),1)= -2
                  TrueEdgesG(PlaneS(I),2)= -2
                End Do

                n_Planes(Sphere)= -1
                NumPlanes= -1
                Return

              Else If( INTERIOR ) Then

                INTERIOR= .False.

C               < Planes intersect inside sphere
C               < Increment the apex pair or edge count for each plane

                n_Edges(Ith)= n_Edges(Ith) + 1
                n_Edges(Jth)= n_Edges(Jth) + 1

C               < Increment the number of the apices
C               < Apex= ApexB= ApexA + 1

                ApexA= ApexB + 1
                ApexB= ApexB + 2

C               < Store coordinates of first apex, apex A

                Apexx(ApexA)= xA(1)
                Apexy(ApexA)= xA(2)
                Apexz(ApexA)= xA(3)

C               < Store coordinates of second apex, apex B

                Apexx(ApexB)= xB(1)
                Apexy(ApexB)= xB(2)
                Apexz(ApexB)= xB(3)

C               < For each apex, flag that it exists as part of
C               < an apex pair (i.e., edge)

                ApexF(ApexA,1)= +2
                ApexF(ApexA,2)= -1

                ApexF(ApexB,1)= +2
                ApexF(ApexB,2)= -1

C               < Initialize the Twin list.

                Twin(ApexA)= 0
                Twin(ApexB)= 0

C               < Store (in order) apices associated with plane i
C               < (Pi):(A,B) stored anti-clockwise looking from
C               < outside. This is the edge formed by the intersection
C               < of plane i with plane j within the sphere.

                Edge_List(Ith,n_Edges(Ith),1)= ApexA
                Edge_List(Ith,n_Edges(Ith),2)= ApexB
                Edge_List(Ith,n_Edges(Ith),3)= Jth

C               < Store (in order) apices associated with plane j
C               < (Pj):(B,A) stored anti-clockwise looking from
C               < outside. This is the edge formed by the intersection
C               < of plane i with plane j within the sphere.

                Edge_List(Jth,n_Edges(Jth),1)= ApexB
                Edge_List(Jth,n_Edges(Jth),2)= ApexA
                Edge_List(Jth,n_Edges(Jth),3)= Ith

C               < Store (in order) the planes around apex A
C               < (A):(Pj,Pi) stored anti-clockwise from outside

                Plane_Order(ApexA,1)= Jth
                Plane_Order(ApexA,2)= Ith
                Plane_Order(ApexA,3)= 0

C               < Store (in order) the planes around apex B
C               < (B):(Pi,Pj) stored anti-clockwise from outside

                Plane_Order(ApexB,1)= Ith
                Plane_Order(ApexB,2)= Jth
                Plane_Order(ApexB,3)= 0

              Else If( REMOVE ) Then

                REMOVE= .False.

C               < Planes intersect outside sphere such
C               < that we have a removal of a plane

C               < REMOVE all vestiges of dead plane
C               < and leave inner loop

C               < Turn-off all of dead plane's apices (i.e., edges)

                Do I= 1, n_Edges(DeadPlane)
                  ApexF(Edge_List(DeadPlane,I,1),1)= -2
                  ApexF(Edge_List(DeadPlane,I,2),1)= -2
                End Do

C               < Turn-off dead plane itself

                n_Edges(DeadPlane)= -1
                NumArcsG(PlaneS(DeadPlane),IndexS(DeadPlane))= -2
                TrueEdgesG(PlaneS(DeadPlane),IndexS(DeadPlane))= -2

C               < Get out of inner loop if necessary.  That is, if the
C               < plane that has been killed is the Jth plane then do
C               < nothing to the inner loop counter.  However, if the
C               < dead plane is Ith then alter the inner loop counter
C               < to be greater than NumPlanes

                Jth= Jth + NumPlanes*(Jth - DeadPlane)

              End If

            End If

          End Do

        End If

      End Do

C     < Store the number of apices.

      ApexCnt= ApexB

      Return
      End


C---------------------------------------------------------------------C
C     Plane_Pair                                                      C
C---------------------------------------------------------------------C

      Subroutine Plane_Pair

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

C     < LOCAL

        Real*8 nixnj(3)

        Real*8 Magsq

C     < Find cross-product of n and n', ni x nj.

      nixnj(1)= nyS(Ith)*nzS(Jth) - nzS(Ith)*nyS(Jth)
      nixnj(2)= nzS(Ith)*nxS(Jth) - nxS(Ith)*nzS(Jth)
      nixnj(3)= nxS(Ith)*nyS(Jth) - nyS(Ith)*nxS(Jth)

C     < Find the squared-magnitude of cross-product.

      Magsq= nixnj(1)*nixnj(1) + nixnj(2)*nixnj(2) + nixnj(3)*nixnj(3)

C     < Are the two planes parallel?  Note- If the cross product of
C     < the normal vectors is the zero vector then the planes are
C     < parallel.

      If( Magsq .Lt. Diminutive ) Then

C       < Planes Are Parallel

        Call Parallel

      Else

C       < Planes Intersect

        Call Intersect( nixnj )

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Parallel Planes                                                 C
C---------------------------------------------------------------------C

      Subroutine Parallel

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

C     < LOCAL

        Integer*4 PiPj, ni_dollar_nj

C     < Need to find dot-product of ni and nj, ni.nj := ni_dollar_nj.
C     < Note that we represent the dot product as an integer where
C     <   ni.nj= +1 means the planes normal vectors are parallel
C     <   ni.nj= -1 means the planes normal vectors are anti-parallel.

      ni_dollar_nj= NInt( nxS(Ith)*nxS(Jth) + nyS(Ith)*nyS(Jth)
     >                               + nzS(Ith)*nzS(Jth) )


C     ### Find The Product Of The Points.
C     -----------------------------------

      PiPj= PointS(Ith)*PointS(Jth)


C     ### Do Either Of The Planes Go Through The Sphere Center?
C     ---------------------------------------------------------

      If( PiPj .Eq. 0 ) Then

C       < One or both planes go through the center

        Call Center( ni_dollar_nj )

      Else

C       < Neither of the planes go through the sphere center

        If( ni_dollar_nj*PiPj .Eq. +1 ) Then

C         < Planes are on the same side of sphere center

          Call Same_Side( ni_dollar_nj )

        Else

C         < Planes are on different sides of sphere center

          Call Diff_Side( ni_dollar_nj )

        End If

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Parallel Planes Through Center                                  C
C---------------------------------------------------------------------C

      Subroutine Center( ni_dollar_nj )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Integer*4 ni_dollar_nj

      Common /Conclusion/ xA, xB, DeadPlane, WIPEOUT, INTERIOR, REMOVE

        Integer*4 DeadPlane

        Logical WIPEOUT, INTERIOR, REMOVE

        Real*8 xA(3), xB(3)

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

      If( ni_dollar_nj .Eq. -1 ) Then

C       < (ni_dollar_nj = -1) ==>
C       < Possible WIPEOUT

        If( Max( PointS(Ith), PointS(Jth) ) .Eq. 0 ) Then

C         < (Pointi,Pointj)= (-1,0), (0,-1) or (0,0)
C         < WIPEOUT of sphere

          WIPEOUT= .True.

        End If

      Else

C       < (ni_dollar_nj = +1) ==>
C       < REMOVAL of a plane

        If( PointS(Ith) .Eq. 0 ) Then

C         < REMOVE one plane

          If( PointS(Jth) .Eq. +1 ) Then

C           < (Pointi,Pointj)= (0,+1)
C           < REMOVE plane not going through center, planej

            REMOVE= .True.
            DeadPlane= Jth

          Else

C           < (Pointi,Pointj)= (0,-1)
C           < REMOVE plane through center, planei
C           <   or
C           < (Pointi,Pointj)= (0,0)
C           < REMOVE either plane (e.g., planei)

            REMOVE= .True.
            DeadPlane= Ith

          End If

        Else

C         < Pointj= 0
C         < REMOVE one plane

          If( PointS(Ith) .Eq. +1 ) Then

C           < (Pointi,Pointj)= (+1,0)
C           < REMOVE plane through center, planei

            REMOVE= .True.
            DeadPlane= Ith

          Else

C           < (Pointi,Pointj)= (-1,0)
C           < REMOVE plane not going through center, planej

            REMOVE= .True.
            DeadPlane= Jth

          End If

        End If

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Parallel Planes On The Same Side Of Center                      C
C---------------------------------------------------------------------C

      Subroutine Same_Side( ni_dollar_nj )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Integer*4 ni_dollar_nj

      Common /Conclusion/ xA, xB, DeadPlane, WIPEOUT, INTERIOR, REMOVE

        Integer*4 DeadPlane

        Logical WIPEOUT, INTERIOR, REMOVE

        Real*8 xA(3), xB(3)

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

C     < We know (ni.nj)*(Pi.Pj)= +1

      If( ni_dollar_nj .Eq. -1 ) Then

C       < (ni_dollar_nj = -1) ==> Points are of opposite sign
C       < Possible WIPEOUT

C       < Does the plane with larger |v| have Point= -1?

        If( Max( v_MagS(Ith), v_MagS(Jth) ) .Eq.
     >    -Min(PointS(Ith)*v_MagS(Ith), PointS(Jth)*v_MagS(Jth)) ) Then

C         < Point= -1 for plane with Max(|v|) or |vi|=|vj|
C         < WIPEOUT of sphere

          WIPEOUT= .True.

        End If

      Else

C       < (ni_dollar_nj = +1) ==> Points are of same sign
C       < REMOVE one plane

        If( (v_MagS(Ith) - v_MagS(Jth))*PointS(Ith) .Gt. 0.0D0 ) Then

C         < REMOVE plane i

          REMOVE= .True.
          DeadPlane= Ith

        Else

C         < REMOVE plane j

          REMOVE= .True.
          DeadPlane= Jth

        End If

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Parallel Planes Different Side Of Center                        C
C---------------------------------------------------------------------C

      Subroutine Diff_Side( ni_dollar_nj )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Integer*4 ni_dollar_nj

      Common /Conclusion/ xA, xB, DeadPlane, WIPEOUT, INTERIOR, REMOVE

        Integer*4 DeadPlane

        Logical WIPEOUT, INTERIOR, REMOVE

        Real*8 xA(3), xB(3)

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

C     < We know (ni.nj)*(Pi.Pj)= -1

      If( ni_dollar_nj .Eq. -1 ) Then

C       < (ni_dollar_nj= -1) ==> Points are of same sign
C       < Possible WIPEOUT

        If( PointS(Ith) .Eq. -1 ) Then

C         < WIPEOUT of sphere

          WIPEOUT= .True.

        End If

      Else

C       < (ni_dollar_nj= +1) ==> Points are of opposite sign
C       < REMOVE one plane

        If( PointS(Ith) .Eq. +1 ) Then

C         < REMOVE plane i

          REMOVE= .True.
          DeadPlane= Ith

        Else

C         < REMOVE plane j

          REMOVE= .True.
          DeadPlane= Jth

        End If

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Intersecting Planes                                             C
C---------------------------------------------------------------------C

      Subroutine Intersect( nixnj )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

         Real*8 nixnj(3)

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

C     < LOCAL

        Real*8 ni_dollar_nj, Mag_Fsq, Alpha, Beta

        Real*8 b(2), F(3)


C     ### Find Point On Line Of Intersection Of Two Planes.
C     -----------------------------------------------------
C
C         F is a vector from the sphere center (r0) to point rH on line
C     l(Pi,Pj) formed by intersection of planes ("F" stands for
C     "Foot").  This vector is perpendicular the line of intersection
C     and hence the shortest possible vector connecting r0 to l(Pi,Pj).
C     Therefore we have
C
C               ni.(F - vi)= 0
C               nj.(F - vj)= 0
C                 (nixnj).F= 0
C
C     with ni.vi= Pointi*|vi|, nj.vj= Pointj*|vj| we have
C
C               ni.F= ni.vi = Pointi*|vi|
C               nj.F= nj.vj = Pointj*|vj|
C          (nixnj).F= 0
C
C     we can solve this system for vector F using Cramer's Rule
C
C               |  ni   |       | Pointi*|vi| |    | bi |
C               |  nj   | . F = | Pointj*|vj| | := | bj |
C               | nixnj |       |      0      |    | 0  |
C
C     We propose
C
C               F= Alpha*ni + Beta*nj + Gamma*(nixnj)
C
C     therefore
C
C            ni.F= Alpha + Beta*(ni.nj) + Gamma*0 := bi
C            nj.F= Alpha*(ni.nj) + Beta + Gamma*0 := bj
C       (nixnj).F= Alpha*0 + Beta*0 + Gamma*(nixnj).(nixnj) := 0
C
C     with a solution of
C
C           Alpha= [bi - Beta*(ni.nj)]
C            Beta= [bj - bi*(ni.nj)]/[1 - (ni.nj)^2]
C           Gamma= 0

      b(1)= PointS(Ith)*v_MagS(Ith)
      b(2)= PointS(Jth)*v_MagS(Jth)


C     ### Need To Find Dot-Product of ni and nj, ni.nj := ni_dollar_nj
C     ---------------------------------------------------------

      ni_dollar_nj= nxS(Ith)*nxS(Jth) + nyS(Ith)*nyS(Jth) 
     >              + nzS(Ith)*nzS(Jth)


C     ### Find The Constants.
C     -----------------------

      Beta= (b(2) - b(1)*ni_dollar_nj)/(1.0D0-ni_dollar_nj*ni_dollar_nj)
      Alpha= b(1) - Beta*ni_dollar_nj


C     ### Find The Vector F.
C     ----------------------

      F(1)= Alpha*nxS(Ith) + Beta*nxS(Jth)
      F(2)= Alpha*nyS(Ith) + Beta*nyS(Jth)
      F(3)= Alpha*nzS(Ith) + Beta*nzS(Jth)


C     ### Find The Square-Magnitude Of F, |F|^2
C     -----------------------------------------

      Mag_Fsq= F(1)*F(1) + F(2)*F(2) + F(3)*F(3)


C     ### Do the planes intersect outside the sphere?
C     ### This will be true if |F|^2 is greater then
C     ### the square of the radius of the sphere.
C     -----------------------------------------------

      If( Mag_Fsq .Gt. Rsq ) Then

C       < Planes Intersect Outside Sphere

        Call Outside( b, ni_dollar_nj, Mag_Fsq )

      Else

C       < Planes Intersect Inside Sphere

        Call Inside( nixnj, b, ni_dollar_nj )

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Planes Intersecting Outside Sphere                              C
C---------------------------------------------------------------------C

      Subroutine Outside( b, ni_dollar_nj, Mag_Fsq )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Real*8 ni_dollar_nj, Mag_Fsq

        Real*8 b(2)

      Common /Conclusion/ xA, xB, DeadPlane, WIPEOUT, INTERIOR, REMOVE

        Integer*4 DeadPlane

        Logical WIPEOUT, INTERIOR, REMOVE

        Real*8 xA(3), xB(3)

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

C     < LOCAL

        Integer*4 PiPj

        Real*8 Scalar


C     ### Find The Product Of The Points.
C     -----------------------------------

      PiPj= PointS(Ith)*PointS(Jth)


C     ### Do either of the planes go through the sphere center?
C     ---------------------------------------------------------

      If( PiPj .Eq. 0 ) Then

C       < One of the planes goes through the center

        Call CenterOut( ni_dollar_nj )

      Else

C       < Neither plane goes through center
C       < Form  [vj x F].[vi x F]/(|vi|*|vj|)

C       Note:
C       There are three vectors eminenting from the sphere center:
C       F= (rH - r0), vi= (rki - r0) and vj= (rkj - r0).
C       The vectors (vi x F) and (vj x F) will be parallel or
C       anti parallel (i.e., their dot product will be positive
C       or negative) depending on the order of the vectors. If F is
C       found between vi and vj then Scalar (below) will be negative
C       and this means r0 is between the planes (Pi) and (Pj).
C
C       Topology Note -
C
C               rki  [vi]                       rki  [vi]
C              /                               /
C             /                               /
C            /                               /
C           /                               /
C       r0 <---------->rH  [F]          r0 <---------->rkj [vj]
C           \                               \
C            \           Scalar < 0          \           Scalar > 0
C             \                               \
C              \                               \
C               rkj  [vj]                       rH  [F]

        Scalar= Mag_Fsq*ni_dollar_nj*PiPj - (v_MagS(Ith)*v_MagS(Jth))

C       < Is Sphere Center Between The Planes?

        If( Scalar .Lt. 0.0D0 ) Then

C         < Sphere center is between planes

          If( PiPj .Eq. +1 ) Then

C           < Possible WIPEOUT

            If( PointS(Ith) .Eq. -1 ) Then

C             < WIPEOUT of sphere

              WIPEOUT= .True.

            End If

          Else

C           < REMOVE one plane

            If( PointS(Jth) .Eq. +1 ) Then

C             < REMOVE plane j

              REMOVE= .True.
              DeadPlane= Jth

            Else

C             < REMOVE plane i

              REMOVE= .True.
              DeadPlane= Ith

            End If

          End If

        Else

C         < Sphere center is not between planes

          If( PiPj .Eq. +1 ) Then

C           < REMOVE one plane
C           < The plane  with the larger |v| if its Point is +1
C           < or the plane with the smaller |v| if its Point is -1.

            If( (v_MagS(Ith)-v_MagS(Jth))*PointS(Ith) .Gt. 0.0D0 ) Then

C             < REMOVE Plane i

              REMOVE= .True.
              DeadPlane= Ith

            Else

C             < REMOVE Plane j

              REMOVE= .True.
              DeadPlane= Jth

            End If

          Else

C           < Possible WIPEOUT
C           < Does the plane with larger |v| have Point= -1?

            If( Max(v_MagS(Ith),v_MagS(Jth)) .Eq. -Min(b(1),b(2))) Then

C             < Point= -1 for plane with Max(|v|) or |vi|=|vj|
C             < WIPEOUT of sphere

              WIPEOUT= .True.

            End If

          End If

        End If

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Planes Intersecting Outside Sphere One Going Through Center     C
C---------------------------------------------------------------------C

      Subroutine CenterOut( ni_dollar_nj )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Real*8 ni_dollar_nj

      Common /Conclusion/ xA, xB, DeadPlane, WIPEOUT, INTERIOR, REMOVE

        Integer*4 DeadPlane

        Logical WIPEOUT, INTERIOR, REMOVE

        Real*8 xA(3), xB(3)

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

C     ### One Of The Planes Goes Through The Center Of The Sphere.
C     ------------------------------------------------------------

      If( ni_dollar_nj .Gt. 0.0D0 ) Then

        If( (v_MagS(Ith) - v_MagS(Jth))*PointS(Ith) .Gt. 0.0D0 ) Then

C         < REMOVE plane i

          REMOVE= .True.
          DeadPlane= Ith

        Else

C         < REMOVE plane j

          REMOVE= .True.
          DeadPlane= Jth

        End If

      Else If( Min(PointS(Ith), PointS(Jth)) .Eq. -1 ) Then

C       < WIPEOUT of sphere

        WIPEOUT= .True.

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Planes Intersecting Inside Sphere                               C
C---------------------------------------------------------------------C

      Subroutine Inside( nixnj, b, ni_dollar_nj )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Real*8 ni_dollar_nj

        Real*8 nixnj(3), b(2)

      Common /Conclusion/ xA, xB, DeadPlane, WIPEOUT, INTERIOR, REMOVE

        Integer*4 DeadPlane

        Logical WIPEOUT, INTERIOR, REMOVE

        Real*8 xA(3), xB(3)

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

C     < LOCAL

        Real*8 Alpha, Beta, GammaA, Group


C     ### Find Where The Line Of Intersection Cuts The Sphere Surface.
C     ----------------------------------------------------------------
C     Note -
C
C     Solve the system
C
C           ni.(x - vi)= 0
C           ni.(x - vj)= 0
C
C     for x= (Point on Line - r0) such that |x|= Sphere Radius or
C
C           ni.x= Pointi*|vi|
C           nj.x= Pointj*|vj|
C            x.x= R^2
C
C     System is quadratic; there are two solutions, xA and xB
C
C             xA= Alpha*ni + Beta*nj + GammaA*(ni x nj)
C             xB= Alpha*ni + Beta*nj + GammaB*(ni x nj)

      INTERIOR= .True.

      Group= 1.0D0/(1.0D0 - ni_dollar_nj*ni_dollar_nj)

      Beta= (b(2) - b(1)*ni_dollar_nj)*Group
      Alpha= b(1) - Beta*ni_dollar_nj
      GammaA= -Sqrt( Group*(Rsq - (Alpha*b(1) + Beta*b(2))) )

      Group= GammaA*nixnj(1)
      xA(1)= Alpha*nxS(Ith) + Beta*nxS(Jth)
      xB(1)= xA(1) - Group
      xA(1)= xA(1) + Group

      Group= GammaA*nixnj(2)
      xA(2)= Alpha*nyS(Ith) + Beta*nyS(Jth)
      xB(2)= xA(2) - Group
      xA(2)= xA(2) + Group

      Group= GammaA*nixnj(3)
      xA(3)= Alpha*nzS(Ith) + Beta*nzS(Jth)
      xB(3)= xA(3) - Group
      xA(3)= xA(3) + Group

C     < Order the points A and B for each plane such that when looking
C     < at each plane from the outside of the cutout (i.e., the part of
C     < the sphere that remains) the points are read in an
C     < anti-clockwise direction.  Thus, the vector connecting A and B
C     < points in an anti-clockwise direction relative to each plane.
C     < This is accomplished by looking at (nixnj).(xA-xB).  If this is
C     < negative then A and B are stored as (A,B) for plane (Pi) and
C     < (B,A) for plane (Pj).  Otherwise if it is positive then A and B
C     < are stored as (B,A) for plane (Pi) and (A,B) for plane (Pj)
C     <
C     < The line of intersection is
C     <      (xA - xB)= (GammaA - GammaB)*(ni x nj)
C     < Therefore the dot product is
C     <       (ni x nj).(xA - xB)= (GammaA - GammaB)*|nixnj|^2
C     < Sign[(ni x nj).(xA - xB)]= Sign(GammaA - GammaB)
C     <
C     < But we know GammaB= - GammaA therefore
C     <
C     < Sign[(ni x nj).(xA - xB)]= Sign(GammaA*2)
C     <                          = Sign(GammaA)
C     <
C     < But GammaA is defined as (-Sqrt) so (nixnj).(xA-xB) is by
C     < definition less than zero.  Therefore we store
C     <
C     <     (Pi):(A,B) and (Pj):(B,A)
C     <     (A):(Pj,Pi) and (B):(Pi,Pj)

      Return
      End


C---------------------------------------------------------------------C
C     Arc-Polygons                                                    C
C---------------------------------------------------------------------C

      Subroutine Arc_Polys( Volume_CP, Area_D, Area_SphSeg , 
     .                      Error_Encountered)

      Implicit None

      Integer*4 n_Max, Max_Plane, Max_Apex, Max_Arc

        Parameter ( n_Max= 300, Max_Plane= n_Max*n_Max/8
     >            , Max_Apex= 20*n_Max
     >            , Max_Arc= n_Max )

      Real*8 Pi, Two_Pi, O_8Pi, O_2Pi

        Parameter ( Pi= 3.14159265358979323846D0, Two_Pi= Pi + Pi
     >            , O_8Pi= 0.125D0/Pi, O_2Pi= 1.0D0/Two_Pi )

C     < PASSED

        Real*8 Volume_CP, Area_D, Area_SphSeg
        logical Error_Encountered

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Spheres Planes/ nxS, nyS, nzS, v_MagS, rho_sqS
     >            , IndexS, PlaneS, PointS, Ith, Jth

        Integer*4 Ith, Jth

        Integer*4 IndexS(Max_Apex), PlaneS(Max_Apex), PointS(Max_Apex)

        Real*8 nxS(Max_Apex), nyS(Max_Apex), nzS(Max_Apex)
     >       , v_MagS(Max_Apex), rho_sqS(Max_Apex)

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Arc Info/ Vertexx, Vertexy, Vertexz
     >       , nGx, nGy, nGz, TanGeo
     >       , Arc_List, Arc_Map, ArcCnt, VertexCnt

        Integer*4 ArcCnt, VertexCnt

        Integer*4 Arc_List(Max_Arc,2), Arc_Map(Max_Arc)

        Real*8 Vertexx(Max_Arc), Vertexy(Max_Arc), Vertexz(Max_Arc)
     >       , nGx(Max_Arc), nGy(Max_Arc), nGz(Max_Arc)
     >       , TanGeo(Max_Arc)

      Common /Diangle Info/ EtaD

        Real*8 EtaD(2)

      Common /Global/ AreaG, NumArcsG, TrueEdgesG

        Real*8 AreaG(Max_Plane,2)

        Integer*4 NumArcsG(Max_Plane,2), TrueEdgesG(Max_Plane,2)

C     < LOCAL

        Integer*4 I, J, EdgeI, N_SphSeg

        Real*8 O_rhoi, Eta_AP, Pvi_Mag

        Error_Encountered = .false.

C     ### Initialize Some Things For This Sphere.
C     -------------------------------------------

C     < Initialize the triple plane intersection counter,
C     < the arc counter and the vertex counter.

      TripCnt= 0
      ArcCnt= 0
      VertexCnt= 0

C     < Initialize the Cone-Pyramid volume, the diangle area, and the
C     < area of the spherical segments.

      Volume_CP= 0.0D0
      Area_D= 0.0D0
      Area_SphSeg= 0.0D0
      N_SphSeg= 0


C     ### Single Loop Over All Planes For This Sphere.
C     ------------------------------------------------

      Do I= 1, NumPlanes

C       < Find the number of apex pairs (i.e., edges) for plane I
C       <    -1 ==> plane does not exist
C       <     0 ==> plane exists as a spherical cap only
C       <    >0 ==> plane contains one or more edges

        NumEdges= n_Edges(I)

        If( NumEdges .Eq. 0 ) Then

C         < This plane does not intersect with other planes.  We have
C         < a spherical cap or SPHERICAL SEGMENT for this plane.
C         < Eta_AP (arc-polygon area divided by rhoi_sq/2) is 2Pi and
C         < we add this to the reduced volume of the cone-pyramid.

          Pvi_Mag= PointS(I)*v_MagS(I)

          Volume_CP= Volume_CP + Two_Pi*Pvi_Mag*rho_sqS(I)

C         < Find the reduced area bordered by this spherical segment.
C         < This will be added to the area of the final spherical
C         < polygon later.

          N_SphSeg= N_SphSeg + 1
          Area_SphSeg= Area_SphSeg + Pvi_Mag

          Indexi= IndexS(I)
          Planei= PlaneS(I)

          NumArcsG(Planei,Indexi)= -3
          TrueEdgesG(Planei,Indexi)= -3

        Else If( NumEdges .Gt. 0 ) Then

C         < We have an ARC-POLYGON
C         < Pull out this plane's information

          Indexi= IndexS(I)
          Planei= PlaneS(I)
          Pointi= PointS(I)

C         < Pull out its unit normal vector.

          ni(1)= nxS(I)
          ni(2)= nyS(I)
          ni(3)= nzS(I)

C         < Pull out its distance from the sphere center and inscribed
C         < radius.

          vi_Mag= v_MagS(I)
          rhoi_sq= rho_sqS(I)

C         < Find the sphere radius and |vi| reduced by the inscribed
C         < radius of the plane.

          O_rhoi_sq= 1.0D0/rhoi_sq
          O_rhoi= Sqrt(O_rhoi_sq)

          R_rhoi= Radius*O_rhoi
          vi_rhoi= vi_Mag*O_rhoi

C         < Precalculate one-half the square curvature of the inscribed
C         < circle.

          Curvei= 0.50D0*O_rhoi_sq

C         < Store all the edges for this plane

          J= 0
          Do EdgeI= 1, NumEdges
            If( ApexF(Edge_List(I,EdgeI,1),1) .Ne. -2 ) Then
              J= J + 1
              Edgesi(J,1)= Edge_List(I,EdgeI,1)
              Edgesi(J,2)= Edge_List(I,EdgeI,2)
            End If
          End Do
          NumEdges= J

C         < Store the local plane number.  This is the number of this
C         < plane for this sphere.

          PlaneNum= I

C         < Find this Arc-Polygon.  Intersect all the edges on this
C         < plane face clipping back or discarding the unneeded parts.

          Call Find_AP

C         < Has this plane been killed?

          If( NumEdges .Le. 0 ) Then

            NumArcsG(Planei,Indexi)= -1
            TrueEdgesG(Planei,Indexi)= -1

          Else

C           < Form the Arc-Polygon. Connect all the edges properly and
C           < find its area and the area of angle of its diangles.
C           < Initialize the diangle accumulators for this plane.

            EtaD(1)= 0.0D0
            EtaD(2)= 0.0D0

            Call Build_AP( Eta_AP, Error_Encountered)
            if (Error_Encountered) return

C           < Accumulate the reduced volume of the cone-pyramid.  The
C           < sign convention is to add the cone-pyramid volume for a
C           < plane if the sphere center is contained in that plane's
C           < cutout and subtract it if not.  Therefore, we
C           < premultiply by Pointi for this plane.  Note that if the
C           < plane goes through the sphere center then Pointi and
C           < vi_Mag is zero the plane makes no contribution to the
C           < cone-pyramid volume.

            Volume_CP= Volume_CP + Pointi*Eta_AP*vi_Mag*rhoi_sq

C           < Find the total diangle area for this plane and add it to
C           < the total for this sphere. The sign convention for the
C           < diangles is to add them positively if the sphere center
C           < is contained in the plane's cutout and subtract them if
C           < not.  Therefore, we premultiply by Pointi for this
C           < plane.  Later this area will be subtracted to find the
C           < lateral area.

            Area_D= Area_D + Pointi*(EtaD(1) - vi_Mag*O_Radius*EtaD(2))

          End If

        End If

      End Do

C     < Find the cone-pyramid volume reduced by 4/3 Pi R^3.

      Volume_CP= Volume_CP*O_Rcu*O_8Pi

C     < Find the diangle area for this sphere reduced by 4 Pi R^2.

      Area_D= Area_D*O_2Pi

C     < Find the area of the spherical segments reduced by 4 Pi R^2.

      Area_SphSeg= 0.50D0*(Dble(N_SphSeg) + O_Radius*Area_SphSeg)

      Return
      End


C---------------------------------------------------------------------C
C     Find The Arc-Polygon                                            C
C---------------------------------------------------------------------C

      Subroutine Find_AP

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 0.10D-8 )

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Debug/ Plane_Order

        Integer*4 Plane_Order(Max_Apex,3)

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Edge Pair/ eAB, eCD, ABxCD, rA, rC, RMagAB, RMagCD
     >       , ApexA, ApexB, ApexC, ApexD
     >       , EdgeI, EdgeJ
     >       , Ith, Jth, Kth

        Integer*4 ApexA, ApexB, ApexC, ApexD
     >     , EdgeI, EdgeJ
     >     , Ith, Jth, Kth

        Real*8 RMagAB, RMagCD

        Real*8 eAB(3), eCD(3), ABxCD(3), rA(3), rC(3)

C     < LOCAL

        Real*8 Magsq, MaxABxCD


C     ### Go Over All Pairs Of Edges (Line Segments).
C     -----------------------------------------------

      Do EdgeI= 1, NumEdges - 1

C       < Pull out apices A and B of EdgeI

        ApexA= Edgesi(EdgeI,1)
        ApexB= Edgesi(EdgeI,2)

        Do EdgeJ= EdgeI + 1, NumEdges

C         < Pull out and store apex A's coordinates.  Note that
C         < this needs to remain in the inner loop since the
C         < coordinates of ApexA or ApexB may be changed within
C         < the inner loop.

          rA(1)= Apexx(ApexA)
          rA(2)= Apexy(ApexA)
          rA(3)= Apexz(ApexA)

C         < Find vector connecting A to B (need to pull out rB)
C         < This is the line segment formed by another plane on
C         < the plane face of plane PlaneNum.

          eAB(1)= Apexx(ApexB) - rA(1)
          eAB(2)= Apexy(ApexB) - rA(2)
          eAB(3)= Apexz(ApexB) - rA(3)

C         < Find the magnitude of this vector.

          RMagAB= 1.0D0/Sqrt( eAB(1)*eAB(1)
     >                      + eAB(2)*eAB(2) + eAB(3)*eAB(3) )

C         < Find the unit vector AB.

          eAB(1)= eAB(1)*RMagAB
          eAB(2)= eAB(2)*RMagAB
          eAB(3)= eAB(3)*RMagAB

C         < Pull out apices C and D of EdgeJ

          ApexC= Edgesi(EdgeJ,1)
          ApexD= Edgesi(EdgeJ,2)

C         < Pull out and store apex C's coordinates

          rC(1)= Apexx(ApexC)
          rC(2)= Apexy(ApexC)
          rC(3)= Apexz(ApexC)

C         < Find vector connecting C to D (need to pull out rD)
C         < This is the line segment formed by another plane on
C         < the plane face of plane PlaneNum.

          eCD(1)= Apexx(ApexD) - rC(1)
          eCD(2)= Apexy(ApexD) - rC(2)
          eCD(3)= Apexz(ApexD) - rC(3)

C         < Find the reciprocal magnitude of this vector.

          RMagCD= 1.0D0/Sqrt( eCD(1)*eCD(1)
     >                      + eCD(2)*eCD(2) + eCD(3)*eCD(3) )

C         < Find the unit vector CD.

          eCD(1)= eCD(1)*RMagCD
          eCD(2)= eCD(2)*RMagCD
          eCD(3)= eCD(3)*RMagCD

C         < Find the cross product of the two vectors.

          Call Cross( eAB, eCD, ABxCD )

C         < Find the squared-magnitude of the cross-product

          Magsq= ABxCD(1)*ABxCD(1) + ABxCD(2)*ABxCD(2)
     >                             + ABxCD(3)*ABxCD(3)

C         < We must determine how these two edges lie relative to one
C         < another and relative to the cutout.  First, we must
C         < determine if the edges are parallel (or anti-parallel).
C         < We do this by examining the squared-magnitude of the the
C         < cross-product of the two edge vectors.  If this is less
C         < than some Epsilon then we find the dot-product of the two
C         < edge vectors to determine if they are parallel or
C         < anti-parallel. If the squared-magnitude is greater than
C         < some Epsilon then the vectors intersect.  We must then
C         < largest, and by definition, non-zero element of the
C         < cross-product for this will be the determinate of the
C         < linear system we will solve to find the point of
C         < intersection of the two edge vectors.

          If( Magsq .Lt. Diminutive ) Then

            If( eAB(1)*eCD(1) + eAB(2)*eCD(2)
     >                        + eAB(3)*eCD(3) .Lt. 0.0D0 ) Then

C             < The two edges are anti-parallel. We should check to
C             < see if they coincide indicating that the arc-polygon
C             < is a line and we should kill this plane face.  If they
C             < do not coincide we then need to see if they point in
C             < an anti-clockwise fashion.  If not then the plane
C             < needs to be killed.

              Call AParallel_Edges

C             < Has this plane been killed?

              If( NumEdges .Eq. -1 ) Return

            Else

C             < The two edges are parallel. The edge further from the
C             < convex body should be removed. The plane can not be
C             < killed in this case so we do not test for it.

              Call Parallel_Edges

            End If

          Else

C           < The edges are not parallel (or anti-parallel for that
C           < matter) and therefore they intersect.  We must find a
C           < non-zero element of the cross-product.  This will be the
C           < determinate of the linear system we need to solve to
C           < find the point of intersection.  Find the largest
C           < element of the cross-product.

            MaxABxCD= Max(Abs(ABxCD(1)), Abs(ABxCD(2)), Abs(ABxCD(3)))

            If( Abs(ABxCD(1)) .Eq. MaxABxCD ) Then
              Kth= 1
              Ith= 2
              Jth= 3
            Else If( Abs(ABxCD(2)) .Eq. MaxABxCD ) Then
              Kth= 2
              Ith= 3
              Jth= 1
            Else
              Kth= 3
              Ith= 1
              Jth= 2
            End If

            Call Three_Plane

C           < Has this plane been killed?

            If( NumEdges .Eq. -1 ) Return

          End If

        End Do

      End Do

      Return
      End


C---------------------------------------------------------------------C
C     Parallel Edges                                                  C
C---------------------------------------------------------------------C

      Subroutine Parallel_Edges

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Edge Pair/ eAB, eCD, ABxCD, rA, rC, RMagAB, RMagCD
     >       , ApexA, ApexB, ApexC, ApexD
     >       , EdgeI, EdgeJ
     >       , Ith, Jth, Kth

        Integer*4 ApexA, ApexB, ApexC, ApexD
     >     , EdgeI, EdgeJ
     >     , Ith, Jth, Kth

        Real*8 RMagAB, RMagCD

        Real*8 eAB(3), eCD(3), ABxCD(3), rA(3), rC(3)

C     < LOCAL

        Integer*4 ELIM

        Real*8 Magsq, ACxAB_dollar_ni

        Real*8 rAC(3), ACxAB(3)

C     < Form the vector connecting the start of the AB vector
C     < with the start of the CD vector.

      rAC(1)= rC(1) - rA(1)
      rAC(2)= rC(2) - rA(2)
      rAC(3)= rC(3) - rA(3)

C     < Find the cross product of this vector with the AB vector

      Call Cross( rAC, eAB, ACxAB )

      Magsq= ACxAB(1)*ACxAB(1) + ACxAB(2)*ACxAB(2) + ACxAB(3)*ACxAB(3)

C     < Find out if these two edges coincide.

      If( Magsq .Lt. Diminutive ) Then

C       < They coincide.

        Call Coincidental

      Else

C       < Find the dot product of the cross product ACxAB with the
C       < normal vector of plane PlaneNum, ni

        ACxAB_dollar_ni= ACxAB(1)*ni(1)+ACxAB(2)*ni(2)+ACxAB(3)*ni(3)

C       < If (ACxAB).ni is positive then edge AB is on the "inside"
C       < (closer to the cutout) and edge CD should be removed.
C       < Otherwise, (ACxAB).ni negative, edge AB is on the "outside"
C       < and should be removed.
C       <   AB Inside --> ACxAB.ni > 0 --> ELIM= EdgeJ
C       <   CD Inside --> ACxAB.ni < 0 --> ELIM= EdgeI

        ELIM= EdgeI + 
     >   (EdgeJ-EdgeI)*Int(1.0D0 + Sign(0.50D0,ACxAB_dollar_ni))

C       < Remove the edge ELIM.

        ApexF(Edgesi(ELIM,1),1)= -1
        ApexF(Edgesi(ELIM,2),1)= -1

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Parallel Coincidental Edges                                     C
C---------------------------------------------------------------------C

      Subroutine Coincidental

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Edge Pair/ eAB, eCD, ABxCD, rA, rC, RMagAB, RMagCD
     >       , ApexA, ApexB, ApexC, ApexD
     >       , EdgeI, EdgeJ
     >       , Ith, Jth, Kth

        Integer*4 ApexA, ApexB, ApexC, ApexD
     >     , EdgeI, EdgeJ
     >     , Ith, Jth, Kth

        Real*8 RMagAB, RMagCD

        Real*8 eAB(3), eCD(3), ABxCD(3), rA(3), rC(3)

      Common /Twin List/ Twin

        Integer*4 Twin(Max_Apex)

C     < We have two edges on a plane face that are parallel and
C     < coincidental.  This means that three planes are coming
C     < together to form a three plane line.  On two plane faces we
C     < will have coincidental parallel edges and on a third plane
C     < faces we will have coincidental anti-parallel edges.  One of
C     < the parallel coincidental edges will need to be perserved on
C     < the plane faces.  For the present calculations we should
C     < arbitrarily kill one of the edges if both are still alive
C     < since they are numerically identical. However, for purposes of
C     < bookkeeping they are different. This will need to be rectified
C     < for later work dealing with festoon formation using twin
C     < information.

      If( ApexF(ApexA,1) .Ne. -1 ) Then

        If( ApexF(ApexC,1) .Ne. -1 ) Then

C         < Both the AB edge and the CD edge are alive.  Have either
C         < of these edges already been involved in a parallel
C         < coincidental plane face?

          If( Twin(ApexA) .Ne. 0 ) Then

C           < The AB edge was previously involved in a parallel
C           < coincidence on this or another plane face. This means
C           < that the AB edge should continue to live. Kill the CD
C           < edge letting it adopt AB as its twin.

            ApexF(ApexC,1)= -1
            ApexF(ApexD,1)= -1

            Twin(ApexC)= ApexA
            Twin(ApexD)= ApexB

          Else If( Twin(ApexC) .Ne. 0 ) Then

C           < The CD edge was previously involved in a parallel
C           < coincidence on this or another plane face. This means
C           < that the CD edge should continue to live. Kill the AB
C           < edge letting it adopt CD as its twin.

            ApexF(ApexA,1)= -1
            ApexF(ApexB,1)= -1

            Twin(ApexA)=  ApexC
            Twin(ApexB)=  ApexD

          Else

C           < Neither of the edges was previously involved in a
C           < parallel coincidence. We must mark each apex with its
C           < twin and then arbitrarily kill one of the edges (i.e.,
C           < edge AB) . The array Twin will contain these twin
C           < mappings.  For instance, ApexA is identical to ApexC. We
C           < give Edge CD the negative apex numbers of AB to indicate
C           < that AB was killed in this confrontation and CD was kept
C           < alive.

            Twin(ApexA)=  ApexC
            Twin(ApexB)=  ApexD
            Twin(ApexC)= -ApexA
            Twin(ApexD)= -ApexB

            ApexF(ApexA,1)= -1
            ApexF(ApexB,1)= -1

          End If

        Else If( Twin(ApexC) .Gt. 0 ) Then

C         < Edge AB is still alive while edge CD is dead.  We also
C         < know that Edge CD was previously involved in a parallel
C         < coincidence but was improperly and arbitrarily killed. If
C         < either the start or end of Edge AB is still a two-plane
C         < intersection and if the corresponding start or end of CD's
C         < twin is still a two-plane intersection then the vertices
C         < should be transferred.  These are assigned in RenumArc for
C         < use in forming the vertices.

          If( ApexF(ApexA,1) .Eq. +2 .And.
     >         ApexF(Twin(ApexC),1) .Eq. +2 )
     >         ApexF(ApexA,2)= ApexF(Twin(ApexC),2)

          If( ApexF(ApexB,1) .Eq. +2 .And.
     >         ApexF(Twin(ApexD),1) .Eq. +2 )
     >         ApexF(ApexB,2)= ApexF(Twin(ApexD),2)

        End If

      Else If( ApexF(ApexC,1) .Ne. -1 .And. Twin(ApexA) .Gt. 0 ) Then

C       < Edge CD is still alive while edge AB is dead.  We also know
C       < that Edge AB was previously involved in a parallel
C       < coincidence but was improperly and arbitrarily killed. If
C       < either the start or end of Edge CD is still a two-plane
C       < intersection and if the corresponding start or end of AB's
C       < twin is still a two-plane intersection then the vertices
C       < should be transferred.  These are assigned in RenumArc for
C       < use in forming the vertices.

        If( ApexF(ApexC,1) .Eq. +2 .And.
     >       ApexF(Twin(ApexA),1) .Eq. +2 )
     >       ApexF(ApexC,2)= ApexF(Twin(ApexA),2)

        If( ApexF(ApexD,1) .Eq. +2 .And.
     >       ApexF(Twin(ApexB),1) .Eq. +2 )
     >       ApexF(ApexD,2)= ApexF(Twin(ApexB),2)

      End If

      Return
      End


C---------------------------------------------------------------------C
C     AParallel Edges                                                 C
C---------------------------------------------------------------------C

      Subroutine AParallel_Edges

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Debug/ Plane_Order

        Integer*4 Plane_Order(Max_Apex,3)

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Edge Pair/ eAB, eCD, ABxCD, rA, rC, RMagAB, RMagCD
     >       , ApexA, ApexB, ApexC, ApexD
     >       , EdgeI, EdgeJ
     >       , Ith, Jth, Kth

        Integer*4 ApexA, ApexB, ApexC, ApexD
     >     , EdgeI, EdgeJ
     >     , Ith, Jth, Kth

        Real*8 RMagAB, RMagCD

        Real*8 eAB(3), eCD(3), ABxCD(3), rA(3), rC(3)

C     < LOCAL

        Real*8 Magsq

        Real*8 rAC(3), ACxAB(3)

C     < Form the vector connecting the start of the AB vector
C     < with the start of the CD vector.

      rAC(1)= rC(1) - rA(1)
      rAC(2)= rC(2) - rA(2)
      rAC(3)= rC(3) - rA(3)

C     < Find the cross product of this vector with the AB vector

      Call Cross( rAC, eAB, ACxAB )

      Magsq= ACxAB(1)*ACxAB(1) + ACxAB(2)*ACxAB(2) + ACxAB(3)*ACxAB(3)

C     < If the magnitude of the cross-product is zero then we know
C     < that AB and CD are anti-parallel and coincide.  The
C     < arc-polygon is contained between these two vectors and since
C     < AB and CD coincide it means that the arc-polygon does not
C     < exist at all.  Turn off this plane face.

      If( Magsq .Lt. Diminutive ) Then

C       < Turn-off plane

        NumEdges= -1
        n_Edges(PlaneNum)= -1

      Else If( ACxAB(1)*ni(1) + ACxAB(2)*ni(2)
     >                        + ACxAB(3)*ni(3) .Gt. 0.0D0 ) Then

C       < The anti-parallel edges do not coincide but they do violate
C       < the laws of convexity. The triple product with the normal
C       < vector of this plane is positive. The anti-parallel edges
C       < are positioned clockwise on this plane face and so there is
C       < no arc-polygon and this plane must be killed.

        NumEdges= -1
        n_Edges(PlaneNum)= -1

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Three Plane Intersection                                        C
C---------------------------------------------------------------------C

      Subroutine Three_Plane

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 50.0D-8 )

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Edge Pair/ eAB, eCD, ABxCD, rA, rC, RMagAB, RMagCD
     >       , ApexA, ApexB, ApexC, ApexD
     >       , EdgeI, EdgeJ
     >       , Ith, Jth, Kth

        Integer*4 ApexA, ApexB, ApexC, ApexD
     >     , EdgeI, EdgeJ
     >     , Ith, Jth, Kth

        Real*8 RMagAB, RMagCD

        Real*8 eAB(3), eCD(3), ABxCD(3), rA(3), rC(3)

      Common /For Move/ Lambda, LambdaF, LambdaS, Nu1
     >                 , First, Second, Movement

        Integer*4 First, Second, Movement

        Real*8 LambdaF, LambdaS, Nu1

        Real*8 Lambda(2)

C     < LOCAL

        Integer*4 K, L, Leading, Lagging

        Logical Elimination

        Logical EdgeAlive(2)

        Real*8 ABxCD_dollar_ni

C     < Find the Lambdas for the point of intersection.

      Call Find_Lambdas( Lambda, Nu1 )

C     < Now that we have the Lambda's we can proceed onto the moving of
C     < apices and the elimination of edges.  This whole algorithm is
C     < based on two planes intersecting on another plane face and the
C     < ordering of the edges formed in a right-hand sense.

C     < The first step is to find out which of the two edges is the
C     < leading edge vector in a right-handed sense.  That is, which of
C     < the two edge vector when crossed with the other and dotted with
C     < the normal vector of the plane yields a positive number.  Form
C     < this triple product.

      ABxCD_dollar_ni= ABxCD(1)*ni(1) + ABxCD(2)*ni(2) + ABxCD(3)*ni(3)

C     < If ABxCD.ni is positive then K will be zero while if ABxCD.ni
C     < is negative then K will be unity:
C     <   AB Leading --> ABxCD.ni > 0 --> K= 0 --> L= 0
C     <   CD Leading --> ABxCD.ni < 0 --> K= 1 --> L= EdgeJ - EdgeI
C     < Thus, below First will equal EdgeI in the former case and
C     < EdgeJ in the later case.

      K= Int(1.0D0 - Sign(0.50D0, ABxCD_dollar_ni))
      L= K*(EdgeJ - EdgeI)
      First= EdgeI + L
      Second= EdgeJ - L

      Leading= 1 + K
      Lagging= 2 - K

C     < Retrieve the proper Lambda's based on the value of K.

      LambdaF= Lambda(Leading)
      LambdaS= Lambda(Lagging)

C     < Find out who is alive and who is dead.

      EdgeAlive(1)= (ApexF(ApexA,1) .Ne. -1)
      EdgeAlive(2)= (ApexF(ApexC,1) .Ne. -1)

C     < Now that we know the values of Lambda corresponding to the
C     < leading and lagging edges we can proceed.
C     <
C     < When two planes intersect they form a line or an edge that
C     < would appear on both plane faces.  When three planes intersect
C     < they form a point and so the face of each plane involved would
C     < have two edges and a point appearing on it.  So our interest is
C     < in the intersection of two lines or edges.  More generally our
C     < interest is in the intersection of two rays eminating from the
C     < edge vectors.  The picture for the intersection of the two edge
C     < rays is the following (this is the "template")
C     <
C     <                        o           o
C     <                         \         /
C     <             +--------->  \       /
C     <             |             \     /
C     <        These sections      \   /
C     <        need to be           \ /       Region of
C     <        discarded             X rT     desired
C     <             |               / \       cutout
C     <             +--------->    /   \
C     <                           /     \
C     <                          /       \
C     <                     mu1 /         \ mu2
C     <                        V_         _V
C     <                   Leading      Lagging
C     <                    Edge         Edge
C     <                    Ray          Ray
C     <
C     < The point rT is the triple plane intersection.  The parts of
C     < the real edge vectors that lie on the sections to the left will
C     < need to be discarded and if the entire edge vector lies in
C     < these sections then the entire edge should be eliminated.  So
C     < the matter of interest is where the edge vectors lie on this
C     < coordinate system relative to the point rT.
C     <
C     < Each apex of the edge vectors (i.e., the start and terminus of
C     < each edge vector) will be transformed from r-space or Cartesian
C     < space to mu-space where the two rays eminating from the edge
C     < vectors will be the mu's with rT as the origin.  We do this by
C     < defining
C     <
C     < (mu1,mu2)= (   -LambdaF,0) for the start of the leading edge
C     < (mu1,mu2)= (1 - LambdaF,0) for the terminus of the leading edge
C     < (mu1,mu2)= (0,   -LambdaS) for the start of the lagging edge
C     < (mu1,mu2)= (0,1 - LambdaS) for the terminus of the lagging edge
C     <
C     < The rules for dealing with these edges are then
C     <
C     <    o  any apex with mu1 > +Epsilon should be moved to mu1= 0
C     <       this is a clipping of the edge that apex is on and is a
C     <       result of a three-plane intersection.
C     <    o  any apex with mu2 < -Epsilon should be moved to mu2= 0
C     <       this is a clipping of the edge that apex is on and is a
C     <       result of a three-plane intersection.
C     <    o  any edge with (mu1,mu2)= (0,0) for both apices should be
C     <       eliminated this is a requirement for convexity of the
C     <       final arc-polygon.
C     <
C     < Epsilon (a very small positive number) is included in the two
C     < inequalities because when two edges intersect to form a
C     < three-plane intersection, rT, they will be both be clipped back
C     < to form two edges of the final arc-polygon for the present
C     < plane face.  Each of these clipped edges (call them 1 and 2)
C     < will appear separately on subsequent plane faces and will
C     < terminate or originate from the triple-plane intersection rT.
C     < One of the edges on each of these subsequent plane faces will
C     < be the previously unknown third edge that forms the
C     < triple-plane intersection.  Thus, it will be found that clipped
C     < edge vector 1 will terminate or originate at the point rT
C     < somewhere on the third edge (thus its edge vector 1's Lambda
C     < will be ~1 or ~0 while the third edge vector's Lambda will be
C     < in (0,1)) and we will need to clip this third edge at the
C     < triple point intersection.  Finally, when we arrive at the
C     < third plane face of the three planes involved in the three
C     < plane intersection, the clipped edge vector 2 will terminate or
C     < originate at rT on the start or terminus of the now clipped
C     < third edge vector.  In this case the one Lambda will be ~0 and
C     < another will be ~1.  From above then the three planes i, j and
C     < k will form a triad with edges 1, 2 and 3 meeting at point rT
C     <
C     <                                (1)
C     <                      Plane     /
C     <                      Face j   /
C     <                              /    Plane
C     <                  (3)--------< rT  Face i
C     <                              \
C     <                      Plane    \
C     <                      Face k    \
C     <                                (2)
C     <
C     < Try to imagine point rT above the plane of the screen (or paper
C     < or whatever media you are using to read this program) and the
C     < three plane faces i, j and k all tilting up converging at it.
C     <
C     < Thus edge vectors 1 and 2 clipped each other on plane face i
C     < forming the point rT but edge vector 3 was unknown on that
C     < plane face since it was formed by the intersection of planes j
C     < and k and when executing plane face j edge vector 1 will be
C     < found to originate on edge vector 3 (Lambda1~0 and Lambda3 in
C     < (0,1)).  We will clip edge vector 3 and then when examining
C     < plane face k we will find that edge vector 2 terminates at the
C     < start of edge vector 3 (Lambda2~1 and Lambda3~0).  To make a
C     < long story come to an end, the Epsilon is needed to handle the
C     < cases of ~0 and ~1.
C     <
C     < The other possible values of Lambda cover the cases where one
C     < edge is farther from the arc-polygon than the other vector and
C     < is either eliminated entirely or partially clipped. It also
C     < covers the cases where nothing should be done.  That is, the
C     < two edges are entirely contained in the rays closest to the
C     < cutout.  Another interesting case that is covered is when we
C     < get a four plane intersection (yes, they can occur) in which
C     < case the two of the edge vectors either have their starts or
C     < their termini touching.  One of these two edge vectors will
C     < reside on one of the rays closest to the cutout and will be
C     < kept.  In this case, the other edge vector will be eliminated.
C     < It is also possible for neither of the edge vectors to appear
C     < on the rays closest to the cutout and thus both will be
C     < eliminated.

C     < We are ready to test the values of the Lambdas.  Initialize the
C     < flags that will let us determine whether or not we have a
C     < double clipping, that is, a true intersection.

      Elimination= .False.
      Movement= 0

C     < ----------------------------
C     < Start with the LEADING EDGE.
C     < ----------------------------

      If( -LambdaF .Gt. -Diminutive ) Then

C       < Since (-LambdaF > -Eps) we can add one to both sides and
C       < write the following
C       <
C       <        1 - LambdaF > 1 - Eps  > Eps
C       <
C       < This means that the leading edge begins and ends in the
C       < undesired region and thus should be eliminated.
C       <
C       <   Lagging
C       <     o
C       <      \
C       <       \
C       <        \ o
C       <         %   <---  rT
C       <        / \
C       <       /   \
C       <      /     \
C       <     V_     _V
C       <   Leading
C       <
C       < Kill both of the apices of the First edge.
C       < ELIMINATION

        ApexF(Edgesi(First,1),1)= -1
        ApexF(Edgesi(First,2),1)= -1

        Elimination= .True.

      Else If( 1.0D0 - LambdaF .Gt. Diminutive
     >                        .And. EdgeAlive(Leading) ) Then

C       < We know that (-LambdaF < -Eps) and now we also know that
C       < (1 - LambdaF > Eps) so we have
C       <
C       <   Lagging
C       <     o
C       <      \
C       <       \   o
C       <        \ /
C       <         %   <---  rT
C       <        / \
C       <       /   \
C       <      /     \
C       <     V_     _V
C       <   Leading
C       <
C       < This means that the leading edge starts at least Eps into
C       < the desired region but its terminus is in the undesired
C       < region for the leading edge.  We have to move (clip) the end
C       < of the first edge back to rT.

        Movement= 1

      End If

C     < ------------------------
C     < Now do the LAGGING EDGE.
C     < ------------------------

      If( 1.0D0 - LambdaS .Lt. Diminutive ) Then

C       < Since (1 - LambdaS < Eps) we can subtract one from both
C       < sides and write the following
C       <
C       <        - LambdaS < Eps - 1 < -Eps
C       <
C       < This means that the lagging edge should be eliminated.
C       <
C       <   Lagging
C       <     o       o
C       <      \     /
C       <       \   /
C       <        \ /
C       <         %   <---  rT
C       <        /_V
C       <       /
C       <      /
C       <     V_
C       <   Leading
C       <
C       < Kill both apices of the Second edge.
C       < ELIMINATION

        If( Elimination ) Then

C         < Both edges have been entirely eliminated from this plane
C         < face indicating that they were clockwise to each other and
C         < so the cutout was not between them.  This means that this
C         < plane should be killed.

          NumEdges= -1
          n_Edges(PlaneNum)= -1

        Else

C         < The leading edge was not killed entirely.  Just kill the
C         < lagging edge.

          ApexF(Edgesi(Second,1),1)= -1
          ApexF(Edgesi(Second,2),1)= -1

        End If

      Else If( -LambdaS .Lt. -Diminutive .And. EdgeAlive(Lagging)) Then

C       < We know that (1 - LambdaS > Eps) and now we also know that
C       < (-LambdaS < -Eps) so we have
C       <
C       <   Lagging
C       <     o       o
C       <      \     /
C       <       \   /
C       <        \ /
C       <         %   <---  rT
C       <        / \
C       <       /  _V
C       <      /
C       <     V_
C       <   Leading
C       <
C       < This means that the terminus of the lagging edge is at least
C       < Eps into the desired region but the start of the lagging
C       < edge is in the undesired region.  We have to move (clip) the
C       < start of the lagging edge up to rT.

        Movement= Movement + 2

      End If

C     < Do we have a movement of an apex?

      If( Movement .Ne. 0 ) Call Move

      Return
      End


C---------------------------------------------------------------------C
C     Move                                                            C
C---------------------------------------------------------------------C

      Subroutine Move

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Debug/ Plane_Order

        Integer*4 Plane_Order(Max_Apex,3)

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Edge Pair/ eAB, eCD, ABxCD, rA, rC, RMagAB, RMagCD
     >       , ApexA, ApexB, ApexC, ApexD
     >       , EdgeI, EdgeJ
     >       , Ith, Jth, Kth

        Integer*4 ApexA, ApexB, ApexC, ApexD
     >     , EdgeI, EdgeJ
     >     , Ith, Jth, Kth

        Real*8 RMagAB, RMagCD

        Real*8 eAB(3), eCD(3), ABxCD(3), rA(3), rC(3)

      Common /For Move/ Lambda, LambdaF, LambdaS, Nu1
     >                 , First, Second, Movement

        Integer*4 First, Second, Movement

        Real*8 LambdaF, LambdaS, Nu1

        Real*8 Lambda(2)

C     < LOCAL

        Integer*4 Apex1, Apex2, P2, P3

C     < Pull out the location of second apex of leading edge and the
C     < location of first apex of lagging edge.

      Apex1= Edgesi(First,2)
      Apex2= Edgesi(Second,1)

      If( Movement .Eq. 3 ) Then

C       < Both apices have to be moved.  We have an intersection of
C       < the two edges on the plane face.
C       <
C       <   Apex2
C       <     o       o
C       <      \     /
C       <       \   /
C       <        \ /
C       <         X   <---  rT
C       <        / \
C       <       /   \
C       <      /     \
C       <     V_     _V
C       <   Apex1
C       <
C       < Replace coordinates of Apex1 with those of rT.

        Apexx(Apex1)= rA(1) + Nu1*eAB(1)
        Apexy(Apex1)= rA(2) + Nu1*eAB(2)
        Apexz(Apex1)= rA(3) + Nu1*eAB(3)

C       < Move Apex2 to the new Apex1.

        Apexx(Apex2)= Apexx(Apex1)
        Apexy(Apex2)= Apexy(Apex1)
        Apexz(Apex2)= Apexz(Apex1)

C       < Increment the triple plane intersection counter.

        TripCnt= TripCnt + 1

C       < Flag this apex as resulting from 3-plane intersection then
C       < store the number of the 3-plane intersection for future
C       < reference.

        ApexF(Apex1,1)= +3
        ApexF(Apex1,2)= TripCnt

        ApexF(Apex2,1)= +3
        ApexF(Apex2,2)= TripCnt

        Call Ordering( Apex1, Apex2, PlaneNum, P2, P3 )

        Plane_Order(Apex1,1)= PlaneNum
        Plane_Order(Apex1,2)= P2
        Plane_Order(Apex1,3)= P3

        Plane_Order(Apex2,1)= PlaneNum
        Plane_Order(Apex2,2)= P2
        Plane_Order(Apex2,3)= P3

      Else If( Movement .Eq. 1 ) Then

        If( Abs(LambdaS) .Lt. Diminutive ) Then

C         < Only the terminus of the leading edge needs to be moved.
C         < Since Movement is +1 this means that the start of the
C         < leading edge is at least Eps into the desired region and
C         < we have this picture.
C         <
C         <             o
C         <            /
C         <           /
C         <          /
C         <         %   <---  rT AND Apex2
C         <        / \
C         <       /   \
C         <      /     \
C         <     V_     _V
C         <   Apex1
C         <
C         < That is, we are sure that the beginning of the leading
C         < edge is in the desired section and therefore it needs to
C         < be clipped back and that the entire leading edge need not
C         < be eliminated.  So only Apex1 has to be moved and the
C         < lagging edge starts on the leading edge.  This is a
C         < degenerate intersection.  Move Apex1 to start of lagging
C         < edge and adopt all of the characteristics of Apex2.

          Apexx(Apex1)= Apexx(Apex2)
          Apexy(Apex1)= Apexy(Apex2)
          Apexz(Apex1)= Apexz(Apex2)

          ApexF(Apex1,1)= +3
          ApexF(Apex1,2)= ApexF(Apex2,2)

          Plane_Order(Apex1,1)= Plane_Order(Apex2,1)
          Plane_Order(Apex1,2)= Plane_Order(Apex2,2)
          Plane_Order(Apex1,3)= Plane_Order(Apex2,3)

        Else

C         < Only Apex1 has to be moved but the start of the lagging
C         < does not make actual contact with the leading edge.
C         <
C         <             o
C         <            /
C         <           /
C         <          /
C         <  rT --> %
C         <        /
C         <       /    o Apex2
C         <      /      \
C         <     V_       \
C         <   Apex1       \
C         <               _V
C         <
C         < Replace coordinates of Apex1 with those of rT.

          Apexx(Apex1)= rA(1) + Nu1*eAB(1)
          Apexy(Apex1)= rA(2) + Nu1*eAB(2)
          Apexz(Apex1)= rA(3) + Nu1*eAB(3)

C         < Increment the triple plane intersection counter.

          TripCnt= TripCnt + 1

C         < Let Apex1 adopt all of the characteristics of Apex2 except
C         < for its TripCnt number.

          ApexF(Apex1,1)= +3
          ApexF(Apex1,2)= TripCnt

C         < Assign plane ordering to the new end of the leading edge as
C         < if the lagging edge truly did pass through the leading
C         < edge and we had real contact.  We will not modify the plane
C         < ordering of Apex2 however.

          Call Ordering( Apex1, Apex2, PlaneNum, P2, P3 )

          Plane_Order(Apex1,1)= PlaneNum
          Plane_Order(Apex1,2)= P2
          Plane_Order(Apex1,3)= P3

        End If

        If( ApexF(Apex2,1) .Eq. +2 ) Then

C         < RARE CASE: The start of the lagging edge is still a
C         < two-plane intersection even though it lies on the leading
C         < edge (more or less). This means that three planes are
C         < intersecting on the sphere surface.  Act as if these two
C         < edges intersected normally forming a three-plane point.

          Apexx(Apex2)= Apexx(Apex1)
          Apexy(Apex2)= Apexy(Apex1)
          Apexz(Apex2)= Apexz(Apex1)

          ApexF(Apex2,1)= +3
          ApexF(Apex2,2)= ApexF(Apex1,2)

          Plane_Order(Apex2,1)= Plane_Order(Apex1,1)
          Plane_Order(Apex2,2)= Plane_Order(Apex1,2)
          Plane_Order(Apex2,3)= Plane_Order(Apex1,3)

        End If

      Else

C       < Movement must equal 2. Only Apex2 has to be moved.
C       < Does the end of the leading edge sit on the lagging edge?

        If( Abs(1.0D0 - LambdaF) .Lt. Diminutive ) Then

C         < Only the start of the lagging needs to be moved.  Since
C         < Movement is +2 this means that the terminus of the lagging
C         < edge is at least Eps into the desired region and we have
C         < this picture.
C         <
C         <   Apex2
C         <     o       o
C         <      \     /
C         <       \   /
C         <        \ /
C         <         V_  <---  rT AND Apex1
C         <          \
C         <           \
C         <            \
C         <            _V
C         <
C         < That is, we are sure that the end of the lagging edge is
C         < in the desired section and therefore it needs to be
C         < clipped back and that the entire lagging edge need not be
C         < eliminated.  So only Apex2 has to be moved and the leading
C         < edge ends on the lagging edge.  This is a degenerate
C         < intersection.  Move Apex2 to end of the leading edge and
C         < have it adopt all of the characteristics of Apex1.

          Apexx(Apex2)= Apexx(Apex1)
          Apexy(Apex2)= Apexy(Apex1)
          Apexz(Apex2)= Apexz(Apex1)

          ApexF(Apex2,1)= +3
          ApexF(Apex2,2)= ApexF(Apex1,2)

          Plane_Order(Apex2,1)= Plane_Order(Apex1,1)
          Plane_Order(Apex2,2)= Plane_Order(Apex1,2)
          Plane_Order(Apex2,3)= Plane_Order(Apex1,3)

        Else

C         < The end of the leading edge cuts the lagging edge but does
C         < not make actual contact with it.
C         <
C         <               o
C         <   Apex2      /
C         <     o       /
C         <      \     /
C         <       \   V_ Apex1
C         <        \
C         <         %  <---  rT
C         <          \
C         <           \
C         <            \
C         <            _V
C         <
C         < Replace coordinates of Apex2 with those of rT.

          Apexx(Apex2)= rA(1) + Nu1*eAB(1)
          Apexy(Apex2)= rA(2) + Nu1*eAB(2)
          Apexz(Apex2)= rA(3) + Nu1*eAB(3)

C         < Increment the triple plane intersection counter.

          TripCnt= TripCnt + 1

C         < Let Apex2 adopt all of the characteristics of Apex1 except
C         < for the TripCnt of Apex1.

          ApexF(Apex2,1)= +3
          ApexF(Apex2,2)= TripCnt

C         < Assign plane ordering to the new start of the lagging edge
C         < as if the leading edge truly did pass through the lagging
C         < edge and we had real contact.  We will not modify the plane
C         < ordering of Apex1 however.

          Call Ordering( Apex1, Apex2, PlaneNum, P2, P3 )

          Plane_Order(Apex2,1)= PlaneNum
          Plane_Order(Apex2,2)= P2
          Plane_Order(Apex2,3)= P3

        End If

        If( ApexF(Apex1,1) .Eq. +2 ) Then

C         < RARE CASE: The end of the leading edge is still a
C         < two-plane intersection even though it lies on the lagging
C         < edge (more or less). This means that three planes are
C         < intersecting on the sphere surface. Act as if these two
C         < edges intersected normally forming a three-plane point.

          Apexx(Apex1)= Apexx(Apex2)
          Apexy(Apex1)= Apexy(Apex2)
          Apexz(Apex1)= Apexz(Apex2)

          ApexF(Apex1,1)= +3
          ApexF(Apex1,2)= ApexF(Apex2,2)

          Plane_Order(Apex1,1)= Plane_Order(Apex2,1)
          Plane_Order(Apex1,2)= Plane_Order(Apex2,2)
          Plane_Order(Apex1,3)= Plane_Order(Apex2,3)

        End If

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Find The Lambdas For The Point Of Intersection                  C
C---------------------------------------------------------------------C

      Subroutine Find_Lambdas( Lambda, Nu1 )

      Implicit None

C     < PASSED

        Real*8 Nu1

        Real*8 Lambda(2)

      Common /Edge Pair/ eAB, eCD, ABxCD, rA, rC, RMagAB, RMagCD
     >       , ApexA, ApexB, ApexC, ApexD
     >       , EdgeI, EdgeJ
     >       , Ith, Jth, Kth

        Integer*4 ApexA, ApexB, ApexC, ApexD
     >     , EdgeI, EdgeJ
     >     , Ith, Jth, Kth

        Real*8 RMagAB, RMagCD

        Real*8 eAB(3), eCD(3), ABxCD(3), rA(3), rC(3)

C     < LOCAL

        Real*8 RDeter, rACi, rACj

C     < The point of intersection between the two vectors
C     < rAB= rB-rA  and rCD= rD-rC is given by point rT
C     <
C     <   rT= rA + Lambda1*rAB
C     <   rT= rC + Lambda2*rCD
C     <
C     < where the Lambda's are constants. Thus we write
C     <
C     <   Lambda1*rAB - Lambda2*rCD= (rC - rA)
C     <                            = rAC
C     <
C     < therefore
C     <
C     < [1]  Lambda1*rAB(1) - Lambda2*rCD(1)= (rC(1) - rA(1))= rACx
C     < [2]  Lambda1*rAB(2) - Lambda2*rCD(2)= (rC(2) - rA(2))= rACy
C     < [3]  Lambda1*rAB(3) - Lambda2*rCD(3)= (rC(3) - rA(3))= rACz
C     <
C     < Three equations in two unknowns.  We express this as the
C     < general two by two system Ax= b where
C     <
C     <         [ rAB(i) -rCD(i) ]     [Lambda1]    [rAC(i)]
C     <      A= [ rAB(j) -rCD(j) ], x= [Lambda2], b=[rAC(j)]
C     <
C     < with (i,j) an element of {(1,2), (3,1), (2,3)}.  This system
C     < is non-singular iff the determinate of A is non-zero, that is
C     < A is non-singular iff
C     <
C     <   -[rABxrCD](k)= -[rAB(i)*rCD(j) - rAB(j)*rCD(i)] < > 0
C     <
C     < where k is an element of {3,2,1}.  Solution of this system
C     < using Cramer's rule is given by
C     <
C     <    Lambda1= Deter(1)(k)/Deter(k)
C     <    Lambda2= Deter(2)(k)/Deter(k)
C     <
C     < where Deter(k) is the determinate of A and Deter(m)(k) is the
C     < determinate of A with the mth column replaced with b, all for
C     < a given value of k.
C     <
C     <    We will use the unit vectors so the expressions for
C     < Lambda1, Lambda2 are
C     <
C     <               1   rAC(i)*eCD(j) - rAC(j)*eCD(i)
C     <    Lambda1= ----- -----------------------------
C     <             |rAB|       [eAB x eCD](k)
C     <
C     <               1   rAC(i)*eAB(j) - rAC(j)*eAB(i)
C     <    Lambda2= ----- -----------------------------
C     <             |rCD|       [eAB x eCD](k)
C     <
C     < Note that the permissible (Kth,Ith,Jth) triplets are cyclic
C     < permutations of (1,2,3), that is, (1,2,3), (2,3,1), and
C     < (3,1,2). First find the reciprocal of Kth determinate

      RDeter= 1.0D0/ABxCD(Kth)

C     < Find the Ith and Jth right-hand side.

      rACi= rC(Ith) - rA(Ith)
      rACj= rC(Jth) - rA(Jth)

C     < Form the Lambda's.  We will also store Lambda1*|rAB| for use
C     < in the expression
C     <
C     <   rT= rA + Lambda1*rAB
C     <     = rA + Lambda1*eAB*|rAB|
C     <     = rA + (Lambda1*|rAB|)*eAB
C     <     = rA + Nu1*eAB

      Nu1= (rACi*eCD(Jth) - rACj*eCD(Ith))*RDeter
      Lambda(1)= Nu1*RMagAB
      Lambda(2)= (rACi*eAB(Jth) - rACj*eAB(Ith))*RDeter*RMagCD

      Return
      End


C---------------------------------------------------------------------C
C     Find The Proper Plane Ordering For The Point Of Intersection    C
C---------------------------------------------------------------------C

      Subroutine Ordering( Apex1, Apex2, PlaneNum, P2, P3 )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Integer*4 Apex1, Apex2, PlaneNum, P2, P3

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Debug/ Plane_Order

        Integer*4 Plane_Order(Max_Apex,3)

C     < The ordering of the planes involved in the three-plane
C     < intersection is defined via this algorithm
C     <
C     <   o Locate the present plane in the ordered plane list of
C     <     the end of the leading edge.
C     <   o Pull out the plane immediately following this plane in
C     <     the cyclic permutation.
C     <   o Locate the present plane in the ordered plane list of
C     <     the start of the lagging edge.
C     <   o Pull out the plane immediately preceding this plane
C     <     in the cyclic permutation.
C     <   o Form the triple plane order by combining the present
C     <     plane with the two other planes located.

      If( Plane_Order(Apex1,1) .Eq. PlaneNum ) Then
        P2= Plane_Order(Apex1,2)
      Else If( Plane_Order(Apex1,2) .Eq. PlaneNum ) Then
        P2= Plane_Order(Apex1,3)
        If( P2 .Eq. 0 ) P2= Plane_Order(Apex1,1)
      Else
        P2= Plane_Order(Apex1,1)
      End If

      If( Plane_Order(Apex2,1) .Eq. PlaneNum ) Then
        P3= Plane_Order(Apex2,3)
        If( P3 .Eq. 0 ) P3= Plane_Order(Apex2,2)
      Else If( Plane_Order(Apex2,2) .Eq. PlaneNum ) Then
        P3= Plane_Order(Apex2,1)
      Else
        P3= Plane_Order(Apex2,2)
      End If

      Return
      End


C---------------------------------------------------------------------C
C     Build The Arc-Polygon                                           C
C---------------------------------------------------------------------C

      Subroutine Build_AP( Eta_AP , Error_Encountered)

      Implicit None

      Integer*4 n_Max, Max_Plane, Max_Apex

        Parameter ( n_Max= 300, Max_Plane= n_Max*n_Max/8
     >            , Max_Apex= 20*n_Max )

      Real*8 Pi, O_2Pi

        Parameter ( Pi= 3.14159265358979323846D0, O_2Pi= 0.50D0/Pi )

C     < PASSED

        Real*8 Eta_AP
        logical Error_Encountered

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Sorted Out/ NumStart, NumEnd, LocStart, LocEnd, Edge_Num

        Integer*4 NumStart(-1:3), NumEnd(-1:3)
     >     , LocStart(Max_Apex,-1:3), LocEnd(Max_Apex,-1:3)
     >     , Edge_Num(0:Max_Apex,-1:1)

      Common /Local Arc Info/ Local_List, NumArcs

        Integer*4 NumArcs

        Integer*4 Local_List(Max_Apex,2)

      Common /Global/ AreaG, NumArcsG, TrueEdgesG

        Real*8 AreaG(Max_Plane,2)

        Integer*4 NumArcsG(Max_Plane,2), TrueEdgesG(Max_Plane,2)

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

C     < LOCAL

        Integer*4 I, Loc_I, Flag_I, TrueEdges

        Logical TwoPlane

        Real*8 ArcAreas, PolyArea

        external getpid
        integer getpid


        Error_Encountered = .false.

C     < Initialize the arc-polygon areas.

      ArcAreas= 0.0D0
      PolyArea= 0.0D0

C     < Sort the two-plane and the three-plane apices and the ends
C     < and starts. At the conclusion of this section we will have
C     <
C     <   NumStart(-1)= number of dead apex starts
C     <   NumStart(+2)= number of 2-plane apex starts
C     <   NumStart(+3)= number of 3-plane apex starts
C     <
C     <   LocStart(1:NumStart(-1),-1)= location of dead apex starts
C     <   LocStart(1:NumStart(+2),+2)= location of 2-plane apex starts
C     <   LocStart(1:NumStart(+3),+3)= location of 3-plane apex starts
C     <
C     < and likewise for the apex ends. Note that number of edge starts
C     < must equal the number of edge ends. We also will have an
C     < updated list of edges containing only living (non-dead) edges.
C     <
C     <   Edge_Num(*,-1)= edge numbers of all dead edges
C     <   Edge_Num(*,+1)= edge numbers of all alive (non-dead) edges
C     <
C     < where '*' represents '1:(NumStart(+2) + NumStart(+3))'

C     < Initialize the counters

      NumStart(-1)= 0
      NumStart(+2)= 0
      NumStart(+3)= 0

      NumEnd(-1)= 0
      NumEnd(+2)= 0
      NumEnd(+3)= 0

C     < Loop through the edges

      Do I= 1, NumEdges

C       < START
C       << Pull out the Edge I's start and its flag
C       << Flag_I will be equal to -1, +2 or +3 depending on
C       << whether the edge is dead, a 2-plane or a 3-plane

        Loc_I= Edgesi(I,1)
        Flag_I= ApexF(Loc_I,1)

C       << Increment the proper flag counter

        NumStart(Flag_I)= NumStart(Flag_I) + 1

C       << Store the location of this type of apex

        LocStart(NumStart(Flag_I),Flag_I)= Loc_I

C       < END
C       << Pull out the Edge I's end and its flag
C       << Flag_I will be equal to -1, +2 or +3 depending on
C       << whether the edge is dead, a 2-plane or a 3-plane

        Loc_I= Edgesi(I,2)
        Flag_I= ApexF(Loc_I,1)

C       << Increment the proper flag counter

        NumEnd(Flag_I)= NumEnd(Flag_I) + 1

C       << Store the location of this type of apex

        LocEnd(NumEnd(Flag_I),Flag_I)= Loc_I

C       < EDGE
C       << Store the number of this type of edge
C       << (NB, the types here are -1 and +1)

        Edge_Num(NumStart(+2) + NumStart(+3),Sign(1,Flag_I))= I

      End Do

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     < NOTE: These lines test whether or not the number of two-plane
C     < START apices match the number of two-plane END apices and if
C     < the number of three-plane START apices match the number of
C     < three-plane END apices.  This was used in debugging the code.
C     < However, it is possible that because of numerical precision,
C     < two or more planes may just graze a sphere and cause one or
C     < more two-plane apices to be mistaken for three-plane apices or
C     < vise-versa.  This happens rarely but when it does the program
C     < should be stopped.  The solution is to perturb the coordinates
C     < or radii of the spheres involved by a tiny amount (< 0.01
C     < Angstroms).  This is a patch but it works and the change in
C     < radii has almost no effect on the final answers.  One could
C     < replace the STOP statement with something like DEAD= .True.
C     < RETURN where the variable DEAD is tested in the calling
C     < routine and all subsequent calling routines like: If( DEAD )
C     < RETURN In this way one could "bubble" up to the top routine,
C     < perturb the radii, and restart the calculation.  I never
C     < automated this but it was effectively what I did.
C     < 
C     < Finally, for the same reason of numerical imprecision, a
C     < similar error can occur in the routine Bin_Arcs.  Perturbing
C     < coordinates slightly solves that problem too.

C     < LRD, Thu Dec 12 08:07:11 PST 1991.

      If( NumStart(+2) .Ne. NumEnd(+2) ) Then
        Write(0,*) ' (+2) <> ', NumStart(+2), NumEnd(+2)
        Write(0,*) NumEdges
        Write(0,*) NumStart
        Write(0,*) NumEnd
        Write(0,*) Sphere

        Error_Encountered = .true.

        return
      End if

      If( NumStart(+3) .Ne. NumEnd(+3) ) Then
        Write(0,*) ' (+3) <> ', NumStart(+3), NumEnd(+3)
        Write(0,*) NumEdges
        Write(0,*) NumStart
        Write(0,*) NumEnd
        Write(0,*) Sphere

        Error_Encountered = .true.

        return
      End If
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     < If we have any 2-plane starts at all we should form the arcs

      NumArcs= NumStart(+2)
      TrueEdges= 0
      TwoPlane= ( NumArcs .Gt. 0 )

      If( TwoPlane ) Then

        Call Form_Arcs( ArcAreas )

      End If

C     < If we have any 3-plane starts at all we should form the
C     < inscribed polygon on the plane face

      If( NumStart(+2) + NumStart(+3) .Gt. 1 ) Then

C       < The number of living (non-dead) edges is equal to sum of the
C       < number of 2-plane edge starts and 3-plane edge starts.  We
C       < will use only this subset of edges as ennumerated in the
C       < Edge_Num(*,1) for all future work on this plane face.

        TrueEdges= NumStart(+2) + NumStart(+3)

C       < Form the polygon

        Call Form_Polygon( PolyArea, TrueEdges, TwoPlane )

      End If

C     < Find the reduced area for this arc-polygon. Eta_AP is defined
C     < as the area of the arc-polygon divided by rhoi_sq/2, where
C     < rhoi_sq is the radius of the inscribed circle due to this
C     < plane cutting the sphere.

      Eta_AP= ArcAreas + PolyArea

C     < Store these globally for purposes of debugging.

      AreaG(Planei,Indexi)= Eta_AP*O_2Pi
      NumArcsG(Planei,Indexi)= NumArcs
      TrueEdgesG(Planei,Indexi)= TrueEdges

      Return
      End


C---------------------------------------------------------------------C
C     Form The Arcs Of The Arc-Polygon                                C
C---------------------------------------------------------------------C

      Subroutine Form_Arcs( ArcAreas )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Real*8 Pi, Two_Pi

        Parameter ( Pi= 3.14159265358979323846D0, Two_Pi= Pi + Pi )

C     < PASSED

        Real*8 ArcAreas

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Debug/ Plane_Order

        Integer*4 Plane_Order(Max_Apex,3)

      Common /Local Arc Info/ Local_List, NumArcs

        Integer*4 NumArcs

        Integer*4 Local_List(Max_Apex,2)

      Common /Sorted Out/ NumStart, NumEnd, LocStart, LocEnd, Edge_Num

        Integer*4 NumStart(-1:3), NumEnd(-1:3)
     >     , LocStart(Max_Apex,-1:3), LocEnd(Max_Apex,-1:3)
     >     , Edge_Num(0:Max_Apex,-1:1)

C     < LOCAL

        Integer*4 I, J, K, Loc_A, Loc_B

        Integer*4 IStart(Max_Apex), IEnd(Max_Apex)

        Real*8 Theta

        Real*8 rA(3), rB(3)
     >       , ThetaStart(Max_Apex), ThetaEnd(Max_Apex)

      Data  ThetaStart(1)/Two_Pi/

C     < Loop through all the 2-plane apices (starts and ends) and
C     < find the anti-clockwise angle between it and the first
C     < 2-plane apex start. This angle will be between 0 and 2 Pi
C     < radians.

C     < Pull out the reference 2-plane apex start, call it A
C     < Apex A's angle is defined as 2Pi

      Loc_A= LocStart(1,+2)

C     < Pull out A's coordinates. Note that these are vectors
C     < connecting the sphere center to the point A on the edge.

      rA(1)= Apexx(Loc_A)
      rA(2)= Apexy(Loc_A)
      rA(3)= Apexz(Loc_A)

C     < Pull out the first 2-plane apex end, call it B

      Loc_B= LocEnd(1,+2)

C     < Pull out B's coordinates. Note that these are vectors
C     < connecting the sphere center to the point B on the edge.

      rB(1)= Apexx(Loc_B)
      rB(2)= Apexy(Loc_B)
      rB(3)= Apexz(Loc_B)

C     < Find the angle between A and B
C     < HERE IS THE PROBLEM
      Call Angle( rA, rB, ThetaEnd(1) )

C     < Loop through all the remaining apices and find all the angles
C     < relative to the first 2-plane apex start, apex A.  We will do
C     < the starts and the ends at the same time.  Note that we may
C     < have only one pair of 2-plane apices.  There may in fact be
C     < just one edge on this entire plane face so that these two
C     < apices are connected via a chord.  Whichever the case may be,
C     < the following code still applies (the next do-loop would not
C     < be executed and the sorting routine would return immediately).

      Do I= 2, NumArcs

C       < Pull out the next 2-plane apex start

        Loc_B= LocStart(I,+2)

C       < Pull out B's coordinates. Note that these are vectors
C       < connecting the sphere center to the point B on the edge.

        rB(1)= Apexx(Loc_B)
        rB(2)= Apexy(Loc_B)
        rB(3)= Apexz(Loc_B)

C       < Find the angle between A and B

        Call Angle( rA, rB, ThetaStart(I) )

C       < Pull out the next 2-plane apex end

        Loc_B= LocEnd(I,+2)

C       < Pull out B's coordinates. Note that these are vectors
C       < connecting the sphere center to the point B on the edge.

        rB(1)= Apexx(Loc_B)
        rB(2)= Apexy(Loc_B)
        rB(3)= Apexz(Loc_B)

C       < Find the angle between A and B

        Call Angle( rA, rB, ThetaEnd(I) )

      End Do

C     < We must form and order the arcs of this plane face. The method
C     < we will use to form these arcs is as follows
C     <
C     <      o sort the theta-starts into ascending order (the angles
C     <        between each 2-plane apex start and the reference
C     <        2-plane apex start, A). Note that the reference
C     <        2-plane apex start will be the first element of this
C     <        ordered list since its theta-start is defined as zero.
C     <
C     <      o sort the theta-ends into ascending order (the angles
C     <        between each 2-plane apex end and the reference
C     <        2-plane apex start, A)
C     <
C     <      o match the two ordered lists; the largest theta-end
C     <        should form an arc with the largest theta-start
C     <        (i.e., the reference apex defined as having an
C     <        angle of 2Pi) and in general the Kth largest
C     <        theta-end should pair with the Kth largest
C     <        theta-start.

C     < Sort the angles for the 2-plane starts into ascending order.

      Call Sort( ThetaStart, IStart, NumArcs )

C     < Sort the angles for the 2-plane ends into ascending order.

      Call Sort( ThetaEnd, IEnd, NumArcs )

C     < Initialize the accumulator for the arc areas.

      ArcAreas= 0.0D0

C     < The arrays IStart and IEnd contain the ordered indices for the
C     < 2-plane apex starts and ends, respectively.  We may now sweep
C     < through the indices, matching the ascending starts and ends.

      Do K= 1, NumArcs

C       < Retrieve the Index of the arc pair. The Kth 2-plane
C       < apex end forms an arc with the Kth 2-plane apex start.

        I= IStart(K)
        J= IEnd(K)

C       < Store this arc as a 2-plane end and 2-plane start

        Local_List(K,1)= LocEnd(J,+2)
        Local_List(K,2)= LocStart(I,+2)

C       < Find the true anti-clockwise angle for the arc in radians

        Theta= ThetaStart(I) - ThetaEnd(J)

C       < Find the area of the circle segment formed by the arc and the
C       < chord connecting the ends of the arc.  Note that since we are
C       < dealing with the true arc angles the circular segment area
C       < formula gives the true circular segment area.  We will sum
C       < these arc areas and multiply by the squared radius of the
C       < inscribed circle's radius at the end.
C       < There might be a way to get Sin of Theta from (rAxrB).n - LRD

        ArcAreas= ArcAreas + ( Theta - Sin(Theta) )

C       < Find all the arc information for use later with the festoons.

        Call RenumArc( LocEnd(J,+2), LocStart(I,+2) )

      End Do

      Return
      End


C---------------------------------------------------------------------C
C     Angle Between Two 2-Plane Apices                                C
C---------------------------------------------------------------------C

      Subroutine Angle( rA, rB, Theta )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Real*8 rA(3), rB(3)

        Real*8 Theta

      Real*8 Pi

        Parameter ( Pi= 3.14159265358979323846D0 )

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

C     < LOCAL

        Real*8 rAB(3), AxB(3)

        Real*8 rAB_dollar_rAB, Cos_w, w

C     < Find the difference point B and point A

      rAB(1)= rB(1) - rA(1)
      rAB(2)= rB(2) - rA(2)
      rAB(3)= rB(3) - rA(3)

C     < Find the squared magnitude of this vector

      rAB_dollar_rAB= rAB(1)*rAB(1) + rAB(2)*rAB(2) + rAB(3)*rAB(3)

C     < Find the cosine of the angle between rA and rB

      Cos_w= 1.0D0 - Curvei*rAB_dollar_rAB

C     < Find the angle. Make sure |Cos_w| < 1

      If( Abs( Cos_w ) .Gt. 1.0D0 ) Then
        w= Pi*(0.50D0 - Sign(0.50D0,Cos_w))
      Else
        w= ACos( Cos_w )
      End If

C     < Determine whether w is the real anti-clockwise angle or not.
C     < Find the cross product of rA and rB, rA x rB

      Call Cross( rA, rB, AxB )

C     < Find the dot product of this vector with the
C     < normal vector of the plane ni, (rA x rB).ni and then
C     < find the true anti-clockwise angle between A and B
C     <
C     <   (rA x rB).ni > 0 ==> Theta= w
C     <   (rA x rB).ni = 0 ==> Theta= w
C     <   (rA x rB).ni < 0 ==> Theta= 2Pi - w

      Theta= Pi + (w - Pi)*
     >     Sign( 1.0D0, AxB(1)*ni(1) + AxB(2)*ni(2) + AxB(3)*ni(3) )

      Return
      End


C---------------------------------------------------------------------C
C     Renumber Arcs                                                   C
C---------------------------------------------------------------------C

      Subroutine RenumArc( Start, Terminus )

      Implicit None

      Integer*4 n_Max, Max_Apex, Max_Arc

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max, Max_Arc= n_Max )

C     < PASSED

        Integer*4 Start, Terminus

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Arc Info/ Vertexx, Vertexy, Vertexz
     >       , nGx, nGy, nGz, TanGeo
     >       , Arc_List, Arc_Map, ArcCnt, VertexCnt

        Integer*4 ArcCnt, VertexCnt

        Integer*4 Arc_List(Max_Arc,2), Arc_Map(Max_Arc)

        Real*8 Vertexx(Max_Arc), Vertexy(Max_Arc), Vertexz(Max_Arc)
     >       , nGx(Max_Arc), nGy(Max_Arc), nGz(Max_Arc)
     >       , TanGeo(Max_Arc)

C     < Increment the arc counter

      ArcCnt= ArcCnt + 1

C     < Every arc for a sphere is made up of two two-plane apices both
C     < are not necessarily created by the same pair of planes.  Each
C     < of these two-plane apices will appear on exactly two plane
C     < faces. We will renumber all of this apices defining them as
C     < vertices of the ultimate spherical polygon.  We must perserve
C     < the non-uniqueness of this apices, that is, a two-plane apex
C     < that is an arc start on one plane face is the end of another
C     < arc on another plane face.  We wish to perserve this
C     < connection.  Therefore, as we assign new vertex numbers to
C     < these arc starts and termini we must check to see that they do
C     < not already have a new vertex number.

C     < START
C     < Have we done this vertex yet?

      If( ApexF(Start,2) .Eq. -1 ) Then

C       < We have NOT done this vertex before.
C       < Assign the apex a new vertex number.

        VertexCnt= VertexCnt + 1

C       < Copy the coordinates.

        Vertexx(VertexCnt)= Apexx(Start)
        Vertexy(VertexCnt)= Apexy(Start)
        Vertexz(VertexCnt)= Apexz(Start)

C       < Mark this apex as done for next time.

        ApexF(Start,2)= VertexCnt

C       < Define this arc as starting with this vertex.

        Arc_List(ArcCnt,1)= VertexCnt

      Else

C       < We have come across this apex before.  Pull out its new
C       < vertex number and define this arc as starting with it.

        Arc_List(ArcCnt,1)= ApexF(Start,2)

      End If

C     < Later we will need to be able to establish the connectiveness
C     < of the arcs for this sphere so it will be incredibly useful to
C     < have a map from the vertex numbers of the arc starts to the
C     < arc numbers.

      Arc_Map(Arc_List(ArcCnt,1))= ArcCnt

C     < TERMINUS
C     < Now repeat for the arc terminus.

      If( ApexF(Terminus,2) .Eq. -1 ) Then

C       < We have NOT done this vertex before.
C       < Assign the apex a new vertex number.

        VertexCnt= VertexCnt + 1

C       < Copy the coordinates.

        Vertexx(VertexCnt)= Apexx(Terminus)
        Vertexy(VertexCnt)= Apexy(Terminus)
        Vertexz(VertexCnt)= Apexz(Terminus)

C       < Mark this apex as done for next time.

        ApexF(Terminus,2)= VertexCnt

C       < Define this arc as terminating with this vertex.

        Arc_List(ArcCnt,2)= VertexCnt

      Else

C       < We have come across this apex before.  Pull out its new
C       < vertex number and define this arc as terminating with it.

        Arc_List(ArcCnt,2)= ApexF(Terminus,2)

      End If

C     < Now that we have the arc start and terminus we need to find
C     < the geodesic quarter-angle, the diangle area, and the geodesic
C     < unit normal for this arc.

      Call Diangle( Arc_List(ArcCnt,1), Arc_List(ArcCnt,2) )

      Return
      End


C---------------------------------------------------------------------C
C     Diangle Information                                             C
C---------------------------------------------------------------------C

      Subroutine Diangle( Start, Terminus )

      Implicit None

      Integer*4 n_Max, Max_Apex, Max_Arc

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max, Max_Arc= n_Max )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Real*8 Sqrt2

        Parameter ( Sqrt2= 1.41421356237309504880D0 )

C     < PASSED

        Integer*4 Start, Terminus

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Arc Info/ Vertexx, Vertexy, Vertexz
     >       , nGx, nGy, nGz, TanGeo
     >       , Arc_List, Arc_Map, ArcCnt, VertexCnt

        Integer*4 ArcCnt, VertexCnt

        Integer*4 Arc_List(Max_Arc,2), Arc_Map(Max_Arc)

        Real*8 Vertexx(Max_Arc), Vertexy(Max_Arc), Vertexz(Max_Arc)
     >       , nGx(Max_Arc), nGy(Max_Arc), nGz(Max_Arc)
     >       , TanGeo(Max_Arc)

      Common /Diangle Info/ EtaD

        Real*8 EtaD(2)

C     < LOCAL

        Real*8 RMag, Cos_Theta, Cos_Chi, Cot_Chi, Denom, Cosx

        Real*8 xS(3), xT(3), SxT(3)

C     < Pull out the coordinates of the vector connecting the sphere
C     < center to the start of this arc, ArcCnt.

      xS(1)= Vertexx(Start)
      xS(2)= Vertexy(Start)
      xS(3)= Vertexz(Start)

C     < Pull out the coordinates of the vector connecting the sphere
C     < center to the terminus of this arc, ArcCnt.

      xT(1)= Vertexx(Terminus)
      xT(2)= Vertexy(Terminus)
      xT(3)= Vertexz(Terminus)

C     < Find the cosine of the angle between these vectors.  Both
C     < vectors are of magnitude R, the radius of the sphere.  The dot
C     < product is defined as,
C     <
C     <        xA.xB= |xA||xB| Cos(Theta)
C     <             = (R*R) Cos(Theta)
C     <
C     < Therefore, Cos(Theta)= xA.xB/(R*R)

      Cos_Theta= (xS(1)*xT(1) + xS(2)*xT(2) + xS(3)*xT(3))*O_Rsq

C     < Find out if Cos(Theta) is close is -1.

      If( Abs(1.0D0 + Cos_Theta) .Lt. Diminutive ) Then

C       < The geodesic angle between the arc start and end is Pi
C       < radians.  This means that the tangent of the geodesic
C       < quarter-angle is unity, the normal vector of the plane is
C       < the geodesic unit normal vector, and there is no diangle.

        TanGeo(ArcCnt)= 1.0D0
        nGx(ArcCnt)= ni(1)
        nGy(ArcCnt)= ni(2)
        nGz(ArcCnt)= ni(3)

C       < Return now.

        Return

      End If

C     < Find the square of the tangent of the geodesic quarter-angle.
C     < The square tangent of the quarter angle is related to the
C     < cosine and sine half-angle via,
C     <
C     <   Tan^2(Theta/4)= Sin^2(Theta/2)/(1 + Cos(Theta/2))^2
C     <
C     < and the half-angles are related the whole angles via,
C     <
C     <   Sin^2(Theta/2)= (1 - Cos(Theta))/2
C     <     Cos(Theta/2)= Sqrt(1 + Cos(Theta))/Sqrt(2)
C     <
C     < thus, after multiplying numerator and denominator by 2, we have
C     <
C     <                           1 - Cos(Theta)
C     <   Tan^2(Theta/4)= ---------------------------------
C     <                   [Sqrt(2) + Sqrt(1 + Cos(Theta)]^2
C     <

      Denom= Sqrt2 + Sqrt(1.0D0 + Cos_Theta)
      TanGeo(ArcCnt)= (1.0D0 - Cos_Theta)/(Denom*Denom)

C     < Find the cross-product of these two vectors.

      Call Cross( xS, xT, SxT )

C     < Find the reciprocal magnitude of the cross product.

      RMag= 1.0D0/Sqrt( SxT(1)*SxT(1) + SxT(2)*SxT(2) + SxT(3)*SxT(3) )

C     < Find the geodesic unit normal vector, nG, defined as the unit
C     < vector normal the the plane that xA and xB lie on.

      nGx(ArcCnt)= SxT(1)*RMag
      nGy(ArcCnt)= SxT(2)*RMag
      nGz(ArcCnt)= SxT(3)*RMag

C     < Find the cosine of the diangle angle.  This is defined by the
C     < dot product of the geodesic unit normal, nG, with the normal
C     < of the plane, ni,
C     <
C     < nG.ni= |nG|*|ni| Cos(Chi)= Cos(Chi)

      Cos_Chi= nGx(ArcCnt)*ni(1) + nGy(ArcCnt)*ni(2)
     >                           + nGz(ArcCnt)*ni(3)

C     < Find out if Cos(Chi) is close to +/-1.

      If( Abs(1.0D0 - Abs(Cos_Chi)) .Lt. Diminutive ) Then

C       < The great plane passing through points A and B is the plane
C       < itself. Note that the formula for the diangle would have
C       < |vi| times the cotangent of Chi. Chi in this case is zero
C       < (or 2Pi) and so Cot(Chi) is infinity.  Therefore, the
C       < diangle is zero in this case.  Return now.

        Return

      End If

C     < Find the cotangent of this angle.  Note that Cot(Chi)'s sign
C     < over 0 to Pi is do strictly to Cos(Chi) thus the sign is
C     < perserved in the formula below.

      Cot_Chi= Cos_Chi/Sqrt(1.0D0 - Cos_Chi*Cos_Chi)

C     < Find the area of the diangle via Gibson's formulae [eqn. A10
C     < of appendix of Gibson and Scheraga, Molecular Physics, 62
C     < (1987) 1247-1265].  We have to consider only type one and type
C     < two diangles (z is always positive with Chi being both acute
C     < or obtuse). AD_2Rsq is the area of diangle over twice the
C     < sphere radius squared.  The quantities are defined as follows
C     <
C     <    R_rhoi := R/rhoi= R/Sqrt(R^2 - z^2)
C     <      vi_R := |vi|/R
C     <   vi_rhoi := |vi|/rhoi= |vi|/Sqrt(R^2 - z^2)
C     <
C     < where z is Gibson's z defined as the distance between the
C     < plane and the geodesic plane and in our terminology is |vi|
C     < and hence Sqrt(R^2 - z^2) is really rhoi the radius of the
C     < inscribed circle.  We break up the formula for the diangle
C     < area into two contributions.  We sum these for all the arcs on
C     < this plane and then we will multiply the second part by
C     < -|vk|/R and add it to the first part.  Eta is equal to the
C     < area of the diangle over twice the sphere radius squared.

      Cosx= R_rhoi*Cos_Chi
      Cosx= Sign( Min( 1.0D0, Abs(Cosx) ), Cosx )
      EtaD(1)= EtaD(1) + ACos( Cosx )

      Cosx= vi_rhoi*Cot_Chi
      Cosx= Sign( Min( 1.0D0, Abs(Cosx) ), Cosx )
      EtaD(2)= EtaD(2) + ACos( Cosx )

      Return
      End


C---------------------------------------------------------------------C
C     Form Polygon                                                    C
C---------------------------------------------------------------------C

      Subroutine Form_Polygon( Area, TrueEdges, TwoPlane )

      Implicit None

C     < PASSED

        Integer*4 TrueEdges

        Logical TwoPlane

        Real*8 Area

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

      Common /This Plane/ ni, vi_Mag, rhoi_sq
     >       , Curvei, O_rhoi_sq, R_rhoi, vi_rhoi
     >       , Edgesi, Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Indexi, Planei, Pointi, PlaneNum, NumEdges

        Integer*4 Edgesi(Max_Apex,2)

        Real*8 vi_Mag, rhoi_sq, Curvei, O_rhoi_sq, R_rhoi, vi_rhoi

        Real*8 ni(3)

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /Sorted Out/ NumStart, NumEnd, LocStart, LocEnd, Edge_Num

        Integer*4 NumStart(-1:3), NumEnd(-1:3)
     >     , LocStart(Max_Apex,-1:3), LocEnd(Max_Apex,-1:3)
     >     , Edge_Num(0:Max_Apex,-1:1)

      Common /Local Arc Info/ Local_List, NumArcs

        Integer*4 NumArcs

        Integer*4 Local_List(Max_Apex,2)

      Common /For Accum/ rA1, axb

        Real*8 rA1(3), axb(3)

C     < LOCAL

        Integer*4 A1, I, EdgeI, Eend, Estart, A1TripNum

C     < The basic algorithm used to find the area of the inscribed
C     < polygon will be the following
C     <
C     < If there are no two-plane apices then
C     <
C     <     o pick a reference 3-plane apex START, call it A1
C     <     o loop through all the REMAINING edges on this plane face
C     <     o find the area of the triangle formed by the reference
C     <       apex A1 and the start and end of the edge
C     <       (A1 Estart Eend)
C     <     o sum all of these areas
C     <
C     <   Since there are no arcs and we simply have an inscribed
C     <   polygon on this plane face. We loop through only the
C     <   REMAINING edges of the plane face comparing the triple plane
C     <   number of the end of each edge with A1's triple plane number.
C     <   If they match then we skip this edge.
C     <

C     < Else if there are two-plane apices then
C     <
C     <     o pick a reference 2-plane apex END, call it A1
C     <     o loop through ALL the edges on this plane face
C     <     o find the area of the triangle formed by the reference
C     <       apex A1 and the start and end of the edge
C     <       (A1 Estart Eend)
C     <     o sum all of these areas
C     <
C     <   We loop through ALL edges, however, we check the end of each
C     <   edge to see whether the edge end is a 2-plane apex and if so
C     <   whether it is A1 or not.  If it is A1 we skip this edge.
C     <
C     <   At the end of the sweep through all edges the total area will
C     <   not contain the areas of the triangles formed by the pseudo
C     <   edges connecting the starts and ends of the arcs.  It will,
C     <   however, contain the area for the arc that A1 starts.  We
C     <   will need to do the following
C     <
C     <     o loop through all the remaining arcs on this plane face
C     <     o find the area of the triangle formed by the reference
C     <       apex A1 and the start and end of each arc
C     <       (A1 Astart Aend)
C     <     o sum all of these areas

C     < Initialize the a x b vector

      axb(1)= 0.0D0
      axb(2)= 0.0D0
      axb(3)= 0.0D0

C     < Find out if our reference apex will be 2-plane or 3-plane

      If( TwoPlane ) Then

C       < Pull out the 2-plane apex end that forms the start of the
C       < first arc we found above.

        A1= Local_List(1,1)

C       < Pull out the coordinates of our reference apex

        rA1(1)= Apexx(A1)
        rA1(2)= Apexy(A1)
        rA1(3)= Apexz(A1)

C       < Loop through ALL the edges

        Do I= 1, TrueEdges

C         < Pull out living (not dead) edge number

          EdgeI= Edge_Num(I,1)

C         < Pull out the end of edge I

          Eend= Edgesi(EdgeI,2)

C         < Is it our reference apex?
C         < Note by definition this edge is alive.

          If( Eend .Ne. A1 ) Then

C           < The end of edge I, Eend is not a 2-plane and is not our
C           < reference apex A1. We must form the area of the triangle
C           < (A1 Estart Eend)

C           < Pull out the start of edge I

            Estart= Edgesi(EdgeI,1)

C           < Accumulate the cross-product of the vector from A1 to the
C           < edge start and the vector from A1 to the edge end.

            Call Accumulate( Estart, Eend )

          End If

        End Do

C       < Loop through the remaining arcs of this plane face.
C       < These will be our pseudo-edges.

        Do I= 2, NumArcs

C         < Accumulate the cross-product of the vector from A1 to the
C         < arc start, Local_List(I,1), and the vector from A1
C         < to the arc end, Local_List(I,2).

          Call Accumulate( Local_List(I,1), Local_List(I,2) )

        End Do

      Else

C       < Pull out the 3-plane apex start that forms the start of the
C       < first living (non-dead) edge in the edge list.  Note that we
C       < know that its not dead because we are selecting it from the
C       < Edge_Num list.

        A1= Edgesi(Edge_Num(1,1),1)

C       < Pull out the coordinates of our reference apex

        rA1(1)= Apexx(A1)
        rA1(2)= Apexy(A1)
        rA1(3)= Apexz(A1)

C       < Pull out the triple plane number for our reference apex
C       < We know that ApexF(A1,1)= +3 because we have no 2-plane
C       < apices and this edge is the first living (non-dead) edge
C       < (i.e., we know that ApexF(A1,1) is not -1).

        A1TripNum= ApexF(A1,2)

C       < Loop through the REMAINING living (non-dead) edges

        Do I= 2, TrueEdges

C         < Pull out living (not dead) edge number

          EdgeI= Edge_Num(I,1)

C         < Pull out the end of edge I

          Eend= Edgesi(EdgeI,2)

C         < Is it our reference apex?
C         < Note we know that this edge is alive (non-dead)

          If(  ApexF(Eend,2) .Ne. A1TripNum ) Then

C           < The end of edge I, Eend does not have the same triple
C           < plane number as the triple plane number of our reference
C           < apex A1.  We must form the area of the triangle
C           < (A1 Estart Eend)

C           < Pull out the start of edge I

            Estart= Edgesi(EdgeI,1)

C           < Accumulate the cross-product of the vector from A1 to the
C           < edge start and the vector from A1 to the edge end.

            Call Accumulate( Estart, Eend )

          End If

        End Do

      End If

C     < Dot the accumulator vector, (a x b), with this plane's unit
C     < normal vector ni.  This will give the sum of the areas of all
C     < the parallelograms on this plane face.  The total area of the
C     < inscribed polygon is the sum of the areas of all the triangles
C     < formed.  Thus the total area is one-half the area of all the
C     < parallelograms formed.

      Area= O_rhoi_sq*( axb(1)*ni(1) + axb(2)*ni(2) + axb(3)*ni(3) )

      Return
      End


C---------------------------------------------------------------------C
C     Accumulate Cross Products                                       C
C---------------------------------------------------------------------C

      Subroutine Accumulate( Estart, Eend )

      Implicit None

      Integer*4 n_Max, Max_Apex

        Parameter ( n_Max= 300, Max_Apex= 20*n_Max )

C     < PASSED

        Integer*4 Estart, Eend

      Common /Apex Info/ Apexx, Apexy, Apexz, ApexF
     >               , n_Edges, Edge_List, ApexCnt

        Integer*4 ApexCnt

        Integer*4 ApexF(Max_Apex,2), n_Edges(n_Max)
     >     , Edge_List(n_Max,Max_Apex,3)

        Real*8 Apexx(Max_Apex), Apexy(Max_Apex), Apexz(Max_Apex)

      Common /For Accum/ rA1, axb

        Real*8 rA1(3), axb(3)

C     < LOCAL

        Real*8 rA(3), rB(3), sxe(3)

C     < Find the vector connecting A1 to Estart, call it rA

      rA(1)= Apexx(Estart) - rA1(1)
      rA(2)= Apexy(Estart) - rA1(2)
      rA(3)= Apexz(Estart) - rA1(3)

C     < Find the vector connecting A1 to Eend, call it rB

      rB(1)= Apexx(Eend) - rA1(1)
      rB(2)= Apexy(Eend) - rA1(2)
      rB(3)= Apexz(Eend) - rA1(3)

C     < Find the cross-product of these two vectors, call it sxe
C     < (for 'start-vector cross-product end-vector')

      Call Cross( rA, rB, sxe )

C     < The area of the parallelogram formed by these two vectors, rA
C     < and rB, is given by the dot of the cross-product with a unit
C     < normal vector pointing in the same direction as the
C     < cross-product.  We have formed the cross-product in a
C     < right-handed sense and crossed the start vector into the end
C     < vector and so the cross-product is normal to the plane face
C     < pointing away from the cutout.  Therefore the plane's unit
C     < normal vector ni points in the same direction as the cross
C     < product, sxe.  Furthermore, this same unit normal vector can be
C     < used for all the cross-products of this plane face and we will
C     < save time by performing the dot product at the end of the loop
C     < on the sum of the cross-products.  Mathematically this means
C     <
C     <    Area of Parallelograms= Sum_{j=1,n} [(sj x ej).ni]
C     <                          = [Sum_{j=1,n} (sj x ej)].ni
C     <                          = (a x b).ni
C     <
C     < where (a x b)= Sum_{j=1,n} (sj x ej).  Note that axb was
C     < initialized earlier.

      axb(1)= axb(1) + sxe(1)
      axb(2)= axb(2) + sxe(2)
      axb(3)= axb(3) + sxe(3)

      Return
      End


C---------------------------------------------------------------------C
C     Festoonery                                                      C
C---------------------------------------------------------------------C

      Subroutine Festoonery( Area_SP , Error_Encountered)

      Implicit None

      Integer*4 n_Max, Max_Arc

        Parameter ( n_Max= 300, Max_Arc= n_Max )

C     < PASSED

        Real*8 Area_SP
        logical Error_Encountered

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Arc Info/ Vertexx, Vertexy, Vertexz
     >       , nGx, nGy, nGz, TanGeo
     >       , Arc_List, Arc_Map, ArcCnt, VertexCnt

        Integer*4 ArcCnt, VertexCnt

        Integer*4 Arc_List(Max_Arc,2), Arc_Map(Max_Arc)

        Real*8 Vertexx(Max_Arc), Vertexy(Max_Arc), Vertexz(Max_Arc)
     >       , nGx(Max_Arc), nGy(Max_Arc), nGz(Max_Arc)
     >       , TanGeo(Max_Arc)

      Error_Encountered = .false.


C     < How many arcs do we have for this sphere?  Note that none of
C     < them will be spherical caps so that an arc is necessariy
C     < connected to another arc.

      If( ArcCnt .Le. 2 ) Then

C       < Initialize the area of the spherical polygon.

        Area_SP= 0.0D0

      Else

C       < There are three or more arcs and so we must bin them into
C       < festoons.  The spherical polygon may be just a simple
C       < spherical triangle.  Note that four arcs could form two
C       < linear spherical polygons and not necessarily a spherical
C       < rectangle.

        Call Bin_Arcs (Error_Encountered)

        if (Error_Encountered) return

C       < Find the area bordered by each spherical polygon.

        Call SpherePoly( Area_SP )

      End If

      Return
      End


C---------------------------------------------------------------------C
C     Bin Arcs Into Festoons                                          C
C---------------------------------------------------------------------C

      Subroutine Bin_Arcs (Error_Encountered)

      Implicit None

      logical Error_Encountered

      Integer*4 n_Max, Max_Arc, Max_Fest

        Parameter ( n_Max= 300, Max_Arc= n_Max, Max_Fest= n_Max/10 )

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Arc Info/ Vertexx, Vertexy, Vertexz
     >       , nGx, nGy, nGz, TanGeo
     >       , Arc_List, Arc_Map, ArcCnt, VertexCnt

        Integer*4 ArcCnt, VertexCnt

        Integer*4 Arc_List(Max_Arc,2), Arc_Map(Max_Arc)

        Real*8 Vertexx(Max_Arc), Vertexy(Max_Arc), Vertexz(Max_Arc)
     >       , nGx(Max_Arc), nGy(Max_Arc), nGz(Max_Arc)
     >       , TanGeo(Max_Arc)

      Common /Festoon Info/ Festoon, FestArcs, FestoonCnt

        Integer*4 FestoonCnt

        Integer*4 Festoon(Max_Arc, Max_Fest), FestArcs(Max_Fest)

C     < LOCAL

        Integer*4 UsedArcs, ArcRef, Arc, Cnt

        Logical NotUsed(Max_Arc)

      external getpid
      integer getpid

      Error_Encountered = .false.

C     < Initialize some variables.

      UsedArcs= 0
      ArcRef= 0
      FestoonCnt= 0

      Do Arc= 1, ArcCnt
        NotUsed(Arc)= .True.
      End Do

C     < We have to put the arcs into festoon bins.  It takes three or
C     < more arcs to make a festoon.  We have to analyze the
C     < connectivity of the arcs for this sphere to determine the
C     < festoons.

      Do While( UsedArcs .Lt. ArcCnt )

C       < Define the reference arc starting this festoon.

        ArcRef= ArcRef + 1

C       < Has this been used?

        If( NotUsed(ArcRef) ) Then

C         < Increment the festoon count.

          FestoonCnt= FestoonCnt + 1

C         < Initialize the arc count for this festoon.

          Cnt= 1

C         < Store the first arc.

          Festoon(Cnt,FestoonCnt)= ArcRef

C         < Find the arc that reference arc's terminus starts.

          Arc= Arc_Map( Arc_List(ArcRef,2) )

C         < Tag this arc as used.

          NotUsed(Arc)= .False.

C         < Connect the arcs until the reference arc appears again.

          Do While( Arc .Ne. ArcRef )

C           < Increment the arc count for this festoon.

            Cnt= Cnt + 1

C           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C           < Test to see if we are endlessly looping.

C           < See note in routine Build_AP concerning two-plane versus
C           < three-plane apices.
C           < LRD - Saturday, February 08, 1992, 10:14:54 PST

            If( Cnt .Gt. Max_Arc ) Then
              Write(0,*) ' Perturb some coordinates '
              Error_Encountered = .true.
              return
            End If
C           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C           < Store the arc.

            Festoon(Cnt,FestoonCnt)= Arc

C           < Find the arc that the last arc's terminus starts.

            Arc= Arc_Map( Arc_List(Arc,2) )

C           < Tag this arc as used.

            NotUsed(Arc)= .False.

          End Do

C         < Increment the number of arcs used.

          UsedArcs= UsedArcs + Cnt

C         < Store the number of arcs belonging to this festoon.

          FestArcs(FestoonCnt)= Cnt

        End If

      End Do

      Return
      End


C---------------------------------------------------------------------C
C     Find The Area Of All The Spherical Polygons                     C
C---------------------------------------------------------------------C

      Subroutine SpherePoly( Area_SP )

      Implicit None

      Integer*4 n_Max, Max_Arc, Max_Fest

        Parameter ( n_Max= 300, Max_Arc= n_Max, Max_Fest= n_Max/10 )

      Real*8 Diminutive

        Parameter ( Diminutive= 1.0D-8 )

      Real*8 Pi, O_Pi

        Parameter ( Pi= 3.14159265358979323846D0, O_Pi= 1.0D0/Pi )

      Real*8 Sqrt2

        Parameter ( Sqrt2= 1.41421356237309504880D0 )

C     < PASSED

        Real*8 Area_SP

      Common /Sphere Info/ Radius, Rsq, O_Radius, O_Rsq, O_Rcu
     >                   , Sphere, NumPlanes
     >                   , TripCnt

        Integer*4 Sphere, NumPlanes, TripCnt

        Real*8 Radius, Rsq, O_Radius, O_Rsq, O_Rcu

      Common /Arc Info/ Vertexx, Vertexy, Vertexz
     >       , nGx, nGy, nGz, TanGeo
     >       , Arc_List, Arc_Map, ArcCnt, VertexCnt

        Integer*4 ArcCnt, VertexCnt

        Integer*4 Arc_List(Max_Arc,2), Arc_Map(Max_Arc)

        Real*8 Vertexx(Max_Arc), Vertexy(Max_Arc), Vertexz(Max_Arc)
     >       , nGx(Max_Arc), nGy(Max_Arc), nGz(Max_Arc)
     >       , TanGeo(Max_Arc)

      Common /Festoon Info/ Festoon, FestArcs, FestoonCnt

        Integer*4 FestoonCnt

        Integer*4 Festoon(Max_Arc, Max_Fest), FestArcs(Max_Fest)

C     < LOCAL

        Integer*4 Fest, ArcRef, Cnt, Arc, Vert, ArcLast, NumArcs

        Logical NotDone

        Real*8 TanOld, TanArc, Denom, Cos_Theta, TanDia, AreaET, ALune
     >       , nG_dollar_x, TanLast, RMag, Cos_Omega

        Real*8 x(3), xT(3), nG(3)

C     < The method for finding the area of the spherical polygon is
C     < analogous to finding the area of a polygon on a plane. There
C     < we pick a vertex of the polygon and find the area of triangles
C     < formed by each adjacent pair of vertices of the polygon and
C     < the reference vertex.  So for a polygon with n sides the area
C     < is defined as a sum of n-2 triangles.  For a spherical polygon
C     < we do the same thing except now the triangles are spherical
C     < triangles and in particular Euler triangles.

C     < Initialize the area of the spherical polygon.

      Area_SP= 0.0D0

C     < Single Loop over all festoons (i.e., a single loop over all
C     < the spherical polygons on this sphere).

      Do Fest= 1, FestoonCnt

C       < Pull out the number of arcs for this festoon.

        NumArcs= FestArcs(Fest)

C       < The festoon has to contain at least three arcs in order to
C       < be non-trivial. If it contains two arcs then the spherical
C       < polygon has degenerated down to a line. It can not contain
C       < only one arc because this would be a spherical segment and
C       < would have be considered earlier.

        If( NumArcs .Gt. 2 ) Then

C         < Pull out the arc number of the first arc for this festoon.
C         < This first arc will be our reference arc.

          ArcRef= Festoon(1,Fest)

C         < Pull out the square tangent of the geodesic quarter-angle
C         < for this reference arc. We will disguise it as diagonal's.

          TanOld= TanGeo(ArcRef)

C         < Pull out the vertex number of the vertex starting this
C         < reference arc.

          Vert= Arc_List(ArcRef,1)

C         < Pull out the vector connecting the sphere center to the
C         < vertex starting this reference arc, call this vector x.

          x(1)= Vertexx(Vert)
          x(2)= Vertexy(Vert)
          x(3)= Vertexz(Vert)

C         < Do a single loop over all the remaining arcs of this
C         < festoon but two. Note that if the total number of arcs for
C         < this festoon is three then this loop will not be executed.

          Cnt= 1
          NotDone= .True.

          Do While( Cnt .Lt. NumArcs - 2 )

C           < Increment counter.

            Cnt= Cnt + 1

C           < Pull out the arc number of arc Cnt.

            Arc= Festoon(Cnt,Fest)

C           < Pull out the square tangent of the geodesic quarter-angle
C           < for the present arc Cnt.

            TanArc= TanGeo(Arc)

C           < Pull out the number of the terminating vertex of arc Cnt.

            Vert= Arc_List(Arc,2)

C           < Pull out the vector connecting the sphere center to the
C           < end of arc Arc.  Call it xT; T is for terminus.

            xT(1)= Vertexx(Vert)
            xT(2)= Vertexy(Vert)
            xT(3)= Vertexz(Vert)

C           < Find the square tangent of the geodesic quarter-angle
C           < between the vector connecting the sphere center to the
C           < terminus of this arc (vertex Vert) and the vector for
C           < the reference arc, x.  This is the diagonal of the Euler
C           < triangle to be considered.  First, find the cosine of
C           < the angle between these vectors.  Both vectors are of
C           < magnitude R, the radius of the sphere. Therefore,
C           < Cos(Theta)= x.xT/(R*R)

            Cos_Theta= (x(1)*xT(1) + x(2)*xT(2) + x(3)*xT(3))*O_Rsq

C           < Is Cos(Theta) close to -1?

            If( Abs(1.0D0 + Cos_Theta) .Lt. Diminutive ) Then

C             < The diagonal of this triangle is exactly Pi radians
C             < long.  This means that the lengths of the other two
C             < sides (one or both being arcs of the spherical
C             < polygon) sum to Pi radians also.  Furthermore, the
C             < next triangle will also contain this diagonal so the
C             < these two triangles have reduced to a spherical lune.
C             < Pull out the unit normal geodesic vector for this arc
C             < and dot it with the unit normal geodesic vector for
C             < the next arc that will share this diagonal.

              Call Lune( xT, nG_dollar_x, ALune, Cnt, Fest, Arc )

C             < Sum the area adding or subtracting the spherical
C             < lune's area when appropriate.

              Area_SP= Area_SP - Sign(ALune,nG_dollar_x)

C             < The variable Cnt has been incremented in Lune because
C             < we have used two arcs and the reference vertex to form
C             < the lune.  If Cnt is still less or equal to than
C             < (NumArcs - 2) then we need to get the geodesic angle
C             < for the diagonal.  Otherwise we are done with this
C             < spherical polygon.

              If( Cnt .Le. NumArcs - 2 ) Then

C               < In order to find the area of the next spherical
C               < triangle, which may be the last, we need the square
C               < of the tangent of the geodesic quarter-angle for the
C               < diagonal formed by the reference vertex and the
C               < terminus of the second arc.  Actually, since the
C               < diagonal is of length Pi, the second arc and this
C               < diagonal have to sum to a length of Pi.  This could
C               < be used as a short cut. Find the Cos(Theta).

                Cos_Theta= (x(1)*xT(1) + x(2)*xT(2) + x(3)*xT(3))*O_Rsq

C               < The square tangent of the quarter angle is related to
C               < the cosine angle, assign it directly to TanOld.

                Denom= Sqrt2 + Sqrt(1.0D0 + Cos_Theta)
                TanOld= (1.0D0 - Cos_Theta)/(Denom*Denom)

              Else

C               < This spherical polygon is complete.  We will need to
C               < skip over the last triangle normally done outside
C               < this loop because we have already incorporated it
C               < into the lune.

                NotDone= .False.

              End If

            Else

C             < The square tangent of the quarter angle is related to
C             < the cosine angle via,
C             <
C             <                            1 - Cos(Theta)
C             <   Tan^2(Theta/4)= ----------------------------------
C             <                   [Sqrt(2) + Sqrt(1 + Cos(Theta))]^2
C             <

              Denom= Sqrt2 + Sqrt(1.0D0 + Cos_Theta)
              TanDia= (1.0D0 - Cos_Theta)/(Denom*Denom)

C             < Find out if the square of the tangent of the geodesic
C             < quarter-angle of the arc is unity.

              If( Abs(1.0D0 - TanArc) .Lt. Diminutive ) Then

C               < If so then the arc is of length Pi and the other two
C               < sides sum to Pi also. This means that the triangle
C               < reduces to a lune with the arc forming one side. The
C               < other side is formed by the diagonals going from the
C               < arc end to the arc start through the reference
C               < vertex. This is also a geodesic plane. We must find
C               < the angle between these two planes. First find the
C               < geodesic unit normal vector for the geodesic going
C               < from the end of the arc Arc through the reference
C               < vertex to the start of the arc Arc.
C               < Find the cross product of xT into x the reference
C               < vertex. Call it nG= (xT) x (x).

                Call Cross( xT, x, nG )

C               < Find the magnitude of nG

                RMag= 1.0D0/Sqrt(nG(1)*nG(1) + nG(2)*nG(2)
     >                                       + nG(3)*nG(3))

C               < Find the dot product between this geodesic normal
C               < and the geodesic unit normal for arc Arc.  Multiply
C               < by the magnitude RMag. This then is the cosine of
C               < the complement of the angle between the two geodesic
C               < planes.  The angle we wish, call it alpha, is Pi/2
C               < minus the angle between the two normal vectors,
C               < cos(w)= -cos(alpha)

                Cos_Omega= (nGx(Arc)*nG(1) + nGy(Arc)*nG(2)
     >                                     + nGz(Arc)*nG(3))*RMag

C               < The area of a lune is given as A= 2R^2 Alpha. Thus,
C               < the area of the lune divided by 4R^2 is Alpha/2. We
C               < will assign this to the variable AreaET since a lune
C               < is simply a special case of an Euler triangle.

                AreaET= ACos(-Cos_Omega)*0.50D0

                Write(0,*) ' Lune: Side '

              Else

C               < Find the area of this Euler triangle. The area is
C               < reduced by 4 Pi R^2.

                Call Euler( TanOld, TanArc, TanDia, AreaET )

              End If

C             < Find the directionality of this triangle.  The
C             < directionality of the triangle represents its location
C             < relative to the cutout.  If it is anti-clockwise then
C             < the triangle is out of the cutout and so the area of
C             < this area should be subtracted from that of a sphere.
C             < Otherwise the area should be added to find the area of
C             < the cutout.  We will establish directionality by
C             < considering the three points on the surface of the
C             < sphere making up the triangle and the sphere center as
C             < a tetrahedron.  The directionality of the base of this
C             < tetrahedron is found by finding the triple product
C             < (rAxrB).rC where rA is the vector connecting the
C             < sphere center to the point A on the sphere surface.
C             < If B and C are both on an arc then the unit geodesic
C             < normal, nG, is equal to rBxrC/|rBxrC| and we need
C             < merely dot this with the third point, rA.  In this
C             < case rA is the reference vertex of the reference arc.
C             < So find nG(Arc).x and if this is positive then it
C             < means that the triangle is directed in an
C             < anti-clockwise manner and we must subtract this area.

              nG_dollar_x= nGx(Arc)*x(1) + nGy(Arc)*x(2) + nGz(Arc)*x(3)

C             < Accumulate the area of the triangle to the area of the
C             < spherical polygon adding or subtracting the when
C             < appropriate.

              Area_SP= Area_SP - Sign(AreaET,nG_dollar_x)

C             < Define the square tangent of the diagonal as the old
C             < square tangent.

              TanOld= TanDia

            End If

          End Do

C         < Check if a lune has already used the last two arcs of the
C         < spherical polygon.

          If( NotDone ) Then

C           < At this point only the last two arcs of the spherical
C           < polygon have not been used yet.  Pull out the number of
C           < the penultimate arc.

            Arc= Festoon(NumArcs - 1,Fest)

C           < Pull out the square tangent of the geodesic quarter-angle
C           < for the penultimate arc.

            TanArc= TanGeo(Arc)

C           < Pull out the number of the final arc.

            ArcLast= Festoon(NumArcs,Fest)

C           < Pull out the square tangent of the geodesic quarter-angle
C           < for the final arc.

            TanLast= TanGeo(ArcLast)

C           < Find out if the square of the tangent of the geodesic
C           < quarter-angle of the penultimate arc is unity.

            If( Abs(1.0D0 - TanArc) .Lt. Diminutive .Or.
     >          Abs(1.0D0 - TanLast) .Lt. Diminutive ) Then

C             < Either the last or the penultimate arc is of length
C             < Pi.  This means that the reference vertex is on the
C             < geodesic plane of the other arc.  Thus the Euler
C             < triangle has reduced to a lune. The angle of the lune
C             < is found by taking the dot product of the unit normal
C             < geodesic vectors of each plane and then adding or
C             < subtracting Pi/2.

              Cos_Omega= nGx(Arc)*nGx(ArcLast) + nGy(Arc)*nGy(ArcLast)
     >                                         + nGz(Arc)*nGz(ArcLast)

C             < The area of a lune is given as A= 2R^2 Alpha. Thus,
C             < the area of the lune divided by 4R^2 is Alpha/2. We
C             < will assign this to the variable AreaET since a lune
C             < is simply a special case of an Euler triangle.

              AreaET= ACos(-Cos_Omega)*0.50D0

              Write(0,*) ' Lune: Side '

            Else

C             < Find the area of this last Euler triangle for this
C             < spherical polygon.

              Call Euler( TanOld, TanArc, TanLast, AreaET )

            End If

C           < Find the directionality of this triangle.

            nG_dollar_x= nGx(Arc)*x(1) +  nGy(Arc)*x(2) + nGz(Arc)*x(3)

C           < Sum the area adding or subtracting the Euler triangle's
C           < area when appropriate.

            Area_SP= Area_SP - Sign(AreaET,nG_dollar_x)

          End If

        End If

      End Do

C     < Divide the reduced spherical polygon area by Pi making it
C     < reduced by 4 Pi R^2 the surface area of a sphere.

      Area_SP= Area_SP*O_Pi

      Return
      End


C---------------------------------------------------------------------C
C     Find The Area Of An Euler Triangle                              C
C---------------------------------------------------------------------C

      Subroutine Euler( Asq, Bsq, Csq, Area )

      Implicit None

C     < PASSED

        Real*8 Asq, Bsq, Csq, Area

C     < LOCAL

        Real*8 Sum, ABsq, Group, x1sq, x2sq, x3sq, SqrtArg

C     < Define some useful quantities.

      Sum= Asq + Bsq
      ABsq= Asq*Bsq
      Group= (1.0D0 + ABsq)

C     < Find the three squares.

      x1sq= Sum - Csq*Group
      x1sq= x1sq*x1sq

      x2sq= Group - Csq*Sum
      x2sq= x2sq*x2sq

      x3sq= (1.0D0 + Csq)
      x3sq= 4.0D0*ABsq*x3sq*x3sq

C     < The square of the tangent of the quarter-angle for the
C     < spherical excess is given by Tan^2(Epsilon/4). Epsilon is the
C     < spherical excess defined as the sum of the three angles of the
C     < Euler triangle minus Pi radians.  The sum of the angles of an
C     < Euler triangle must be greater than Pi but less than 3Pi thus
C     < the spherical excess must be between zero and 2Pi. The reduced
C     < area is given as A/4R^2= Epsilon/4, which is less than Pi/2.

      SqrtArg= Max( 0.0D0, (x3sq - x1sq)/(x2sq - x3sq) )
      Area= ATan(Sqrt(SqrtArg))

      Return
      End


C---------------------------------------------------------------------C
C     Lune                                                            C
C---------------------------------------------------------------------C

      Subroutine Lune( xT, nG_dollar_xT, Area, Cnt, Fest, Arc )

      Implicit None

      Integer*4 n_Max, Max_Arc, Max_Fest

        Parameter ( n_Max= 300, Max_Arc= n_Max, Max_Fest= n_Max/10 )

C     < PASSED

        Integer*4 Cnt, Fest, Arc

        Real*8 nG_dollar_xT, Area

        Real*8 xT(3)

      Common /Arc Info/ Vertexx, Vertexy, Vertexz
     >       , nGx, nGy, nGz, TanGeo
     >       , Arc_List, Arc_Map, ArcCnt, VertexCnt

        Integer*4 ArcCnt, VertexCnt

        Integer*4 Arc_List(Max_Arc,2), Arc_Map(Max_Arc)

        Real*8 Vertexx(Max_Arc), Vertexy(Max_Arc), Vertexz(Max_Arc)
     >       , nGx(Max_Arc), nGy(Max_Arc), nGz(Max_Arc)
     >       , TanGeo(Max_Arc)

      Common /Festoon Info/ Festoon, FestArcs, FestoonCnt

        Integer*4 FestoonCnt

        Integer*4 Festoon(Max_Arc, Max_Fest), FestArcs(Max_Fest)

C     < LOCAL

        Integer*4 Vert

        Real*8 Cos_Omega

        Real*8 nG(3)

      Write(0,*) ' Lune: Diagonal '

C     < Pull out the unit normal geodesic vector

      nG(1)= nGx(Arc)
      nG(2)= nGy(Arc)
      nG(3)= nGz(Arc)

C     < Increment the arc number.

      Cnt= Cnt + 1

C     < Pull out the arc number of arc Cnt.

      Arc= Festoon(Cnt,Fest)

C     < Find the cosine of the angle between nG and the unit
C     < geodesic norm for the second arc, cos(w)= nG.nG(Arc).

      Cos_Omega= nG(1)*nGx(Arc) + nG(2)*nGy(Arc) + nG(3)*nGz(Arc)

C     < This angle is the complement of the angle between the two
C     < geodesic planes: one forming the first arc and the other
C     < forming the second arc. Find the angle between the planes.
C     < The area of a lune is defined as A= 2R^2 Alpha so the area
C     < reduced by 4R^2 is Alpha/2.

      Area= ACos(-Cos_Omega)*0.50D0

C     < Pull out the number of the terminating vertex of this
C     < second arc Cnt.

      Vert= Arc_List(Arc,2)

C     < Pull out the vector connecting the sphere center to the
C     < end of arc Arc.  Call it xT; T is for terminus.

      xT(1)= Vertexx(Vert)
      xT(2)= Vertexy(Vert)
      xT(3)= Vertexz(Vert)

C     < The area of this lune should be subtracted from total
C     < spherical polygon area if it directed in an anti-clockwise
C     < manner on the sphere surface and added if directed in a
C     < clockwise manner.  The directionality is determined by dot
C     < product of geodesic unit normal of the first arc with the
C     < vector connecting the sphere center with the terminus of the
C     < second arc, Vert. Note that we can not use the vector x
C     < connecting the sphere center to the reference vertex as we
C     < normally do since this lies on the geodesic plane of the first
C     < arc and the second arc.

      nG_dollar_xT= nG(1)*xT(1) + nG(2)*xT(2) + nG(3)*xT(3)

      Return
      End


C=====================================================================C
C                          UTILITY ROUTINES                           C
C=====================================================================C

C---------------------------------------------------------------------C
C     Cross Product Of Two 3x1 Vectors                                C
C---------------------------------------------------------------------C

      Subroutine Cross( a, b, axb )

      Implicit None

C     < PASSED

        Real*8 a(3), b(3), axb(3)

      axb(1)= a(2)*b(3) - a(3)*b(2)
      axb(2)= a(3)*b(1) - a(1)*b(3)
      axb(3)= a(1)*b(2) - a(2)*b(1)

      Return
      End


C---------------------------------------------------------------------C
C     Sort                                                            C
C---------------------------------------------------------------------C

      Subroutine Sort( A, Index, N )

      Implicit None

C     < PASSED

        Integer*4 N

        Integer*4 Index(N)

        Real*8 A(N)

C     < LOCAL

        Integer*4 I, J, HoldI

        Real*8 Hold

C     < Initialize the index list

      Do I= 1, N
        Index(I)= I
      End Do

C     < Sort A(Index) into ascending order

      Do J= 2, N

C       < Pull out Jth value

        Hold= A(Index(J))
        HoldI= Index(J)

        I= J - 1

        Do While( I .Gt. 0 .And. A(Index(I)) .Gt. Hold )
          Index(I+1)= Index(I)
          I= I - 1
        End Do

        Index(I+1)= HoldI

      End Do

      Return
      End
