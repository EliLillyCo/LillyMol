      INTEGER FUNCTION IDAMAX(N,SX,INCX)
      implicit none
C***BEGIN PROLOGUE  ISAMAX
C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C***END PROLOGUE  ISAMAX
C
      double precision SX(*),SMAX,XMAG
      integer n
      integer incx

      integer i, ii, ns

C***FIRST EXECUTABLE STATEMENT  ISAMAX
      IDAMAX = 0
      IF(N.LE.0) RETURN
      IDAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      SMAX = ABS(SX(1))
      NS = N*INCX
      II = 1
          DO 10 I=1,NS,INCX
          XMAG = ABS(SX(I))
          IF(XMAG.LE.SMAX) GO TO 5
          IDAMAX = II
          SMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N 
         XMAG = ABS(SX(I))
         IF(XMAG.LE.SMAX) GO TO 30
         IDAMAX = I 
         SMAX = XMAG
   30 CONTINUE
      RETURN
      END 
