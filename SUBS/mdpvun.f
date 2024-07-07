C
C mdpvun
C Make Double Precision Vector Unit Norm.
C
C We have a double precision vector DV with N elements
C
C We want it such that the norm calculated with the
C BLAS routine DNRM2 is 1.0d0
C Returns leaving vector unchanged if it is too small
C
      SUBROUTINE MDPVUN( N, DV )
      IMPLICIT NONE
C
      INTEGER     N
      REAL*8      DV( N )
C
      INTEGER     I
      REAL*8      DSCALE
      REAL*8      DLOW
      PARAMETER ( DLOW = 1.0d-9 )
      INTEGER     INCX
      PARAMETER ( INCX = 1      )
C
      REAL*8      DNRM2
C
      DSCALE = DNRM2( N, DV, INCX )
      IF ( DSCALE.LT.DLOW ) RETURN
C
      DSCALE = 1.0d0/DSCALE
C
      DO I = 1, N
        DV( I ) = DV( I )*DSCALE
      ENDDO
C
      RETURN
      END
C
