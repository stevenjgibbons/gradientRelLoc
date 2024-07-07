C
C Gradient Relative Location Solve
C Steve Gibbons 2024/07/07
C
C We have a maximum of NRELVM - observations of
C DTIME1() - an epoch time for event 1
C DTIME2() - an epoch time for event 2
C DOMPC2() - a correction term for event 2 such that
C            the actual time used is DTIME2(i) + DOMPC2(i).
C            i.e. DOMPC2 negative makes the time earlier.
C NRELVL is the actual number of observations.
C DSXVLS() - The x-slowness vector for relative time i (s/km)
C DSYVLS() - The y-slowness vector for relative time i (s/km)
C DCVARR() - The weight for relative time i (e.g. correlation coeff.)
C
C DISTMX   - is a time in seconds.
C            We calculate the median time, DMEDIA, of all the
C              DTIME2 + DOMPC2 and any i with
C               DABS( DTIME2 + DOMPC2 ).gt.DISTMX
C            has DCVARR (weight) set to almost zero.
C
C We return DELXKM and DELYKM and the vector DRESVL()
C which should contain the time residual for each of the
C observations - and should hopefully be able to be used as
C input in the array DOMPC2 for a new calculation to give a minimal
C residual.
C DRESNM is a scalar giving the norm.
C
C The following arrays are all work arrays.
C
C  DDVEC( NRELVM )
C  DWORK1( 2*NRELVM )
C  DAMAT( NRELVM, 3 )
C  DTEMP1( NRELVM )
C  DOLDRV( NRELVM )
C
      SUBROUTINE GRDRLS( IERR, NRELVM, NRELVL, DTIME1, DTIME2,
     1                   DOMPC2, DSXVLS, DSYVLS, DCVARR, DRESVL,
     2                   DISTMX, DELXKM, DELYKM, DRESNM,
     3                   DDVEC, DWORK1, DAMAT, DTEMP1, DOLDRV )
      IMPLICIT NONE
C
      INTEGER IERR
      INTEGER NRELVM
      INTEGER NRELVL
      REAL*8  DTIME1( NRELVM )
      REAL*8  DTIME2( NRELVM )
      REAL*8  DOMPC2( NRELVM )
      REAL*8  DSXVLS( NRELVM )
      REAL*8  DSYVLS( NRELVM )
      REAL*8  DCVARR( NRELVM )
      REAL*8  DRESVL( NRELVM )
C
      REAL*8  DISTMX
      REAL*8  DELXKM
      REAL*8  DELYKM
      REAL*8  DRESNM
C
      REAL*8  DDVEC( NRELVM )
      REAL*8  DWORK1( 2*NRELVM )
      REAL*8  DAMAT( NRELVM, 3 )
      REAL*8  DTEMP1( NRELVM )
      REAL*8  DOLDRV( NRELVM )
C
C Working variables
C
      INTEGER     LWORK
      INTEGER     LWOPT
      INTEGER     NMMAX
      PARAMETER ( NMMAX = 3 )
      INTEGER     NM
      PARAMETER ( NM    = 3 )
      INTEGER     INCX
      PARAMETER ( INCX  = 1 )
      REAL*8      DWORK2( NMMAX, NM )
      REAL*8      DWORK3( NM )
      REAL*8      DMVEC( NM )
      INTEGER     IWORK1( NM )
      INTEGER     IC
      INTEGER     IRELVL
      INTEGER     MXITER
      PARAMETER ( MXITER = 10000 )
      REAL*8      DTOL
      PARAMETER ( DTOL = 1.0d-6 )
      REAL*8      DMEDIA
C external routine
      REAL*8      DNRM2
C
      IERR   = 0
C
      LWORK  = 2*NRELVM
C
      IF ( NRELVL.LT.4 .OR. NRELVL.GT.NRELVM ) THEN
        WRITE (6,*) 'Subroutine GRDRLS: error '
        WRITE (6,*) 'NRELVL  = ', NRELVL
        WRITE (6,*) 'NRELVM  = ', NRELVM
        IERR   = 1
        RETURN
      ENDIF
C
C  First zero the DAMAT array
C
      DO IC = 1, 3
        DO IRELVL = 1, NRELVL
          DAMAT( IRELVL, IC ) = 0.0d0
        ENDDO
      ENDDO
C
C Now enter the slowness values into DAMAT
C
      DO IRELVL = 1, NRELVL
        DAMAT( IRELVL, 1 ) = DSXVLS( IRELVL )
        DAMAT( IRELVL, 2 ) = DSYVLS( IRELVL )
        DAMAT( IRELVL, 3 ) =   1.0d0
      ENDDO
C
C First we make a mock-up of the RHS vector in order
C to calculate the median. (The RHS vector will be destroyed
C in the process so we will need to build it again.)
C
      DO IRELVL = 1, NRELVL
        DDVEC( IRELVL ) = DTIME1( IRELVL ) - DTIME2( IRELVL ) -
     1                      DOMPC2( IRELVL )
      ENDDO
C
      CALL DMEDVL( IERR, NRELVL, DDVEC, DMEDIA )
C
C DMEDIA now contains the median value. Reconstruct DDVAL
C
      DO IRELVL = 1, NRELVL
        DDVEC( IRELVL ) = DTIME1( IRELVL ) - DTIME2( IRELVL ) -
     1                      DOMPC2( IRELVL ) - DMEDIA
C       .
C       . down-weight outliers
C       .
        IF ( DABS( DDVEC( IRELVL ) ).GT.DISTMX ) THEN
          DCVARR( IRELVL ) = 0.000001d0
        ENDIF
      ENDDO
C
C Now solve the system.
C
      CALL IRWMPS( IERR, NRELVM, NRELVL, NM, NM,
     1             DDVEC, DAMAT, DCVARR, DMVEC, DRESVL,
     2             LWORK, LWOPT, DWORK1, DWORK2, DWORK3,
     3             IWORK1, MXITER, DTOL, DTEMP1, DOLDRV )
      IF ( IERR.NE.0 ) THEN
        WRITE (6,*) 'Subroutine GRDRLS: error '
        WRITE (6,*) 'IRWMPS returned IERR = ', IERR
        IERR   = 1
        RETURN
      ENDIF
C
      DELXKM = DMVEC( 1 )
      DELYKM = DMVEC( 2 )
      DRESNM = DNRM2( NRELVL, DRESVL, INCX )
C
      RETURN
      END
C
