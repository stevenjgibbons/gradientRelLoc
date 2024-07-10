C
C mcgrrl
C Monte Carlo Gradient Relative Location.
C
C Want a "Monte Carlo" version of the program
C which reads in altenative relative time measurements into
C "banks" of measurements.
C 
C The integer array
C 
C NTMEAS( NSPMAX )
C 
C stores how many measurements we have for each ISP.
C 
C NLSP is the number of "live" station and phase pairs -
C that is to say ISP with NTMEAS( ISP ).gt.0
C 
C NMEASM is the highest possible value of NTMEAS( ISP )
C 
C DT1BNK( NSPMAX, NMEASM ) contains the event 1 time measurements for each ISP
C DT2BNK( NSPMAX, NMEASM ) contains the event 2 time measurements for each ISP
C DCCBNK( NSPMAX, NMEASM ) contains the quality measures (CC)
C                            for each measurement
C 
C Now we want to perform NREAL realizations, for which NREALM is
C the maximum possible.
C 
C For each realization, IREAL, we loop around the NLSP station/phase pairs
C and and for each such ISP we pick at random a value IMEAS between
C 1 and NTMEAS( ISP ) and set
C (Note there is one "irelv" per "live" ISP but some specified
C  station/phase pairs may not have measurements and so are omitted.)
C DTIME1( irelv ) = DT1BNK( ISP, IMEAS )
C DTIME2( irelv ) = DT2BNK( ISP, IMEAS )
C DCVARR( irelv ) = DT2BNK( ISP, IMEAS )
C 
C We then call
C GRDRLS
C for this particular combination and record the output
C DELXKM, DELYKM, DRESNM
C 
C DELXUT( NREALM ) = DELXKM
C DELYUT( NREALM ) = DELYKM
C DRESUT( NREALM ) = DRESNM
C 
C We also store the times and the residuals in the arrays
C 
C DT1USD( NSPMAX, NREALM )
C DT2USD( NSPMAX, NREALM )
C DCCUSD( NSPMAX, NREALM )
C DRESTM( NSPMAX, NREALM )
C 
      SUBROUTINE MCGRRL( IERR, NSPMAX, NMEASM, NSP, NLSP, NTMEAS,
     1                   DSXARR, DSYARR, DSXVLS, DSYVLS,
     2                   DTIME1, DTIME2, DT1BNK, DT2BNK, DCCBNK,
     3                   DOMPC2, DCVARR, DRESVL, DISTMX,
     4                   DDVEC, DWORK1, DAMAT, DTEMP1, DOLDRV,
     5                   NREALM, NREAL, DELXUT, DELYUT, DRESUT,
     6                   DT1USD, DT2USD, DCCUSD, DRESTM )
      IMPLICIT NONE
C
      INTEGER IERR
      INTEGER NSPMAX
      INTEGER NMEASM
      INTEGER NSP
      INTEGER NLSP
      INTEGER NTMEAS( NSPMAX )
      REAL*8  DSXARR( NSPMAX )
      REAL*8  DSYARR( NSPMAX )
      REAL*8  DSXVLS( NSPMAX )
      REAL*8  DSYVLS( NSPMAX )
      REAL*8  DT1BNK( NSPMAX, NMEASM )
      REAL*8  DT2BNK( NSPMAX, NMEASM )
      REAL*8  DCCBNK( NSPMAX, NMEASM )
C note that the following arrays are limited here by NSPMAX
      REAL*8  DTIME1( NSPMAX )
      REAL*8  DTIME2( NSPMAX )
      REAL*8  DOMPC2( NSPMAX )
      REAL*8  DCVARR( NSPMAX )
      REAL*8  DRESVL( NSPMAX )
C
      REAL*8  DISTMX
C
      REAL*8  DDVEC( NSPMAX )
      REAL*8  DWORK1( 2*NSPMAX )
      REAL*8  DAMAT( NSPMAX, 3 )
      REAL*8  DTEMP1( NSPMAX )
      REAL*8  DOLDRV( NSPMAX )
C
      INTEGER NREALM
      INTEGER NREAL
C
      REAL*8  DELXUT( NREALM )
      REAL*8  DELYUT( NREALM )
      REAL*8  DRESUT( NREALM )
C
      REAL*8  DT1USD( NSPMAX, NREALM )
      REAL*8  DT2USD( NSPMAX, NREALM )
      REAL*8  DCCUSD( NSPMAX, NREALM )
      REAL*8  DRESTM( NSPMAX, NREALM )
C
C Working variables
C
      INTEGER NRELVM
      INTEGER NRELVL
      INTEGER NMEAS
      INTEGER IREAL
      INTEGER ISP
      INTEGER ILSP
      INTEGER IMEAS
C
      INTEGER IDUM
      REAL*8  DRAN1N
      REAL*8  DVAL01
C
      REAL*8  DELXKM
      REAL*8  DELYKM
      REAL*8  DRESNM
C
      IERR   = 0
      NRELVM = NSPMAX
      NRELVL = NLSP
C     . random number seed
      IDUM   = NLSP + 3781
C
      DO IREAL = 1, NREAL
C       .
C       . We do a single realization of our relative times.
C       . Need to loop around ISP = 1, NSP
C       . and increment ILSP when we find a station/phase
C       . with measurements.
C       .
        ILSP   = 0
        DO ISP = 1, NSP
          NMEAS  = NTMEAS( ISP )
          IF ( NMEAS.GT.0 ) THEN
            ILSP   = ILSP + 1
            DVAL01 = DRAN1N( IDUM )
            IMEAS  = INT( DBLE( NMEAS )*DVAL01 ) + 1
            IF ( IMEAS.LT.1     ) IMEAS = 1
            IF ( IMEAS.GT.NMEAS ) IMEAS = NMEAS
            DSXVLS( ILSP ) = DSXARR( ISP )
            DSYVLS( ILSP ) = DSYARR( ISP )
            DTIME1( ILSP ) = DT1BNK( ISP, IMEAS )
            DTIME2( ILSP ) = DT2BNK( ISP, IMEAS )
            DOMPC2( ILSP ) = 0.0d0
            DCVARR( ILSP ) = DCCBNK( ISP, IMEAS )
          ENDIF
        ENDDO
C       .
        IF ( ILSP.NE.NLSP ) THEN
          WRITE (6,*) 'Subroutine MCGRRL: Error'
          WRITE (6,*) 'Expected ILSP to be ', NLSP
          WRITE (6,*) 'Actual ILSP =       ', ILSP
          IERR   = 1
          RETURN
        ENDIF
C       .
C       . Now solve for the relative location under this realization
C       .
        CALL GRDRLS( IERR, NRELVM, NRELVL, DTIME1, DTIME2,
     1               DOMPC2, DSXVLS, DSYVLS, DCVARR, DRESVL,
     2               DISTMX, DELXKM, DELYKM, DRESNM,
     3               DDVEC, DWORK1, DAMAT, DTEMP1, DOLDRV )
        IF ( IERR.NE.0 ) THEN
          WRITE (6,*) 'Subroutine MCGRRL: Error'
          WRITE (6,*) 'GRDRLS returned IERR = ', IERR
          IERR   = 1
          RETURN
        ENDIF
C       .
C       . We should have a valid outcome of DELXKM, DELYKM, DRESNM
C       .
        DELXUT( IREAL ) = DELXKM
        DELYUT( IREAL ) = DELYKM
        DRESUT( IREAL ) = DRESNM
C       .
C       . Now we put the actual times used in the logs
C       .
        ILSP   = 0
        DO ISP = 1, NSP
          NMEAS  = NTMEAS( ISP )
          IF ( NMEAS.GT.0 ) THEN
            ILSP   = ILSP + 1
            DT1USD( ISP, IREAL ) = DTIME1( ILSP )
            DT2USD( ISP, IREAL ) = DTIME2( ILSP )
            DCCUSD( ISP, IREAL ) = DCVARR( ILSP )
            DRESTM( ISP, IREAL ) = DRESVL( ILSP )
          ELSE
            DT1USD( ISP, IREAL ) = 0.0d0
            DT2USD( ISP, IREAL ) = 0.0d0
            DCCUSD( ISP, IREAL ) = 0.0d0
            DRESTM( ISP, IREAL ) = 0.0d0
          ENDIF
        ENDDO
C       .
      ENDDO
C     . end of ireal = 1, nreal
      RETURN
      END
C
