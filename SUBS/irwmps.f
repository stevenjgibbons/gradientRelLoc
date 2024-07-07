C
C Steve Gibbons
C 2023-12-07
C
C Iteratively Reweighted Model Parameter Solve
C
C The routine SIWMPS solves for NM model parameters DMVEC
C where the observation vector DDVEC has dimension ND
C and they are linked by DDVEC = DAMAT * DMVEC
C with DAMAT being an ND * NM matrix.
C 
C The vector DWEIGH (ND) contains the weights to the observations
C (high weight indicating good quality)
C
C Here we call SIWMPS repeatedly with a weight
C
C DWORK1( I ) = 1.0d0/( DEPS + DSQRT( DABS( DRESID(I) ) ) )
C 
C where DRESID(I) is the residual from the previous iteration.
C
C We stop when ITER exceeds MXITER or when the
C residual from the previous iteration is identical to the residual
C from the current iteration.
C
C See SIWMPS for details of the inputs.
C
C We have the additional integer variable MXITER
C the additional real*8 variables DTOL
C and 2 work arrays of length ND:
C  DTEMP1, DOLDRV
C
      SUBROUTINE IRWMPS( IERR, LDA, ND, NM, NMMAX,
     1                   DDVEC, DAMAT, DWEIG, DMVEC, DRESV,
     2                   LWORK, LWOPT, DWORK1, DWORK2, DWORK3,
     3                   IWORK1,
     4                   MXITER, DTOL, DTEMP1, DOLDRV )
      IMPLICIT NONE
C
      INTEGER      IERR
      INTEGER      LDA
      INTEGER      ND
      INTEGER      NM
      INTEGER      NMMAX
      REAL*8       DDVEC( ND )
      REAL*8       DAMAT( LDA, NM )
      REAL*8       DWEIG( ND )
      REAL*8       DMVEC( NM )
      REAL*8       DRESV( ND )
C
      INTEGER      LWORK
      INTEGER      LWOPT
      REAL*8       DWORK1( LWORK )
      REAL*8       DWORK2( NMMAX, NM )
      REAL*8       DWORK3( NM )
      INTEGER      IWORK1( NM )
C
C ---------
      INTEGER     MXITER
      REAL*8      DTOL
      REAL*8      DTEMP1( ND )
      REAL*8      DOLDRV( ND )
C ---------
      INTEGER     I
      INTEGER     NITER
C
      INTEGER     INCX
      PARAMETER ( INCX = 1      )
C
      REAL*8      DEPS
      PARAMETER ( DEPS = 1.0d-4 )
C
      REAL*8      DDIFF
C
      REAL*8      DNRM2
C
      IERR   = 0
C
C First simply copy the original weights into DTEMP1
C and then normalize and then solve for the initial set of
C residuals in DRESV
C
      DO I = 1, ND
        DTEMP1( I ) = DWEIG( I )
      ENDDO
      CALL MDPVUN( ND, DTEMP1 )
C
      CALL SIWMPS( IERR, LDA, ND, NM, NMMAX,
     1             DDVEC, DAMAT, DTEMP1, DMVEC, DRESV,
     2             LWORK, LWOPT, DWORK1, DWORK2, DWORK3,
     3             IWORK1 )
C
      IF ( IERR.NE.0 ) THEN
        WRITE (6,*) 'Subroutine IRWMPS'
        WRITE (6,*) 'SIWMPS returned IERR = ', IERR
        RETURN
      ENDIF
C
C Now copy the vector of residuals into DOLDRV
C
      DO I = 1, ND
        DOLDRV( I ) = DRESV( I )
      ENDDO
C
C Now start our loop around the iterations
C
      DO NITER = 1, MXITER
C       .
C       . First construct the vector of weights
C       . from the previous set of residuals
C       .
        DO I = 1, ND
          DTEMP1( I ) = 1.0d0/( DEPS + DSQRT( DABS( DOLDRV(I) ) ) )
        ENDDO
        CALL MDPVUN( ND, DWORK1 )
C       .
        CALL SIWMPS( IERR, LDA, ND, NM, NMMAX,
     1               DDVEC, DAMAT, DTEMP1, DMVEC, DRESV,
     2               LWORK, LWOPT, DWORK1, DWORK2, DWORK3,
     3               IWORK1 )
C       .
        IF ( IERR.NE.0 ) THEN
          WRITE (6,*) 'Subroutine IRWMPS'
          WRITE (6,*) 'SIWMPS returned IERR = ', IERR
          RETURN
        ENDIF
C       .
C       . Put the difference between the old residuals and the
C       . new into DWORK1 as we do not need this array right now
C       .
        DO I = 1, ND
          DTEMP1( I ) = DRESV( I ) - DOLDRV( I )
        ENDDO
        DDIFF = DNRM2( ND, DTEMP1, INCX )

        WRITE (6,81) NITER, DDIFF
 81     FORMAT('Iteration ',I8,' : resid_diff = ', f20.6)
        IF ( DDIFF.LT.DTOL ) RETURN
C       .
C       . Now need to copy the new residuals into the old
C       . for the next iteration.
C       .
        DO I = 1, ND
          DOLDRV( I ) = DRESV( I )
        ENDDO
C       .
      ENDDO
C
      WRITE (6,*) 'IRWMPS: MXITER exceeded without convergence'
C
      RETURN
      END
C
