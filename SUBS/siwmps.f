C
C SIWMPS
C Single Iteration Weighted Model Parameter Solve
C
C Steve Gibbons 2023-12-07
C
C We have ND measurements
C We have NM model parameters
C
C Input
C
C DDVEC is the vector of observations (ND)
C DAMAT is the matrix of relations (ND,NM)
C DWEIG is the vector of weights (ND)
C   such that a high weight is a good measurement.
C
C Output
C
C DMVEC is the vector of model parameters (NM)
C DRESV is the vector of residuals to the observations (ND)
C
C LDA is the leading dimension of DAMAT
C NMMAX is the leading dimension of the DWORK2 array
C LWORK is the leading dimension of the DWORK1 array
C LWOPT is returned with the optimal value of LWORK
C
C
      SUBROUTINE SIWMPS( IERR, LDA, ND, NM, NMMAX,
     1                   DDVEC, DAMAT, DWEIG, DMVEC, DRESV,
     2                   LWORK, LWOPT, DWORK1, DWORK2, DWORK3,
     3                   IWORK1 )
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
      INTEGER      INFO
      INTEGER      IM
      INTEGER      JM
      INTEGER      ID
      REAL*8       DTEMP
C
      CHARACTER*1  TRANS
      REAL*8       ALPHA
      REAL*8       BETA
      INTEGER      INCX
      INTEGER      INCY
C
      IERR  = 0
      IF ( LDA.LT.ND .OR. NM.LT.1 .OR. ND.LT.1 ) THEN
        WRITE (6,*) 'Subroutine SIWMPS.'
        WRITE (6,*) 'NM       = ', NM
        WRITE (6,*) 'ND       = ', ND
        WRITE (6,*) 'LDA      = ', LDA
        IERR  = 9
        RETURN
      ENDIF
C
C (W is an ND*ND matrix with the weights along the diagonal)
C Want to form A^T W D in the DWORK3 vector
C
      DO IM = 1, NM
        DWORK3( IM ) = 0.0d0
        DO ID = 1, ND
          DWORK3( IM ) = DWORK3( IM ) + 
     1      DAMAT( ID, IM )*DWEIG( ID )*DDVEC( ID )
        ENDDO
      ENDDO
C
C Now fill DWORK2 with A^{T} W A
C
      DO JM = 1, NM
        DO IM = 1, NM
          DTEMP  = 0.0d0
          DO ID = 1, ND
            DTEMP  = DTEMP + 
     1                DAMAT( ID, IM )*DWEIG( ID )*DAMAT( ID, JM )
          ENDDO
          DWORK2( IM, JM ) = DTEMP
        ENDDO
      ENDDO
C
C Invert DWORK2 using DGETRF and DGETRI
C
      CALL DGETRF( NM, NM, DWORK2, NMMAX, IWORK1, INFO )
      IF ( INFO.NE.0 ) THEN
        WRITE (6,*) 'Subroutine SIWMPS'
        WRITE (6,*) 'DGETRF returned INFO = ', INFO
        IERR  = 9
        RETURN
      ENDIF      
      CALL DGETRI( NM, DWORK2, NMMAX, IWORK1, DWORK1, LWORK, INFO )
      LWOPT  = INT( DWORK1(1)+0.0001d0 )
      IF ( INFO.NE.0 ) THEN
        WRITE (6,*) 'Subroutine SIWMPS'
        WRITE (6,*) 'DGETRI returned INFO = ', INFO
        IERR  = 9
        RETURN
      ENDIF
C
C Our vector of model parameters (DMVEC, dim NM)
C is then just DWORK2( NM, NM ) * DWORK3
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      TRANS  = 'N'
      INCX   = 1
      INCY   = 1
      CALL DGEMV( TRANS, NM, NM, ALPHA, DWORK2, NMMAX,
     1            DWORK3, INCX, BETA, DMVEC, INCY )
C
C Finally, we want to calculate our residual vector DRESV (dim ND)
C To get OBS - PRED we need +ddvec - damat*dmvec
C
      DO ID = 1, ND
        DRESV( ID ) = DDVEC( ID )
      ENDDO
      ALPHA  = -1.0d0
      BETA   =  1.0d0
      TRANS  = 'N'
      INCX   = 1
      INCY   = 1
      CALL DGEMV( TRANS, ND, NM, ALPHA, DAMAT, LDA, DMVEC, INCX,
     1            BETA, DRESV, INCY )
C
      RETURN
      END
C
