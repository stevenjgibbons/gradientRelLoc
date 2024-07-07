C
C StG dmedvl 
C Written by ChatGPT.
C Edited and checked by me.
C
C IERR   - error flag - zero for good output
C N      - number of points in array
C DARR   - real*8 array (N)
C DMEDN  - median value output
C
      SUBROUTINE DMEDVL( IERR, N, DARR, DMEDN )
      IMPLICIT NONE
C
      INTEGER   IERR
      INTEGER   N
      REAL*8    DARR( N )
      REAL*8    DMEDN
C
      INTEGER   I
      INTEGER   J
      INTEGER   MID
      REAL*8    DTMP
C
      IERR   = 0
      IF ( N.LT.1 ) THEN
        IERR   = 1
        RETURN
      ENDIF
C
C Quick escapes for 1 and 2
C
      IF ( N.EQ.1 ) THEN
        DMEDN = DARR( 1 )
        RETURN
      ENDIF
C
      IF ( N.EQ.2 ) THEN
        DMEDN = 0.5d0*( DARR( 1 ) + DARR( 2 ) )
        RETURN
      ENDIF
C
C Sort the array using bubble sort
      DO I = 1, N-1
        DO J = 1, N-i
          IF ( DARR(J) > DARR(J+1) ) THEN
            DTMP      = DARR(J)
            DARR(J)   = DARR(J+1)
            DARR(J+1) = DTMP
          ENDIF
        ENDDO
      ENDDO
C
      MID = N / 2
      IF ( MOD(N, 2) == 1) THEN
        DMEDN = DARR(MID + 1)
      ELSE
        DMEDN = 0.5d0*( DARR(MID) + DARR(MID+1) )
      ENDIF
C
      RETURN
      END
C
