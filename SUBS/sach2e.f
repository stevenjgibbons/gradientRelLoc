C
C Steve Gibbons
C NGI
C
C SAC human time 2 epoch time
C
C ISACTM( 1 ) = year
C ISACTM( 2 ) = Julian day
C ISACTM( 3 ) = hour
C ISACTM( 4 ) = minute
C ISACTM( 5 ) = second
C ISACTM( 6 ) = millisecond
C
      SUBROUTINE SACH2E( ISACTM, SACEPT )
      IMPLICIT NONE
C
      INTEGER ISACTM( 6 )
      REAL*8  SACEPT
C
      INTEGER NDAYS
      INTEGER IY
      INTEGER LEAPYR
C
      NDAYS  = 0
C
      IF ( ISACTM(1).LT.1970 ) THEN
        DO IY = ISACTM(1), 1969
          NDAYS   = NDAYS - 365
          CALL ISLEAP( IY, LEAPYR )
          NDAYS   = NDAYS - LEAPYR
        ENDDO
      ENDIF
C
      IF ( ISACTM(1).GT.1970 ) THEN
        DO IY = 1970, ISACTM(1) - 1
          NDAYS   = NDAYS + 365
          CALL ISLEAP( IY, LEAPYR )
          NDAYS   = NDAYS + LEAPYR
        ENDDO
      ENDIF
C
      NDAYS   = NDAYS + ISACTM( 2 ) - 1
      SACEPT  = 86400.0d0*DBLE( NDAYS )
      SACEPT  = SACEPT + 3600.0d0*DBLE( ISACTM(3) )
      SACEPT  = SACEPT +   60.0d0*DBLE( ISACTM(4) )
      SACEPT  = SACEPT +          DBLE( ISACTM(5) )
      SACEPT  = SACEPT +  0.001d0*DBLE( ISACTM(6) )
C
      RETURN
      END
C
