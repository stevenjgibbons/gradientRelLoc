C
C Steve Gibbons
C NGI 2023-09-18
C
C      isleap.f
C
C Inputs a year IYEAR and outputs LEAPYR = 0 if not leap year
C and LEAPYR = 1 if it is.
C
      SUBROUTINE ISLEAP( IYEAR, LEAPYR )
      IMPLICIT NONE
C
      INTEGER IYEAR
      INTEGER LEAPYR
C
      INTEGER ITERM1
      INTEGER ITERM2
      INTEGER ITERM3
C
      LEAPYR = 0
      ITERM1 = IYEAR - (IYEAR/4)*4
      ITERM2 = IYEAR - (IYEAR/100)*100
      ITERM3 = IYEAR - (IYEAR/400)*400
      IF ( ITERM1.EQ.0 .AND. ITERM2.NE.0 .OR. ITERM3.EQ.0 ) LEAPYR = 1
C
      RETURN
      END
C
