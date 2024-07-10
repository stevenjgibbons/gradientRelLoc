C
C StG test the Numerical Recipes dran1n function
C
      PROGRAM testdran1n
      IMPLICIT NONE
C
      INTEGER I
      REAL*8  DRAN1N
      REAL*8  DVAL
      REAL*4  RAND
      REAL*8  DRAND
C
      CHARACTER*(8)  CHDATE
      CHARACTER*(10) CHTIME
      CHARACTER*(5)  CHZONE
      INTEGER        IVALUE(8)
C
c     CALL DATE_AND_TIME( CHDATE, CHTIME, CHZONE, IVALUE )
c     I      =      1000*IVALUE(7) + IVALUE(8)
c     I      =                       IVALUE(8)
c     I      = INT( 1000.0d0*DCOS( DBLE(I) ) )
      READ (5,*,ERR=60,END=60) I
 50   CONTINUE
c     I      = mod ( I , 257 )
c     I      = -I + 13
c     I      = INT( 1000.0d0*DCOS( DBLE(I) ) )
      I      = I + 31
      DVAL   = DRAN1N( I )
      CALL SRAND( I )
      DRAND  = DBLE( RAND( I ) )
      WRITE (6,89) I, DVAL, DRAND
 89   FORMAT(I10,1X,f20.7,1X,f20.7)
      GOTO 50
 60   CONTINUE
      END
 
