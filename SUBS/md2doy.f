C                                                                 
C Steven J Gibbons/NGI
C Inspired by SAC routines
C                                                                 
C Takes 3 integer input parameters,
C IYEAR (year),
C IMON (1=Jan, 2=Feb, 3=Mar etc.)
C and day of month, IDOM.                     
C Returns the day-of-year (IJUL)
C and error flag IERR.            
C
      SUBROUTINE MD2DOY( IYEAR, IMON, IDOM, IJUL, IERR )
      IMPLICIT NONE
C
C Parameters 
      INTEGER IYEAR, IMON, IDOM, IJUL, IERR
C
C Working variables 
      INTEGER IM, NUMDIM( 12 )
      LOGICAL OLEAP
C
      IERR        = 0
C
      NUMDIM(  1 ) = 31
      NUMDIM(  2 ) = 28
      NUMDIM(  3 ) = 31
      NUMDIM(  4 ) = 30
      NUMDIM(  5 ) = 31
      NUMDIM(  6 ) = 30
      NUMDIM(  7 ) = 31
      NUMDIM(  8 ) = 31
      NUMDIM(  9 ) = 30
      NUMDIM( 10 ) = 31
      NUMDIM( 11 ) = 30
      NUMDIM( 12 ) = 31
C
      OLEAP = .FALSE.
C
      IF (     (IYEAR - (IYEAR/4)*4).EQ.0    .AND.
     1         (IYEAR - (IYEAR/100)*100).NE.0      ) OLEAP = .TRUE.
C
      IF (     (IYEAR - (IYEAR/400)*400).EQ.0      ) OLEAP = .TRUE.
C
      IF ( OLEAP ) NUMDIM( 2 ) = 29
C
      IF ( IMON.LT.1 .OR. IMON.GT.12 ) THEN
        PRINT *,' Subroutine MD2DOY'
        PRINT *,' IMON = ', IMON
        IERR   = 1
        RETURN
      ENDIF
C
      IF ( IDOM.LT.1 .OR. IDOM.GT.NUMDIM( IMON ) ) THEN
        PRINT *,' Subroutine MD2DOY'
        PRINT *,' IDOM = ', IDOM
        PRINT *,' Month ', IMON,' has only ',NUMDIM( IMON ),' days.'
        IERR   = 2
        RETURN
      ENDIF
C
      IJUL   = IDOM
      IF ( IMON.EQ.1 ) RETURN
      DO IM = 1, IMON - 1
        IJUL   = IJUL   + NUMDIM( IM )
      ENDDO
C
      RETURN
      END
C
