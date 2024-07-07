C*********************************************************************
C Subroutine Character String ARGument Map ***************************
C            -         -      ---      -   ***************************
C Steve Gibbons August 2001 (University of Leeds)                    C
C____________________________________________________________________C
C                                                                    C
C Takes in a character string, LINE, a number of arguments, NARGS,   C
C and a maximum line length. Assuming the arguments to be separated  C
C by spaces, CSARGM will fill the integer array IARGL with the       C
C character positions of the start and end of the argument.          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input Variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     LINE      : String (*) containing input arguments.             C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NARGS     : Number of arguments to be searched for.            C
C     NARGM     : Maximum no. of arguments (i.e. dimension for IARGL C
C     CSLEN     : Length of character string.                        C
C                                                                    C
C Output Variables :-                                                C
C ================                                                   C
C                                                                    C
C     IARGL     : Array dim ( NARGM, 2 ).                            C
C                 On exit, IARGL( I, 1 ) contains the character      C
C                 position where argument I begins: IARGL( I, 2 )    C
C                 contains the character position where arg I ends.  C
C                                                                    C
C     IERR      : =0 all NARGS are found and located.                C
C                 >0 CSLEN was exceeded with only IERR arguments     C
C                 being located.                                     C
C                 IERR = -1 means that no arguments were found.      C
C                 IERR = -2 means that the last argument was not     C
C                 terminated, indicating a possible string longer    C
C                 than was anticipated.                              C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CSARGM( LINE, NARGS, NARGM, CSLEN, IARGL, IERR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
C
      CHARACTER *(*) LINE
      INTEGER        NARGS, NARGM, CSLEN, IARGL( NARGM, 2 ), IERR
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER        I, IARGF
      LOGICAL        OWORDP, OWORDC
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NARGS.LT.0 .OR. NARGS.GT.NARGM ) THEN
        PRINT *,' Subroutine CSARGM.'
        PRINT *,' NARGS = ', NARGS
        PRINT *,' NARGM = ', NARGM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IARGF = 1, NARGM
        IARGL( IARGF, 1 ) = -1
        IARGL( IARGF, 2 ) = -1
      ENDDO
C
      IARGF  = 0
      OWORDP = .FALSE.
      OWORDC = .FALSE.
      DO I = 1, CSLEN
        IF ( LINE(I:I).EQ.' ' ) THEN
          OWORDC = .FALSE.
        ELSE
          OWORDC = .TRUE.
        ENDIF
        IF ( I.EQ.1 ) OWORDP = .FALSE.
C
C Treat the case of us starting a new word
C
        IF ( .NOT. OWORDP .AND. OWORDC ) THEN
          IARGF = IARGF + 1
          IARGL( IARGF, 1 ) = I
        ENDIF
C
C Treat the case of us ending a word
C
        IF ( OWORDP .AND. .NOT. OWORDC ) IARGL( IARGF, 2 ) = I - 1
C
        IF ( IARGF.EQ.NARGS .AND. IARGL( IARGF, 2 ).GT.0 ) THEN
          IERR = 0
          RETURN
        ENDIF
C
        OWORDP = OWORDC
C
      ENDDO
C
      IF ( IARGF.EQ.0 ) THEN
        IERR = -1
        RETURN
      ENDIF
c     print *,' iargf = ', iargf
      IF ( IARGF.LT.NARGS ) IERR = IARGF
      IF ( IARGF.EQ.NARGS .AND. IARGL( IARGF, 2 ).LT.0 ) THEN
        IERR = -2
        IARGL( IARGF, 2 ) = CSLEN
      ENDIF
C
      RETURN
      END
C*********************************************************************
