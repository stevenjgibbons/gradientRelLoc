C gradientRelLoc
C 
C You specify a number of lines of a file statphase.txt containing the following items
C 
C station phase   statlat   statlon    reflat    reflon       sx            sy
C 
C YSDH  Pn        35.16340  132.85580  41.29500  129.08000    0.05634205   -0.11011804
C YUZH  Pn        39.19130  140.47100  41.29500  129.08000    0.12188289   -0.02109431
C YWTH  Pn        38.97010  140.03330  41.29500  129.08000    0.12104773   -0.02545301
C 
C say that statphase.txt has NSP lines
C 
C and then you read in all the lines of a file CC_times.txt with the format
C 
C event1   event2    time_ev1                 time_ev2                stat    phase   CC
C 
C DPRK5    DPRK6    2016-09-09T00:32:09.6900 2017-09-03T03:32:10.0850 YZWH     Pn   0.6710
C DPRK6    DPRK1    2017-09-03T03:32:10.0100 2006-10-09T01:37:35.9247 YZWH     Pn   0.3624
C DPRK6    DPRK2    2017-09-03T03:32:10.0100 2009-05-25T00:56:51.3870 YZWH     Pn   0.6056
C 
C 
C So do
C 
C event1=DPRK3
C event2=DPRK4
C NSP=`wc statphase.txt | awk '{print $1}'`
C cp statphase.txt gradientRelLoc.input
C cat CC_times.txt >> gradientRelLoc.input
C gradientRelLoc ${NSP} ${event1} ${event2} < gradientRelLoc.input
C 
C It will then read the NSP lines reading (at least)
C 
C STATION PHASE SX  SY
C 
C and then read in the remaining lines and ignoring those for which 
C event1 is not EVENT1 and event2 is not EVENT2
C 
C We read in time_ev1, time_ev2, stat, phase, CC
C 
C We identify which SX and SY our stat and phase corresponds to and abort
C with an error if we do not find a suitable value.
C 
C We will borrow a lot of code from onestatphaserel2abs.f
C 
      PROGRAM gradientRelLoc
      IMPLICIT NONE
C
      INTEGER     NIARGS
      INTEGER     IARGC
      INTEGER     IARG
C
c     INTEGER     INCX
c     REAL*8      DNRM2
      REAL*8      DRESNM
c     INTEGER     LWOPT
      INTEGER     IERR
      INTEGER     ISP
c     INTEGER     MXITER
c     PARAMETER ( MXITER = 10000 )
c     REAL*8           DTOL
c     PARAMETER      ( DTOL = 1.0d-6 )
C
c     CHARACTER*(24)  CHUTMI
c     CHARACTER*(24)  CHUTMJ
      CHARACTER*(200) CHARG
      CHARACTER*(80)  CMESS
C for the Finland dataset we need 55 events and over 30000 relative readings.
      INTEGER     NRELVL
      INTEGER     NRELVM
      PARAMETER ( NRELVM = 50000 )
      INTEGER     NSP
      INTEGER     NSPMAX
      PARAMETER ( NSPMAX = 500 )
C
      CHARACTER*(8) CSTARR( NSPMAX )
      INTEGER       LSTARR( NSPMAX )
      CHARACTER*(8) CSTAT
      CHARACTER*(8) CPHARR( NSPMAX )
      INTEGER       LPHARR( NSPMAX )
      CHARACTER*(8) CPHAS
      REAL*8        DSXARR( NSPMAX )
      REAL*8        DSX
      REAL*8        DSYARR( NSPMAX )
      REAL*8        DSY
      INTEGER       LSTAT   
      INTEGER       LPHAS   
C
      INTEGER          LWORK
      PARAMETER      ( LWORK  = 2*NRELVM )
      INTEGER       NM
      PARAMETER   ( NM = 3 )
c     INTEGER          IWORK1( NM )
      REAL*8        DWORK1( LWORK )
c     REAL*8        DWORK2( NRELVM,  NM )
c     REAL*8        DWORK3( NM )
      REAL*8        DAMAT( NRELVM, NM )
      REAL*8        DDVEC( NRELVM )
c     REAL*8        DMVEC( NM )
      REAL*8        DRESVL( NRELVM )
      REAL*8        DOLDRV( NRELVM )
      REAL*8        DTEMP1( NRELVM )
c     CHARACTER*(90) CRLINE( NRELVM )
C
c     INTEGER       ISPARR( NRELVM )
      REAL*8        DCVARR( NRELVM )
      REAL*8        DEPOT1( NRELVM )
      REAL*8        DEPOT2( NRELVM )
      REAL*8        DOMPC2( NRELVM )
c     REAL*8        DIFARR( NRELVM )
c     REAL*8        DWKARR( NRELVM )
      REAL*8        DSXVAR( NRELVM )
      REAL*8        DSYVAR( NRELVM )
      REAL*8         DTIME1
      REAL*8         DTIME2
      REAL*8         DCCVAL
c     REAL*8         DDIFF
c     REAL*8         DVAL
      REAL*8         DISTMX
      REAL*8         DELXKM
      REAL*8         DELYKM
C 
C Variables for SAC routines
C
      INTEGER     ISACTM( 6 )
      INTEGER     IMON
      INTEGER     IDOM
      INTEGER     IJUL
      INTEGER     IYYYY
      REAL*8      DSECS
C
C     Variables for the routine CSARGM
C
      INTEGER     ILEN
      INTEGER     I0
      INTEGER     I1
      INTEGER     I2
      INTEGER     NARGM
      PARAMETER ( NARGM = 8 )
      INTEGER     NARGS, CSLEN, IARGL( NARGM, 2 )
C
C CTARE1 is the name of event one
C LTARE1 is the length of the character string CTARE1
      CHARACTER*(24) CTARE1
      INTEGER        LTARE1
C CTARE2 is the name of event two
C LTARE2 is the length of the character string CTARE2
      CHARACTER*(24) CTARE2
      INTEGER        LTARE2
C
      INTEGER        I
C
      CMESS   = ' '
      NIARGS  = IARGC()
      IF ( NIARGS.NE.3 ) THEN
        WRITE (6,*) 'Usage:  NSP       EVENT1    EVENT2     '
        WRITE (6,*) '        111       DPRK3      DPRK4      '
        CALL EXIT(1)
      ENDIF
C
      CHARG  = ' '
      IARG   = 1
      CALL GETARG( IARG, CHARG )
      CMESS  = 'Error reading integer NSP'
      READ ( CHARG, *, ERR=99, END=99 ) NSP
      IF ( NSP.LT.1 .OR. NSP.GT.NSPMAX ) THEN
        WRITE (6,*) 'NSP    = ', NSP
        WRITE (6,*) 'NSPMAX = ', NSPMAX
        CMESS  = 'NSPMAX exceeded'
        GOTO 99
      ENDIF
C
      CHARG  = ' '
      IARG   = 2
      CALL GETARG( IARG, CHARG )
      LTARE1 = 24
      DO I = 2, 24
        IF ( CHARG(I:I).EQ.' ' .AND. LTARE1.EQ.24 ) LTARE1 = I-1
      ENDDO
      CTARE1           = '                        '
      CTARE1(1:LTARE1) = CHARG(1:LTARE1)
C
      CHARG  = ' '
      IARG   = 3
      CALL GETARG( IARG, CHARG )
      LTARE2 = 24
      DO I = 2, 24
        IF ( CHARG(I:I).EQ.' ' .AND. LTARE2.EQ.24 ) LTARE2 = I-1
      ENDDO
      CTARE2           = '                        '
      CTARE2(1:LTARE2) = CHARG(1:LTARE2)
C
C Successfully read in NSP, CTARE1, CTARE2
c     print *,' NSP = ', NSP, CTARE1(1:LTARE1), CTARE2(1:LTARE2)
C Now need to read in NSP lines containing STAT PHASE lat lon lat lon SX SY
C (ignore the lats and lons for the time being)
C
      ISP    = 0
 30   CONTINUE
      CHARG  = ' '
      CMESS  = 'Error reading in slowness line'
      READ ( 5, '(A)', ERR=99, END=99 ) CHARG
      IF ( CHARG(1:1).EQ.'#' ) GOTO 30
      IF ( CHARG(1:1).EQ.'*' ) GOTO 30
      CSLEN  = 200
      NARGS  = 8
      CALL CSARGM( CHARG, NARGS, NARGM, CSLEN, IARGL, IERR )
      IF ( IERR.NE.0 ) THEN
        CMESS  = 'Error from CSARGM '
        WRITE (6,*) 'CSARGM returned IERR = ', IERR
        GOTO 99
      ENDIF
C
C Need to find STAT
C
      IARG   = 1
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      LSTAT  = I2 - I1 + 1
      IF ( LSTAT.GT.8 ) THEN
        WRITE (6,*) 'Invalid length for station ', CHARG(I1:I2)
        GOTO 99
      ENDIF
      CSTAT  = ' '
      CSTAT  = CHARG(I1:I2)
C
C Need to find PHAS
C
      IARG   = 2
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      LPHAS  = I2 - I1 + 1
      IF ( LPHAS.GT.8 ) THEN
        WRITE (6,*) 'Invalid length for station ', CHARG(I1:I2)
        GOTO 99
      ENDIF
      CPHAS  = ' '
      CPHAS  = CHARG(I1:I2)
C
C Need to find SX
C
      IARG   = 7
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      CMESS  = 'Error reading SX'
      READ ( CHARG(I1:I2), *, ERR=99, END=99 ) DSX
C
C Need to find SY
C
      IARG   = 8
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      CMESS  = 'Error reading SY'
      READ ( CHARG(I1:I2), *, ERR=99, END=99 ) DSY
C
      ISP    = ISP + 1
      CSTARR( ISP ) = CSTAT
      LSTARR( ISP ) = LSTAT
      CPHARR( ISP ) = CPHAS
      LPHARR( ISP ) = LPHAS
      DSXARR( ISP ) = DSX
      DSYARR( ISP ) = DSY
      IF ( ISP.EQ.NSP ) GOTO 40
      GOTO 30
 40   CONTINUE

c     do i = 1, nsp
c       write (6,*) CSTARR( i ), CPHARR( i ),
c    1             DSXARR( i ), DSYARR( i )
c     enddo
C
      NRELVL = 0
 50   CONTINUE
      CHARG = ' '
      READ (5,'(A)',ERR=50,END=60) CHARG
      IF ( CHARG(1:1).EQ.'*' ) GOTO 50
      IF ( CHARG(1:1).EQ.'#' ) GOTO 50
C
C This line should contain 7 terms
C CEV1
C CEV2
C CHTIM1
C CHTIM2
C CSTAT
C CPHASE
C CCV
C
      NARGS  = 7
      CALL CSARGM( CHARG, NARGS, NARGM, CSLEN, IARGL, IERR )
      IF ( IERR.NE.0 ) THEN
        WRITE (6,*) 'CSARGM returned IERR = ', IERR
        GOTO 99
      ENDIF
C
C We read in all items!
C So make an early exit of this line if CEV1 is not CTARE1(1:LTARE1)
C or if CEV2 is not CTARE2(1:LTARE2)
C
C Need to find CEV1
C
      IARG   = 1
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      ILEN   = I2 - I1 + 1
      IF ( ILEN.NE.LTARE1 ) GOTO 50
      IF ( CHARG(I1:I2).NE.CTARE1(1:LTARE1)  )   GOTO 50
C
C Need to find CEV2
C
      IARG   = 2
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      ILEN   = I2 - I1 + 1
      IF ( ILEN.NE.LTARE2 ) GOTO 50
      IF ( CHARG(I1:I2).NE.CTARE2(1:LTARE2)  )   GOTO 50
C
C So we know this line corresponds to the right event pair.
C Need to find ISP - the station and phase combination the
C line belongs to ...
C
      IARG   = 5
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      LSTAT  = I2 - I1 + 1
      IF ( LSTAT.GT.8 ) THEN
        CMESS = 'Invalid length of station name'
        GOTO 99
      ENDIF
      CSTAT  = ' '
      CSTAT  = CHARG(I1:I2)
C
      IARG   = 6
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      LPHAS  = I2 - I1 + 1
      IF ( LPHAS.GT.8 ) THEN
        CMESS = 'Invalid length of phase name'
        GOTO 99
      ENDIF
      CPHAS  = ' '
      CPHAS  = CHARG(I1:I2)
C
      ISP    = 0
      DO I = 1, NSP
        IF ( LPHAS.EQ.LPHARR( I ) .AND. LSTAT.EQ.LSTARR( I ) ) THEN
          IF (    CSTARR(I)(1:LSTAT).EQ.CSTAT(1:LSTAT)  .AND.
     1            CPHARR(I)(1:LPHAS).EQ.CPHAS(1:LPHAS)  ) ISP = I
        ENDIF
      ENDDO
      IF ( ISP.EQ.0 ) THEN
        WRITE ( 6,*) 'Station = ', CSTAT
        WRITE ( 6,*) 'Phase   = ', CPHAS
        CMESS = 'Station/phase not found'
        GOTO 99
      ENDIF
C     .
C     . Now need to try to read the epoch time for event 1
C     .
      IARG   = 3
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      ILEN   = I2 - I1 + 1
      IF ( ILEN.LT.23 ) THEN
        WRITE (6,*) 'Invalid string for time ', CHARG(I1:I2)
        GOTO 99
      ENDIF
      I0     = I1 - 1
      READ ( CHARG(I0+ 1:I0+ 4), '(I4)', END=99, ERR=99 ) ISACTM( 1 )
      READ ( CHARG(I0+ 6:I0+ 7), '(I2)', END=99, ERR=99 ) IMON
      READ ( CHARG(I0+ 9:I0+10), '(I2)', END=99, ERR=99 ) IDOM
      IYYYY  = ISACTM( 1 )
      CALL MD2DOY( IYYYY, IMON, IDOM, IJUL, IERR )
      IF ( IERR.NE.0 ) THEN
        WRITE (6,*) 'Error from MD2DOY: ', IYYYY, IMON, IDOM, IJUL
        GOTO 99
      ENDIF
      ISACTM( 2 ) = IJUL
      READ ( CHARG(I0+12:I0+13), '(I2)', END=99, ERR=99 ) ISACTM( 3 )
      READ ( CHARG(I0+15:I0+16), '(I2)', END=99, ERR=99 ) ISACTM( 4 )
      READ ( CHARG(I0+18:I0+24), *,      END=99, ERR=99 ) DSECS
      ISACTM( 5 ) = INT( DSECS )
      ISACTM( 6 ) = INT( (DSECS - DBLE( INT( DSECS ) ))*1000.0d0 )
      CALL SACH2E( ISACTM, DTIME1 )
c     WRITE (6,'(f20.4)') DTIME1
C
C     . Now need to try to read the epoch time for event 2
C     .
      IARG   = 4
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      ILEN   = I2 - I1 + 1
      IF ( ILEN.LT.23 ) THEN
        WRITE (6,*) 'Invalid string for time ', CHARG(I1:I2)
        GOTO 99
      ENDIF
      I0     = I1 - 1
      READ ( CHARG(I0+ 1:I0+ 4), '(I4)', END=99, ERR=99 ) ISACTM( 1 )
      READ ( CHARG(I0+ 6:I0+ 7), '(I2)', END=99, ERR=99 ) IMON
      READ ( CHARG(I0+ 9:I0+10), '(I2)', END=99, ERR=99 ) IDOM
      IYYYY  = ISACTM( 1 )
      CALL MD2DOY( IYYYY, IMON, IDOM, IJUL, IERR )
      IF ( IERR.NE.0 ) THEN
        WRITE (6,*) 'Error from MD2DOY: ', IYYYY, IMON, IDOM, IJUL
        GOTO 99
      ENDIF
      ISACTM( 2 ) = IJUL
      READ ( CHARG(I0+12:I0+13), '(I2)', END=99, ERR=99 ) ISACTM( 3 )
      READ ( CHARG(I0+15:I0+16), '(I2)', END=99, ERR=99 ) ISACTM( 4 )
      READ ( CHARG(I0+18:I0+24), *,      END=99, ERR=99 ) DSECS
      ISACTM( 5 ) = INT( DSECS )
      ISACTM( 6 ) = INT( (DSECS - DBLE( INT( DSECS ) ))*1000.0d0 )
      CALL SACH2E( ISACTM, DTIME2 )
c     WRITE (6,'(f20.4)') DTIME2
C
C Finally, read the weight in the 7th column
C
      IARG   = 7
      I1     = IARGL( IARG, 1 )
      I2     = IARGL( IARG, 2 )
      READ ( CHARG(I1:I2), *, ERR=99, END=99 ) DCCVAL
C
C OK - so we have all of the values we need!
C So we can increment NRELVL
C
      IF ( NRELVL.EQ.NRELVM ) THEN
        CMESS  = 'NRELVL about to be exceed ... '
        GOTO 99
      ENDIF
      NRELVL = NRELVL + 1
c     DDIFF  = DTIME2 - DTIME1
c     ISPARR( NRELVL ) = ISP
      DEPOT1( NRELVL ) = DTIME1
      DEPOT2( NRELVL ) = DTIME2
      DCVARR( NRELVL ) = DCCVAL
c     DIFARR( NRELVL ) = DDIFF
c     DWKARR( NRELVL ) = DDIFF
      DSXVAR( NRELVL ) = DSXARR( ISP )
      DSYVAR( NRELVL ) = DSYARR( ISP )
      DAMAT( NRELVL, 1 ) = DSXVAR( NRELVL )
      DAMAT( NRELVL, 2 ) = DSYVAR( NRELVL )
      DAMAT( NRELVL, 3 ) = 1.0d0
c     CRLINE( NRELVL ) = ' '
c     CRLINE( NRELVL ) = CHARG(1:90)
c     PRINT *, ISP, DDIFF, DCCVAL
c     print *, ISP, CSTAT, CPHAS
C
      GOTO 50
 60   CONTINUE
C
C So we have a total of NRELVL measurements.
C The first thing we want to do is to subtract the median
C difference from each element of DIFARR( )
C
      DISTMX = 0.5d0
      CALL GRDRLS( IERR, NRELVM, NRELVL, DEPOT1, DEPOT2,
     1             DOMPC2, DSXVAR, DSYVAR, DCVARR, DRESVL,
     2             DISTMX, DELXKM, DELYKM, DRESNM,
     3             DDVEC, DWORK1, DAMAT, DTEMP1, DOLDRV )
      IF ( IERR.NE.0 ) THEN
        WRITE (6,*) 'Error from GRDRLS '
        CALL EXIT(1)
      ENDIF
C
      DO I = 1, NRELVL
        WRITE (6,83) I, DSXVAR(I), DSYVAR(I), DDVEC( I ),
     1                  DDVEC( I ) - DRESVL( I ), DRESVL( I ),
     2                  CTARE1(1:LTARE1), CTARE2(1:LTARE2), 'res'
      ENDDO
 83   FORMAT(I5,1X,f10.6,1X,f10.6,1X,f10.4,1X,f10.4,1X,f6.3,
     1       1X,A,1X,A,1X,A)
C
      WRITE (6,82) CTARE1(1:LTARE1),
     1             CTARE2(1:LTARE2),
     2             DELXKM, DELYKM,
     3             DRESNM
 82   FORMAT('Result: ',A,1X,A,1X,f10.5,1X,f10.5,1X,f12.6)
c     PRINT *,'  delx, dely = ', DMVEC( 1 ), DMVEC( 2 )
C
C Now check the residuals by repeating the procedure
C with DRESVL put into DOMPC2
C
      DO I = 1, NRELVL
        DOMPC2( I ) = DRESVL( I )
      ENDDO
C
      DISTMX = 0.5d0
      CALL GRDRLS( IERR, NRELVM, NRELVL, DEPOT1, DEPOT2,
     1             DOMPC2, DSXVAR, DSYVAR, DCVARR, DRESVL,
     2             DISTMX, DELXKM, DELYKM, DRESNM,
     3             DDVEC, DWORK1, DAMAT, DTEMP1, DOLDRV )
      IF ( IERR.NE.0 ) THEN
        WRITE (6,*) 'Error from GRDRLS '
        CALL EXIT(1)
      ENDIF
C
      DO I = 1, NRELVL
        WRITE (6,83) I, DSXVAR(I), DSYVAR(I), DDVEC( I ),
     1                  DDVEC( I ) - DRESVL( I ), DRESVL( I ),
     2                  CTARE1(1:LTARE1), CTARE2(1:LTARE2), 'res'
      ENDDO
c83   FORMAT(I5,1X,f10.6,1X,f10.6,1X,f10.4,1X,f10.4,1X,f6.3,
c    1       1X,A,1X,A,1X,A)
C
      WRITE (6,82) CTARE1(1:LTARE1),
     1             CTARE2(1:LTARE2),
     2             DELXKM, DELYKM,
     3             DRESNM
c82   FORMAT('Result: ',A,1X,A,1X,f10.5,1X,f10.5,1X,f12.6)
C
      CALL EXIT(0)
 99   CONTINUE
      WRITE (6,'(A)') CMESS
      CALL EXIT(1)
      END
C
