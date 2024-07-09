C
C Steve Gibbons - Numerical Recipes
C University of Leeds.
C function RAN1 ... slightly modified with IMPLICIT NONE,
C DOUBLE PRECISION
C
      DOUBLE PRECISION FUNCTION DRAN1N( IDUM )
      IMPLICIT NONE
C
      INTEGER IDUM
C
      INTEGER          IA, IM, IQ, IR, NTAB, NDIV, IDUM2
      DOUBLE PRECISION AM, EPS, RNMX
      PARAMETER      ( IA = 16807, IM = 2147483647,
     1                 AM = 1.0d0/IM, IQ = 127773,
     2                 IR = 2836, NTAB = 32,
     3                 NDIV = 1+(IM-1)/NTAB,
     4                 EPS = 1.2d-7,
     5                 RNMX = 1.0d0-EPS )
      INTEGER          j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
C
C  Copy IDUM to IDUM2 so that we can control the value of IDUM
C  on subsequent calls.
C
      IDUM2  = IDUM
C
      if (IDUM2.le.0.or.iy.eq.0) then
        IDUM2=max(-IDUM2,1)
        do j=NTAB+8,1,-1
          k=IDUM2/IQ
          IDUM2=IA*(IDUM2-k*IQ)-IR*k
          if (IDUM2.lt.0) IDUM2=IDUM2+IM
          if (j.le.NTAB) iv(j)=IDUM2
        enddo
        iy=iv(1)
      endif
      k=IDUM2/IQ
      IDUM2=IA*(IDUM2-k*IQ)-IR*k
      if (IDUM2.lt.0) IDUM2=IDUM2+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=IDUM2
      DRAN1N=min(AM*iy,RNMX)
C
C Set input argument IDUM back to the value IDUM2.
C We want to return the modified value.
C
      IDUM  = IDUM2
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.
C
