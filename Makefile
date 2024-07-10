#
# Makefile for gradientRelLoc
# Steven J Gibbons 2024/07/07 (Oslo)
# Before trying to compile, please ensure that the following
# lines have the correct locations of LAPACK, and BLAS
# BINDIR must be set to the directory in which you want the executable to reside
#
LAPACK= /usr/lib/x86_64-linux-gnu/liblapack.a
BLAS= /usr/lib/x86_64-linux-gnu/libblas.a
TOPDIR=   .
BINDIR=  $(TOPDIR)/bin
#
# PLEASE CHECK ALL THE ABOVE LINES FOR YOUR SYSTEM ----
#
PROGNAME1= gradientRelLoc
PROGNAME2= testdran1n
PROGNAME3= mc_gradientRelLoc
PROGNAME4= fixedSlovecsEventSolve
PROGNAME5= fixedEventsSlovecsSolve
PROGNAME6= checkcircvc
#
ALLSOURCECODE=  \
   $(PROGNAME1).f   \
   $(PROGNAME2).f   \
   $(PROGNAME3).f   \
   $(PROGNAME4).f   \
   $(PROGNAME5).f   \
   $(PROGNAME6).f  
#
SOURCES= \
        $(ALLSOURCECODE)
#
OPTIM=	  -O3
EXEFILE1= $(BINDIR)/$(PROGNAME1)
EXEFILE2= $(BINDIR)/$(PROGNAME2)
EXEFILE3= $(BINDIR)/$(PROGNAME3)
EXEFILE4= $(BINDIR)/$(PROGNAME4)
EXEFILE5= $(BINDIR)/$(PROGNAME5)
EXEFILE6= $(BINDIR)/$(PROGNAME6)
FORTRAN= gfortran
#
LIBDIR=  .
SUBSLIB= $(LIBDIR)/subslib.a
LIBS=    $(SUBSLIB) $(LAPACK) $(BLAS) $(SACLIB) 
#
backup:
	cp -ip $(ALLSOURCECODE) ./BACKUP ; \
	cd ./BACKUP ; \
	\rm -f *.gz ; \
	gzip $(ALLSOURCECODE) ; \
	cd ../
#
$(PROGNAME1):	$(PROGNAME1).f $(LIBS)
	$(FORTRAN) -o $(EXEFILE1) $(PROGNAME1).f $(LIBS) $(OPTIM)
$(PROGNAME2):	$(PROGNAME2).f $(LIBS)
	$(FORTRAN) -o $(EXEFILE2) $(PROGNAME2).f $(LIBS) $(OPTIM)
$(PROGNAME3):	$(PROGNAME3).f $(LIBS)
	$(FORTRAN) -o $(EXEFILE3) $(PROGNAME3).f $(LIBS) $(OPTIM)
$(PROGNAME4):	$(PROGNAME4).f $(LIBS)
	$(FORTRAN) -o $(EXEFILE4) $(PROGNAME4).f $(LIBS) $(OPTIM)
$(PROGNAME5):	$(PROGNAME5).f $(LIBS)
	$(FORTRAN) -o $(EXEFILE5) $(PROGNAME5).f $(LIBS) $(OPTIM)
$(PROGNAME6):	$(PROGNAME6).f $(LIBS)
	$(FORTRAN) -o $(EXEFILE6) $(PROGNAME6).f $(LIBS) $(OPTIM)
