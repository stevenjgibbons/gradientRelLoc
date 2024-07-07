#
ALLSOURCECODE=  \
	csargm.f \
	dmedvl.f \
	grdrls.f \
	irwmps.f \
	isleap.f \
	md2doy.f \
	mdpvun.f \
	sach2e.f \
	siwmps.f
#
VALUABLE= \
        $(ALLSOURCECODE)
#
SOURCES= \
        $(ALLSOURCECODE)
#
OPTIM=	  -O3
LIBDIR=  ..
SUBSLIB= $(LIBDIR)/subslib.a
FORTRAN=  gfortran
.PRECIOUS: 	$(SUBSLIB)
#
$(SUBSLIB):	$(SOURCES)
	$(FORTRAN) $(OPTIM) $(STDLIBS) $? -c
	ar rv $(SUBSLIB) *.o
	ranlib $(SUBSLIB)
	chmod 644 $(SUBSLIB)
	\rm *.o
#
backup:
	cp -ip $(VALUABLE) ./BACKUP ; \
	cd ./BACKUP ; \
	\rm -f *.gz ; \
	gzip $(VALUABLE) ; \
	cd ../
#
