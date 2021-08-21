FORTRAN = gfortran
#FFLAGS = -c -O3 -fPIC -fopenmp -fallow-argument-mismatch # for old gcc/gfortran, replace -fallow-argument-mismatch by -Wno-argument-mismatch
FFLAGS = -c -O3 -fPIC -fopenmp -Wno-argument-mismatch # for old gcc/gfortran, replace -fallow-argument-mismatch by -Wno-argument-mismatch


OBJS=\
bhcommon_oldtree.o\
bhfmm2dpart.o\
bhfmm2drouts.o\
d2tstrcr_omp.o\
d2mtreeplot.o\
l2dterms.o\
laprouts2d.o\
bhrouts.o\
prini.o\
hkrand.o\
dlaran.o\

libfmm: $(OBJS)
ifeq ($(shell uname),Darwin)
	$(FORTRAN) -install_name "`pwd`/libfmm.so" -fopenmp -shared -fPIC $(OBJS) -o libfmm.so
else
	$(FORTRAN) -fopenmp -shared -fPIC $(OBJS) -o libfmm.so
endif

%.o:  %.f
	$(FORTRAN) $(FFLAGS) -c -o $@ $<	

clean:
	-rm *.o
	-rm *.so
