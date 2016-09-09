FORTRAN = gfortran
FFLAGS = -c -O3

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
	$(FORTRAN) -shared -fPIC $(OBJS) -o libfmm.so

%.o:  %.f
	$(FORTRAN) $(FFLAGS) -c -o $@ $<


clean:
	rm *.o
	rm *.so
