FORTRAN = gfortran
#FFLAGS = -c -O3 -fPIC -fopenmp -fallow-argument-mismatch # for old gcc/gfortran, replace -fallow-argument-mismatch by -Wno-argument-mismatch
FFLAGS = -c -O3 -fPIC -fopenmp -Wno-argument-mismatch # for old gcc/gfortran, replace -fallow-argument-mismatch by -Wno-argument-mismatch


OBJS=\
stokesSolver.o\
gmres_routines.o\
dcfft.o

libstokes: $(OBJS) fmm/libfmm.so
ifeq ($(shell uname),Darwin)
	$(FORTRAN) -shared -fopenmp -fPIC -install_name '$$ORIGIN/libstokes.so.dylib' -Lfmm -lfmm $(OBJS) -o libstokes.so.dylib
else
	$(FORTRAN) -shared -fopenmp -fPIC -Lfmm -Wl,-rpath,'$$ORIGIN'/fmm -lfmm $(OBJS) -o libstokes.so
endif

fmm/libfmm.so:
	make -C fmm

%.o:  %.f
	$(FORTRAN) $(FFLAGS) -c -o $@ $<	

clean:
	make -C fmm clean
	-rm *.o
	-rm *.so
	-rm *.dylib
