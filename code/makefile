FORTRAN = gfortran
FFLAGS = -c -O3

OBJS=\
stokesSolver.o\
gmres_routines.o\
dcfft.o\
libfmm.so

libstokes: $(OBJS)
	$(FORTRAN) -shared -fPIC $(OBJS) -o libstokes.so

%.o:  %.f
	$(FORTRAN) $(FFLAGS) -c -o $@ $<	

clean: 
	rm *.o
	rm *.so
