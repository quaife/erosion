FORTRAN = gfortran-mp-5
FFLAGS = -c -O3

OBJS=\
stokesSolver.o\
gmres_routines.o\
dcfft.o

libstokes: $(OBJS)
	$(FORTRAN) -shared -fPIC $(OBJS) -o libstokes.so

%.o:  %.f
	$(FORTRAN) $(FFLAGS) -c -o $@ $<	

clean: 
	rm *.o
	rm *.so
