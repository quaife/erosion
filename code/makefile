FORTRAN = gfortran-mp-5
# gfortran for Bryan gfortran-mp-5 for Nick
FFLAGS = -c -O3

OBJS=\
stokesDriver.o\
stokesSolver.o\
gmres_routines.o\
dcfft.o

libstokes: $(OBJS)
	$(FORTRAN) -shared -fPIC $(OBJS) -o libstokes.so

%.o:  %.f
	$(FORTRAN) $(FFLAGS) -c -o $@ $<	

clean: 
	rm stokesDriver
	rm *.o
