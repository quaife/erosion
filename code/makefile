FORTRAN = gfortran
FFLAGS = -c -O3

OBJS=\
stokesDriver.o\
stokesSolver.o\
gmres_routines.o\
dcfft.o

stokesDriver: $(OBJS)
	$(FORTRAN) $(OBJS) -o stokesDriver

%.o:  %.f 
	$(FORTRAN) $(FFLAGS) -c -o $@ $<	

clean: 
	rm stokesDriver
	rm *.o
