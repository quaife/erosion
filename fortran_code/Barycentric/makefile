FORTRAN = gfortran
FFLAGS = -c -O3 -fPIC

OBJS=\
stokesSolver_Barycentric.o\
gmres_routines.o\
dcfft.o\
subroutine_Barycentric.o\
libfmm.so


libstokes: $(OBJS)
	$(FORTRAN) -shared -fPIC $(OBJS) -o libstokes.so

%.o:  %.f
	$(FORTRAN) $(FFLAGS) -c -o $@ $<	

clean: 
	rm *.o
	rm *.so
