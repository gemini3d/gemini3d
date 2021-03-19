FC=gfortran -O3

all:
	$(FC) -c newton.f90
	$(FC) -c dipole.f90

clean:
	rm *.mod
	rm *.smod
	rm *.o
