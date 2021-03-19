FC=gfortran -O3
FL=gfortran

all:
	$(FC) -c newton.f90
	$(FC) -c dipole.f90
	$(FC) -c dipole_testdriver.f90
	$(FL) -o dipole_testdriver newton.o dipole.o dipole_testdriver.o

clean:
	rm *.mod
	rm *.smod
	rm *.o
	rm dipole_testdriver
