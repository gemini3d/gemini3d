FC=gfortran -O3
FL=gfortran

all:
	$(FC) -c newton.f90
	$(FC) -c dipole.f90
	$(FC) -c newton_testdriver.f90
	$(FC) -c dipole_testdriver.f90
	$(FC) -c fullgrid_testdriver.f90
	$(FL) -o newton_testdriver newton.o dipole.o newton_testdriver.o
	$(FL) -o dipole_testdriver newton.o dipole.o dipole_testdriver.o
	$(FL) -o fullgrid_testdriver newton.o dipole.o fullgrid_testdriver.o

clean:
	rm *.mod
	rm *.smod
	rm *.o
	rm newton_testdriver
	rm dipole_testdriver
	rm fullgrid_testdriver
