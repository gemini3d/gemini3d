FC=gfortran -O3
FL=gfortran

all: 
	$(FC) -c interpolation.f90 testinterp1.f90 testinterp2.f90	
	$(FL) -o testinterp1 interpolation.o testinterp1.o
	$(FL) -o testinterp2 interpolation.o testinterp2.o
