# Requires that one first build the gemini model so that we have libraries to link to
FC=gfortran -O3
FL=gfortran
INCDIR=../../objects/
LINKDIR=../../objects/numerical/

all:
	$(FC) -c testinterp1.f90 -I$(INCDIR)
	$(FL) -o testinterp1 testinterp1.o -L$(LINKDIR) -linterp
	$(FC) -c testinterp2.f90 -I$(INCDIR)
	$(FL) -o testinterp2 testinterp2.o -L$(LINKDIR) -linterp
	$(FC) -c testinterp3.f90 -I$(INCDIR)
	$(FL) -o testinterp3 testinterp3.o -L$(LINKDIR) -linterp
