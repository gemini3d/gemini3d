FC=mpifort
FL=mpifort
INCDIR=/opt/local/include/
INCDIR2=../../build/
MUMPSDIR=/opt/local/lib/
SCALDIR=/opt/local/lib/
BLASDIR=/usr/lib/
PDEDIR=../../build/numerical/potential/
CONSTDIR=../../build/numerical/

all:
#	$(FC) -I$(INCDIR) -c test_potential2D.f90 -o test_potential2D.o
#	$(FL) test_potential2D.o -o test_potential2D -L$(MUMPSDIR) -ldmumps -lmumps_common -L$(SCALDIR) -lscalapack -L$(BLASDIR) -lblas
#	$(FC) -I$(INCDIR) -c test_potential3D.f90 -o test_potential3D.o
#	$(FL) test_potential3D.o -o test_potential3D -L$(MUMPSDIR) -ldmumps -lmumps_common -L$(SCALDIR) -lscalapack -L$(BLASDIR) -lblas
	$(FC) -I$(INCDIR) -I$(INCDIR2) -c test_potential2D.f90 -o test_potential2D.o
	$(FL) test_potential2D.o -o test_potential2D -L$(PDEDIR) -lPDEelliptic -L$(MUMPSDIR) -ldmumps -lmumps_common -L$(SCALDIR) -lscalapack -L$(BLASDIR) -lblas -L$(CONSTDIR) -lconst
