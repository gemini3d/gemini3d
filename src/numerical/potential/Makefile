FC=mpifort
FL=mpifort
INCDIR=/Users/zettergm/Projects/gemini_libs_gcc_clang/include
INCDIR2=../../build/
INCDIR3=../../../build/include/
MUMPSDIR=/Users/zettergm/Projects/gemini_libs_gcc_clang/lib
SCALDIR=/Users/zettergm/Projects/gemini_libs_gcc_clang/lib
BLASDIR=/Users/zettergm/Projects/gemini_libs_gcc_clang/lib
BLACSDIR=/Users/zettergm/Projects/gemini_libs_gcc_clang/lib
PDEDIR=../../build/numerical/potential/
CONSTDIR=../../build/numerical/

all:
#	$(FC) -I$(INCDIR) -c test_potential2D.f90 -o test_potential2D.o
#	$(FL) test_potential2D.o -o test_potential2D -L$(MUMPSDIR) -ldmumps -lmumps_common -L$(SCALDIR) -lscalapack -L$(BLASDIR) -lblas
	$(FC) -I$(INCDIR) -I$(INCDIR3) -c test_potential3D.f90 -o test_potential3D.o
	$(FL) test_potential3D.o -o test_potential3D -L$(MUMPSDIR) -ldmumps -lmumps_common -L$(SCALDIR) -lscalapack -L$(BLASDIR) -lblas -L$(BLACSDIR) -lblacs -lpord -llapack
#	$(FC) -I$(INCDIR) -I$(INCDIR2) -c test_potential2D.f90 -o test_potential2D.o
#	$(FL) test_potential2D.o -o test_potential2D -L$(PDEDIR) -lPDEelliptic -L$(MUMPSDIR) -ldmumps -lmumps_common -L$(SCALDIR) -lscalapack -L$(BLASDIR) -lblas -L$(CONSTDIR) -lconst

