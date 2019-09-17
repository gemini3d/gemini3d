# Requires that one first build the gemini model so that we have libraries to link to
FC=gfortran -O3
FL=mpifort
INCDIR=../../build/
LINKDIR1=../../build/numerical/diffusion
LINKDIR2=../../build/numerical/
LINKDIR3=/usr/lib/


#This is a bit awkward because I need to link to all this other crap that tehcnically isn't necessary (mpi,grid,etc.) just to compile.  Possibly we need to separate basic numerical routines from specific calls that involve using the GEMINI grid...
all:
	$(FC) -c test_diffusion1D.f90 -I$(INCDIR)
#	$(FL) -o test_diffusion1D test_diffusion1D.o -L$(LINKDIR1) -ldiffusion -L$(LINKDIR2) -lgrid -lmpimod -lconst -L$(LINKDIR3) -llapack
	$(FL) -o test_diffusion1D test_diffusion1D.o -L$(LINKDIR1) -ldiffusion -lPDEparabolic -L$(LINKDIR3) -llapack
