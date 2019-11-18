#BY DEFAULT ONLY MAKE THE MAIN GEMINI PROGRAM (MOST PEOPLE WILL NOT CARE ABOUT COMPUTING MAGNETIC FIELDS
.DEFAULT_GOAL := gemini_mumps

FC=mpifort
FL=mpifort

CURRDIR = $(shell pwd)
HOMEDIR = $(HOME)
LAPACK77 = $(HOMEDIR)/lib/lapack-3.5.0/liblapack.a
BLAS     = $(HOMEDIR)/lib/BLAS/libblas.a
LAPACK95 = $(HOMEDIR)/lib/LAPACK95/lapack95.a
TMG = $(HOMEDIR)/lib/lapack-3.5.0/libtmglib.a
MUMPS = $(MUMPSDIR)/lib/libdmumps.a $(MUMPSDIR)/lib/libmumps_common.a  -L$(MUMPSDIR)/PORD/lib/ -lpord

OBJDIR = ./objects/
MUMPSDIR=$(HOMEDIR)/lib/MUMPS_4.10.0/
SCALADIR=$(HOMEDIR)/lib/scalapack-2.0.2/
BLACSDIR=$(HOMEDIR)/lib/BLACS/
BLASDIR=$(HOMEDIR)/lib/BLAS/
ILUDIR=$(HOMEDIR)/lib/ilupackV2.4_GNU64_long_MUMPS/

LIBBLACS=$(BLACSDIR)/LIB/blacs_MPI-LINUX-0.a $(BLACSDIR)/LIB/blacsF77init_MPI-LINUX-0.a

#ONLY WORKS WITHOUT OPTIMIZATION RIGHT NOW
#OPTIM=-O0 -fmax-errors=10
OPTIM=-O3 -fmax-errors=10

#-03 on everything but phys_consts.f90 appears to work
OPTIM2=-O3 -fmax-errors=10

#TO KILL OFF WARNINGS FROM MSIS CODE
OPTIONS=-w

#all: gemini_mumps gemini_ilupack
all: gemini_mumps magcalc eq

gemini_mumps: $(OBJDIR)advec_mpi.o $(OBJDIR)diffusion.o $(OBJDIR)ionization.o $(OBJDIR)phys_consts.o $(OBJDIR)potential_mumps.o \
	$(OBJDIR)calculus.o $(OBJDIR)gemini.o $(OBJDIR)msis00_gfortran.o $(OBJDIR)potentialBCs_mumps.o $(OBJDIR)precipBCs_mod.o \
	$(OBJDIR)temporal.o $(OBJDIR)collisions.o $(OBJDIR)io.o $(OBJDIR)neutral.o $(OBJDIR)potential_comm_mumps.o $(OBJDIR)sources.o  \
	$(OBJDIR)mpimod.o $(OBJDIR)multifluid.o $(OBJDIR)grid.o $(OBJDIR)interpolation.o $(OBJDIR)expanduser.o
	$(FL) -o gemini_mumps $(OBJDIR)expanduser.o $(OBJDIR)interpolation.o $(OBJDIR)grid.o $(OBJDIR)multifluid.o $(OBJDIR)mpimod.o $(OBJDIR)advec_mpi.o $(OBJDIR)diffusion.o $(OBJDIR)ionization.o $(OBJDIR)phys_consts.o $(OBJDIR)potential_mumps.o $(OBJDIR)calculus.o $(OBJDIR)gemini.o $(OBJDIR)msis00_gfortran.o $(OBJDIR)potentialBCs_mumps.o $(OBJDIR)precipBCs_mod.o $(OBJDIR)temporal.o $(OBJDIR)collisions.o $(OBJDIR)io.o $(OBJDIR)neutral.o $(OBJDIR)potential_comm_mumps.o $(OBJDIR)sources.o $(MUMPS)  $(SCALADIR)/libscalapack.a $(LIBBLACS) $(LAPACK95) $(TMG) $(LAPACK77) $(BLAS) -L/usr/local/lib/ -lmpi -lpthread -L/usr/local/lib/ -L/usr/lib/x86_64-linux-gnu/


magcalc: $(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o $(OBJDIR)grid.o $(OBJDIR)io.o $(OBJDIR)magcalc.o $(OBJDIR)temporal.o \
	$(OBJDIR)expanduser.o
	$(FL) -o magcalc $(OBJDIR)expanduser.o $(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o $(OBJDIR)grid.o $(OBJDIR)io.o $(OBJDIR)magcalc.o $(OBJDIR)temporal.o

eq: $(OBJDIR)msis00_gfortran.o $(OBJDIR)call_msis_gfortran.o
	$(FL) -o ./setup/msis $(OBJDIR)msis00_gfortran.o $(OBJDIR)call_msis_gfortran.o

clean:
	$(RM) -r $(OBJDIR)/*

$(OBJDIR)call_msis_gfortran.o: ./setup/MSIS00/call_msis_gfortran.f90
	$(FC) -c $(OPTIM2) ./setup/MSIS00/call_msis_gfortran.f90  -o $(OBJDIR)/call_msis_gfortran.o -J$(OBJDIR)

$(OBJDIR)phys_consts.o: ./numerical/constants/phys_consts.f90
	$(FC) -c $(OPTIM) $^ -o $@ -J$(OBJDIR)

$(OBJDIR)grid.o: ./numerical/grid/grid.f90 $(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o
	$(FC) -c $(OPTIM2) ./numerical/grid/grid.f90  -o $(OBJDIR)/grid.o -J$(OBJDIR)

$(OBJDIR)calculus.o: ./numerical/calculus/calculus.f90 $(OBJDIR)grid.o
	$(FC) -c $(OPTIM2) ./numerical/calculus/calculus.f90 -o $(OBJDIR)/calculus.o -J$(OBJDIR)

$(OBJDIR)interpolation.o: ./numerical/interpolation/interpolation.f90 $(OBJDIR)phys_consts.o
	$(FC) -c $(OPTIM2) ./numerical/interpolation/interpolation.f90 -o $(OBJDIR)/interpolation.o -J$(OBJDIR)

$(OBJDIR)msis00_gfortran.o: ./neutral/neutral.f90
	$(FC) -c $(OPTIM2) $(OPTIONS) ./vendor/msis00/msis00_gfortran.f -o $(OBJDIR)/msis00_gfortran.o -J$(OBJDIR)

$(OBJDIR)neutral.o: ./neutral/neutral.f90 $(OBJDIR)phys_consts.o $(OBJDIR)msis00_gfortran.o $(OBJDIR)temporal.o $(OBJDIR)interpolation.o $(OBJDIR)mpimod.o $(OBJDIR)io.o
	$(FC) -c $(OPTIM2) ./neutral/neutral.f90 -o $(OBJDIR)/neutral.o -J$(OBJDIR)

$(OBJDIR)precipBCs_mod.o: ./ionization/boundary_conditions/precipBCs_mod.f90 $(OBJDIR)grid.o $(OBJDIR)interpolation.o $(OBJDIR)io.o $(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o $(OBJDIR)temporal.o
	$(FC) -c $(OPTIM2) ./ionization/boundary_conditions/precipBCs_mod.f90 -o $(OBJDIR)/precipBCs_mod.o -J$(OBJDIR)

$(OBJDIR)ionization.o: ./ionization/ionization.f90 $(OBJDIR)phys_consts.o $(OBJDIR)calculus.o $(OBJDIR)grid.o $(OBJDIR)neutral.o
	$(FC) -c $(OPTIM2) ./ionization/ionization.f90 -o $(OBJDIR)/ionization.o -J$(OBJDIR)

$(OBJDIR)collisions.o: ./collisions/collisions.f90 $(OBJDIR)phys_consts.o
	$(FC) -c $(OPTIM2) ./collisions/collisions.f90 -o $(OBJDIR)/collisions.o -J$(OBJDIR)

$(OBJDIR)sources.o: ./sources/sources.f90 $(OBJDIR)phys_consts.o $(OBJDIR)collisions.o $(OBJDIR)calculus.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM2) ./sources/sources.f90  -o $(OBJDIR)/sources.o -J$(OBJDIR)

$(OBJDIR)temporal.o: ./temporal/temporal.f90 $(OBJDIR)phys_consts.o $(OBJDIR)mpimod.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM2) ./temporal/temporal.f90 -o $(OBJDIR)/temporal.o -J$(OBJDIR)

$(OBJDIR)diffusion.o: ./numerical/diffusion/diffusion.f90 $(OBJDIR)phys_consts.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM2) ./numerical/diffusion/diffusion.f90 -I$(HOMEDIR)/lib/LAPACK95/lapack95_modules/ -o $(OBJDIR)/diffusion.o -J$(OBJDIR)

$(OBJDIR)mpimod.o: ./numerical/mpimod/mpimod.f90 $(OBJDIR)phys_consts.o
	$(FC) -c $(OPTIM2) $^ -o $@ -J$(OBJDIR)

$(OBJDIR)advec_mpi.o: ./numerical/advection/advec_mpi.f90 $(OBJDIR)phys_consts.o $(OBJDIR)mpimod.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM2) ./numerical/advection/advec_mpi.f90  -o $(OBJDIR)/advec_mpi.o -J$(OBJDIR)

$(OBJDIR)io.o: ./io/io.f90 $(OBJDIR)phys_consts.o $(OBJDIR)mpimod.o $(OBJDIR)calculus.o $(OBJDIR)grid.o $(OBJDIR)expanduser.o
	$(FC) -c $(OPTIM2) ./io/io.f90  -o $(OBJDIR)/io.o -J$(OBJDIR)

$(OBJDIR)expanduser.o: ./io/expanduser.f90
	$(FC) -c $(OPTIM2) ./io/expanduser.f90  -o $(OBJDIR)/expanduser.o -J$(OBJDIR)

$(OBJDIR)potential_mumps.o:  ./numerical/potential/potential_mumps.f90 $(OBJDIR)calculus.o
	$(FC) -c $(OPTIM2) ./numerical/potential/potential_mumps.f90 -o $(OBJDIR)/potential_mumps.o -I$(MUMPSDIR) -I$(MUMPSDIR)/include -I/usr/local/include -J$(OBJDIR)

$(OBJDIR)multifluid.o:  ./multifluid/multifluid.f90 $(OBJDIR)phys_consts.o $(OBJDIR)mpimod.o $(OBJDIR)calculus.o $(OBJDIR)diffusion.o \
	$(OBJDIR)advec_mpi.o $(OBJDIR)ionization.o $(OBJDIR)sources.o $(OBJDIR)grid.o $(OBJDIR)precipBCs_mod.o $(OBJDIR)temporal.o
	$(FC) -c $(OPTIM2) ./multifluid/multifluid.f90 -o $(OBJDIR)/multifluid.o -J$(OBJDIR)

$(OBJDIR)potentialBCs_mumps.o: ./numerical/potential/boundary_conditions/potentialBCs_mumps.f90 $(OBJDIR)grid.o $(OBJDIR)phys_consts.o $(OBJDIR)calculus.o $(OBJDIR)temporal.o $(OBJDIR)interpolation.o $(OBJDIR)io.o
	$(FC) -c $(OPTIM2) ./numerical/potential/boundary_conditions/potentialBCs_mumps.f90  -o $(OBJDIR)/potentialBCs_mumps.o -J$(OBJDIR)

$(OBJDIR)potential_comm_mumps.o: ./numerical/potential/potential_comm_mumps.f90 $(OBJDIR)potential_mumps.o $(OBJDIR)phys_consts.o $(OBJDIR)calculus.o $(OBJDIR)mpimod.o $(OBJDIR)collisions.o $(OBJDIR)grid.o $(OBJDIR)potentialBCs_mumps.o
	$(FC) -c $(OPTIM2) ./numerical/potential/potential_comm_mumps.f90  -o $(OBJDIR)/potential_comm_mumps.o -J$(OBJDIR)

$(OBJDIR)gemini.o: gemini.f90 $(OBJDIR)phys_consts.o $(OBJDIR)potential_comm_mumps.o $(OBJDIR)temporal.o $(OBJDIR)neutral.o \
        $(OBJDIR)io.o $(OBJDIR)mpimod.o $(OBJDIR)multifluid.o $(OBJDIR)grid.o $(OBJDIR)precipBCs_mod.o
	$(FC) -c $(OPTIM2) gemini.f90  -o $(OBJDIR)/gemini.o -J$(OBJDIR)

$(OBJDIR)magcalc.o: ./magcalc.f90 $(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o $(OBJDIR)grid.o $(OBJDIR)io.o $(OBJDIR)temporal.o
	$(FC) -c $(OPTIM2) magcalc.f90 -o $(OBJDIR)/magcalc.o -J$(OBJDIR)
