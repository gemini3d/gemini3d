.DEFAULT_GOAL := gemini_mumps

FC=mpifort #compiler
FL=mpifort #linker

#OPTS = -g -Wall -Wextra -Warray-temporaries -Wconversion -fmax-errors=10 -fimplicit-none -ffree-line-length-0 -fcheck=all -finit-real=nan -fbacktrace -std=f2008

LIBDIR = $(HOME)/source/
MUMPSDIR = $(LIBDIR)/MUMPS_5.1.1/
LAPACKDIR = $(LIBDIR)/lapack-3.8.0/
SCALADIR = $(LIBDIR)/scalapack-2.0.2/
BLASDIR = $(LIBDIR)/lapack-3.8.0/

#Linear Algebra Libraries
LLAPACK = -L$(LAPACKDIR) -llapack
LBLAS = -L$(BLASDIR) -lrefblas
LSCALA = -L$(SCALADIR) -lscalapack

#Matrix ordering libraries
#LSEQMETIS = -lmetis
#DSEQMETIS = -Dmetis

#LPARMETIS = -lparmetis
#DPARMETIS = -Dparmetis

#LSCOTCH = -lesmumps -lscotch -lscotcherr -lptesmumps -lptscotch -lptscotcherr
#DSCOTCH = -Dscotch -Dptscotch

IPORD = -I$(MUMPSDIR)/PORD/include/
LPORD = -L$(MUMPSDIR)/PORD/lib/ -lpord
DPORD = -Dpord

IORDERINGS = $(IPARMETIS) $(ISEQMETIS) $(ISCOTCH) $(IPORD)
LORDERINGS = $(LPARMETIS) $(LSEQMETIS) $(LSCOTCH) $(LPORD)
DORDERINGS = $(DPARMETIS) $(DSEQMETIS) $(DSCOTCH) $(DPORD)

#Parallel solver libraries
LMUMPS = -L$(MUMPSDIR)/lib -ldmumps -lmumps_common
IMUMPS = -I$(MUMPSDIR)/include

#Gemini object output
OBJDIR = ./objects/

OPTIM = -O3 $(OPTS)

OBJS = $(OBJDIR)advec_mpi.o $(OBJDIR)diffusion.o $(OBJDIR)ionization.o $(OBJDIR)phys_consts.o $(OBJDIR)potential_mumps.o \
	$(OBJDIR)calculus.o $(OBJDIR)gemini.o $(OBJDIR)msis00_gfortran.o $(OBJDIR)potentialBCs_mumps.o $(OBJDIR)precipBCs_mod.o \
	$(OBJDIR)temporal.o $(OBJDIR)collisions.o $(OBJDIR)io.o $(OBJDIR)neutral.o $(OBJDIR)potential_comm_mumps.o $(OBJDIR)sources.o  \
	$(OBJDIR)mpimod.o $(OBJDIR)multifluid.o $(OBJDIR)grid.o $(OBJDIR)interpolation.o

all: gemini_mumps magcalc eq

gemini_mumps: $(OBJS)
	$(FL) $(OPTIM) -o gemini_mumps $(OBJDIR)interpolation.o $(OBJDIR)grid.o $(OBJDIR)multifluid.o $(OBJDIR)mpimod.o $(OBJDIR)advec_mpi.o $(OBJDIR)diffusion.o $(OBJDIR)ionization.o $(OBJDIR)phys_consts.o $(OBJDIR)potential_mumps.o $(OBJDIR)calculus.o $(OBJDIR)gemini.o $(OBJDIR)msis00_gfortran.o $(OBJDIR)potentialBCs_mumps.o $(OBJDIR)precipBCs_mod.o $(OBJDIR)temporal.o $(OBJDIR)collisions.o $(OBJDIR)io.o $(OBJDIR)neutral.o $(OBJDIR)potential_comm_mumps.o $(OBJDIR)sources.o $(LMUMPS) $(LORDERINGS) $(DORDERINGS) $(LSCALA) $(LLAPACK) $(LBLAS)

magcalc: $(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o $(OBJDIR)grid.o $(OBJDIR)io.o $(OBJDIR)magcalc.o $(OBJDIR)temporal.o
	$(FL) -o magcalc $(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o $(OBJDIR)grid.o $(OBJDIR)io.o $(OBJDIR)magcalc.o $(OBJDIR)temporal.o

eq: $(OBJDIR)msis00_gfortran.o $(OBJDIR)call_msis_gfortran.o
	$(FL) -o ./setup/msis $(OBJDIR)msis00_gfortran.o $(OBJDIR)call_msis_gfortran.o

clean:
	$(RM) -r $(OBJDIR)/*

$(OBJDIR)call_msis_gfortran.o: ./setup/MSIS00/call_msis_gfortran.f90
	$(FC) -c $(OPTIM) ./setup/MSIS00/call_msis_gfortran.f90  -o $(OBJDIR)/call_msis_gfortran.o -J$(OBJDIR)

$(OBJDIR)phys_consts.o: ./numerical/constants/phys_consts.f90
	$(FC) -c $(OPTIM) ./numerical/constants/phys_consts.f90 -o $(OBJDIR)/phys_consts.o -J$(OBJDIR)

$(OBJDIR)grid.o: ./numerical/grid/grid.f90 $(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o
	$(FC) -c $(OPTIM) ./numerical/grid/grid.f90 -o $(OBJDIR)/grid.o -J$(OBJDIR)

$(OBJDIR)calculus.o: ./numerical/calculus/calculus.f90 $(OBJDIR)phys_consts.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM) ./numerical/calculus/calculus.f90 -o $(OBJDIR)/calculus.o -J$(OBJDIR)

$(OBJDIR)msis00_gfortran.o: ./neutral/msis00/msis00_gfortran.f
	$(FC) -c $(OPTIM) -w ./neutral/msis00/msis00_gfortran.f -o $(OBJDIR)/msis00_gfortran.o -J$(OBJDIR)

$(OBJDIR)interpolation.o: ./numerical/interpolation/interpolation.f90 $(OBJDIR)phys_consts.o
	$(FC) -c $(OPTIM) ./numerical/interpolation/interpolation.f90 -o $(OBJDIR)/interpolation.o -J$(OBJDIR)

$(OBJDIR)neutral.o: ./neutral/neutral.f90 $(OBJDIR)phys_consts.o $(OBJDIR)msis00_gfortran.o $(OBJDIR)temporal.o $(OBJDIR)interpolation.o \
	$(OBJDIR)mpimod.o $(OBJDIR)io.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM) ./neutral/neutral.f90 -o $(OBJDIR)/neutral.o -J$(OBJDIR)

$(OBJDIR)precipBCs_mod.o: ./ionization/boundary_conditions/precipBCs_mod.f90 $(OBJDIR)grid.o $(OBJDIR)interpolation.o $(OBJDIR)io.o \
	$(OBJDIR)mpimod.o $(OBJDIR)phys_consts.o $(OBJDIR)temporal.o
	$(FC) -c $(OPTIM) ./ionization/boundary_conditions/precipBCs_mod.f90 -o $(OBJDIR)/precipBCs_mod.o -J$(OBJDIR)

$(OBJDIR)ionization.o: ./ionization/ionization.f90 $(OBJDIR)phys_consts.o $(OBJDIR)calculus.o $(OBJDIR)grid.o $(OBJDIR)neutral.o \
	$(OBJDIR)temporal.o
	$(FC) -c $(OPTIM) ./ionization/ionization.f90 -o $(OBJDIR)/ionization.o -J$(OBJDIR)

$(OBJDIR)collisions.o: ./collisions/collisions.f90 $(OBJDIR)phys_consts.o
	$(FC) -c $(OPTIM) ./collisions/collisions.f90 -o $(OBJDIR)/collisions.o -J$(OBJDIR)

$(OBJDIR)sources.o: ./sources/sources.f90 $(OBJDIR)phys_consts.o $(OBJDIR)grid.o $(OBJDIR)collisions.o $(OBJDIR)calculus.o $(OBJDIR)mpimod.o
	$(FC) -c $(OPTIM) ./sources/sources.f90 -o $(OBJDIR)/sources.o -J$(OBJDIR)

$(OBJDIR)temporal.o: ./temporal/temporal.f90 $(OBJDIR)phys_consts.o $(OBJDIR)mpimod.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM) ./temporal/temporal.f90 -o $(OBJDIR)/temporal.o -J$(OBJDIR)

$(OBJDIR)diffusion.o: ./numerical/diffusion/diffusion.f90 $(OBJDIR)phys_consts.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM) ./numerical/diffusion/diffusion.f90 $(LLAPACK) -o $(OBJDIR)/diffusion.o -J$(OBJDIR)

$(OBJDIR)mpimod.o: ./numerical/mpimod/mpimod.f90 $(OBJDIR)phys_consts.o
	$(FC) -c $(OPTIM) ./numerical/mpimod/mpimod.f90 -o $(OBJDIR)/mpimod.o -J$(OBJDIR)

$(OBJDIR)advec_mpi.o: ./numerical/advection/advec_mpi.f90 $(OBJDIR)phys_consts.o $(OBJDIR)mpimod.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM) ./numerical/advection/advec_mpi.f90 -o $(OBJDIR)/advec_mpi.o -J$(OBJDIR)

$(OBJDIR)io.o: ./io/io.f90 $(OBJDIR)phys_consts.o $(OBJDIR)mpimod.o $(OBJDIR)calculus.o $(OBJDIR)grid.o
	$(FC) -c $(OPTIM) ./io/io.f90 -o $(OBJDIR)/io.o -J$(OBJDIR)

$(OBJDIR)potential_mumps.o:  ./numerical/potential/potential_mumps.f90 $(OBJDIR)calculus.o
	$(FC) -c $(OPTIM) ./numerical/potential/potential_mumps.f90 $(IMUMPS) -o $(OBJDIR)/potential_mumps.o -J$(OBJDIR)

$(OBJDIR)multifluid.o:  ./multifluid/multifluid.f90 $(OBJDIR)phys_consts.o $(OBJDIR)mpimod.o $(OBJDIR)calculus.o $(OBJDIR)diffusion.o \
	$(OBJDIR)advec_mpi.o $(OBJDIR)ionization.o $(OBJDIR)sources.o $(OBJDIR)grid.o $(OBJDIR)precipBCs_mod.o $(OBJDIR)temporal.o
	$(FC) -c $(OPTIM) ./multifluid/multifluid.f90 -o $(OBJDIR)/multifluid.o -J$(OBJDIR)

$(OBJDIR)potentialBCs_mumps.o: ./numerical/potential/boundary_conditions/potentialBCs_mumps.f90 $(OBJDIR)grid.o $(OBJDIR)phys_consts.o $(OBJDIR)calculus.o \
	$(OBJDIR)io.o $(OBJDIR)temporal.o $(OBJDIR)interpolation.o
	$(FC) -c $(OPTIM) ./numerical/potential/boundary_conditions/potentialBCs_mumps.f90 -o $(OBJDIR)/potentialBCs_mumps.o -J$(OBJDIR)

$(OBJDIR)potential_comm_mumps.o: ./numerical/potential/potential_comm_mumps.f90 $(OBJDIR)potential_mumps.o $(OBJDIR)phys_consts.o $(OBJDIR)calculus.o \
	$(OBJDIR)mpimod.o $(OBJDIR)collisions.o $(OBJDIR)grid.o $(OBJDIR)potentialBCs_mumps.o
	$(FC) -c $(OPTIM) ./numerical/potential/potential_comm_mumps.f90 -o $(OBJDIR)/potential_comm_mumps.o -J$(OBJDIR)

$(OBJDIR)gemini.o: gemini.f90 $(OBJDIR)phys_consts.o $(OBJDIR)potential_comm_mumps.o $(OBJDIR)temporal.o $(OBJDIR)neutral.o \
        $(OBJDIR)io.o $(OBJDIR)mpimod.o $(OBJDIR)multifluid.o $(OBJDIR)grid.o $(OBJDIR)precipBCs_mod.o $(OBJDIR)potentialBCs_mumps.o \
	$(OBJDIR)glow_run.o
	$(FC) -c $(OPTIM) gemini.f90 -o $(OBJDIR)/gemini.o -J$(OBJDIR)

$(OBJDIR)magcalc.o: ./magcalc.f90 $(OBJDIR)phys_consts.o $(OBJDIR)grid.o $(OBJDIR)temporal.o $(OBJDIR)io.o $(OBJDIR)mpimod.o
	$(FC) -c $(OPTIM) magcalc.f90 -o $(OBJDIR)/magcalc.o -J$(OBJDIR)

$(OBJDIR)glow_run.o: ./ionization/glow_run.f90
	$(FC) -c $(OPTIM) ./ionization/glow_run.f90 -o $(OBJDIR)/glow_run.o -J$(OBJDIR)
