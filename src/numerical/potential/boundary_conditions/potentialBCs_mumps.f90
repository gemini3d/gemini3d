module potentialBCs_mumps

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

use mpimod, only : mpi_cfg
use phys_consts, only: wp, pi, Re, debug
use grid, only: lx1, lx2, lx2all, lx3all, gridflag
use meshobj, only: curvmesh
use interpolation, only : interp1,interp2
use timeutils, only : dateinc, date_filename, find_lastdate
use reader, only : get_grid2, get_simsize2, get_Efield
use gemini3d_config, only: gemini_cfg
use efielddataobj, only: efielddata

implicit none (type, external)
private
public :: potentialbcs2D, potentialbcs2D_fileinput, init_Efieldinput, &
            compute_rootBGEfields
contains
  subroutine init_Efieldinput(dt,cfg,ymd,UTsec,x,efield)
    !> Initialize variables to hold electric field input file data, can be called by any worker but only root does anything
    !    We need some alternate code to deal with the situation where we are running without mpi...
    real(wp), intent(in) :: dt
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    type(efielddata), intent(inout) :: efield

    !> initializes the auroral electric field/current and particle inputs to read in a file corresponding to the first time step
    if (mpi_cfg%myid==0 .and. cfg%flagE0file==1) then    !only root needs these...
      call efield%init(cfg,cfg%E0dir,x,dt,cfg%dtE0,ymd,UTsec)
    end if
  end subroutine init_Efieldinput


  subroutine potentialBCs2D_fileinput(dtmodel,t,ymd,UTsec,cfg,x,efield,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                    Vmaxx3,E01all,E02all,E03all,flagdirich)
    !! A FILE INPUT BASED BOUNDARY CONDITIONS FOR ELECTRIC POTENTIAL OR
    !! FIELD-ALIGNED CURRENT.
    !! NOTE: THIS IS ONLY CALLED BY THE ROOT PROCESS
    real(wp), intent(in) :: dtmodel
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    type(efielddata), intent(inout) :: efield
    real(wp), dimension(:,:), intent(inout), target :: Vminx1,Vmaxx1
    !! intent(out)
    real(wp), dimension(:,:), intent(inout) :: Vminx2,Vmaxx2
    !! intent(out)
    real(wp), dimension(:,:), intent(inout) :: Vminx3,Vmaxx3
    !! intent(out)
    real(wp), dimension(:,:,:), intent(inout) :: E01all,E02all,E03all
    !! intent(out)
    integer, intent(out) :: flagdirich
    integer :: ix1,ix2,ix3


    !> COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
    E01all = 0
    !! do not allow a background parallel field

    !! tell the efield data object to update
    call efield%update(cfg,dtmodel,t,x,ymd,UTsec)

    !! now we need to take data in the efield object and map along the field lines
    call compute_rootBGEfields(x,E02all,E03all,efield)

    !! object double flag to int
    flagdirich=nint(efield%flagdirich)

    !! set boundary condition output arguments based on object data
    !!   Note that we are effectively making copies of object properties; however, additional processing is
    !!   being done here that is specific to the potential solver so that may be justified extra memory use...
    do ix3=1,lx3all
      do ix2=1,lx2all
        Vminx1(ix2,ix3)=efield%Vminx1inow(ix2,ix3)
        Vmaxx1(ix2,ix3)=efield%Vmaxx1inow(ix2,ix3)
      end do
    end do

    !! This forces certain patterns to the boundary conditions to make sure solvers don't get garbage data...
    if (lx2all/=1 .and. lx3all/=1) then
      ! full 3D grid
      do ix3=1,lx3all
        do ix1=1,lx1
          Vminx2(ix1,ix3)=efield%Vminx2isnow(ix3)
          Vmaxx2(ix1,ix3)=efield%Vmaxx2isnow(ix3)
        end do
      end do

      do ix2=1,lx2all
        do ix1=1,lx1
          Vminx3(ix1,ix2)=efield%Vminx3isnow(ix2)
          Vmaxx3(ix1,ix2)=efield%Vmaxx3isnow(ix2)
        end do
      end do
    else
      ! some type of 2D grid, lateral boundary will be overwritten
      Vminx2 = 0
      Vmaxx2 = 0
      if (flagdirich==1) then
        ! Dirichlet:  needs to be the same as the physical top corner grid points
        if (gridflag/=1) then   !non-inverted so logical end is the max alt.
          if (lx2all==1) then
            do ix1=1,lx1
              Vminx3(ix1,:)=Vmaxx1(:,1)
              Vmaxx3(ix1,:)=Vmaxx1(:,lx3all)
            end do
            Vminx2=0
            Vmaxx2=0
          else
            do ix1=1,lx1
              Vminx2(ix1,:)=Vmaxx1(1,:)
              Vmaxx2(ix1,:)=Vmaxx1(lx2all,:)
            end do
            Vminx3=0
            Vmaxx3=0
          end if
        else                    !inverted so logical beginning is max alt.
          if (lx2all==1) then
            do ix1=1,lx1
              Vminx3(ix1,:)=Vminx1(:,1)
              Vmaxx3(ix1,:)=Vminx1(:,lx3all)
            end do
            Vminx2=0
            Vmaxx2=0
          else
            do ix1=1,lx1
              Vminx2(ix1,:)=Vminx1(1,:)
              Vmaxx2(ix1,:)=Vminx1(lx2all,:)
            end do
            Vminx3=0
            Vmaxx3=0
          end if
        end if
      else
        ! Neumann in x1:  sides are grounded...
        if (lx2all==1) then
          Vminx3 = 0
          Vmaxx3 = 0
        else
          Vminx2=0
          Vmaxx2=0
        end if
      end if
    end if
  end subroutine potentialBCs2D_fileinput


  subroutine compute_rootBGEfields(x,E02all,E03all,efield)
    !> Returns a background electric field calculation for use by external program units.
    !   This requires that all necessary files, etc. have already been loaded into module
    !   variables.  This is only to be called by a root process as it deals with fullgrid
    !   data.  An interface for workers and root is in the top-level potential module.  This
    !   particular bit of code is needed both when setting boundary conditions and also when
    !   initializing background electric field; hence it is a subroutine as opposed to block of code
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:), intent(inout) :: E02all,E03all
    !! intent(out)
    type(efielddata), intent(inout) :: efield
    integer :: ix1,ix2,ix3
    real(wp) :: h2ref,h3ref
    integer :: ix1ref,ix2ref,ix3ref    ! reference locations for field line mapping

    !! the only danger here is that this routine could be called before any module data are loaded
    !   so check just to make sure it isn't being misused in this way
    !!    FIXME: does this accomplish anything???
    if (.not. associated(efield%E0xinow)) error stop  &
          'potentialBCs:compute_rootBGEfields is trying to access unallocated module data'

    !! recompute reference locations here (also computed in object)
    if (lx2all > 1 .and. lx3all>1) then ! 3D sim
      ix2ref = lx2all/2      !note integer division
      ix3ref = lx3all/2
    else if (lx2all==1 .and. lx3all>1) then
      ix2ref = 1
      ix3ref=lx3all/2
    else if (lx2all>1 .and. lx3all==1) then
      ix2ref=lx2all/2
      ix3ref=1
    else
      error stop 'Unable to orient boundary conditions for electric potential'
    endif

    !! by default the code uses 300km altitude as a reference location, using the center x2,x3 point
    !! These are the coordinates for inputs varying along axes 2,3
    ix1ref = minloc(abs(x%rall(:,ix2ref,ix3ref) - Re - 300e3_wp), dim=1)

    !! scale electric fields at some reference point into the full grid
    do ix3=1,lx3all
      do ix2=1,lx2all
        h2ref=x%h2all(ix1ref,ix2,ix3)
        !! define a reference metric factor for a given field line
        h3ref=x%h3all(ix1ref,ix2,ix3)
        do ix1=1,lx1
          E02all(ix1,ix2,ix3)=efield%E0xinow(ix2,ix3)*h2ref/x%h2all(ix1,ix2,ix3)
          E03all(ix1,ix2,ix3)=efield%E0yinow(ix2,ix3)*h3ref/x%h3all(ix1,ix2,ix3)
        end do
      end do
    end do
  end subroutine compute_rootBGEfields


  subroutine potentialBCs2D(UTsec,cfg,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                        Vmaxx3,E01all,E02all,E03all,flagdirich)
    ! This is a default routine for setting electromagnetic boundary conditions in cases where user file input is not specified.  It also computes the equatorial vertical drift for the EIA if requested by the user. This routine *could* be modified to hard-code specific conditions in if needed but we really recommend using file input for that.
    real(wp), intent(in) :: UTsec
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:), intent(inout), target :: Vminx1,Vmaxx1
    !! intent(out)
    real(wp), dimension(:,:), intent(inout) :: Vminx2,Vmaxx2
    !! intent(out)
    real(wp), dimension(:,:), intent(inout) :: Vminx3,Vmaxx3
    !! intent(out)
    real(wp), dimension(:,:,:), intent(inout) :: E01all,E02all,E03all
    !! intent(out)
    integer, intent(out) :: flagdirich
    real(wp) :: Phipk
    integer :: ix1,ix2,ix3    !grid sizes are borrowed from grid module
    !    integer, parameter :: lmodes=8; this type of thing done from input scripts now...
    real(wp) :: meanx2,sigx2,meanx3,sigx3 !for setting background field
    real(wp), dimension(:,:), pointer :: Vtopalt,Vbotalt
    real(wp) :: vamp,LThrs,veltime,z,glonmer
    integer :: ix1eq=-1                ! index for the equatorial location in terms of index into the x%x1 array; used by default boundary conditions

    !CALCULATE/SET TOP BOUNDARY CONDITIONS
    sigx2=1/20._wp*(x%x2all(lx2all)-x%x2all(1))
    meanx2=0.5_wp*(x%x2all(1)+x%x2all(lx2all))
    sigx3=1/20._wp*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
    meanx3=0.5_wp*(x%x3all(1)+x%x3all(lx3all))


    ! FIXME: the pointer swapping to deal with top vs. bottom here is confusing/superfluous; it may be better simply to have the input preparation scripts assign things accordingly, which is what we will do from now on.  For this routine right now it doesn't matter since both just zeroed out anyway...
    if (gridflag/=2) then
      Vtopalt=>Vminx1
      Vbotalt=>Vmaxx1
    else
      Vtopalt=>Vmaxx1
      Vbotalt=>Vminx1
    end if

    Phipk = 0      !pk current density
    flagdirich = 0    !Neumann conditions
    do ix3=1,lx3all
      do ix2=1,lx2all
        Vtopalt(ix2,ix3) = 0
      end do
    end do


    !SOME USER INFO
    if (debug) print *, 'At time (UT seconds):  ',UTsec,'  Max FAC set to be:  ',maxval(abs(Vtopalt))


    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
    Vbotalt = 0   !since we need to have no current through bottom boundary
    Vminx2 = 0
    Vmaxx2 = 0
    Vminx3 = 0
    Vmaxx3 = 0


    !PI's EIA code COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
    if (cfg%flagEIA) then
      if (ix1eq<=0) then   !recompute the position of the equator in terms of the x1 variable
        ix1eq = minloc(abs(x%x1), dim=1)    !equator location is that closest to zero in the x1 (q) variable
        if (debug) print*, 'equator ix1:  ',ix1eq,x%x1(ix1eq)
      end if

      vamp=cfg%v0equator    !amplitude of vertical drift at equator from input config.nml file

      E01all=0
      E02all=0

      do ix2=1,lx2all    !for a swapped grid this is longitude
        !for each meridional slice define a local time
        glonmer=x%glonall(ix1eq,ix2,lx3all/2)     !just use halfway up in altitude at the magnetic equator
        do while (glonmer<0)
          glonmer=glonmer+360
        end do

        LThrs=UTsec/3600+glonmer/15                 !Local time of center of meridian
        veltime = sin(2*pi*(LThrs-7)/24)    ! Huba's formulate for velocity amplitude vs. time

        do ix3=1,lx3all     !here this is L-shell
          z = x%altall(ix1eq,ix2,ix3)  !Current altitude of center of this flux tube
          do ix1=1,lx1
            if (z<=150e3_wp) then
              E03all(ix1,ix2,ix3) = 0
            elseif ((z>=150e3_wp) .and. (z<=300e3_wp)) then
              E03all(ix1,ix2,ix3) = (veltime*vamp*(z-150e3_wp)/150e3_wp)*x%Bmagall(ix1eq,ix2,ix3)    !minus sign to deal with permuted dimensions
            elseif (z>300e3_wp) then
              E03all(ix1,ix2,ix3) = veltime*vamp*x%Bmagall(ix1eq,ix2,ix3)
            end if
          end do
        end do
      end do

    !  print*, '  Applied EIA perturbation to background electric field...'
    !  print*, '    ',minval(E03all),maxval(E03all)
    else
      E01all = 0
      E02all = 0
      E03all = 0
    end if
  end subroutine potentialBCs2D
end module potentialBCs_mumps
