module potentialBCs_nompi

use phys_consts, only: wp, pi, Re, debug
use gemini3d_config, only: gemini_cfg
use grid, only: lx1, lx2, lx3, lx2all, lx3all, gridflag
use efielddataobj, only: efielddata
use meshobj, only: curvmesh

implicit none (type, external)
private
public :: init_Efieldinput_nompi, potentialBCs2D_fileinput_nompi, compute_BGEfields_nompi

contains
  subroutine init_Efieldinput_nompi(dt,cfg,ymd,UTsec,x,efield)
    !> Initialize variables to hold electric field input file data, can be called by any worker but only root does anything
    !    We need some alternate code to deal with the situation where we are running without mpi...
    real(wp), intent(in) :: dt
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    type(efielddata), intent(inout) :: efield

    !> initializes the auroral electric field/current and particle inputs to read in a file corresponding to the first time step
    !    In this case don't check root v. all; assuming calling function knows already
    !if ( (mpi_cfg%myid==0 .and. cfg%flagE0file==1) .or. (.not. efield%flagrootonly) ) then    !only root needs these...
      call efield%init(cfg,cfg%E0dir,x,dt,cfg%dtE0,ymd,UTsec)
    !end if
  end subroutine init_Efieldinput_nompi


  !> sidewall etc. boundaries do not need to be set in the nompi case, particularly if we are doing an actual solve
  subroutine potentialBCs2D_fileinput_nompi(dtmodel,t,ymd,UTsec,cfg,x,efield,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                    Vmaxx3,E01,E02,E03,flagdirich)
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
    real(wp), dimension(:,:,:), intent(inout) :: E01,E02,E03
    !! intent(out)
    integer, intent(out) :: flagdirich
    integer :: ix1,ix2,ix3


    !> COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
    E01 = 0
    !! do not allow a background parallel field

    !! tell the efield data object to update
    call efield%update(cfg,dtmodel,t,x,ymd,UTsec)

    !! now we need to take data in the efield object and map along the field lines
    call compute_BGEfields_nompi(x,E02,E03,efield)

    !! object double flag to int
    flagdirich=nint(efield%flagdirich)

    !! set boundary condition output arguments based on object data
    !!   Note that we are effectively making copies of object properties; however, additional processing is
    !!   being done here that is specific to the potential solver so that may be justified extra memory use...
    do ix3=1,lx3
      do ix2=1,lx2
        Vminx1(ix2,ix3)=efield%Vminx1inow(ix2,ix3)
        Vmaxx1(ix2,ix3)=efield%Vmaxx1inow(ix2,ix3)
      end do
    end do

    !! This forces certain patterns to the boundary conditions to make sure solvers don't get garbage data...
    if (lx2/=1 .and. lx3/=1) then
      ! full 3D grid
      do ix3=1,lx3
        do ix1=1,lx1
          Vminx2(ix1,ix3)=efield%Vminx2isnow(ix3)
          Vmaxx2(ix1,ix3)=efield%Vmaxx2isnow(ix3)
        end do
      end do

      do ix2=1,lx2
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
          if (lx2==1) then
            do ix1=1,lx1
              Vminx3(ix1,:)=Vmaxx1(:,1)
              Vmaxx3(ix1,:)=Vmaxx1(:,lx3)
            end do
            Vminx2=0
            Vmaxx2=0
          else
            do ix1=1,lx1
              Vminx2(ix1,:)=Vmaxx1(1,:)
              Vmaxx2(ix1,:)=Vmaxx1(lx2,:)
            end do
            Vminx3=0
            Vmaxx3=0
          end if
        else                    !inverted so logical beginning is max alt.
          if (lx2==1) then
            do ix1=1,lx1
              Vminx3(ix1,:)=Vminx1(:,1)
              Vmaxx3(ix1,:)=Vminx1(:,lx3)
            end do
            Vminx2=0
            Vmaxx2=0
          else
            do ix1=1,lx1
              Vminx2(ix1,:)=Vminx1(1,:)
              Vmaxx2(ix1,:)=Vminx1(lx2,:)
            end do
            Vminx3=0
            Vmaxx3=0
          end if
        end if
      else
        ! Neumann in x1:  sides are grounded...
        if (lx2==1) then
          Vminx3 = 0
          Vmaxx3 = 0
        else
          Vminx2=0
          Vmaxx2=0
        end if
      end if
    end if
  end subroutine potentialBCs2D_fileinput_nompi


  subroutine compute_BGEfields_nompi(x,E02,E03,efield)
    !> Returns a background electric field calculation for use by external program units.
    !   This requires that all necessary files, etc. have already been loaded into module
    !   variables.  This is only to be called by a root process as it deals with fullgrid
    !   data.  An interface for workers and root is in the top-level potential module.  This
    !   particular bit of code is needed both when setting boundary conditions and also when
    !   initializing background electric field; hence it is a subroutine as opposed to block of code
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:), intent(inout) :: E02,E03
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
    if (lx2 > 1 .and. lx3>1) then ! 3D sim
      ix2ref = lx2/2      !note integer division
      ix3ref = lx3/2
    else if (lx2==1 .and. lx3>1) then
      ix2ref = 1
      ix3ref=lx3/2
    else if (lx2>1 .and. lx3==1) then
      ix2ref=lx2/2
      ix3ref=1
    else
      error stop 'Unable to orient boundary conditions for electric potential'
    endif

    !! by default the code uses 300km altitude as a reference location, using the center x2,x3 point
    !! These are the coordinates for inputs varying along axes 2,3
    ix1ref = minloc(abs(x%r(:,ix2ref,ix3ref) - Re - 300e3_wp), dim=1)

    !! scale electric fields at some reference point into the full grid
    do ix3=1,lx3
      do ix2=1,lx2
        h2ref=x%h2(ix1ref,ix2,ix3)
        !! define a reference metric factor for a given field line
        h3ref=x%h3(ix1ref,ix2,ix3)
        do ix1=1,lx1
          E02(ix1,ix2,ix3)=efield%E0xinow(ix2,ix3)*h2ref/x%h2(ix1,ix2,ix3)
          E03(ix1,ix2,ix3)=efield%E0yinow(ix2,ix3)*h3ref/x%h3(ix1,ix2,ix3)
        end do
      end do
    end do
  end subroutine compute_BGEfields_nompi
end module
