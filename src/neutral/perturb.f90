submodule (neutral) perturb

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use grid, only : gridflag
use reader, only : get_neutral2
use timeutils, only : dateinc, date_filename
use mpimod, only: mpi_realprec, mpi_integer, mpi_comm_world, mpi_status_ignore, &
tag=>gemini_mpi
use h5fortran, only: hdf5_file
use pathlib, only: get_suffix,get_filename

implicit none (type, external)
external :: mpi_send, mpi_recv

interface ! proj.f90 submodule
  module subroutine gridproj_dneu2D(cfg,flagcart,x)
    type(gemini_cfg), intent(in) :: cfg
    logical, intent(in) :: flagcart
    class(curvmesh), intent(inout) :: x
  end subroutine gridproj_dneu2D

  module subroutine gridproj_dneu3D(cfg,x)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(inout) :: x
  end subroutine gridproj_dneu3D
end interface


interface ! interp.f90 submodule
  module subroutine spaceinterp_dneu2D(flagcart)
    logical, intent(in) :: flagcart
  end subroutine spaceinterp_dneu2D

  module subroutine spaceinterp_dneu3D()
  end subroutine spaceinterp_dneu3D

  module subroutine timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)
    real(wp), intent(in) :: t,dt
    real(wp), dimension(:,:,:), intent(inout) :: dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow
    !! intent(out)
  end subroutine timeinterp_dneu
end interface

contains
  !>  top-level procedure for computing neutral perturbations
  module procedure neutral_perturb
    ! advance object state
    call atmosperturb%update(cfg,dt,t,x,ymd,UTsec)
  
    ! copy object data into module variables used in neutral_update()
    ! FIXME: these need to be pointers to avoid wasting space by duplicating
    dnOinow=atmosperturb%dnOinow
    dnN2inow=atmosperturb%dnN2inow
    dnO2inow=atmosperturb%dnO2inow
    dvn1inow=atmosperturb%dvn1inow
    dvn2inow=atmosperturb%dvn2inow
    dvn3inow=atmosperturb%dvn3inow
    dTninow=atmosperturb%dTninow
    
    !Add interpolated perturbations to module reference atmosphere arrays
    call neutral_update(nn,Tn,vn1,vn2,vn3,v2grid,v3grid)
  end procedure neutral_perturb
end submodule perturb
