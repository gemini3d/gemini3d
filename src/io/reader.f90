module reader
!! simple file reading procedures

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

use phys_consts, only: wp, debug
use filesystem, only : is_file

implicit none (type, external)
private
public :: get_simsize3, get_simsize2, get_grid2, get_Efield, get_precip, get_neutral2, get_neutral3

interface !< reader_{raw,hdf5,nc4}.f90
module subroutine get_simsize2_hdf5(path, llon, llat)
character(*), intent(in) :: path
integer, intent(out) :: llon, llat
end subroutine

module subroutine get_simsize3_hdf5(path, lx1, lx2all, lx3all)
character(*), intent(in) :: path
integer, intent(out) :: lx1, lx2all
integer, intent(out), optional :: lx3all
end subroutine

module subroutine get_grid2_hdf5(path, mlonp, mlatp)
character(*), intent(in) :: path
real(wp), dimension(:), intent(inout) :: mlonp, mlatp
!! intent(out)
end subroutine

module subroutine get_Efield_hdf5(path, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
character(*), intent(in) :: path
integer, intent(out) :: flagdirich
real(wp), dimension(:,:), intent(inout) :: E0xp,E0yp,Vminx1p,Vmaxx1p
!! intent(out)
real(wp), dimension(:), intent(inout) :: Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice
!! intent(out)
end subroutine

module subroutine get_precip_hdf5(path, Qp, E0p)
character(*), intent(in) :: path
real(wp), dimension(:,:), intent(inout) :: Qp, E0p
!! intent(out)
end subroutine

module subroutine get_neutral2_hdf5(path, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(inout) :: dnO,dnN2,dnO2,dvnrho,dvnz,dTn
!! intent(out)
end subroutine

module subroutine get_neutral3_hdf5(path, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(inout) :: dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall
!! intent(out)
end subroutine

end interface

contains

subroutine get_simsize2(path, llon, llat)
character(*), intent(in) :: path
integer, intent(out) :: llon, llat

call get_simsize2_hdf5(path, llon, llat)

end subroutine get_simsize2


subroutine get_simsize3(path, lx1, lx2all, lx3all)
character(*), intent(in) :: path
integer, intent(out) :: lx1, lx2all
integer, intent(out), optional :: lx3all

call get_simsize3_hdf5(path, lx1, lx2all, lx3all)

end subroutine get_simsize3


subroutine get_grid2(path, mlonp, mlatp)
character(*), intent(in) :: path
real(wp), dimension(:), intent(inout) :: mlonp, mlatp
!! intent(out)

call get_grid2_hdf5(path, mlonp, mlatp)

end subroutine get_grid2


subroutine get_Efield(path, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
character(*), intent(in) :: path
integer, intent(out) :: flagdirich
real(wp), dimension(:,:), intent(inout) :: E0xp,E0yp,Vminx1p,Vmaxx1p
!! intent(out)
real(wp), dimension(:), intent(inout) :: Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice
!! intent(out)

call get_Efield_hdf5(path, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)

end subroutine get_Efield


subroutine get_precip(path, Qp, E0p)
!! Qp, E0p are (llon, llat)
character(*), intent(in) :: path
real(wp), dimension(:,:), intent(inout) :: Qp, E0p
!! intent(out)

integer :: i
character(:), allocatable :: fn

fn = path
i = len_trim(path)

if (.not.is_file(fn)) fn(i:i) = '1'
!! workaround for old files like *.000001.h5
if (.not. is_file(fn)) error stop "reader:precip No file found on " // path

call get_precip_hdf5(fn, Qp, E0p)

!> sanity check precipitation input
if(.not. all(ieee_is_finite(Qp))) error stop 'precipBCs_mod.f90:precipBCs_fileinput: Qp must be finite'

if(any(Qp < 0)) error stop 'precipBCs_mod.f90:precipBCs_fileinput: Qp must be non-negative'

if(any(Qp >= 1e-6_wp)) then
  !! if flux \equiv 0 for this time step, we don't care what E0 is
  !! We use E0 = NaN to indicate this is a time where particle flux is not specified
  if(.not. all(ieee_is_finite(E0p))) error stop 'precipBCs_mod.f90:precipBCs_fileinput: E0p must be finite'

  if(any(E0p < 0)) error stop 'precipBCs_mod.f90:precipBCs_fileinput: E0p must be non-negative'
endif

end subroutine get_precip


subroutine get_neutral2(path, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(inout) :: dnO,dnN2,dnO2,dvnrho,dvnz,dTn
!! intent(out)

call get_neutral2_hdf5(path, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)

end subroutine get_neutral2


subroutine get_neutral3(path, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(inout) :: dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall
!! intent(out)

call get_neutral3_hdf5(path, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)

end subroutine get_neutral3


end module reader
