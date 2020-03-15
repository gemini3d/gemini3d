module reader
!! simple file reading procedures

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp, debug
use pathlib, only : get_suffix

implicit none
private
public :: get_simsize3, get_simsize2, get_grid2, get_Efield, get_precip, get_neutral2, get_neutral3

interface ! reader_{raw,hdf5,nc4}.f90
module subroutine get_simsize2_raw(path, llon, llat)
character(*), intent(in) :: path
integer, intent(out) :: llon, llat
end subroutine get_simsize2_raw

module subroutine get_simsize2_hdf5(path, llon, llat)
character(*), intent(in) :: path
integer, intent(out) :: llon, llat
end subroutine get_simsize2_hdf5

module subroutine get_simsize2_nc4(path, llon, llat)
character(*), intent(in) :: path
integer, intent(out) :: llon, llat
end subroutine get_simsize2_nc4


module subroutine get_simsize3_raw(path, lx1, lx2all, lx3all)
character(*), intent(in) :: path
integer, intent(out) :: lx1, lx2all
integer, intent(out), optional :: lx3all
end subroutine get_simsize3_raw

module subroutine get_simsize3_hdf5(path, lx1, lx2all, lx3all)
character(*), intent(in) :: path
integer, intent(out) :: lx1, lx2all
integer, intent(out), optional :: lx3all
end subroutine get_simsize3_hdf5

module subroutine get_simsize3_nc4(path, lx1, lx2all, lx3all)
character(*), intent(in) :: path
integer, intent(out) :: lx1, lx2all
integer, intent(out), optional :: lx3all
end subroutine get_simsize3_nc4


module subroutine get_grid2_raw(path, mlonp, mlatp)
character(*), intent(in) :: path
real(wp), dimension(:), intent(out) :: mlonp, mlatp
end subroutine get_grid2_raw

module subroutine get_grid2_hdf5(path, mlonp, mlatp)
character(*), intent(in) :: path
real(wp), dimension(:), intent(out) :: mlonp, mlatp
end subroutine get_grid2_hdf5

module subroutine get_grid2_nc4(path, mlonp, mlatp)
character(*), intent(in) :: path
real(wp), dimension(:), intent(out) :: mlonp, mlatp
end subroutine get_grid2_nc4


module subroutine get_Efield_raw(path, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
character(*), intent(in) :: path
integer, intent(out) :: flagdirich
real(wp), dimension(:,:), intent(out) :: E0xp,E0yp,Vminx1p,Vmaxx1p
real(wp), dimension(:), intent(out) :: Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice
end subroutine get_Efield_raw

module subroutine get_Efield_hdf5(path, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
character(*), intent(in) :: path
integer, intent(out) :: flagdirich
real(wp), dimension(:,:), intent(out) :: E0xp,E0yp,Vminx1p,Vmaxx1p
real(wp), dimension(:), intent(out) :: Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice
end subroutine get_Efield_hdf5

module subroutine get_Efield_nc4(path, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
character(*), intent(in) :: path
integer, intent(out) :: flagdirich
real(wp), dimension(:,:), intent(out) :: E0xp,E0yp,Vminx1p,Vmaxx1p
real(wp), dimension(:), intent(out) :: Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice
end subroutine get_Efield_nc4


module subroutine get_precip_raw(path, Qp, E0p)
character(*), intent(in) :: path
real(wp), dimension(:,:), intent(out) :: Qp, E0p
end subroutine get_precip_raw

module subroutine get_precip_hdf5(path, Qp, E0p)
character(*), intent(in) :: path
real(wp), dimension(:,:), intent(out) :: Qp, E0p
end subroutine get_precip_hdf5

module subroutine get_precip_nc4(path, Qp, E0p)
character(*), intent(in) :: path
real(wp), dimension(:,:), intent(out) :: Qp, E0p
end subroutine get_precip_nc4


module subroutine get_neutral2_raw(path, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnO,dnN2,dnO2,dvnrho,dvnz,dTn
end subroutine get_neutral2_raw

module subroutine get_neutral2_hdf5(path, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnO,dnN2,dnO2,dvnrho,dvnz,dTn
end subroutine get_neutral2_hdf5

module subroutine get_neutral2_nc4(path, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnO,dnN2,dnO2,dvnrho,dvnz,dTn
end subroutine get_neutral2_nc4


module subroutine get_neutral3_raw(path, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall
end subroutine get_neutral3_raw

module subroutine get_neutral3_hdf5(path, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall
end subroutine get_neutral3_hdf5

module subroutine get_neutral3_nc4(path, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall
end subroutine get_neutral3_nc4

end interface

contains

function get_filename(path, stem) result(fn)
!! given a path and stem, find the full filename
!! assuming suffix is
character(*), intent(in) :: path
character(*), intent(in), optional :: stem

character(:), allocatable :: fn, stem1
integer :: i
logical :: exists
character(*), parameter :: suffix(3) = [character(4) :: '.h5', '.nc', '.dat' ]

if(present(stem)) then
  stem1 = '/' // stem
  !> find suffix location, if it exists
  i = index(path, '.', back=.true.)

  !> find filename if path specified
  if (i > 8 .and. path(i-7:i-1) == stem) then
    fn = path  !< assuming it's the full path
    return
  endif
else
  stem1 = ''
endif

do i = 1, size(suffix)
  fn = path // stem1 // trim(suffix(i))
  inquire(file=fn, exist=exists)
  ! print *, 'reader:get_filename: ', fn, exists
  if (exists) return
enddo

write(stderr,*) 'ERROR: filename not found in ', path
fn = ''

end function get_filename


subroutine get_simsize2(path, llon, llat)
character(*), intent(in) :: path
integer, intent(out) :: llon, llat

character(:), allocatable :: fn

fn = get_filename(path, 'simsize')

select case (get_suffix(fn))
case ('.h5')
  call get_simsize2_hdf5(fn, llon, llat)
case ('.nc')
  call get_simsize2_nc4(fn, llon, llat)
case ('.dat')
  call get_simsize2_raw(fn, llon, llat)
case default
  write(stderr,*) 'reader:get_simsize2: unknown file suffix on ' // fn
  error stop 6
end select
end subroutine get_simsize2


subroutine get_simsize3(path, lx1, lx2all, lx3all)
character(*), intent(in) :: path
integer, intent(out) :: lx1, lx2all
integer, intent(out), optional :: lx3all

select case (get_suffix(path))
case ('.h5')
  call get_simsize3_hdf5(path, lx1, lx2all, lx3all)
case ('.nc')
  call get_simsize3_nc4(path, lx1, lx2all, lx3all)
case ('.dat')
  call get_simsize3_raw(path, lx1, lx2all, lx3all)
case default
  write(stderr,*) 'reader:get_simsize3: unknown file suffix' // get_suffix(path)
  error stop 6
end select
end subroutine get_simsize3


subroutine get_grid2(path, mlonp, mlatp)
character(*), intent(in) :: path
real(wp), dimension(:), intent(out) :: mlonp, mlatp

character(:), allocatable :: fn

fn = get_filename(path, 'simgrid')

select case (get_suffix(fn))
case ('.h5')
  call get_grid2_hdf5(fn, mlonp, mlatp)
case ('.nc')
  call get_grid2_nc4(fn, mlonp, mlatp)
case ('.dat')
  call get_grid2_raw(fn, mlonp, mlatp)
case default
  write(stderr,*) 'reader:get_grid2: unknown file suffix on ' // fn
  error stop 6
end select
end subroutine get_grid2


subroutine get_Efield(path, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
character(*), intent(in) :: path
integer, intent(out) :: flagdirich
real(wp), dimension(:,:), intent(out) :: E0xp,E0yp,Vminx1p,Vmaxx1p
real(wp), dimension(:), intent(out) :: Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice

character(:), allocatable :: fn

fn = get_filename(path)

select case (get_suffix(fn))
case ('.h5')
  call get_Efield_hdf5(fn, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
case ('.nc')
  call get_Efield_nc4(fn, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
case ('.dat')
  call get_Efield_raw(fn, flagdirich,E0xp,E0yp,Vminx1p,Vmaxx1p,Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)
case default
  write(stderr,*) 'reader:Efield: unknown file suffix on ' // fn
  error stop 6
end select

end subroutine get_Efield


subroutine get_precip(path, Qp, E0p)
character(*), intent(in) :: path
real(wp), dimension(:,:), intent(out) :: Qp, E0p

character(:), allocatable :: fn

fn = get_filename(path)

select case (get_suffix(fn))
case ('.h5')
  call get_precip_hdf5(fn, Qp, E0p)
case ('.nc')
  call get_precip_nc4(fn, Qp, E0p)
case ('.dat')
  call get_precip_raw(fn, Qp, E0p)
case default
  write(stderr,*) 'reader:get_precip: unknown file suffix on ' // fn
  error stop 6
end select
end subroutine get_precip


subroutine get_neutral2(path, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnO,dnN2,dnO2,dvnrho,dvnz,dTn

character(:), allocatable :: fn

fn = get_filename(path)

select case (get_suffix(fn))
case ('.h5')
  call get_neutral2_hdf5(fn, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
case ('.nc')
  call get_neutral2_nc4(fn, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
case ('.dat')
  call get_neutral2_raw(fn, dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
case default
  write(stderr,*) 'reader:get_neutral2: unknown file suffix' // get_suffix(fn)
  error stop 6
end select
end subroutine get_neutral2


subroutine get_neutral3(path, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
character(*), intent(in) :: path
real(wp), dimension(:,:,:), intent(out) :: dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall

character(:), allocatable :: fn

fn = get_filename(path)

select case (get_suffix(fn))
case ('.h5')
  call get_neutral3_hdf5(fn, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
case ('.nc')
  call get_neutral3_nc4(fn, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
case ('.dat')
  call get_neutral3_raw(fn, dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
case default
  write(stderr,*) 'reader:get_neutral3: unknown file suffix' // get_suffix(fn)
  error stop 6
end select
end subroutine get_neutral3


end module reader