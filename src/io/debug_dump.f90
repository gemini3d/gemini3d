module debug_dump
!! dump to disk up to 9 real(wp) variables of rank 0..7
!! requires Fortran 2018 assumed rank, available in compilers including:
!! GCC >= 10, Intel Fortran >= 20.0, or Intel oneAPI

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use h5fortran, only : hdf5_file, h5write
use phys_consts, only : wp
implicit none (type, external)

private
public :: h5dump

contains


subroutine h5dump(filename, n1, v1, n2, v2, n3, v3, n4, v4, n5, v5, n6, v6, n7, v7, n8, v8, n9, v9)
!! dump up to nine: name, variable pairs.
character(*), intent(in) :: filename
character(*), intent(in), optional :: n1,n2,n3,n4,n5,n6,n7,n8,n9
real(wp), intent(in), optional :: v1(..), v2(..), v3(..), v4(..), v5(..), v6(..), v7(..), v8(..), v9(..)

if(present(n1) .and. present(v1)) call dumper(filename, n1, v1)
if(present(n2) .and. present(v2)) call dumper(filename, n2, v2)
if(present(n3) .and. present(v3)) call dumper(filename, n3, v3)
if(present(n4) .and. present(v4)) call dumper(filename, n4, v4)
if(present(n5) .and. present(v5)) call dumper(filename, n5, v5)
if(present(n6) .and. present(v6)) call dumper(filename, n6, v6)
if(present(n7) .and. present(v7)) call dumper(filename, n7, v7)
if(present(n8) .and. present(v8)) call dumper(filename, n8, v8)
if(present(n9) .and. present(v9)) call dumper(filename, n9, v9)

end subroutine h5dump


subroutine dumper(filename, name, value)
!! it's more concise to unpack rank with a separate procedure

character(*), intent(in) :: filename, name
real(wp), intent(in) :: value(..)

select rank (value)
rank (0)
  call h5write(filename, name, value)
rank (1)
  call h5write(filename, name, value)
rank (2)
  call h5write(filename, name, value)
rank (3)
  call h5write(filename, name, value)
rank (4)
  call h5write(filename, name, value)
rank (5)
  call h5write(filename, name, value)
rank (6)
  call h5write(filename, name, value)
rank (7)
  call h5write(filename, name, value)
rank default
  write(stderr, '(A,I0)') "Error: debug_dump:h5dump: h5fortran does not handle " // name // " of rank: ", rank(value)
end select

end subroutine dumper

end module debug_dump
