! SPDX-License-Identifier: Apache-2.0
module atomic_file
! Same-directory publication on the supported POSIX research platform.
! Atomic visibility does not imply durability across power loss.
use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
implicit none (type, external)
private
public :: stage_file, publish_file
interface
  function c_mkstemp(template) bind(C,name='gemini_stage_create') result(fd)
    import c_int,c_char
    character(c_char), intent(inout) :: template(*)
    integer(c_int) :: fd
  end function
  function c_close(fd) bind(C,name='gemini_stage_close') result(status)
    import c_int
    integer(c_int), value :: fd
    integer(c_int) :: status
  end function
  function c_rename(old,new) bind(C,name='gemini_stage_publish') result(status)
    import c_int,c_char
    character(c_char), intent(in) :: old(*),new(*)
    integer(c_int) :: status
  end function
end interface
contains
function stage_file(final) result(temporary)
  character(*), intent(in) :: final
  character(:), allocatable :: temporary,template
  integer(c_int) :: fd,status
  template=final//'.partial.XXXXXX'//c_null_char
  fd=c_mkstemp(template)
  if(fd<0) error stop 'Cannot reserve checkpoint staging file: '//final
  status=c_close(fd)
  if(status/=0) error stop 'Cannot close checkpoint staging descriptor'
  temporary=template(:len(template)-1)
end function
subroutine publish_file(temporary,final)
  character(*), intent(in) :: temporary,final
  integer(c_int) :: status
  status=c_rename(temporary//c_null_char,final//c_null_char)
  if(status/=0) error stop 'Atomic checkpoint publication failed: '//final
end subroutine
end module
