!! Tests the file dump on error capability.
program test_error_write

use errors, only : error_stop
use phys_consts, only : wp

implicit none (type, external)

character(1000) :: buf
character(:), allocatable :: mode, filename

integer :: flagoutput

integer, dimension(3) :: ymd
real(wp) :: UTsec
real(wp), dimension(-1:5,-1:4,-1:3,7) :: nsall,vs1all,Tsall
real(wp), dimension(5,4,3) :: Phiall, v2avgall,v3avgall,J1all,J2all,J3all,neall,v1avgall,Tavgall,Teall

real(wp), dimension(-1:5,-1:4,-1:3, 7) :: vs2,vs3,ns,vs1,Ts
real(wp), dimension(5,4,3) :: J1,J2,J3

if(command_argument_count() /= 2) error stop "need test_type and dump_filename"

call get_command_argument(1, buf)
mode = trim(buf)

call get_command_argument(2, buf)
filename = trim(buf)

select case (mode)
case ("root")
  call error_stop(filename,"testing:root",flagoutput,ymd,UTsec,v2avgall,v3avgall,nsall,vs1all,Tsall, &
    Phiall,J1all,J2all,J3all,neall,v1avgall,Tavgall,Teall)
case ("worker")
  call error_stop(filename, "testing:worker", vs2,vs3,ns,vs1,Ts,J1,J2,J3)
case ("input")
  call error_stop(filename, "testing:input", ns, vs1, Ts)
end select

error stop "Unknown test dump case: " // mode

end program
