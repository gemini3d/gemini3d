program test_errors

use errors, only : error_stop
use phys_consts, only : wp
implicit none (type, external)

integer :: i, argc
character(1000) :: buf
character(*), parameter :: file = 'errdump.h5'

integer :: flagoutput

integer, dimension(3) :: ymd
real(wp) :: UTsec
real(wp), dimension(-1:5,-1:4,-1:3,7) :: nsall,vs1all,Tsall
real(wp), dimension(5,4,3) :: Phiall, v2avgall,v3avgall,J1all,J2all,J3all,neall,v1avgall,Tavgall,Teall

real(wp), dimension(-1:5,-1:4,-1:3, 7) :: vs2,vs3,ns,vs1,Ts
real(wp), dimension(5,4,3) :: J1,J2,J3

call get_command_argument(1, buf)
read(buf,'(I1)') i

select case (buf)
case ("root")
  call error_stop(file,"testing:root",flagoutput,ymd,UTsec,v2avgall,v3avgall,nsall,vs1all,Tsall, &
    Phiall,J1all,J2all,J3all,neall,v1avgall,Tavgall,Teall)
case ("worker")
  call error_stop(file, "testing:worker", vs2,vs3,ns,vs1,Ts,J1,J2,J3)
case ("input")
  call error_stop(file, "testing:input", ns, vs1, Ts)
case default
  error stop "Error: unknown case: " // trim(buf)
end select


end program
