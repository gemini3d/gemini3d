program audit_time
use phys_consts, only: wp
use timeutils, only: find_lastdate,find_time_elapsed,shift_datetime
implicit none
character(40) :: mode,buf
integer :: a(3),b(3),out(3),i
real(wp) :: ua,ub,step,uout
call get_command_argument(1,mode)
do i=1,3
 call get_command_argument(i+1,buf);read(buf,*)a(i)
enddo
call get_command_argument(5,buf);read(buf,*)ua
if(trim(mode)=='shift')then
 call get_command_argument(6,buf);read(buf,*)step
 call shift_datetime(step,a,ua)
 print *,a,ua
else
 do i=1,3
  call get_command_argument(i+5,buf);read(buf,*)b(i)
 enddo
 call get_command_argument(9,buf);read(buf,*)ub
 call get_command_argument(10,buf);read(buf,*)step
 if(trim(mode)=='last')then
  call find_lastdate(a,ua,b,ub,step,out,uout)
  print *,out,uout
 else
  print *,find_time_elapsed(a,ua,b,ub,step)
 endif
endif
end program
