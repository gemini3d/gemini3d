program atomic_probe
use atomic_file, only: stage_file,publish_file
implicit none(type,external)
character(4096) :: directory,mode
character(:), allocatable :: temporary,final
integer :: u
call get_command_argument(1,directory)
call get_command_argument(2,mode)
final=trim(directory)//'/checkpoint.h5'
if(trim(mode)=='missing_parent') final=trim(directory)//'/missing/checkpoint.h5'
temporary=stage_file(final)
open(newunit=u,file=temporary,status='old',action='write')
write(u,'(A)') 'complete new checkpoint'
close(u)
if(trim(mode)=='leave') stop
if(trim(mode)=='fail_publish') final=trim(directory)//'/directory'
call publish_file(temporary,final)
end program
