! SPDX-License-Identifier: Apache-2.0
program restart_probe
use mpi_f08, only: MPI_Init,MPI_Finalize,MPI_Barrier,MPI_COMM_WORLD
use mpimod, only: mpisetup,mpi_manualgrid,mpi_cfg
use phys_consts, only: wp
use gemini3d_config, only: gemini_cfg
use restart_runtime, only: runtime_select,runtime_read,runtime_save,runtime_restore, &
  runtime_output_header,runtime_output_path,runtime_times,runtime_dt,runtime_ready
use atomic_file, only: publish_file
use h5fortran, only: hdf5_file
use timeutils, only: date_filename
use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
implicit none
type(gemini_cfg) :: cfg
type(hdf5_file) :: f
real(wp) :: fluid(6,6,6,35),electro(6,6,6,7),vi1(3,2,2,7),vi2(2,3,2,7),vi3(2,2,3,7)
real(wp) :: dt,tout,tglow,neutral,base
integer :: iteration,layout(2)
character(4096) :: arg
character(32) :: mode
character(:), allocatable :: final,staged
call get_command_argument(1,mode)
call get_command_argument(2,arg)
cfg%outdir=trim(arg)
call get_command_argument(3,arg)
read(arg,*) layout(1)
call get_command_argument(4,arg)
read(arg,*) layout(2)
call MPI_Init()
call mpisetup()
call mpi_manualgrid(4,4,layout(1),layout(2))
cfg%dtout=60._wp
call runtime_select(cfg)
final=date_filename(cfg%outdir,[2013,2,20],60._wp)//'.h5'
base=real(mpi_cfg%myid+1,wp)
if(trim(mode)=='write'.or.trim(mode)=='write_nan') then
  fluid=base;electro=base+10;vi1=base+20;vi2=base+30;vi3=base+40
  if(trim(mode)=='write_nan'.and.mpi_cfg%myid==mpi_cfg%lid-1) &
    vi3(2,2,3,7)=ieee_value(0._wp,ieee_quiet_nan)
  dt=0.25_wp
  call runtime_dt(dt,.false.)
  call runtime_save(cfg,[2013,2,20],60._wp,fluid,electro,vi1,vi2,vi3,17,0._wp,120._wp,0._wp)
  if(mpi_cfg%myid==0) then
    staged=runtime_output_path(final)
    call f%open(staged,action='w')
    call runtime_output_header(f)
    call f%close()
    call publish_file(staged,final)
  endif
  call MPI_Barrier(MPI_COMM_WORLD)
elseif(trim(mode)=='read') then
  fluid=-1;electro=-1;vi1=-1;vi2=-1;vi3=-1
  call runtime_read(final,cfg,60._wp)
  if(runtime_ready) error stop 'Runtime became ready before payload validation'
  tout=-1;tglow=-1
  call runtime_times(tout,tglow)
  if(tout/=120._wp.or.tglow/=0._wp) error stop 'Restored deadlines mismatch'
  call runtime_restore(fluid,electro,vi1,vi2,vi3,iteration,neutral)
  if(.not.runtime_ready) error stop 'Runtime checkpoint not restored'
  if(any(fluid/=base).or.any(electro/=base+10).or.any(vi1/=base+20).or. &
     any(vi2/=base+30).or.any(vi3/=base+40)) error stop 'Restored payload mismatch'
  if(iteration/=17.or.neutral/=0._wp) error stop 'Restored runtime metadata mismatch'
  dt=-1
  call runtime_dt(dt,.true.)
  if(dt/=0.25_wp.or.runtime_ready) error stop 'Restored timestep mismatch'
else
  error stop 'Unknown restart probe mode'
endif
call MPI_Finalize()
end program
