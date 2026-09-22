! SPDX-License-Identifier: Apache-2.0
module restart_runtime
! Opt-in full-step state for fixed-grid, electrostatic FANG/frozen-neutral runs.
! Rank-local records are published before the final root checkpoint name.
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp
use gemini3d_config, only: gemini_cfg
use h5fortran, only: hdf5_file
use hdf5, only: H5T_NATIVE_DOUBLE, h5tequal_f
use atomic_file, only: stage_file, publish_file
use timeutils, only: date_filename
use mpimod, only: mpi_cfg
use mpi_f08, only: MPI_COMM_WORLD,MPI_CHARACTER,MPI_Bcast,MPI_Barrier
implicit none (type,external)
private
public :: runtime_select,runtime_read,runtime_save,runtime_restore,runtime_output_header, &
          runtime_output_path,runtime_times,runtime_dt,runtime_ready,runtime_exclude_model
logical :: enabled=.false., pending=.false., runtime_ready=.false.
integer, parameter :: runtime_schema=2
real(wp) :: previous_dt=0, next_output=0, next_glow=0, saved_neutral=0
integer :: saved_iteration=1
character(:), allocatable :: input_digest, executable_digest, restore_file, staged_output, state_prefix
character(:), allocatable :: checkpoint_name, restore_checkpoint, restore_generation
contains
subroutine runtime_exclude_model(model)
  character(*), intent(in) :: model
  character(32) :: value
  integer :: status
  call get_environment_variable('GEMINI_EXACT_RESTART',value,status=status)
  if(status==0.and.trim(value)=='1') &
    error stop 'Exact restart research profile excludes model: '//model
end subroutine
subroutine check_mode(cfg)
  type(gemini_cfg), intent(in) :: cfg
  if(cfg%potsolve/=1.or.cfg%flagcap/=0.or.cfg%flagglow/=0.or.cfg%flagneuBG.or. &
     cfg%flagdneu/=0.or.cfg%flagneutralBGfile/=0.or.cfg%flaglagrangian.or.cfg%flagoutput/=1.or. &
     cfg%flagsolfluxfile/=0.or.cfg%allow_missing_spatial) &
    error stop 'Exact restart profile requires full-output fixed-grid electrostatic FANG and frozen empirical neutrals'
  if(storage_size(1._wp)/=64) error stop 'Exact restart profile requires real64'
end subroutine
subroutine runtime_select(cfg)
  type(gemini_cfg), intent(in) :: cfg
  character(256) :: value
  integer :: status,length,i
  enabled=.false.;pending=.false.;runtime_ready=.false.;previous_dt=0
  if(allocated(staged_output)) deallocate(staged_output)
  if(allocated(state_prefix)) deallocate(state_prefix)
  if(allocated(checkpoint_name)) deallocate(checkpoint_name)
  call get_environment_variable('GEMINI_EXACT_RESTART',value,length,status)
  if(status==0.and.trim(value)=='1') enabled=.true.
  if(.not.enabled) return
  call check_mode(cfg)
  call get_environment_variable('GEMINI_INPUT_SHA256',value,length,status)
  if(status/=0.or.length/=64) error stop 'Exact restart requires the verified-input launcher and GEMINI_INPUT_SHA256'
  do i=1,64
    if(index('0123456789abcdef',value(i:i))==0) error stop 'Invalid input SHA256 token'
  enddo
  input_digest=value(:64)
  call get_environment_variable('GEMINI_EXECUTABLE_SHA256',value,length,status)
  if(status/=0.or.length/=64) error stop 'Exact restart requires an executable SHA256 token'
  do i=1,64
    if(index('0123456789abcdef',value(i:i))==0) error stop 'Invalid executable SHA256 token'
  enddo
  executable_digest=value(:64)
end subroutine
subroutine runtime_read(core,cfg,t)
  real(wp), intent(in) :: t
  logical :: exists
  character(*), intent(in) :: core
  type(gemini_cfg), intent(in) :: cfg
  type(hdf5_file) :: f
  character(4096) :: prefix,checkpoint
  character(64) :: digest
  character(16) :: rank
  integer :: schema,layout(2),iteration
  call f%open(core,action='r')
  if(.not.f%exist('/restart_runtime')) then
    call f%close()
    if(enabled) error stop 'Exact restart requested but full-step runtime state is missing'
    return
  endif
  if(.not.enabled) error stop 'This checkpoint requires the verified exact-restart profile launcher'
  call check_mode(cfg)
  call f%read('/restart_runtime/schema',schema)
  if(schema/=runtime_schema) error stop 'Unsupported runtime checkpoint schema: exact restart requires schema 2'
  call f%read('/restart_runtime/layout',layout)
  call f%read('/restart_runtime/prefix',prefix)
  call f%read('/restart_runtime/checkpoint',checkpoint)
  call f%read('/restart_runtime/executable_sha256',digest)
  if(digest/=executable_digest) error stop 'Changed restart executable: SHA256 mismatch'
  call f%read('/restart_runtime/input_sha256',digest)
  call f%close()
  if(any(layout/=[mpi_cfg%lid2,mpi_cfg%lid3])) error stop 'Incompatible runtime checkpoint layout'
  if(digest/=input_digest) error stop 'Changed restart inputs: input SHA256 mismatch'
  if(trim(checkpoint)/=checkpoint_basename(core)) error stop 'Runtime checkpoint root identity mismatch'
  if(len_trim(checkpoint)==0.or.len_trim(prefix)<=len_trim(checkpoint)+9.or. &
     index(prefix,'/')/=0.or.index(prefix,'\')/=0.or.index(prefix,'..')/=0) &
    error stop 'Invalid runtime checkpoint prefix'
  if(index(trim(prefix),trim(checkpoint)//'.partial.')/=1) error stop 'Runtime checkpoint generation mismatch'
  restore_checkpoint=trim(checkpoint)
  restore_generation=trim(prefix)
  write(rank,'(I8.8)') mpi_cfg%myid
  restore_file=cfg%outdir//'/'//trim(prefix)//'.r'//trim(rank)//'.h5'
  inquire(file=restore_file,exist=exists)
  if(.not.exists) error stop 'Missing runtime state file: '//restore_file
  call f%open(restore_file,action='r')
  call check_identity(f)
  call f%read('/iteration',iteration)
  call f%read('/dt',previous_dt)
  call f%read('/next_output',next_output)
  call f%read('/next_glow',next_glow)
  call f%read('/neutral_time',saved_neutral)
  if(iteration<1) error stop 'Invalid runtime state metadata'
  if(.not.all(ieee_is_finite([previous_dt,next_output,next_glow,saved_neutral])).or.previous_dt<=0) &
    error stop 'Invalid runtime checkpoint clock'
  if(abs(next_output-(t+cfg%dtout))>1e-5_wp.or.previous_dt>cfg%dtout) &
    error stop 'Inconsistent runtime checkpoint clock'
  saved_iteration=iteration
  call f%close()
  pending=.true.
end subroutine
subroutine check_identity(f)
  type(hdf5_file), intent(inout) :: f
  integer :: schema,rank,layout(2),complete
  character(64) :: digest
  character(4096) :: checkpoint,generation
  call f%read('/schema',schema)
  if(schema/=runtime_schema) error stop 'Unsupported runtime state schema: exact restart requires schema 2'
  call f%read('/complete',complete)
  if(complete/=1) error stop 'Incomplete runtime checkpoint'
  call f%read('/rank',rank)
  if(rank/=mpi_cfg%myid) error stop 'Runtime state rank mismatch'
  call f%read('/layout',layout)
  if(any(layout/=[mpi_cfg%lid2,mpi_cfg%lid3])) error stop 'Runtime state layout mismatch'
  call f%read('/checkpoint',checkpoint)
  if(trim(checkpoint)/=restore_checkpoint) error stop 'Runtime state root identity mismatch'
  call f%read('/generation',generation)
  if(trim(generation)/=restore_generation) error stop 'Runtime state generation mismatch'
  call f%read('/input_sha256',digest)
  if(digest/=input_digest) error stop 'Changed restart inputs: runtime state input SHA256 mismatch'
  call f%read('/executable_sha256',digest)
  if(digest/=executable_digest) error stop 'Changed restart executable: runtime state SHA256 mismatch'
end subroutine
subroutine check_payload(fluid,electro,vi1,vi2,vi3)
  real(wp), intent(in) :: fluid(:,:,:,:),electro(:,:,:,:),vi1(:,:,:,:),vi2(:,:,:,:),vi3(:,:,:,:)
  if(.not.all(ieee_is_finite(fluid))) error stop 'Nonfinite runtime checkpoint payload: fluid'
  if(.not.all(ieee_is_finite(electro))) error stop 'Nonfinite runtime checkpoint payload: electro'
  if(.not.all(ieee_is_finite(vi1))) error stop 'Nonfinite runtime checkpoint payload: vi1'
  if(.not.all(ieee_is_finite(vi2))) error stop 'Nonfinite runtime checkpoint payload: vi2'
  if(.not.all(ieee_is_finite(vi3))) error stop 'Nonfinite runtime checkpoint payload: vi3'
end subroutine
subroutine runtime_times(tout,tglowout)
  real(wp), intent(inout) :: tout,tglowout
  if(.not.pending) return
  tout=next_output;tglowout=next_glow
end subroutine
subroutine runtime_dt(dt,read_saved)
  real(wp), intent(inout) :: dt
  logical, intent(in) :: read_saved
  if(read_saved.and.runtime_ready) then
    dt=previous_dt
    runtime_ready=.false.
  elseif(.not.read_saved) then
    previous_dt=dt
  endif
end subroutine
subroutine runtime_restore(fluid,electro,vi1,vi2,vi3,iteration,neutral_time)
  real(wp), intent(inout) :: fluid(:,:,:,:),electro(:,:,:,:),vi1(:,:,:,:),vi2(:,:,:,:),vi3(:,:,:,:)
  integer, intent(inout) :: iteration
  real(wp), intent(inout) :: neutral_time
  type(hdf5_file) :: f
  integer :: i,type_error
  logical :: is_double
  character(7), parameter :: fields(5)=[character(7)::'fluid','electro','vi1','vi2','vi3']
  if(.not.pending) return
  call f%open(restore_file,action='r')
  call check_identity(f)
  do i=1,size(fields)
    call h5tequal_f(f%dtype('/'//trim(fields(i))),H5T_NATIVE_DOUBLE,is_double,type_error)
    if(type_error/=0) error stop 'Cannot compare runtime checkpoint datatype'
    if(.not.is_double.or.f%ndim('/'//trim(fields(i)))/=4) &
      error stop 'Runtime checkpoint requires float64 rank-four state'
  enddo
  call f%read('/fluid',fluid);call f%read('/electro',electro)
  call f%read('/vi1',vi1);call f%read('/vi2',vi2);call f%read('/vi3',vi3)
  call f%close()
  call check_payload(fluid,electro,vi1,vi2,vi3)
  iteration=saved_iteration;neutral_time=saved_neutral
  pending=.false.;runtime_ready=.true.
end subroutine
subroutine runtime_save(cfg,ymd,ut,fluid,electro,vi1,vi2,vi3,iteration,neutral_time,tout,tglowout)
  type(gemini_cfg), intent(in) :: cfg
  integer, intent(in) :: ymd(3),iteration
  real(wp), intent(in) :: ut,neutral_time,tout,tglowout
  real(wp), intent(in) :: fluid(:,:,:,:),electro(:,:,:,:),vi1(:,:,:,:),vi2(:,:,:,:),vi3(:,:,:,:)
  character(:), allocatable :: path,temporary,final
  character(4096) :: message
  character(16) :: rank
  type(hdf5_file) :: f
  if(.not.enabled) return
  call check_payload(fluid,electro,vi1,vi2,vi3)
  final=date_filename(cfg%outdir,ymd,ut)//'.h5'
  checkpoint_name=final(scan(final,'/\',back=.true.)+1:)
  message=''
  if(mpi_cfg%myid==0) then
    staged_output=stage_file(final)
    if(len(staged_output)>len(message)) error stop 'Checkpoint path too long'
    message=staged_output
  endif
  call MPI_Bcast(message,len(message),MPI_CHARACTER,0,MPI_COMM_WORLD)
  staged_output=trim(message)
  state_prefix=staged_output(len_trim(cfg%outdir)+2:)
  write(rank,'(I8.8)') mpi_cfg%myid
  path=staged_output//'.r'//trim(rank)//'.h5'
  temporary=stage_file(path)
  call f%open(temporary,action='w')
  call f%write('/schema',runtime_schema);call f%write('/input_sha256',input_digest)
  call f%write('/executable_sha256',executable_digest)
  call f%write('/rank',mpi_cfg%myid)
  call f%write('/layout',[mpi_cfg%lid2,mpi_cfg%lid3])
  call f%write('/checkpoint',checkpoint_name)
  call f%write('/generation',state_prefix)
  call f%write('/iteration',iteration);call f%write('/dt',previous_dt)
  call f%write('/next_output',tout);call f%write('/next_glow',tglowout)
  call f%write('/neutral_time',neutral_time)
  call f%write('/fluid',fluid);call f%write('/electro',electro)
  call f%write('/vi1',vi1);call f%write('/vi2',vi2);call f%write('/vi3',vi3)
  call f%write('/complete',1);call f%close()
  call publish_file(temporary,path)
  call MPI_Barrier(MPI_COMM_WORLD)
end subroutine
function runtime_output_path(final) result(path)
  character(*), intent(in) :: final
  character(:), allocatable :: path
  if(enabled.and.allocated(staged_output)) then
    path=staged_output
  else
    path=stage_file(final)
  endif
end function
subroutine runtime_output_header(f)
  type(hdf5_file), intent(inout) :: f
  if(.not.enabled.or..not.allocated(state_prefix)) return
  call f%write('/restart_runtime/schema',runtime_schema)
  call f%write('/restart_runtime/layout',[mpi_cfg%lid2,mpi_cfg%lid3])
  call f%write('/restart_runtime/input_sha256',input_digest)
  call f%write('/restart_runtime/executable_sha256',executable_digest)
  call f%write('/restart_runtime/prefix',state_prefix)
  call f%write('/restart_runtime/checkpoint',checkpoint_name)
end subroutine
function checkpoint_basename(path) result(name)
  character(*), intent(in) :: path
  character(:), allocatable :: name
  integer :: n,last
  n=len_trim(path)
  if(n<1) error stop 'Runtime checkpoint root identity mismatch'
  if(path(n:n)=='/'.or.path(n:n)=='\') error stop 'Runtime checkpoint root identity mismatch'
  last=scan(path(:n),'/\',back=.true.)
  name=path(last+1:n)
  if(len_trim(name)==0) error stop 'Runtime checkpoint root identity mismatch'
end function
end module
