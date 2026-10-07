program audit_mpi_time
use mpi_f08, only: MPI_Init,MPI_Finalize
use mpimod, only: mpisetup,mpi_cfg
use phys_consts, only: wp,lsp
use temporal_mpi, only: dt_comm,enforce_post_update_cfl
use meshobj_cart, only: cartmesh
use gemini3d_config, only: gemini_cfg
implicit none
type(cartmesh) :: x
type(gemini_cfg) :: cfg
real(wp) :: ns(-1:4,-1:3,-1:3,lsp),ts(-1:4,-1:3,-1:3,lsp),v(-1:4,-1:3,-1:3,lsp)
real(wp) :: v2(-1:4,-1:3,-1:3,lsp),b(-1:4,-1:3,-1:3),dt
character(32) :: mode
call get_command_argument(1,mode)
call MPI_Init()
call mpisetup()
allocate(x%dl1i(2,1,1),x%dl2i(2,1,1),x%dl3i(2,1,1))
x%dl1i=1e6_wp;x%dl2i=1;x%dl3i=1
ns=1e11_wp;ts=300;v=0;v2=0;b=0
cfg%tcfl=0.5_wp;cfg%tdur=100;cfg%potsolve=1
if(trim(mode)=='worker_cfl')then
  ! Only a non-root worker violates CFL: the global reduction must detect it.
  if(mpi_cfg%lid<2)error stop 'worker CFL test needs at least two ranks'
  if(mpi_cfg%myid==1)v2(1,1,1,7)=20
  call enforce_post_update_cfl(ts,v,v2,v,x,0.1_wp)
  call MPI_Finalize()
  stop 0  ! WILL_FAIL must reject an undetected violation.
endif
call dt_comm(0._wp,0._wp,0._wp,cfg,ns,ts,v,v2,v,b,b,b,x,dt)
if(abs(dt-1e-6_wp)>1e-15_wp)error stop 'bootstrap cap mismatch'
! A worker-limited step below a microsecond must not be raised by bootstrap.
if(mpi_cfg%myid==1)v2=1e8_wp
call dt_comm(0._wp,0._wp,0._wp,cfg,ns,ts,v,v2,v,b,b,b,x,dt)
if(mpi_cfg%lid>1 .and. dt>5e-9_wp*(1+1e-12_wp))error stop 'unsafe minimum step'
v2=0
call dt_comm(99.9_wp,200._wp,200._wp,cfg,ns,ts,v,v2,v,b,b,b,x,dt)
if(abs(dt-0.1_wp)>1e-12_wp)error stop 'final time overshoot'
cfg%flagE0file=1;cfg%dtE0=10
call dt_comm(9.5_wp,100._wp,100._wp,cfg,ns,ts,v,v2,v,b,b,b,x,dt)
if(abs(dt-0.5_wp)>1e-12_wp)error stop 'driver knot overshoot'
call enforce_post_update_cfl(ts,v,v2,v,x,dt)
deallocate(x%dl1i,x%dl2i,x%dl3i)
call MPI_Finalize()
print *, 'MPI time selection and safe CFL passed'
end program
