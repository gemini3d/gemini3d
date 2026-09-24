program audit_halo_end
use mpi_f08, only: MPI_Init,MPI_Finalize,MPI_Barrier,MPI_COMM_WORLD,MPI_Wtime
use mpimod, only: mpisetup,mpi_manualgrid,mpi_cfg,halo_end
use phys_consts, only: wp
implicit none
integer :: p2,p3,iteration
call MPI_Init()
call mpisetup()
do p2=1,mpi_cfg%lid
  if (mod(mpi_cfg%lid,p2)/=0) cycle
  p3=mpi_cfg%lid/p2
  call mpi_manualgrid(4*p2,4*p3,p2,p3)
  do iteration=1,3
    call check_halo(1,iteration)
    call check_halo(32768,iteration)
  enddo
enddo
call MPI_Finalize()
print *, 'Halo end, top and corner exchanges passed'
contains
subroutine check_halo(n1,iteration)
  integer, intent(in) :: n1,iteration
  integer, parameter :: n2=3,n3=4
  real(wp) :: param(n1,n2,n3),pend(n1,n2),ptop(n1,n3),corner(n1),expected(n1),start
  integer :: i,j,k,g2,g3
  do k=1,n3
    do j=1,n2
      g2=mpi_cfg%myid2*n2+j
      g3=mpi_cfg%myid3*n3+k
      param(:,j,k)=[(value(i,g2,g3,iteration),i=1,n1)]
    enddo
  enddo
  pend=-1; ptop=-1; corner=-1
  call MPI_Barrier(MPI_COMM_WORLD)
  if (mpi_cfg%myid==mpi_cfg%lid-1) then
    start=MPI_Wtime()
    do while (MPI_Wtime()-start<0.01_wp)
    enddo
  endif
  ! The long vector exceeds typical eager limits, exercising receive completion.
  call halo_end(param,pend,ptop,corner,450)
  do j=1,n2
    expected=0
    if (mpi_cfg%myid3<p3-1) expected=[(value(i,mpi_cfg%myid2*n2+j,(mpi_cfg%myid3+1)*n3+1,iteration),i=1,n1)]
    if (any(pend(:,j)/=expected)) error stop 'Incorrect halo end'
  enddo
  do k=1,n3
    expected=0
    if (mpi_cfg%myid2<p2-1) expected=[(value(i,(mpi_cfg%myid2+1)*n2+1,mpi_cfg%myid3*n3+k,iteration),i=1,n1)]
    if (any(ptop(:,k)/=expected)) error stop 'Incorrect halo top'
  enddo
  expected=0
  if (mpi_cfg%myid2<p2-1 .and. mpi_cfg%myid3<p3-1) &
    expected=[(value(i,(mpi_cfg%myid2+1)*n2+1,(mpi_cfg%myid3+1)*n3+1,iteration),i=1,n1)]
  if (any(corner/=expected)) error stop 'Incorrect halo corner'
end subroutine
pure real(wp) function value(i,j,k,iteration)
  integer, intent(in) :: i,j,k,iteration
  value=real(i+100000*j+10000000*k+100000000*iteration,wp)
end function
end program
