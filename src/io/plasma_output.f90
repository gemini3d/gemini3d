submodule (io) plasma_output

implicit none (type, external)

interface ! plasma_output_*.f90

module subroutine output_root_stream_mpi_hdf5(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
character(*), intent(in) :: outdir
integer, intent(in) :: flagoutput

integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts

real(wp), dimension(:,:,:), intent(in) :: Phiall
real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3
end subroutine output_root_stream_mpi_hdf5


module subroutine output_root_stream_mpi_nc4(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
character(*), intent(in) :: outdir
integer, intent(in) :: flagoutput

integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts

real(wp), dimension(:,:,:), intent(in) :: Phiall
real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3
end subroutine output_root_stream_mpi_nc4


module subroutine output_root_stream_mpi_raw(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
character(*), intent(in) :: outdir
integer, intent(in) :: flagoutput

integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts

real(wp), dimension(:,:,:), intent(in) :: Phiall
real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3
end subroutine output_root_stream_mpi_raw

end interface

contains

subroutine output_workers_mpi(vs2,vs3,ns,vs1,Ts,J1,J2,J3)

!------------------------------------------------------------
!-------SEND COMPLETE DATA FROM WORKERS TO ROOT PROCESS FOR OUTPUT.
!-------STATE VARS ARE EXPECTED TO INCLUDE GHOST CELLS
!------------------------------------------------------------

real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3

integer :: lx1,lx2,lx3,lx3all,isp
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: v2avg,v3avg


!SYSTEM SIZES (W/O GHOST CELLS)
lx1=size(ns,1)-4
lx2=size(ns,2)-4
lx3=size(ns,3)-4


!ONLY AVERAGE DRIFTS PERP TO B NEEDED FOR OUTPUT
v2avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs2(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
v2avg=v2avg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
v3avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs3(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
v3avg=v3avg/ns(1:lx1,1:lx2,1:lx3,lsp)


!SEND MY GRID DATA TO THE ROOT PROCESS
call gather_send(v2avg,tag%v2)
call gather_send(v3avg,tag%v3)
call gather_send(ns,tag%ns)
call gather_send(vs1,tag%vs1)
call gather_send(Ts,tag%Ts)


!------- SEND ELECTRODYNAMIC PARAMETERS TO ROOT
call gather_send(J1,tag%J1)
call gather_send(J2,tag%J2)
call gather_send(J3,tag%J3)

end subroutine output_workers_mpi


module procedure output_plasma
! subroutine output_plasma(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
!! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

if (myid == 0) then
  call output_root_stream_mpi(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, out_format)
else
  call output_workers_mpi(vs2,vs3,ns,vs1,Ts,J1,J2,J3)
end if

end procedure output_plasma


subroutine output_root_stream_mpi(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, out_format)
character(*), intent(in) :: outdir
character(*), intent(in) :: out_format
integer, intent(in) :: flagoutput

integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts

real(wp), dimension(:,:,:), intent(in) :: Phiall
real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3

select case (out_format)
case ('h5', 'hdf5')
  call output_root_stream_mpi_hdf5(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
case ('nc', 'nc4')
  call output_root_stream_mpi_nc4(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
case ('raw')
  call output_root_stream_mpi_raw(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
case default
  write(stderr,*) 'plasma_output:output_root_stream_api: unknown format' // out_format
  error stop 6
end select


end subroutine output_root_stream_mpi

end submodule plasma_output