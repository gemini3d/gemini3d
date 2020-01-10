submodule (io:plasma) plasma_output_hdf5

use timeutils, only : date_filename
use h5fortran, only: hdf5_file

contains

module procedure output_root_stream_mpi

!! COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE USING STREAM I/O.
!! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS

integer :: lx1,lx2,lx3,lx2all,lx3all,isp, ierr
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: v2avg,v3avg
real(wp), dimension(-1:size(Phiall,1)+2,-1:size(Phiall,2)+2,-1:size(Phiall,3)+2,1:lsp) :: nsall,vs1all,Tsall
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: v2avgall,v3avgall,v1avgall,Tavgall,neall,Teall
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: J1all,J2all,J3all
character(:), allocatable :: filenamefull

type(hdf5_file) :: hout


!! SYSTEM SIZES
! FIXME: should these be pulled from the grid module???
lx1=size(ns,1)-4
lx2=size(ns,2)-4
lx3=size(ns,3)-4
lx2all=size(Phiall,2)
lx3all=size(Phiall,3)


print *, 'System sizes according to Phiall:  ',lx1,lx2all,lx3all
!ONLY AVERAGE DRIFTS PERP TO B NEEDED FOR OUTPUT
v2avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs2(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
v2avg=v2avg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
v3avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs3(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
v3avg=v3avg/ns(1:lx1,1:lx2,1:lx3,lsp)


!GET THE SUBGRID DATA FORM THE WORKERS
call gather_recv(v2avg,tagv2,v2avgall)
call gather_recv(v3avg,tagv3,v3avgall)
call gather_recv(ns,tagns,nsall)
call gather_recv(vs1,tagvs1,vs1all)
call gather_recv(Ts,tagTs,Tsall)


!> RADD--- NEED TO ALSO GATHER FULL GRID ELECTRODYANMICS PARAMTERS FROM WORKERS
call gather_recv(J1,tagJ1,J1all)
call gather_recv(J2,tagJ2,J2all)
call gather_recv(J3,tagJ3,J3all)


!COMPUTE AVERAGE VALUES FOR ION PLASMA PARAMETERS
v1avgall=sum(nsall(1:lx1,1:lx2all,1:lx3all,1:lsp-1)*vs1all(1:lx1,1:lx2all,1:lx3all,1:lsp-1),4)
v1avgall=v1avgall/nsall(1:lx1,1:lx2all,1:lx3all,lsp)    !compute averages for output.
Tavgall=sum(nsall(1:lx1,1:lx2all,1:lx3all,1:lsp-1)*Tsall(1:lx1,1:lx2all,1:lx3all,1:lsp-1),4)
Tavgall=Tavgall/nsall(1:lx1,1:lx2all,1:lx3all,lsp)    !compute averages for output.
neall=nsall(1:lx1,1:lx2all,1:lx3all,lsp)
Teall=Tsall(1:lx1,1:lx2all,1:lx3all,lsp)


!> FIGURE OUT THE FILENAME
filenamefull = date_filename(outdir,ymd,UTsec) // '.h5'
print *, 'HDF5 Output file name:  ', filenamefull

call hout%initialize(filenamefull, ierr, status='new',action='w',comp_lvl=1)

call hout%write('/time/ymd', ymd, ierr)
call hout%write('/time/UThour',UTsec/3600._wp, ierr)

if (flagswap/=1) then
  select case (flagoutput)
    case (2)    !output ISR-like average parameters
      call hout%write('neall', neall(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('v1avgall', v1avgall(1:lx1,1:lx2all,1:lx3all), ierr)
      !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
      call hout%write('Tavgall', Tavgall(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('TEall', Teall(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('J1all', J1all(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('J2all', J2all(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('J3all', J3all(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('v2avgall', v2avgall(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('v3avgall', v3avgall(1:lx1,1:lx2all,1:lx3all), ierr)
    case (3)     !just electron density
      print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'
      call hout%write('neall', neall(1:lx1,1:lx2all,1:lx3all), ierr)
    case default    !output everything
      print *, '!!!NOTE:  Input file has selected full output, large files may result!'
      call hout%write('nsall', nsall(1:lx1,1:lx2all,1:lx3all,:), ierr)
      call hout%write('vs1all', vs1all(1:lx1,1:lx2all,1:lx3all,:), ierr)
      !this is full output of all parameters in 3D
      call hout%write('Tsall', Tsall(1:lx1,1:lx2all,1:lx3all,:), ierr)

      call hout%write('J1all', J1all(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('J2all', J2all(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('J3all', J3all(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('v2avgall', v2avgall(1:lx1,1:lx2all,1:lx3all), ierr)
      call hout%write('v3avgall', v3avgall(1:lx1,1:lx2all,1:lx3all), ierr)
    end select
else
!! 2D simulation
  select case (flagoutput)
    case (2)    !averaged parameters
      call hout%write('neall', neall, ierr)
      call hout%write('v1avgall', v1avgall, ierr)
      call hout%write('Tavgall', Tavgall, ierr)
      call hout%write('TEall', Teall, ierr)

      call hout%write('J1all', J1all, ierr)

      ! J3,J2 and V3, V2 are swapped
      call hout%write('J2all', J3all, ierr)
      call hout%write('J3all', J2all, ierr)
      call hout%write('v2avgall', v3avgall, ierr)
      call hout%write('v3avgall', v2avgall, ierr)
    case (3)     !electron density only output
      print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'

      call hout%write('neall', neall, ierr)

    case default
      print *, '!!!NOTE:  Input file has selected full output, large files may result!'

      call hout%write('nsall', nsall(1:lx1,1:lx2all,1:lx3all,:), ierr)
      call hout%write('vs1all', vs1all(1:lx1,1:lx2all,1:lx3all,:), ierr)
      call hout%write('Tsall', Tsall(1:lx1,1:lx2all,1:lx3all,:), ierr)

      call hout%write('J1all', J1all, ierr)

      !! NOTE: J3,J2 and V3, V2 are swapped in name like this
      call hout%write('J2all', J3all, ierr)
      call hout%write('J3all', J2all, ierr)
      call hout%write('v2avgall', v3avgall, ierr)
      call hout%write('v3avgall', v2avgall, ierr)
  end select
end if
if (gridflag==1) then
  print *, 'Writing topside boundary conditions for inverted-type grid...'
  call hout%write('Phiall', Phiall(1,:,:), ierr)
else
  print *, 'Writing topside boundary conditions for non-inverted-type grid...'
  call hout%write('Phiall', Phiall(lx1,:,:), ierr)
end if

call hout%finalize(ierr)

!! Check for any NaN before proceeding to next time step

if (.not.all(ieee_is_finite(neall))) error stop 'neall has non-finite value(s)'
if (.not.all(ieee_is_finite(v1avgall))) error stop 'v1avgall has non-finite value(s)'
if (.not.all(ieee_is_finite(Teall))) error stop 'Teall has non-finite value(s)'
if (.not.all(ieee_is_finite(J1all))) error stop 'J1all has non-finite value(s)'
if (.not.all(ieee_is_finite(J2all))) error stop 'J2all has non-finite value(s)'
if (.not.all(ieee_is_finite(J3all))) error stop 'J3all has non-finite value(s)'
if (.not.all(ieee_is_finite(v2avgall))) error stop 'v2avgall has non-finite value(s)'
if (.not.all(ieee_is_finite(v3avgall))) error stop 'v3avgall has non-finite value(s)'
if (.not.all(ieee_is_finite(Phiall))) error stop 'Phiall has non-finite value(s)'

end procedure output_root_stream_mpi


end submodule plasma_output_hdf5
