submodule (io:plasma) plasma_output_hdf5

use hdf5_interface, only: hdf5_file

contains

module procedure output_root_stream_mpi

!! COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE USING STREAM I/O.
!! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS

integer :: lx1,lx2,lx3,lx3all,isp, u
real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: v2avg,v3avg
real(wp), dimension(-1:size(Phiall,1)+2,-1:size(Phiall,2)+2,-1:size(Phiall,3)+2,1:lsp) :: nsall,vs1all,Tsall
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: v2avgall,v3avgall,v1avgall,Tavgall,neall,Teall
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: J1all,J2all,J3all
character(:), allocatable :: filenamefull, h5filenamefull
integer(8) :: recordlength   !can be 8 byte with compiler flag -frecord-marker=8

type(hdf5_file) :: h5f


!! SYSTEM SIZES
! FIXME: should these be pull from the grid module???
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
filenamefull = date_filename(outdir,ymd,UTsec)
h5filenamefull = filenamefull(1:len(filenamefull)-4)//'.h5'
print *, 'HDF5 Output file name:  ',h5filenamefull


!SOME DEBUG OUTPUT ON FILE SIZE
recordlength=int(8,8)+int(8,8)*int(3,8)*int(lx1,8)*int(lx2all,8)*int(lx3all,8)*int(lsp,8)+ &
             int(8,8)*int(5,8)*int(lx1,8)*int(lx2all,8)*int(lx3all,8)+ &
             int(8,8)*int(lx2,8)*int(lx3all,8)
print *, 'Output bit length:  ',recordlength,lx1,lx2all,lx3all,lsp


call h5f%initialize(h5filenamefull,status='new',action='w',comp_lvl=1)

call h5f%add('/time/ymd', ymd)
call h5f%add('/time/UThour',UTsec/3600._wp)

if (flagswap/=1) then
  select case (flagoutput)
    case (2)    !output ISR-like average parameters
      call h5f%add('neall', neall(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('v1avgall', v1avgall(1:lx1,1:lx2all,1:lx3all))
      !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
      call h5f%add('Tavgall', Tavgall(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('TEall', Teall(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('J1all', J1all(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('J2all', J2all(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('J3all', J3all(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('v2avgall', v2avgall(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('v3avgall', v3avgall(1:lx1,1:lx2all,1:lx3all))
    case (3)     !just electron density
      print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'
      call h5f%add('neall', neall(1:lx1,1:lx2all,1:lx3all))
    case default    !output everything
      print *, '!!!NOTE:  Input file has selected full output, large files may result!'
      call h5f%add('nsall', nsall(1:lx1,1:lx2all,1:lx3all,:))
      call h5f%add('vs1all', vs1all(1:lx1,1:lx2all,1:lx3all,:))
      !this is full output of all parameters in 3D
      call h5f%add('Tsall', Tsall(1:lx1,1:lx2all,1:lx3all,:))

      call h5f%add('J1all', J1all(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('J2all', J2all(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('J3all', J3all(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('v2avgall', v2avgall(1:lx1,1:lx2all,1:lx3all))
      call h5f%add('v3avgall', v3avgall(1:lx1,1:lx2all,1:lx3all))
    end select
else
!! 2D simulation
  select case (flagoutput)
    case (2)    !averaged parameters
      call h5f%add('neall', neall)
      call h5f%add('v1avgall', v1avgall)
      call h5f%add('Tavgall', Tavgall)
      call h5f%add('TEall', Teall)

      call h5f%add('J1all', J1all)

      ! J3,J2 and V3, V2 are swapped
      call h5f%add('J2all', J3all)
      call h5f%add('J3all', J2all)
      call h5f%add('v2avgall', v3avgall)
      call h5f%add('v3avgall', v2avgall)
    case (3)     !electron density only output
      print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'

      call h5f%add('neall', neall)

    case default
      print *, '!!!NOTE:  Input file has selected full output, large files may result!'

      call h5f%add('nsall', nsall(1:lx1,1:lx2all,1:lx3all,:))
      call h5f%add('vs1all', vs1all(1:lx1,1:lx2all,1:lx3all,:))
      call h5f%add('Tsall', Tsall(1:lx1,1:lx2all,1:lx3all,:))

      call h5f%add('J1all', J1all)

      !! J3,J2 and V3, V2 are swapped in name like this
      call h5f%add('J2all', J3all)
      call h5f%add('J3all', J2all)
      call h5f%add('v2avgall', v3avgall)
      call h5f%add('v3avgall', v2avgall)
  end select
end if
if (gridflag==1) then
  print *, 'Writing topside boundary conditions for inverted-type grid...'
  call h5f%add('Phiall', Phiall(1,:,:))
else
  print *, 'Writing topside boundary conditions for non-inverted-type grid...'
  call h5f%add('Phiall', Phiall(lx1,:,:))
end if

call h5f%finalize()

end procedure output_root_stream_mpi


end submodule plasma_output_hdf5
