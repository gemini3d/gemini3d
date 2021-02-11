submodule (io:plasma_output) plasma_output_hdf5

use timeutils, only : date_filename
use h5fortran, only: hdf5_file

implicit none (type, external)

contains

module procedure output_root_stream_mpi_hdf5

!! COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE USING STREAM I/O.
!! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS

character(:), allocatable :: filenamefull

integer :: lx1,lx2all,lx3all,isp
type(hdf5_file) :: hout

! single precision work arrays
real, dimension(:,:,:), allocatable :: permarray
real, dimension(:,:,:,:), allocatable :: permarray4D
real, dimension(:,:), allocatable :: permarray2D

!! SYSTEM SIZES
lx1=size(Phiall,1)
lx2all=size(Phiall,2)
lx3all=size(Phiall,3)


!> FIGURE OUT THE FILENAME
filenamefull = date_filename(outdir,ymd,UTsec) // '.h5'
print *, 'HDF5 Output file name:  ', filenamefull

call hout%initialize(filenamefull, status='new',action='w',comp_lvl=comp_lvl)

call hout%write('/time/ymd', ymd)
call hout%write('/time/UThour',   real(UTsec/3600.))

if (flagswap/=1) then
  select case (flagoutput)
    case (2)    !output ISR-like average parameters
      call hout%write('neall',    real(neall(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('v1avgall', real(v1avgall(1:lx1,1:lx2all,1:lx3all)))
      !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
      call hout%write('Tavgall',  real(Tavgall(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('TEall',    real(Teall(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('J1all',    real(J1all(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('J2all',    real(J2all(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('J3all',    real(J3all(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('v2avgall', real(v2avgall(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('v3avgall', real(v3avgall(1:lx1,1:lx2all,1:lx3all)))
    case (3)     !just electron density
      print *, 'INFO:  Input file has selected electron density only output, make sure this is what you really want!'
      call hout%write('neall',    real(neall(1:lx1,1:lx2all,1:lx3all)))
    case default    !output everything
      print *, 'INFO:  Input file has selected full output or milestones, large files may result!'
      call hout%write('nsall',    real(nsall(1:lx1,1:lx2all,1:lx3all,:)))
      call hout%write('vs1all',   real(vs1all(1:lx1,1:lx2all,1:lx3all,:)))
      !this is full output of all parameters in 3D
      call hout%write('Tsall',    real(Tsall(1:lx1,1:lx2all,1:lx3all,:)))

      call hout%write('J1all',    real(J1all(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('J2all',    real(J2all(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('J3all',    real(J3all(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('v2avgall', real(v2avgall(1:lx1,1:lx2all,1:lx3all)))
      call hout%write('v3avgall', real(v3avgall(1:lx1,1:lx2all,1:lx3all)))
    end select
else
!! 2D simulation that has been swapped around
  allocate(permarray(lx1,lx3all,lx2all))    !temporary work array that has been permuted  
  select case (flagoutput)
    case (2)    !averaged parameters
      permarray=reshape(real(neall),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('neall',    permarray)
      permarray=reshape(real(v1avgall),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('v1avgall', permarray)
      permarray=reshape(real(Tavgall),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('Tavgall',  permarray)
      permarray=reshape(real(Teall),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('TEall',    permarray)
      permarray=reshape(real(J1all),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('J1all',    permarray)

      ! J3,J2 and V3, V2 are swapped
      permarray=reshape(real(J3all),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('J2all',    permarray)
      permarray=reshape(real(J2all),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('J3all',    permarray)
      permarray=reshape(real(v3avgall),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('v2avgall', permarray)
      permarray=reshape(real(v2avgall),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('v3avgall', permarray)
    case (3)     !electron density only output
      print *, 'INFO:  Input file has selected electron density only output, make sure this is what you really want!'
      permarray=reshape(real(neall),[lx1,lx3all,lx2all],order=[1,3,2])       
      call hout%write('neall',    permarray)

    case default
      print *, 'INFO:  Input file has selected full output or milestones, large files may result!'

      permarray4D=reshape(real(nsall(1:lx1,1:lx2all,1:lx3all,:)),[lx1,lx3all,lx2all,lsp],order=[1,3,2,4])
      call hout%write('nsall',    permarray4D)
      permarray4D=reshape(real(vs1all(1:lx1,1:lx2all,1:lx3all,:)),[lx1,lx3all,lx2all,lsp],order=[1,3,2,4])
      call hout%write('vs1all',   permarray4D)
      permarray4D=reshape(real(Tsall(1:lx1,1:lx2all,1:lx3all,:)),[lx1,lx3all,lx2all,lsp],order=[1,3,2,4])
      call hout%write('Tsall',    permarray4D)

      permarray=reshape(real(J1all),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('J1all',    permarray)

      !! NOTE: J3,J2 and V3, V2 are swapped in name like this
      permarray=reshape(real(J3all),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('J2all',    permarray)
      permarray=reshape(real(J2all),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('J3all',    permarray)
      permarray=reshape(real(v3avgall),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('v2avgall', permarray)
      permarray=reshape(real(v2avgall),[lx1,lx3all,lx2all],order=[1,3,2])
      call hout%write('v3avgall', permarray)
  end select
end if
if (gridflag==1) then
  print *, 'Writing topside boundary conditions for inverted-type grid...'
  if (flagswap/=1) then
    call hout%write('Phiall',       real(Phiall(1,:,:)))
  else
    permarray2D=reshape(real(Phiall(1,:,:)),[lx3all,lx2all],order=[2,1])
    call hout%write('Phiall',permarray2D)
  end if
else
  print *, 'Writing topside boundary conditions for non-inverted-type grid...'
  if (flagswap/=1) then
    call hout%write('Phiall',       real(Phiall(lx1,:,:)))
  else
    permarray2D=reshape(real(Phiall(lx1,:,:)),[lx3all,lx2all],order=[2,1])
    call hout%write('Phiall',permarray2D)
  end if
end if

call hout%finalize()

end procedure output_root_stream_mpi_hdf5


end submodule plasma_output_hdf5
