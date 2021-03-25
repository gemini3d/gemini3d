submodule (io:plasma_output) plasma_output_nc

use timeutils, only : date_filename
use nc4fortran, only: netcdf_file

implicit none (type, external)

contains

module procedure output_root_stream_mpi_nc4

!! COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE USING STREAM I/O.
!! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS

integer :: lx1,lx2all,lx3all,isp
character(:), allocatable :: filenamefull

character(*), parameter :: dims4(4) = [character(7) :: 'x1', 'x2', 'x3', 'species'], &
  dims3(3) = [character(2) :: 'x1', 'x2', 'x3'], &
  dims23(2) = [character(2) :: 'x2', 'x3']

type(netcdf_file) :: hout


!! SYSTEM SIZES
lx1=size(Phiall,1)
lx2all=size(Phiall,2)
lx3all=size(Phiall,3)


!> FIGURE OUT THE FILENAME
filenamefull = date_filename(outdir,ymd,UTsec) // '.nc'
print *, 'Output file name:  ', filenamefull

call hout%initialize(filenamefull, status='new',action='w',comp_lvl=comp_lvl)

call hout%write("flagoutput", flagoutput)

call hout%write('ymd', ymd)
call hout%write('UThour',UTsec/3600._wp)

if (flagswap/=1) then
  select case (flagoutput)
    case (2)    !output ISR-like average parameters
      call hout%write('neall', neall(1:lx1,1:lx2all,1:lx3all), dims3)
      call hout%write('v1avgall', v1avgall(1:lx1,1:lx2all,1:lx3all), dims3)
      !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
      call hout%write('Tavgall', Tavgall(1:lx1,1:lx2all,1:lx3all), dims3)
      call hout%write('TEall', Teall(1:lx1,1:lx2all,1:lx3all), dims3)
      call hout%write('J1all', J1all(1:lx1,1:lx2all,1:lx3all), dims3)
      call hout%write('J2all', J2all(1:lx1,1:lx2all,1:lx3all), dims3)
      call hout%write('J3all', J3all(1:lx1,1:lx2all,1:lx3all), dims3)
      call hout%write('v2avgall', v2avgall(1:lx1,1:lx2all,1:lx3all), dims3)
      call hout%write('v3avgall', v3avgall(1:lx1,1:lx2all,1:lx3all), dims3)
    case (3)     !just electron density
      print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'
      call hout%write('neall', neall(1:lx1,1:lx2all,1:lx3all), dims3)
    case default    !output everything
      print *, '!!!NOTE:  Input file has selected full output, large files may result!'
      call hout%write('nsall', nsall(1:lx1,1:lx2all,1:lx3all, :), dims4)
      call hout%write('vs1all', vs1all(1:lx1,1:lx2all,1:lx3all, :), dims4)
      !this is full output of all parameters in 3D
      call hout%write('Tsall', Tsall(1:lx1,1:lx2all,1:lx3all, :), dims4)

      call hout%write('J1all', J1all(1:lx1, 1:lx2all, 1:lx3all), dims3)
      call hout%write('J2all', J2all(1:lx1, 1:lx2all, 1:lx3all), dims3)
      call hout%write('J3all', J3all(1:lx1, 1:lx2all, 1:lx3all), dims3)
      call hout%write('v2avgall', v2avgall(1:lx1, 1:lx2all, 1:lx3all), dims3)
      call hout%write('v3avgall', v3avgall(1:lx1, 1:lx2all, 1:lx3all), dims3)
    end select
else
!! 2D simulation
  select case (flagoutput)
    case (2)    !averaged parameters
      call hout%write('neall', neall, dims3)
      call hout%write('v1avgall', v1avgall, dims3)
      call hout%write('Tavgall', Tavgall, dims3)
      call hout%write('TEall', Teall, dims3)

      call hout%write('J1all', J1all, dims3)

      ! J3,J2 and V3, V2 are swapped
      call hout%write('J2all', J3all, dims3)
      call hout%write('J3all', J2all, dims3)
      call hout%write('v2avgall', v3avgall, dims3)
      call hout%write('v3avgall', v2avgall, dims3)
    case (3)     !electron density only output
      print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'

      call hout%write('neall', neall, dims3)

    case default
      print *, '!!!NOTE:  Input file has selected full output, large files may result!'

      call hout%write('nsall', nsall(1:lx1,1:lx2all,1:lx3all, :), dims4)
      call hout%write('vs1all', vs1all(1:lx1,1:lx2all,1:lx3all, :), dims4)
      call hout%write('Tsall', Tsall(1:lx1,1:lx2all,1:lx3all, :), dims4)

      call hout%write('J1all', J1all, dims3)

      !! NOTE: J3,J2 and V3, V2 are swapped in name like this
      call hout%write('J2all', J3all, dims3)
      call hout%write('J3all', J2all, dims3)
      call hout%write('v2avgall', v3avgall, dims3)
      call hout%write('v3avgall', v2avgall, dims3)
  end select
end if
if (gridflag==1) then
  print *, 'Writing topside boundary conditions for inverted-type grid...'
  call hout%write('Phiall', Phiall(1,:,:), dims23)
else
  print *, 'Writing topside boundary conditions for non-inverted-type grid...'
  call hout%write('Phiall', Phiall(lx1,:,:), dims23)
end if

call hout%finalize()

end procedure output_root_stream_mpi_nc4


end submodule plasma_output_nc
