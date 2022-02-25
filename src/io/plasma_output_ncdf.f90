submodule (io:plasma_output) plasma_output_nc

use timeutils, only : date_filename
use nc4fortran, only: netcdf_file

implicit none (type, external)

contains

module procedure output_root_stream_mpi_nc4
  !! COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE USING STREAM I/O.
  !! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS
  integer :: isp
  character(:), allocatable :: filenamefull
  character(*), parameter :: dims4(4) = [character(7) :: 'x1', 'x2', 'x3', 'species'], &
    dims3(3) = [character(2) :: 'x1', 'x2', 'x3'], &
    dims23(2) = [character(2) :: 'x2', 'x3']
  type(netcdf_file) :: hout

  !> FIGURE OUT THE FILENAME
  filenamefull = date_filename(outdir,ymd,UTsec) // '.nc'
  print *, 'Output file name:  ', filenamefull

  call hout%open(filenamefull, action='w',comp_lvl=comp_lvl)

  call hout%write("flagoutput", flagoutput)

  call hout%write('ymd', ymd)
  call hout%write('UThour',UTsec/3600._wp)

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

  if (gridflag==1) then
    print *, 'Writing topside boundary conditions for inverted-type grid...'
    call hout%write('Phiall', Phiall(1,1:lx2all,1:lx3all), dims23)
  else
    print *, 'Writing topside boundary conditions for non-inverted-type grid...'
    call hout%write('Phiall', Phiall(lx1,1:lx2all,1:lx3all), dims23)
  end if

  call hout%close()
end procedure output_root_stream_mpi_nc4

end submodule plasma_output_nc
