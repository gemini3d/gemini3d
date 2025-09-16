submodule (io:plasma_output) plasma_output_hdf5

use timeutils, only : date_filename
use h5fortran, only: hdf5_file, HID_T, H5T_NATIVE_DOUBLE, H5T_NATIVE_REAL

implicit none (type, external)

contains

module procedure output_root_stream_mpi_hdf5
  !! COLLECT OUTPUT FROM WORKERS AND WRITE TO AN HDF5 FILE
  !! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS
  character(:), allocatable :: hdf5_filename
  type(hdf5_file) :: hout
  integer(HID_T) :: dtype

  logical :: hdf5_real64 = .true.
  !! FIXME: Make "hdf5_real64" a global parameter or even better a config file parameter

  if (hdf5_real64) then
    dtype = H5T_NATIVE_DOUBLE
  else
    dtype = H5T_NATIVE_REAL
  end if

  !> FIGURE OUT THE FILENAME
  hdf5_filename = date_filename(outdir,ymd,UTsec) // '.h5'
  print '(a)', 'HDF5 Output file name:  ' // hdf5_filename

  call hout%open(hdf5_filename, action='w', comp_lvl=comp_lvl)

  call hout%write("/flagoutput", flagoutput)
  call hout%write('/time/ymd', ymd)
  call hout%write('/time/UThour', UTsec/3600.)

  select case (flagoutput)
    case (2)    !output ISR-like average parameters
      call hout%write('neall',    neall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('v1avgall', v1avgall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
      call hout%write('Tavgall',  Tavgall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('TEall',    Teall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('J1all',    J1all(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('J2all',    J2all(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('J3all',    J3all(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('v2avgall', v2avgall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('v3avgall', v3avgall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
    case (3)     !just electron density
      print *, 'INFO:  Input file has selected electron density only output, make sure this is what you really want!'
      call hout%write('neall',    neall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
    case default    !output everything
      print *, 'INFO:  Input file has selected full output or milestones, large files may result!'
      call hout%write('nsall',    nsall(1:lx1,1:lx2all,1:lx3all,:), datatype=dtype)
      call hout%write('vs1all',   vs1all(1:lx1,1:lx2all,1:lx3all,:), datatype=dtype)
      !this is full output of all parameters in 3D
      call hout%write('Tsall',    Tsall(1:lx1,1:lx2all,1:lx3all,:), datatype=dtype)
      call hout%write('J1all',    J1all(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('J2all',    J2all(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('J3all',    J3all(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('v2avgall', v2avgall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)
      call hout%write('v3avgall', v3avgall(1:lx1,1:lx2all,1:lx3all), datatype=dtype)

      ! these are user-specified output variables
      if (size(user_outputall,4)>0) then
        call hout%write('user_outputall', user_outputall(1:lx1,1:lx2all,1:lx3all,:), datatype=dtype)
        !print*, 'Min/max user var written:  ',minval(real(user_outputall)),maxval(real(user_outputall))
      end if
  end select

  if (gridflag==1) then
    print *, 'Writing topside boundary conditions for inverted-type grid...'
    call hout%write('Phiall',     Phiall(1,1:lx2all,1:lx3all), datatype=dtype)
  else
    print *, 'Writing topside boundary conditions for non-inverted-type grid...'
    call hout%write('Phiall',     Phiall(lx1,1:lx2all,1:lx3all), datatype=dtype)
  end if

  call hout%close()
end procedure output_root_stream_mpi_hdf5

end submodule plasma_output_hdf5
