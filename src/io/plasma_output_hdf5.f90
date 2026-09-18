submodule (io:plasma_output) plasma_output_hdf5

use timeutils, only : date_filename
use h5fortran, only: hdf5_file
use atomic_file, only: publish_file
use restart_runtime, only: runtime_output_path, runtime_output_header

implicit none (type, external)

contains

module procedure output_root_stream_mpi_hdf5
  !! COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE USING STREAM I/O.
  !! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS
  character(:), allocatable :: filenamefull, temporary
  type(hdf5_file) :: hout

  !> FIGURE OUT THE FILENAME
  filenamefull = date_filename(outdir,ymd,UTsec) // '.h5'
  print *, 'HDF5 Output file name:  ', filenamefull

  temporary=runtime_output_path(filenamefull)
  call hout%open(temporary, action='w',comp_lvl=comp_lvl)

  call hout%write("/flagoutput", flagoutput)
  call hout%write('/time/ymd', ymd)
  call hout%write('/time/UThour', UTsec/3600.)

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

      ! these are user-specified output variables
      if (size(user_outputall,4)>0) then
        call hout%write('user_outputall', real(user_outputall(1:lx1,1:lx2all,1:lx3all,:)))
        !print*, 'Min/max user var written:  ',minval(real(user_outputall)),maxval(real(user_outputall))
      end if

      !if (size(production_rateall,4)>0) then
      !   call hout%write('production_rateall', real(production_rateall(1:lx1,1:lx2all,1:lx3all,:)))
      !end if
  end select

  if (gridflag==1) then
    print *, 'Writing topside boundary conditions for inverted-type grid...'
    call hout%write('Phiall',       real(Phiall(1,1:lx2all,1:lx3all)))
  else
    print *, 'Writing topside boundary conditions for non-inverted-type grid...'
    call hout%write('Phiall',       real(Phiall(lx1,1:lx2all,1:lx3all)))
  end if

  if (flagoutput==1) then
    ! Analysis arrays retain their historical layout and precision. The core
    ! restart record preserves solver precision and the full field potential.
    ! This is not a claim that mode-specific auxiliary state is complete.
    call hout%write('/restart_core/schema',1)
    call hout%write('/restart_core/realbits',storage_size(UTsec))
    call hout%write('/restart_core/ns',nsall(1:lx1,1:lx2all,1:lx3all,:))
    call hout%write('/restart_core/vs1',vs1all(1:lx1,1:lx2all,1:lx3all,:))
    call hout%write('/restart_core/Ts',Tsall(1:lx1,1:lx2all,1:lx3all,:))
    call hout%write('/restart_core/Phi',Phiall(1:lx1,1:lx2all,1:lx3all))
    ! Written last so an interrupted record is rejected on restart.
    call hout%write('/restart_core/complete',1)
  endif
  call runtime_output_header(hout)
  call hout%close()
  call publish_file(temporary,filenamefull)
end procedure output_root_stream_mpi_hdf5

end submodule plasma_output_hdf5
