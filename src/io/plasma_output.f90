submodule (io) plasma_output

use timeutils, only : date_filename
use h5fortran, only : hdf5_file

implicit none (type, external)

contains

module procedure output_plasma

! subroutine output_plasma(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

!! STATE VARS ARE EXPECTED TO INCLUDE GHOST CELLS
!! This is the same regardless of what type of output is being done.

integer :: isp
real(wp), dimension(lx1, lx2, lx3) :: v2avg,v3avg,v1avg,Tavg, ne, Te

character(:), allocatable :: filenamefull
type(hdf5_file) :: h

integer, dimension(4) :: dims4, i0_4, i1_4
integer, dimension(3) :: dims, i0, i1
integer, dimension(2) :: ix1, ix2, ix3
!! HDF5-MPI partitioning


if(mpi_cfg%myid == 0) print '(a,i0,1x,i0,1x,i0)', 'System sizes according to Phiall:  ',lx1,lx2all,lx3all
!ONLY AVERAGE DRIFTS PERP TO B NEEDED FOR OUTPUT
v2avg = sum(ns(1:lx1, 1:lx2, 1:lx3, 1:lsp-1)*vs2(1:lx1, 1:lx2, 1:lx3, 1:lsp-1),4)
v2avg = v2avg / ns(1:lx1, 1:lx2, 1:lx3, lsp)
!! compute averages for output.

v3avg = sum(ns(1:lx1, 1:lx2, 1:lx3, 1:lsp-1)*vs3(1:lx1, 1:lx2, 1:lx3, 1:lsp-1),4)
v3avg=v3avg/ns(1:lx1, 1:lx2, 1:lx3, lsp)

!COMPUTE AVERAGE VALUES FOR ION PLASMA PARAMETERS
!> possible bottleneck; should have workers help?
!> also only compute these if they are actually being output
if (flagoutput==2 .or. flagoutput==3) then
  ne = ns(1:lx1, 1:lx2, 1:lx3, lsp)
end if

if (flagoutput==2) then
  v1avg = sum(ns(1:lx1, 1:lx2, 1:lx3, 1:lsp-1) * vs1(1:lx1, 1:lx2, 1:lx3, 1:lsp-1), dim=4)
  v1avg = v1avg / ns(1:lx1, 1:lx2, 1:lx3, lsp)
  !! compute averages for output.

  Tavg = sum(ns(1:lx1, 1:lx2, 1:lx3, 1:lsp-1) * Ts(1:lx1, 1:lx2, 1:lx3, 1:lsp-1), dim=4)
  Tavg = Tavg / ns(1:lx1, 1:lx2, 1:lx3, lsp)
  !! compute averages for output.
  Te = Ts(1:lx1, 1:lx2, 1:lx3, lsp)
end if


!! COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE
!! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS

dims = [lx1, lx2all, lx3all]
dims4(:3) = dims
dims4(4) = lsp

ix1 = [1, lx1]
ix2 = [mpi_cfg%myid*lx2+1, (mpi_cfg%myid+1)*lx2]
ix3 = [mpi_cfg%myid*lx3+1, (mpi_cfg%myid+1)*lx3]
if(lx2 == 1) then
  ix2 = [1, 1]
elseif(lx3 == 1) then
  ix3 = [1, 1]
endif

i0 = [ix1(1), ix2(1), ix3(1)]
i1 = [ix1(2), ix2(2), ix3(2)]

i0_4(:3) = i0
i0_4(4) = 1
i1_4(:3) = i1
i1_4(4) = lsp

! print '(a,5i4)', "TRACE:gemini3d:output_plasma: lx2all,lx3all,lx2,lx3",lx2all,lx3all,lx2,lx3
! print '(a,i0,a,i0,1x,i0,1x,i0,a,i0,1x,i0,1x,i0)', "TRACE:gemini3d:output_plasma: myid: ", mpi_cfg%myid,&
!  " istart: ", i0, " iend: ", i1
! print *, "ns: shape, lbound, ubound ", shape(ns), lbound(ns), ubound(ns)


filenamefull = date_filename(outdir,ymd,UTsec) // '.h5'
if(mpi_cfg%myid == 0) print *, 'GEMINI3D: plasma file name: ', filenamefull

call h%open(filenamefull, action='w', comp_lvl=comp_lvl, mpi=.true.)

call h%write("/flagoutput", flagoutput)

call h%write('/time/ymd', ymd)
call h%write('/time/UThour',   real(UTsec/3600))

select case (flagoutput)
  case (2)    !output ISR-like average parameters
    call h%write('neall',    real(ne), dset_dims=dims, istart=i0, iend=i1)
    call h%write('v1avgall', real(v1avg), dset_dims=dims, istart=i0, iend=i1)
    !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
    call h%write('Tavgall',  real(Tavg), dset_dims=dims, istart=i0, iend=i1)
    call h%write('TEall',    real(Te), dset_dims=dims, istart=i0, iend=i1)
    call h%write('J1all',    real(J1(1:lx1, 1:lx2, 1:lx3)), dset_dims=dims, istart=i0, iend=i1)
    call h%write('J2all',    real(J2(1:lx1, 1:lx2, 1:lx3)), dset_dims=dims, istart=i0, iend=i1)
    call h%write('J3all',    real(J3(1:lx1, 1:lx2, 1:lx3)), dset_dims=dims, istart=i0, iend=i1)
    call h%write('v2avgall', real(v2avg), dset_dims=dims, istart=i0, iend=i1)
    call h%write('v3avgall', real(v3avg), dset_dims=dims, istart=i0, iend=i1)
  case (3)     !just electron density
    if(mpi_cfg%myid == 0) print *, 'GEMINI3D: electron density only output, make sure this is what you really want!'
    call h%write('neall',    real(ne), dset_dims=dims, istart=i0, iend=i1)
  case default    !output everything
    if(mpi_cfg%myid == 0) print *, 'GEMINI3D: full output or milestones, large files may result!'
    call h%write('nsall',    real(ns(1:lx1, 1:lx2, 1:lx3, :)), dset_dims=dims4, istart=i0_4, iend=i1_4)
    call h%write('vs1all',   real(vs1(1:lx1, 1:lx2, 1:lx3, :)), dset_dims=dims4, istart=i0_4, iend=i1_4)
    !this is full output of all parameters in 3D
    call h%write('Tsall',    real(Ts(1:lx1, 1:lx2, 1:lx3, :)), dset_dims=dims4, istart=i0_4, iend=i1_4)

    call h%write('J1all',    real(J1(1:lx1, 1:lx2, 1:lx3)), dset_dims=dims, istart=i0, iend=i1)
    call h%write('J2all',    real(J2(1:lx1, 1:lx2, 1:lx3)), dset_dims=dims, istart=i0, iend=i1)
    call h%write('J3all',    real(J3(1:lx1, 1:lx2, 1:lx3)), dset_dims=dims, istart=i0, iend=i1)
    call h%write('v2avgall', real(v2avg), dset_dims=dims, istart=i0, iend=i1)
    call h%write('v3avgall', real(v3avg), dset_dims=dims, istart=i0, iend=i1)
end select

call h%close()


if(mpi_cfg%myid == 0) then

  call h%open(filenamefull, action='a', comp_lvl=comp_lvl, mpi=.false.)

  if (gridflag==1) then
    print *, 'Writing topside boundary conditions for inverted-type grid...'
    call h%write('Phiall',   real(Phiall(1, 1:lx2all, 1:lx3all)))
  else
    print *, 'Writing topside boundary conditions for non-inverted-type grid...'
    call h%write('Phiall',   real(Phiall(lx1, 1:lx2all, 1:lx3all)))
  end if

  call h%close()

endif


end procedure output_plasma

end submodule plasma_output
