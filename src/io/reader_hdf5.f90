submodule (reader) reader_hdf5

use h5fortran, only: hdf5_file, hsize_t

implicit none (type, external)

contains


module procedure get_simsize2_hdf5
!! get x2 and x3 dimension sizes
type(hdf5_file) :: hf

if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid size from file:', path

call hf%open(path, action='r')

!> scripts can use variety of variable names
if(hf%exist('/llat')) then
  call hf%read('/llat', llat)
elseif(hf%exist('/Nlat')) then
  call hf%read('/Nlat', llat)
elseif(hf%exist('/lx2')) then
  call hf%read('/lx2', llat)
else
  error stop 'ERROR:gemini3d:reader_hdf5:get_simsize2: llat / lx2'
endif

if(hf%exist('/llon')) then
  call hf%read('/llon', llon)
elseif(hf%exist('/Nlat')) then
  call hf%read('/Nlat', llon)
elseif(hf%exist('/lx3')) then
  call hf%read('/lx3', llon)
else
  error stop 'ERROR:gemini3d:reader_hdf5:get_simsize2: llon / lx3'
endif

call hf%close()

end procedure get_simsize2_hdf5


module procedure get_simsize3_hdf5
!! get x1, x2, x3 dimension sizes
!! sizes include Ghost Cells
type(hdf5_file) :: hf

integer :: lx(3)

if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid  size from file:', path

call hf%open(path, action='r')

if (hf%exist("/lx1")) then
  call hf%read('/lx1', lx1)
  call hf%read('/lx2', lx2all)
  if (present(lx3all)) call hf%read('/lx3', lx3all)
elseif (hf%exist("/lxs")) then
  call hf%read("/lxs", lx)
  lx1 = lx(1)
  lx2all = lx(2)
  if (present(lx3all)) lx3all = lx(3)
elseif (hf%exist("/lx")) then
  call hf%read("/lx", lx)
  lx1 = lx(1)
  lx2all = lx(2)
  if (present(lx3all)) lx3all = lx(3)
endif

call hf%close()

end procedure get_simsize3_hdf5


module procedure get_grid2_hdf5
type(hdf5_file) :: hf

if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid:', path

call hf%open(path, action='r')
call hf%read('/mlon', mlonp)
call hf%read('/mlat', mlatp)
call hf%close()

end procedure get_grid2_hdf5


module procedure get_Efield_hdf5
type(hdf5_file) :: hf

if (debug) print *, 'READ electric field data from file:  ',path

call hf%open(path, action='r')

call hf%read('/flagdirich', flagdirich)

call hf%read('/Exit', E0xp)
call hf%read('/Eyit', E0yp)
call hf%read('/Vminx1it', Vminx1p)
call hf%read('/Vmaxx1it', Vmaxx1p)

!! background fields and top/bottom boundary conditions
call hf%read('/Vminx2ist', Vminx2pslice)
call hf%read('/Vmaxx2ist', Vmaxx2pslice)

!! these only used for 3D simulations
call hf%read('/Vminx3ist', Vminx3pslice)
call hf%read('/Vmaxx3ist', Vmaxx3pslice)

call hf%close()

end procedure get_Efield_hdf5


module procedure get_precip_hdf5
!! Qp, E0p are (llon, llat)
type(hdf5_file) :: hf

if (debug) print *, 'READ precipitation data from file:  ',path

call hf%open(path, action='r')

call hf%read('/Qp', Qp)
call hf%read('/E0p', E0p)

call hf%close()

end procedure get_precip_hdf5


module procedure get_neutral2_hdf5
type(hdf5_file) :: hf

real(wp), dimension(1:size(dnO,1),1:size(dnO,2)) :: buffer
!! FIXME: second dimension is singleton, so read data into a buffer
!! This is something that should be corrected in the Matlab/Python scripts writing the files
!! FIXME (MZ): inputdata refactor uses dim=3 as the singleton so this needed to be patched accordingly
integer :: lz,lx,ly
integer(hsize_t), allocatable :: dims(:)

lz=size(dnO,1)
lx=size(dnO,3)   !should be 1
ly=size(dnO,2)

if (debug) print *, 'READ neutral 2D data from file: ', path

call hf%open(path, action='r')

call hf%shape("/dn0all", dims)
if(size(dims) == 3) then
  call hf%read('/dn0all', dnO)
  call hf%read('/dnN2all', dnN2)
  call hf%read('/dnO2all', dnO2)
  call hf%read('/dvnrhoall', dvnrho)
  call hf%read('/dvnzall', dvnz)
  call hf%read('/dTnall', dTn)
else
  !! FIXME: workaround, the Matlab script should be fixed to write shape [lz,1,ly]
  call hf%read('/dn0all', buffer)
  dnO=reshape(buffer,[lz,ly,1])
  call hf%read('/dnN2all', buffer)
  dnN2=reshape(buffer,[lz,ly,1])
  call hf%read('/dnO2all', buffer)
  dnO2=reshape(buffer,[lz,ly,1])
  call hf%read('/dvnrhoall', buffer)
  dvnrho=reshape(buffer,[lz,ly,1])
  call hf%read('/dvnzall', buffer)
  dvnz=reshape(buffer,[lz,ly,1])
  call hf%read('/dTnall', buffer)
  dTn=reshape(buffer,[lz,ly,1])
endif

call hf%close()

end procedure get_neutral2_hdf5


module procedure get_neutral3_hdf5
type(hdf5_file) :: hf

if (debug) print *, 'READ neutral 3D data from file: ',path

call hf%open(path, action='r')

call hf%read('/dn0all', dnOall)
call hf%read('/dnN2all', dnN2all)
call hf%read('/dnO2all', dnO2all)
call hf%read('/dvnxall', dvnxall)
call hf%read('/dvnrhoall', dvnrhoall)
call hf%read('/dvnzall', dvnzall)
call hf%read('/dTnall', dTnall)

call hf%close()

end procedure get_neutral3_hdf5


end submodule reader_hdf5
