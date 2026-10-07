program size_alias
use reader, only: get_simsize2
use h5fortran, only: hdf5_file
implicit none(type,external)
type(hdf5_file) :: f
integer :: nlat,nlon
character(4096) :: path
call get_command_argument(1,path)
call f%open(trim(path),action='w')
call f%write('Nlat',3);call f%write('Nlon',5)
call f%close()
call get_simsize2(trim(path),nlon,nlat)
if(nlat/=3.or.nlon/=5) error stop 'Nlat/Nlon aliases mapped to wrong dimension'
end program
