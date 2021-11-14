module errors
!! cleanup before shutdown on error, preserving data if possible
!! We use interface error_stop to work with GCC < 10 instead of "select rank".
!! NOTE: To keep procedures Fortran 2003 TKR-distinct requires
!! at least one TKR-distinct non-optional variable

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use h5fortran, only : hdf5_file
use phys_consts, only : wp

implicit none (type, external)

interface error_stop
  !! polymorphic error handling
  procedure dump_worker, dump_root, dump_input, dump_mag, dump_J, dump_step, dump_perturb
end interface error_stop


private
public :: error_stop

contains

subroutine dump_mag(filename, info, Br, Btheta, Bphi)

character(*), intent(in) :: filename, info
real(wp), dimension(:), intent(in) :: Br, Btheta, Bphi

type (hdf5_file) :: h

call h%open(filename, action="w")
call h%write("/Br", Br)
call h%write("/Btheta", Btheta)
call h%write("/Bphi", Bphi)
call h%write("/info", info)
call h%close()

error stop info // " dumped to " // filename

end subroutine dump_mag


subroutine dump_J(filename, info, J1, J2, J3)

character(*), intent(in) :: filename, info
real(wp), dimension(:,:,:), intent(in) :: J1, J2, J3

type (hdf5_file) :: h

call h%open(filename, action="w")
call h%write("/J1", J1)
call h%write("/J2", J2)
call h%write("/J3", J3)
call h%write("/info", info)
call h%close()

error stop info // " dumped to " // filename

end subroutine dump_J


subroutine dump_perturb(filename, info, t_elapsed, worker_id, nn, Tn, vn1, vn2, vn3)

character(*), intent(in) :: filename, info
real(wp), intent(in) :: t_elapsed
integer, intent(in) :: worker_id

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: Tn, vn1, vn2, vn3

character(6) :: wid

type (hdf5_file) :: h

call h%open(filename, action="w")
call h%write("/worker_id", worker_id)
call h%write("/time/elapsed_seconds", t_elapsed)
call h%write("/nn", nn)
call h%write("/Tn", Tn)
call h%write("/vn1", vn1)
call h%write("/vn2", vn2)
call h%write("/vn3", vn3)
call h%write("/info", info)
call h%close()

write(wid, "(I0)") worker_id

error stop info // " worker " // trim(wid) // " dumped to " // filename

end subroutine dump_perturb


subroutine dump_input(filename, info, ns, vs1, Ts)

character(*), intent(in) :: filename, info
real(wp), dimension(:,:,:,:), intent(in) :: ns, vs1, Ts

type (hdf5_file) :: h

call h%open(filename, action="w")
call h%write("/vs1", vs1)
call h%write("/ns", ns)
call h%write("/Ts", Ts)
call h%write("/info", info)
call h%close()

error stop info // " dumped to " // filename

end subroutine dump_input


subroutine dump_step(filename,info, t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

character(*), intent(in) :: filename, info

real(wp), intent(in) :: t_elapsed
integer, intent(in) :: worker_id

real(wp), intent(in), dimension(:,:,:,:) :: vs2, vs3, ns, vs1, Ts
real(wp), intent(in), dimension(:,:,:) :: Phi, J1, J2, J3

type (hdf5_file) :: h

character(6) :: wid

call h%open(filename, action="w")
call h%write("/worker_id", worker_id)
call h%write("/time/elapsed_seconds", t_elapsed)
call h%write("/vs2", vs2)
call h%write("/vs3", vs3)
call h%write("/ns", ns)
call h%write("/vs1", vs1)
call h%write("/Ts", Ts)
call h%write("/Phi", Phi)
call h%write("/J1", J1)
call h%write("/J2", J2)
call h%write("/J3", J3)
call h%write("/info", info)
call h%close()

write(wid, "(I0)") worker_id

error stop info // " worker " // trim(wid) // " dumped to " // filename

end subroutine dump_step


subroutine dump_worker(filename, info, vs2,vs3,ns,vs1,Ts,J1,J2,J3)

character(*), intent(in) :: filename, info
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3

type (hdf5_file) :: h

call h%open(filename, action="w")
call h%write("/vs1", vs1)
call h%write("/vs2", vs2)
call h%write("/vs3", vs3)
call h%write("/ns", ns)
call h%write("/Ts", Ts)
call h%write("/J1", J1)
call h%write("/J2", J2)
call h%write("/J3", J3)
call h%write("/info", info)
call h%close()

error stop info // " dumped to " // filename

end subroutine dump_worker


subroutine dump_root(filename,info,flagoutput,ymd,UTsec,v2avgall,v3avgall,nsall,vs1all,Tsall, &
  Phiall,J1all,J2all,J3all,neall,v1avgall,Tavgall,Teall)

character(*), intent(in) :: filename, info
integer, intent(in) :: flagoutput

integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: v2avgall,v3avgall
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: nsall,vs1all,Tsall

real(wp), dimension(:,:,:), intent(in) :: Phiall
real(wp), dimension(:,:,:), intent(in) :: J1all,J2all,J3all
real(wp), dimension(:,:,:), intent(in) :: neall,v1avgall,Tavgall,Teall

type (hdf5_file) :: h

call h%open(filename, action="w")
call h%write("/flagoutput", flagoutput)
call h%write("/time/ymd", ymd)
call h%write("/time/UTsec", UTsec)
call h%write("/v2avgall", v2avgall)
call h%write("/v3avgall", v3avgall)
call h%write("/nsall", nsall)
call h%write("/vs1all", vs1all)
call h%write("/Tsall", Tsall)
call h%write("/Phiall", Phiall)
call h%write("/J1all", J1all)
call h%write("/J2all", J2all)
call h%write("/J3all", J3all)
call h%write("/neall", neall)
call h%write("/v1avgall", v1avgall)
call h%write("/Tavgall", Tavgall)
call h%write("/Teall", Teall)
call h%write("/info", info)
call h%close()

error stop info // " dumped to " // filename

end subroutine dump_root

end module errors
