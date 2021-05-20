submodule (neutral:perturb) intrp

use interpolation, only : interp2, interp3

implicit none (type, external)

contains

module procedure spaceinterp_dneu2D
  !must take into account the type of interpolation that is being done
  
  integer :: lhorzn
  real(wp), dimension(:,:), allocatable :: tmpinterp
  real(wp), dimension(lx1*lx2*lx3) :: parami    !work array for temp storage of interpolated data, note sizes taken from grid module data
  logical :: flag2
  real(wp), dimension(:), allocatable :: coord2n,coord2i
  
  
  ! Array for packing neutral data
  if (flagcart) then
    lhorzn=lyn
  else
    lhorzn=lrhon
  end if
  allocate(tmpinterp(lzn,lhorzn))
  
  ! find the singleton dimension
!  if (lx2/=1) then
!    flag2=.true.
!  else if (lx3/=1) then
    flag2=.false.    ! note that because of packing below it doesn't matter whether the singleton dim is 2 or 3 so always leave it as 3...
!  else
!    error stop ' spaceinterp_dneu2D:  cannot determine singleton dimension!!!'
!  end if

  ! set up pointer alias for second coordinate to reduce repetition
  if (flagcart) then
    allocate(coord2n(size(yn,1)),coord2i(size(yi,1)))
    coord2n=yn
    coord2i=yi
  else
    allocate(coord2n(size(rhon,1)),coord2i(size(rhoi,1)))
    coord2n=rhon
    coord2i=rhoi
  end if

  ! interpolate input neutral data
  if (flag2) then
    tmpinterp=dnO(:,:,1)
  else
    tmpinterp=dnO(:,1,:)                    !pack into 2D array for interp2
  end if
  parami=interp2(zn,coord2n,tmpinterp,zi,coord2i)         !interp to temp var.
  dnOiprev=dnOinext                       !save new previous
  dnOinext=reshape(parami,[lx1,lx2,lx3])  !overwrite next with new interpolated input

  if (flag2) then
    tmpinterp=dnN2(:,:,1)  
  else
    tmpinterp=dnN2(:,1,:)
  end if
  parami=interp2(zn,coord2n,tmpinterp,zi,coord2i)
  dnN2iprev=dnN2inext
  dnN2inext=reshape(parami,[lx1,lx2,lx3])

  if (flag2) then
    tmpinterp=dnO2(:,:,1)
  else  
    tmpinterp=dnO2(:,1,:)
  end if
  parami=interp2(zn,coord2n,tmpinterp,zi,coord2i)
  dnO2iprev=dnO2inext
  dnO2inext=reshape(parami,[lx1,lx2,lx3])

  if (flag2) then
    tmpinterp=dvnrho(:,:,1)
  else
    tmpinterp=dvnrho(:,1,:)
  end if
  parami=interp2(zn,coord2n,tmpinterp,zi,coord2i)
  dvnrhoiprev=dvnrhoinext    !interpreted as y-component in this (cartesian) function
  dvnrhoinext=reshape(parami,[lx1,lx2,lx3])

  if (flag2) then
    tmpinterp=dvnz(:,:,1)
  else 
    tmpinterp=dvnz(:,1,:)
  end if
  parami=interp2(zn,coord2n,tmpinterp,zi,coord2i)
  dvnziprev=dvnzinext
  dvnzinext=reshape(parami,[lx1,lx2,lx3])
 
  if (flag2) then
    tmpinterp=dTn(:,:,1)
  else 
    tmpinterp=dTn(:,1,:)
  end if
  parami=interp2(zn,coord2n,tmpinterp,zi,coord2i)
  dTniprev=dTninext
  dTninext=reshape(parami,[lx1,lx2,lx3])

  ! diagnostic print for debugging
  if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
    print *, 'Min/max values for dnOi:  ',minval(dnOinext),maxval(dnOinext)
    print *, 'Min/max values for dnN2i:  ',minval(dnN2inext),maxval(dnN2inext)
    print *, 'Min/max values for dnO2i:  ',minval(dnO2inext),maxval(dnO2inext)
    print *, 'Min/max values for dvrhoi:  ',minval(dvnrhoinext),maxval(dvnrhoinext)
    print *, 'Min/max values for dvnzi:  ',minval(dvnzinext),maxval(dvnzinext)
    print *, 'Min/max values for dTni:  ',minval(dTninext),maxval(dTninext)
  end if
  
  !ROTATE VECTORS INTO X1 X2 DIRECTIONS (Need to include unit vectors with grid
  !structure)
  dvn1iprev=dvn1inext   !save the old data
  dvn2iprev=dvn2inext
  dvn3iprev=dvn3inext
  if(flagcart) then
    dvn1inext=dvnrhoinext*proj_eyp_e1+dvnzinext*proj_ezp_e1    !apply projection to complete rotation into dipole coordinates; drhoi interpreted here at teh y component (northward)
    dvn2inext=dvnrhoinext*proj_eyp_e2+dvnzinext*proj_ezp_e2
    dvn3inext=dvnrhoinext*proj_eyp_e3+dvnzinext*proj_ezp_e3
  else
    dvn1inext=dvnrhoinext*proj_erhop_e1+dvnzinext*proj_ezp_e1    !apply projection to complete rotation into dipole coordinates
    dvn2inext=dvnrhoinext*proj_erhop_e2+dvnzinext*proj_ezp_e2
    dvn3inext=dvnrhoinext*proj_erhop_e3+dvnzinext*proj_ezp_e3
  end if
  
  !MORE DIAGNOSTICS
  if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
    print *, 'Min/max values for dnOi:  ',minval(dnOinext),maxval(dnOinext)
    print *, 'Min/max values for dnN2i:  ',minval(dnN2inext),maxval(dnN2inext)
    print *, 'Min/max values for dnO2i:  ',minval(dnO2inext),maxval(dnO2inext)
    print *, 'Min/max values for dvn1i:  ',minval(dvn1inext),maxval(dvn1inext)
    print *, 'Min/max values for dvn2i:  ',minval(dvn2inext),maxval(dvn2inext)
    print *, 'Min/max values for dvn3i:  ',minval(dvn3inext),maxval(dvn3inext)
    print *, 'Min/max values for dTni:  ',minval(dTninext),maxval(dTninext)
  end if
  
  !CLEAR ALLOCATED VARS
  deallocate(tmpinterp)
  deallocate(coord2n,coord2i)
end procedure spaceinterp_dneu2D


module procedure spaceinterp_dneu3D

!performs spatial interpolation for 3D input neutral data from MAGIC or some other source

real(wp), dimension(lx1*lx2*lx3) :: parami    !work array for temp storage of interpolated data, note sizes taken from grid module data


!INTERPOLATE IN THREE DIMENSIONS
parami=interp3(zn,xn,yn,dnO,zi,xi,yi)         !interp to temp var.
dnOiprev=dnOinext                       !save new previous
dnOinext=reshape(parami,[lx1,lx2,lx3])  !overwrite next with new interpolated input

parami=interp3(zn,xn,yn,dnN2,zi,xi,yi)
dnN2iprev=dnN2inext
dnN2inext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dnO2,zi,xi,yi)
dnO2iprev=dnO2inext
dnO2inext=reshape(parami,[lx1,lx2,lx3])

!ZZZ - do we want to make dvnrho-->dvny???
parami=interp3(zn,xn,yn,dvnrho,zi,xi,yi)
dvnrhoiprev=dvnrhoinext    !interpreted as y-component in this (cartesian) function
dvnrhoinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dvnz,zi,xi,yi)
dvnziprev=dvnzinext
dvnzinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dvnx,zi,xi,yi)
dvnxiprev=dvnxinext
dvnxinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dTn,zi,xi,yi)
dTniprev=dTninext
dTninext=reshape(parami,[lx1,lx2,lx3])


!MORE DIAG
if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
  print *, 'Min/max values for dnOi:  ',mpi_cfg%myid,minval(dnOinext),maxval(dnOinext)
  print *, 'Min/max values for dnN2i:  ',mpi_cfg%myid,minval(dnN2inext),maxval(dnN2inext)
  print *, 'Min/max values for dnO2i:  ',mpi_cfg%myid,minval(dnO2inext),maxval(dnO2inext)
  print *, 'Min/max values for dvrhoi:  ',mpi_cfg%myid,minval(dvnrhoinext),maxval(dvnrhoinext)
  print *, 'Min/max values for dvnzi:  ',mpi_cfg%myid,minval(dvnzinext),maxval(dvnzinext)
  print *, 'Min/max values for dvnxi:  ',mpi_cfg%myid,minval(dvnxinext),maxval(dvnxinext)
  print *, 'Min/max values for dTni:  ',mpi_cfg%myid,minval(dTninext),maxval(dTninext)
end if


!ROTATE VECTORS INTO X1 X2 DIRECTIONS (Need to include unit vectors with grid
!structure)
dvn1iprev=dvn1inext   !save the old data for the rotated vectors
dvn2iprev=dvn2inext
dvn3iprev=dvn3inext
dvn1inext=dvnrhoinext*proj_eyp_e1+dvnzinext*proj_ezp_e1+dvnxinext*proj_exp_e1    !apply projection to complete rotation into dipole coordinates; drhoi interpreted here at teh y component (northward)
dvn2inext=dvnrhoinext*proj_eyp_e2+dvnzinext*proj_ezp_e2+dvnxinext*proj_exp_e2
dvn3inext=dvnrhoinext*proj_eyp_e3+dvnzinext*proj_ezp_e3+dvnxinext*proj_exp_e3


!MORE DIAGNOSTICS
if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
  print *, 'Min/max values for dvn1i:  ',mpi_cfg%myid,minval(dvn1inext),maxval(dvn1inext)
  print *, 'Min/max values for dvn2i:  ',mpi_cfg%myid,minval(dvn2inext),maxval(dvn2inext)
  print *, 'Min/max values for dvn3i:  ',mpi_cfg%myid,minval(dvn3inext),maxval(dvn3inext)
end if

end procedure spaceinterp_dneu3D


module procedure timeinterp_dneu

!interpolation in time - no sensitive to dimensionality of the input neutral data so this can be
! the same for 2D vs. 3D

integer :: ix1,ix2,ix3
real(wp) :: slope

do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      slope=(dnOinext(ix1,ix2,ix3)-dnOiprev(ix1,ix2,ix3))/(tnext-tprev)
      dnOinow(ix1,ix2,ix3)=dnOiprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dnN2inext(ix1,ix2,ix3)-dnN2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dnN2inow(ix1,ix2,ix3)=dnN2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dnO2inext(ix1,ix2,ix3)-dnO2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dnO2inow(ix1,ix2,ix3)=dnO2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn1inext(ix1,ix2,ix3)-dvn1iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn1inow(ix1,ix2,ix3)=dvn1iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn2inext(ix1,ix2,ix3)-dvn2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn2inow(ix1,ix2,ix3)=dvn2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn3inext(ix1,ix2,ix3)-dvn3iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn3inow(ix1,ix2,ix3)=dvn3iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dTninext(ix1,ix2,ix3)-dTniprev(ix1,ix2,ix3))/(tnext-tprev)
      dTninow(ix1,ix2,ix3)=dTniprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)
    end do
  end do
end do


!SOME BASIC DIAGNOSTICS
if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
  print *, 'tprev,t,tnext:  ',mpi_cfg%myid,tprev,t+dt/2d0,tnext
  print *, 'Min/max values for dnOinow:  ',mpi_cfg%myid,minval(dnOinow),maxval(dnOinow)
  print *, 'Min/max values for dnN2inow:  ',mpi_cfg%myid,minval(dnN2inow),maxval(dnN2inow)
  print *, 'Min/max values for dnO2inow:  ',mpi_cfg%myid,minval(dnO2inow),maxval(dnO2inow)
  print *, 'Min/max values for dvn1inow:  ',mpi_cfg%myid,minval(dvn1inow),maxval(dvn1inow)
  print *, 'Min/max values for dvn2inow:  ',mpi_cfg%myid,minval(dvn2inow),maxval(dvn2inow)
  print *, 'Min/max values for dvn3inow:  ',mpi_cfg%myid,minval(dvn3inow),maxval(dvn3inow)
  print *, 'Min/max values for dTninow:  ',mpi_cfg%myid,minval(dTninow),maxval(dTninow)
end if

end procedure timeinterp_dneu


end submodule intrp
