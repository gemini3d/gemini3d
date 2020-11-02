submodule (neutral:perturb) proj

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit

use grid, only : clear_unitvecs
use reader, only : get_simsize3

implicit none (type, external)

contains

module procedure gridproj_dneu2D
!! gridproj_dneu2D(cfg,flagcart,x)
!Read in the grid for the neutral data and project unit vectors into the appropriiate directions.
!Also allocate module-scope variables for storing neutral perturbations read in from input files.

!whether or not the input data are to be interpreted as Cartesian
!inout to allow deallocation of unit vectors once we are done with them, should consider exporting this to another functino to be called from main program to avoid having x writeable...

real(wp) :: dhorzn           !neutral grid spacing in horizontal "rho or y" and vertical directions
integer :: lhorzn
real(wp) :: meanyn

real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
real(wp) :: xp,yp
real(wp), dimension(3) :: erhop,ezp,eyp,tmpvec
real(wp) :: tmpsca

integer :: ix1,ix2,ix3,ihorzn,izn,iid,ierr
real(wp), dimension(x%lx1,x%lx2,x%lx3) :: zimat,rhoimat,yimat


!horizontal grid spacing
dhorzn=cfg%drhon

!Establish the size of the grid based on input file and distribute to workers
if (myid==0) then    !root
  print '(A,/,A)', 'Inputting neutral size from:  ',cfg%sourcedir

  ! bit of a tricky issue here; for neutral input, according to makedneuframes.m, the first integer in the size file is
  !  the horizontal grid point count for the input - which get_simsize3 interprets as lx1...
  !call get_simsize3(cfg%sourcedir, lx1=lzn, lx2all=lhorzn)
call get_simsize3(cfg%sourcedir, lx1=lhorzn, lx2all=lzn)

  print *, 'Neutral data has lhorzn,lz size:  ',lhorzn,lzn,' with spacing dhorzn,dz',dhorzn,cfg%dzn
  if (lhorzn < 1 .or. lzn < 1) then
    write(stderr,*) 'ERROR: reading ' // cfg%sourcedir
    error stop 'neutral:gridproj_dneu2D: grid size must be strictly positive'
  endif
  do iid=1,lid-1
    call mpi_send(lhorzn,1,MPI_INTEGER,iid,tag%lrho,MPI_COMM_WORLD,ierr)
    call mpi_send(lzn,1,MPI_INTEGER,iid,tag%lz,MPI_COMM_WORLD,ierr)
  end do
else                 !workers
  call mpi_recv(lhorzn,1,MPI_INTEGER,0,tag%lrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(lzn,1,MPI_INTEGER,0,tag%lz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if


!Everyone must allocate space for the grid of input data
allocate(zn(lzn))    !these are module-scope variables
if (flagcart) then
  allocate(rhon(1))  !not used in Cartesian code so just set to something
  allocate(yn(lhorzn))
  lyn=lhorzn
else
  allocate(rhon(lhorzn))
  allocate(yn(1))    !not used in the axisymmetric code so just initialize to something
  lrhon=lhorzn
end if


!Note that the second dimension ("longitude") is singleton so that we are able to also use these vars for 3D input
allocate(dnO(lzn,1,lhorzn),dnN2(lzn,1,lhorzn),dnO2(lzn,1,lhorzn),dvnrho(lzn,1,lhorzn),dvnz(lzn,1,lhorzn),dTn(lzn,1,lhorzn))


!Define a grid (input data) by assuming that the spacing is constant
if (flagcart) then     !Cartesian neutral simulation
  yn=[ ((real(ihorzn,8)-1._wp)*dhorzn, ihorzn=1,lhorzn) ]
  meanyn=sum(yn,1)/size(yn,1)
  yn=yn-meanyn     !the neutral grid should be centered on zero for a cartesian interpolation
else
  rhon=[ ((real(ihorzn,8)-1._wp)*dhorzn, ihorzn=1,lhorzn) ]
end if
zn=[ ((real(izn,8)-1._wp)*cfg%dzn, izn=1,lzn) ]

if (myid==0) then
  if (flagcart) then
    print *, 'Creating neutral grid with y,z extent:',minval(yn),maxval(yn),minval(zn),maxval(zn)
  else
    print *, 'Creating neutral grid with rho,z extent:  ',minval(rhon),maxval(rhon),minval(zn),maxval(zn)
  end if
end if


!Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
phi1=cfg%sourcemlon*pi/180d0
theta1=pi/2d0-cfg%sourcemlat*pi/180d0


!Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
if (myid==0) then
  print *, 'Computing alt,radial distance values for plasma grid and completing rotations'
end if
zimat=x%alt     !vertical coordinate
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      !INTERPOLATION BASED ON GEOMAGNETIC COORDINATES
      theta2=x%theta(ix1,ix2,ix3)                    !field point zenith angle
      if (lx2/=1) then
        phi2=x%phi(ix1,ix2,ix3)                      !field point azimuth, full 3D calculation
      else
        phi2=phi1                                    !assume the longitude is the samem as the source in 2D, i.e. assume the source epicenter is in the meridian of the grid
      end if


      !COMPUTE DISTANCES
      gammarads=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)     !this is actually cos(gamma)
      if (gammarads>1._wp) then     !handles weird precision issues in 2D
        gammarads=1._wp
      else if (gammarads<-1._wp) then
        gammarads=-1._wp
      end if
      gammarads=acos(gammarads)                     !angle between source location annd field point (in radians)
      rhoimat(ix1,ix2,ix3)=Re*gammarads    !rho here interpreted as the arc-length defined by angle between epicenter and ``field point''

      !we need a phi locationi (not spherical phi, but azimuth angle from epicenter), as well, but not for interpolation - just for doing vector rotations
      theta3=theta2
      phi3=phi1
      gamma1=cos(theta2)*cos(theta3)+sin(theta2)*sin(theta3)*cos(phi2-phi3)
      if (gamma1>1._wp) then     !handles weird precision issues in 2D
        gamma1=1._wp
      else if (gamma1<-1._wp) then
        gamma1=-1._wp
      end if
      gamma1=acos(gamma1)

      gamma2=cos(theta1)*cos(theta3)+sin(theta1)*sin(theta3)*cos(phi1-phi3)
      if (gamma2>1._wp) then     !handles weird precision issues in 2D
        gamma2=1._wp
      else if (gamma2<-1._wp) then
        gamma2=-1._wp
      end if
      gamma2=acos(gamma2)

      xp=Re*gamma1
      yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement


      !COMPUTE COORDIANTES FROM DISTANCES
      if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
        yp = -yp            !do we want an abs here to be safe
      end if
      if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than soure point (phi3=phi1)
        xp = -xp
      end if
      phip=atan2(yp,xp)


      if(flagcart) then
        yimat(ix1,ix2,ix3)=yp
      end if


      !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
      !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections so we don't have to do it repeatedly as sim runs
      ezp=x%er(ix1,ix2,ix3,:)
      tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_ezp_e2(ix1,ix2,ix3)=tmpsca

      tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_ezp_e1(ix1,ix2,ix3)=tmpsca

      tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
      proj_ezp_e3(ix1,ix2,ix3)=tmpsca

      if (flagcart) then
        eyp=-1._wp*x%etheta(ix1,ix2,ix3,:)

        tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_eyp_e1(ix1,ix2,ix3)=tmpsca

        tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_eyp_e2(ix1,ix2,ix3)=tmpsca

        tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_eyp_e3(ix1,ix2,ix3)=tmpsca
      else
        erhop=cos(phip)*x%e3(ix1,ix2,ix3,:)+(-1._wp)*sin(phip)*x%etheta(ix1,ix2,ix3,:)     !unit vector for azimuth (referenced from epicenter - not geocenter!!!) in cartesian geocentric-geomagnetic coords.

        tmpvec=erhop*x%e1(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_erhop_e1(ix1,ix2,ix3)=tmpsca

        tmpvec=erhop*x%e2(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_erhop_e2(ix1,ix2,ix3)=tmpsca

        tmpvec=erhop*x%e3(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_erhop_e3(ix1,ix2,ix3)=tmpsca
      end if
    end do
  end do
end do


!Assign values for flat lists of grid points
zi=pack(zimat,.true.)     !create a flat list of grid points to be used by interpolation ffunctions
if (flagcart) then
  yi=pack(yimat,.true.)
else
  rhoi=pack(rhoimat,.true.)
end if


!GRID UNIT VECTORS NO LONGER NEEDED ONCE PROJECTIONS ARE CALCULATED...
call clear_unitvecs(x)


!PRINT OUT SOME BASIC INFO ABOUT THE GRID THAT WE'VE LOADED
if (myid==0 .and. debug) then
  if (flagcart) then
    print *, 'Min/max yn,zn values',minval(yn),maxval(yn),minval(zn),maxval(zn)
    print *, 'Min/max yi,zi values',minval(yi),maxval(yi),minval(zi),maxval(zi)
  else
    print *, 'Min/max rhon,zn values',minval(rhon),maxval(rhon),minval(zn),maxval(zn)
    print *, 'Min/max rhoi,zi values',minval(rhoi),maxval(rhoi),minval(zi),maxval(zi)
  end if

  print *, 'Source lat/long:  ',cfg%sourcemlat,cfg%sourcemlon
  print *, 'Plasma grid lat range:  ',minval(x%glat(:,:,:)),maxval(x%glat(:,:,:))
  print *, 'Plasma grid lon range:  ',minval(x%glon(:,:,:)),maxval(x%glon(:,:,:))
end if

end procedure gridproj_dneu2D


module procedure gridproj_dneu3D
!! gridproj_dneu3D(cfg,x)

!Read in the grid for the neutral data and project unit vectors into the appropriiate directions.
!Also allocate module-scope variables for storing neutral perturbations read in from input files.

!inout to allow deallocation of unit vectors once we are done with them, should consider exporting this to another functino to be called from main program to avoid having x writeable...

real(wp) :: meanyn
real(wp) :: meanxn

real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
real(wp) :: xp,yp
real(wp), dimension(3) :: erhop,ezp,eyp,tmpvec,exprm
real(wp) :: tmpsca

integer :: ix1,ix2,ix3,iyn,izn,ixn,iid,ierr
real(wp), dimension(x%lx1,x%lx2,x%lx3) :: zimat,rhoimat,yimat,ximat

real(wp) :: maxzn
real(wp), dimension(2) :: xnrange,ynrange
integer, dimension(6) :: indices


!Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
phi1=cfg%sourcemlon*pi/180._wp
theta1=pi/2._wp-cfg%sourcemlat*pi/180._wp


!Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
if (myid==0) then
  print *, 'Computing alt,radial distance values for plasma grid and completing rotations'
end if
zimat=x%alt     !vertical coordinate
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      !INTERPOLATION BASED ON GEOMAGNETIC COORDINATES
      theta2=x%theta(ix1,ix2,ix3)                    !field point zenith angle
      if (lx2/=1) then
        phi2=x%phi(ix1,ix2,ix3)                      !field point azimuth, full 3D calculation
      else
        phi2=phi1                                    !assume the longitude is the samem as the source in 2D, i.e. assume the source epicenter is in the meridian of the grid
      end if


      !!COMPUTE DISTANCES - ZZZ possibly superfluous for 3D case???
      !gammarads=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)     !this is actually cos(gamma)
      !if (gammarads>1._wp) then     !handles weird precision issues in 2D
      !  gammarads=1._wp
      !else if (gammarads<-1._wp) then
      !  gammarads=-1._wp
      !end if
      !gammarads=acos(gammarads)                     !angle between source location annd field point (in radians)
      !rhoimat(ix1,ix2,ix3)=Re*gammarads    !rho here interpreted as the arc-length defined by angle between epicenter and ``field point''
      !! ZZZ end possibly superfluous block of code...

      !we need a phi locationi (not spherical phi, but azimuth angle from epicenter), as well, but not for interpolation - just for doing vector rotations
      theta3=theta2
      phi3=phi1
      gamma1=cos(theta2)*cos(theta3)+sin(theta2)*sin(theta3)*cos(phi2-phi3)
      if (gamma1>1._wp) then     !handles weird precision issues in 2D
        gamma1=1._wp
      else if (gamma1<-1._wp) then
        gamma1=-1._wp
      end if
      gamma1=acos(gamma1)

      gamma2=cos(theta1)*cos(theta3)+sin(theta1)*sin(theta3)*cos(phi1-phi3)
      if (gamma2>1._wp) then     !handles weird precision issues in 2D
        gamma2=1._wp
      else if (gamma2<-1._wp) then
        gamma2=-1._wp
      end if
      gamma2=acos(gamma2)
      xp=Re*gamma1
      yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement


      !COMPUTE COORDIANTES FROM DISTANCES
      if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
        yp=-1._wp*yp            !do we want an abs here to be safe
      end if
      if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than soure point (phi3=phi1)
        xp=-1._wp*xp
      end if
      !phip=atan2(yp,xp)

      ximat(ix1,ix2,ix3)=xp     !eastward distance
      yimat(ix1,ix2,ix3)=yp     !northward distance


      !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
      !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections so we don't have to do it repeatedly as sim runs
      ezp=x%er(ix1,ix2,ix3,:)

      tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_ezp_e2(ix1,ix2,ix3)=tmpsca

      tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_ezp_e1(ix1,ix2,ix3)=tmpsca

      tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
      proj_ezp_e3(ix1,ix2,ix3)=tmpsca

      eyp=-1._wp*x%etheta(ix1,ix2,ix3,:)

      tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_eyp_e1(ix1,ix2,ix3)=tmpsca

      tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_eyp_e2(ix1,ix2,ix3)=tmpsca

      tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_eyp_e3(ix1,ix2,ix3)=tmpsca

      exprm=x%ephi(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)

      tmpvec=exprm*x%e1(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_exp_e1(ix1,ix2,ix3)=tmpsca

      tmpvec=exprm*x%e2(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_exp_e2(ix1,ix2,ix3)=tmpsca

      tmpvec=exprm*x%e3(ix1,ix2,ix3,:)
      tmpsca=sum(tmpvec)
      proj_exp_e3(ix1,ix2,ix3)=tmpsca

    end do
  end do
end do


!Assign values for flat lists of grid points
zi=pack(zimat,.true.)     !create a flat list of grid points to be used by interpolation functions
yi=pack(yimat,.true.)
xi=pack(ximat,.true.)


!GRID UNIT VECTORS NO LONGER NEEDED ONCE PROJECTIONS ARE CALCULATED, so go ahead and free some space
if (myid==0) then
  print*, '...Clearing out unit vectors (after projections)...'
end if
call clear_unitvecs(x)


if(myid==0) then
  print*, 'Projection checking:  ',minval(proj_exp_e1),maxval(proj_exp_e1),minval(proj_exp_e2),maxval(proj_exp_e2), &
                                    minval(proj_exp_e3),maxval(proj_exp_e3)
end if


!Establish the size of the grid based on input file and distribute to workers
if (myid==0) then    !root
  print '(A,/,A)', 'READ neutral size from:', cfg%sourcedir

  call get_simsize3(cfg%sourcedir, lx1=lxnall, lx2all=lynall, lx3all=lzn)

  print *, 'Neutral data has lx,ly,lz size:  ',lxnall,lynall,lzn,' with spacing dx,dy,dz',cfg%dxn,cfg%drhon,cfg%dzn
  if (lxnall < 1 .or. lynall < 1 .or. lzn < 1) then
    write(stderr,*) 'ERROR: reading ' // cfg%sourcedir
    error stop 'neutral:gridproj_dneu3D: grid size must be strictly positive'
  endif

  !root must allocate space for the entire grid of input data - this might be doable one parameter at a time???
  allocate(zn(lzn))        !the z coordinate is never split up in message passing - want to use full altitude range...
  allocate(rhon(1))        !not used in Cartesian or 3D code so just set to something; could likely be left unallocated
  allocate(xnall(lxnall))
  allocate(ynall(lynall))
  allocate(dnOall(lzn,lxnall,lynall),dnN2all(lzn,lxnall,lynall),dnO2all(lzn,lxnall,lynall),dvnrhoall(lzn,lxnall,lynall), &
              dvnzall(lzn,lxnall,lynall),dvnxall(lzn,lxnall,lynall),dTnall(lzn,lxnall,lynall))    !ZZZ - note that these might be deallocated after each read to clean up memory management a bit...


  !calculate the z grid (same for all) and distribute to workers so we can figure out their x-y slabs
  print*, '...creating vertical grid and sending to workers...'
  zn=[ ((real(izn,8)-1._wp)*cfg%dzn, izn=1,lzn) ]    !root calculates and distributes but this is the same for all workers - assmes that the max neutral grid extent in altitude is always less than the plasma grid (should almost always be true)
  maxzn=maxval(zn)
  do iid=1,lid-1
    call mpi_send(lzn,1,MPI_INTEGER,iid,tag%lz,MPI_COMM_WORLD,ierr)
    call mpi_send(zn,lzn,mpi_realprec,iid,tag%zn,MPI_COMM_WORLD,ierr)
  end do


  !Define a global neutral grid (input data) by assuming that the spacing is constant
  ynall=[ ((real(iyn,8)-1._wp)*cfg%drhon, iyn=1,lynall) ]
  meanyn=sum(ynall,1)/size(ynall,1)
  ynall=ynall-meanyn     !the neutral grid should be centered on zero for a cartesian interpolation
  xnall=[ ((real(ixn,8)-1._wp)*cfg%dxn, ixn=1,lxnall) ]
  meanxn=sum(xnall,1)/size(xnall,1)
  xnall=xnall-meanxn     !the neutral grid should be centered on zero for a cartesian interpolation
  print *, 'Created full neutral grid with y,z extent:',minval(xnall),maxval(xnall),minval(ynall), &
                maxval(ynall),minval(zn),maxval(zn)


  !calculate the extent of my piece of the grid using max altitude specified for the neutral grid
  call slabrange(maxzn,ximat,yimat,zimat,cfg%sourcemlat,xnrange,ynrange)
  allocate(extents(0:lid-1,6),indx(0:lid-1,6),slabsizes(0:lid-1,2))
  extents(0,1:6)=[0._wp,maxzn,xnrange(1),xnrange(2),ynrange(1),ynrange(2)]


  !receive extents of each of the other workers: extents(lid,6)
  print*, 'Receiving xn and yn ranges from workers...'
  do iid=1,lid-1
    call mpi_recv(xnrange,2,mpi_realprec,iid,tag%xnrange,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call mpi_recv(ynrange,2,mpi_realprec,iid,tag%ynrange,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    extents(iid,1:6)=[0._wp,maxzn,xnrange(1),xnrange(2),ynrange(1),ynrange(2)]     !need to store values as xnrange overwritten for each worker
    print*, 'Subgrid extents:  ',iid,extents(iid,:)
  end do


  !find index into into neutral arrays for each worker:  indx(lid,6)
  print*, 'Root grid check:  ',ynall(1),ynall(lynall)
  print*, 'Converting ranges to indices...'
  do iid=0,lid-1
    call range2inds(extents(iid,1:6),zn,xnall,ynall,indices)
    indx(iid,1:6)=indices
    print*, 'Subgrid indices',iid,indx(iid,:)
  end do


  !send each worker the sizes for their particular chunk (all different) and send worker that grid chunk
  print*,'Sending sizes and xn,yn subgrids to workers...'
  do iid=1,lid-1
    lxn=indx(iid,4)-indx(iid,3)+1
    lyn=indx(iid,6)-indx(iid,5)+1
    slabsizes(iid,1:2)=[lxn,lyn]
    call mpi_send(lyn,1,MPI_INTEGER,iid,tag%lrho,MPI_COMM_WORLD,ierr)
    call mpi_send(lxn,1,MPI_INTEGER,iid,tag%lx,MPI_COMM_WORLD,ierr)
    allocate(xn(lxn),yn(lyn))
    xn=xnall(indx(iid,3):indx(iid,4))
    yn=ynall(indx(iid,5):indx(iid,6))
    call mpi_send(xn,lxn,mpi_realprec,iid,tag%xn,MPI_COMM_WORLD,ierr)
    call mpi_send(yn,lyn,mpi_realprec,iid,tag%yn,MPI_COMM_WORLD,ierr)
    deallocate(xn,yn)
  end do


  !have root store its part to the full neutral grid
  print*, 'Root is picking out its own subgrid...'
  lxn=indx(0,4)-indx(0,3)+1
  lyn=indx(0,6)-indx(0,5)+1
  slabsizes(0,1:2)=[lxn,lyn]
  allocate(xn(lxn),yn(lyn))
  xn=xnall(indx(0,3):indx(0,4))
  yn=ynall(indx(0,5):indx(0,6))
else                 !workers
  !get teh z-grid from root so we know what the max altitude we have to deal with will be
  call mpi_recv(lzn,1,MPI_INTEGER,0,tag%lz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  allocate(zn(lzn))
  call mpi_recv(zn,lzn,mpi_realprec,0,tag%zn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  maxzn=maxval(zn)


  !calculate the extent of my grid
  call slabrange(maxzn,ximat,yimat,zimat,cfg%sourcemlat,xnrange,ynrange)


  !send ranges to root
  call mpi_send(xnrange,2,mpi_realprec,0,tag%xnrange,MPI_COMM_WORLD,ierr)
  call mpi_send(ynrange,2,mpi_realprec,0,tag%ynrange,MPI_COMM_WORLD,ierr)


  !receive my sizes from root, allocate then receive my pieces of the grid
  call mpi_recv(lxn,1,MPI_INTEGER,0,tag%lx,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(lyn,1,MPI_INTEGER,0,tag%lrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  allocate(xn(lxn),yn(lyn))
  call mpi_recv(xn,lxn,mpi_realprec,0,tag%xn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(yn,lyn,mpi_realprec,0,tag%yn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if


!AT THIS POINT WE CAN ALLOCATE THE SUBGRID SIZES
allocate(dnO(lzn,lxn,lyn),dnN2(lzn,lxn,lyn),dnO2(lzn,lxn,lyn),dvnrho(lzn,lxn,lyn), &
            dvnz(lzn,lxn,lyn),dvnx(lzn,lxn,lyn),dTn(lzn,lxn,lyn))


!PRINT OUT SOME BASIC INFO ABOUT THE GRID THAT WE'VE LOADED
if (debug) then
  print *, 'Min/max zn,xn,yn values',myid,minval(zn),maxval(zn),minval(xn),maxval(xn),minval(yn),maxval(yn)
  print *, 'Min/max zi,xi,yi values',myid,minval(zi),maxval(zi),minval(xi),maxval(xi),minval(yi),maxval(yi)
  !print *, 'Source lat/long:  ',myid,meanlat,meanlong
  !print *, 'Plasma grid lat range:  ',myid,minval(x%glat(:,:,:)),maxval(x%glat(:,:,:))
  !print *, 'Plasma grid lon range:  ',myid,minval(x%glon(:,:,:)),maxval(x%glon(:,:,:))
end if

end procedure gridproj_dneu3D


subroutine slabrange(maxzn,ximat,yimat,zimat,sourcemlat,xnrange,ynrange)

!takes in a subgrid and the max altitude of interest for neutral interpolation and then computes
!what the maximum xn and yn will be for that slab
! ZZZ - also this is specific to dipole grids right now...

real(wp), intent(in) :: maxzn
real(wp), dimension(:,:,:), intent(in) :: ximat,yimat,zimat
real(wp), intent(in) :: sourcemlat
real(wp), dimension(2), intent(out) :: xnrange,ynrange     !for min and max

real(wp), dimension(:,:,:), allocatable :: xitmp,yitmp,zitmp
integer :: lx1tmp
integer, dimension(size(ximat,2),size(ximat,3)) :: ix1stmp
integer :: ix1tmp
logical :: flagSH
integer :: ix1


!in what hemisphere is our source?
if (sourcemlat<=0) then
  flagSH=.true.
else
  flagSH=.false.
end if


!peel the grid in half (source hemisphere if closed dipole)
if (gridflag==0) then    !closed dipole grid

  ix1=maxloc(pack(zimat(:,1,1),.true.),1)    !apex is by definition the highest altitude along a given field line
  if (flagSH) then
    lx1tmp=ix1                  !first piece of arrays
  else
    lx1tmp=lx1-ix1    !second (end) piece of arrays
  end if
  allocate(xitmp(lx1tmp,lx2,lx3),yitmp(lx1tmp,lx2,lx3),zitmp(lx1tmp,lx2,lx3))   !could this be done more less wastefully with pointers???

  if(flagSH) then    !southern hemisphere
    xitmp=ximat(1:ix1,1:lx2,1:lx3)          !select beginning of the array - the southern half
    yitmp=yimat(1:ix1,1:lx2,1:lx3)
    zitmp=zimat(1:ix1,1:lx2,1:lx3)
  else               !northern hemisphere
      xitmp=ximat(ix1+1:lx1,1:lx2,1:lx3)    !select end half of the array
      yitmp=yimat(ix1+1:lx1,1:lx2,1:lx3)
      zitmp=zimat(ix1+1:lx1,1:lx2,1:lx3)
  end if
else     !this is not an interhemispheric grid so our approach is to just use all of the data
  lx1tmp=lx1
  allocate(xitmp(lx1tmp,lx2,lx3),yitmp(lx1tmp,lx2,lx3),zitmp(lx1tmp,lx2,lx3))   !could this be done more less wastefully with pointers?
  xitmp=ximat(1:lx1,1:lx2,1:lx3)
  yitmp=yimat(1:lx1,1:lx2,1:lx3)
  zitmp=zimat(1:lx1,1:lx2,1:lx3)
!  flagSH=.true.    !treat is as southern, doesn't really matter in this case...
end if


!the min and max x are simply determined by longitude...
xnrange(1)=minval(xitmp)
xnrange(2)=maxval(xitmp)


!situation is more complicated for latitude due to dipole grid, need to determine by L-shell
if (flagSH) then
  ix1=minloc(zitmp(:,1,1)-maxzn,1,zitmp(:,1,1)-maxzn>0._wp)    !find the min distance from maxzn subject to constraint that it is > 0, just use the first longitude slice since they will all have the same L-shell-field line relations
  ynrange(2)=yitmp(ix1,1,1)
  ix1=minloc(zitmp(:,lx2,1),1,zitmp(:,lx2,1)<0._wp)
  ynrange(1)=yitmp(ix1,lx2,1)
else    !things are swapped around in NH
  ix1=minloc(zitmp(:,1,1)-maxzn,1,zitmp(:,1,1)-maxzn>0._wp)    !find the min distance from maxzn subject to constraint that it is > 0; this is the southernmost edge of the neutral slab we need
  ynrange(1)=yitmp(ix1,1,1)
  ix1=minloc(zitmp(:,lx2,1),1,zitmp(:,lx2,1)<0._wp)
  ynrange(2)=yitmp(ix1,lx2,1)
end if

deallocate(xitmp,yitmp,zitmp)

end subroutine slabrange


subroutine range2inds(ranges,zn,xnall,ynall,indices)

!! determine where the slab described by ranges falls within the global neutral grid

real(wp), dimension(6), intent(in) :: ranges
real(wp), dimension(:), intent(in) :: zn,xnall,ynall
integer, dimension(6), intent(out) :: indices

real(wp) :: minzn,maxzn,minxn,maxxn,minyn,maxyn
integer :: ixn,iyn


!for clarity
minzn=ranges(1)
maxzn=ranges(2)
minxn=ranges(3)
maxxn=ranges(4)
minyn=ranges(5)
maxyn=ranges(6)


!always use the full z-range
indices(1)=1
indices(2)=lzn


!x-range
ixn=1
do while (ixn<lxnall .and. xnall(ixn)<minxn)
  ixn=ixn+1
end do
indices(3)=max(ixn-1,1)    !just to be sure go back one index so that we cover the min range, don't let index fall below zero
do while (ixn<lxnall .and. xnall(ixn)<maxxn)
  ixn=ixn+1
end do
indices(4)=ixn


!y-range
iyn=1
do while (iyn<lynall .and. ynall(iyn)<minyn)
  iyn=iyn+1
end do
indices(5)=max(iyn-1,1)    !just to be sure go back one index so that we cover the min range
do while (iyn<lynall .and. ynall(iyn)<maxyn)
  iyn=iyn+1
end do
indices(6)=iyn

print*, '!!!!!!!!!!!!!!!!!'
print*, myid
print*, ranges
print*, indices
print*, lxnall,lynall
print*, xnall(indices(3)),xnall(indices(4))
print*, ynall(indices(5)),ynall(indices(6))
print*, '!!!!!!!!!!!!!!!!!'


!corner cases - range is not at all within the neutral grid...  Manifests as both indices being either 1 or lxi, interpolation should zero these out...

end subroutine range2inds


end submodule proj
