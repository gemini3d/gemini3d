submodule (neutraldata3Dobj_mpi) neuslab_mpi

!> this submodule contains utility procedures specifically for computing overlaps between GEMINI and MAGIC grids
!    Because this involving mpi splitting of input data we also have e.g. various custom message passing procedures
!    use to manipulate neutral object data stashed here.

contains
!  !takes in a subgrid and the max altitude of interest for neutral interpolation and then computes
!  !what the maximum xn and yn will be for that slab
!  ! ZZZ - also this is specific to dipole grids right now...
!  module procedure slabrange
!    real(wp), dimension(:,:,:), allocatable :: xitmp,yitmp,zitmp
!    integer :: lx1tmp
!    logical :: flagSH
!    integer :: ix1,ix2,ix3
!    integer :: lx1,lx2,lx3
!    integer,dimension(3) :: ixs
!    integer,dimension(2) :: ixs13
!    real(wp), dimension(:,:), allocatable :: zitmpslice
!
!
!    ! compute sizes from input arrays
!    lx1=size(zimat,1); lx2=size(zimat,2); lx3=size(zimat,3);
!
!    !in what hemisphere is our source?
!    if (sourcemlat<=0) then
!      flagSH=.true.
!    else
!      flagSH=.false.
!    end if
!
!    !peel the grid in half (source hemisphere if closed dipole)
!    if (gridflag==0) then    !closed dipole grid
!
!      ixs = maxloc(zimat(:,:,:))    !apex is by definition the highest altitude along a given field line
!      ix1=ixs(1); ix2=ixs(2); ix3=ixs(3);
!
!      if (flagSH) then
!        lx1tmp=ix1                  !first piece of arrays
!      else
!        lx1tmp=lx1-ix1    !second (end) piece of arrays
!      end if
!      allocate(xitmp(lx1tmp,lx2,lx3), &
!               yitmp(lx1tmp,lx2,lx3), &
!               zitmp(lx1tmp,lx2,lx3))
!      !! could this be done more less wastefully with pointers???
!
!      if(flagSH) then    !southern hemisphere
!        xitmp=ximat(1:ix1,1:lx2,1:lx3)          !select beginning of the array - the southern half
!        yitmp=yimat(1:ix1,1:lx2,1:lx3)
!        zitmp=zimat(1:ix1,1:lx2,1:lx3)
!      else               !northern hemisphere
!          xitmp=ximat(ix1+1:lx1,1:lx2,1:lx3)    !select end half of the array
!          yitmp=yimat(ix1+1:lx1,1:lx2,1:lx3)
!          zitmp=zimat(ix1+1:lx1,1:lx2,1:lx3)
!      end if
!    else     !this is not an interhemispheric grid so our approach is to just use all of the data
!      lx1tmp=lx1
!      allocate(xitmp(lx1tmp,lx2,lx3), &
!               yitmp(lx1tmp,lx2,lx3), &
!               zitmp(lx1tmp,lx2,lx3))
!      !! could this be done more less wastefully with pointers?
!      xitmp=ximat(1:lx1,1:lx2,1:lx3)
!      yitmp=yimat(1:lx1,1:lx2,1:lx3)
!      zitmp=zimat(1:lx1,1:lx2,1:lx3)
!    !  flagSH=.true.    !treat is as southern, doesn't really matter in this case...
!    end if
!
!
!    !the min and max x are simply determined by min/max over all longitudes in the grid.
!    xnrange(1) = minval(xitmp)
!    xnrange(2) = maxval(xitmp)
!
!
!    !situation is more complicated for latitude due to dipole grid, need to determine by L-shell
!    allocate(zitmpslice(lx1tmp,lx3))
!    if (flagSH) then    ! in the southern hemisphere the min Lshell is northernmost latitude
!      zitmpslice=zitmp(:,1,:)
!      if (any(zitmpslice - maxzn > 0)) then
!        ixs13 = minloc(zitmpslice-maxzn, mask=zitmpslice - maxzn > 0)
!        ix1=ixs13(1); ix3=ixs13(2);
!        !! find the min distance from maxzn subject to constraint that it is > 0,
!      else
!        ix1 = lx1
!      end if
!      ynrange(2) = yitmp(ix1,1,ix3)
!
!      zitmpslice=zitmp(:,lx2,:)
!      if (any(zitmpslice < 0)) then
!        ixs13 = minloc(zitmpslice, mask=zitmpslice < 0)
!        ix1=ixs13(1); ix3=ixs13(2);
!      else
!        ix1 = 1
!      end if
!      !ix1=max(ix1,1)
!      ynrange(1)=yitmp(ix1,lx2,ix3)
!    else    !things are swapped around in NH, min L-shell is lowest latitude
!            !  southern edge is dictated by high altitudes in NH; this is ix1=1 in dipole indices
!      zitmpslice=zitmp(:,1,:)          
!      if (any(zitmpslice - maxzn > 0)) then
!        ixs13 = minloc(zitmpslice-maxzn, mask=zitmpslice - maxzn > 0)
!        ix1=ixs13(1); ix3=ixs13(2);
!        ! find the min distance from maxzn subject to constraint that it is > 0; this is the southernmost edge of the neutral slab we need
!      else
!        ix1=1    ! default to first grid point
!      end if
!      ynrange(1)=yitmp(ix1,1,ix3)
!
!      !! an issue here is that the behavior in the case that the mask condition it not met is not well-defined so
!      !!    we really need to check this separately and have the code do something sensible in this case.  I.e. if there is no
!      !!    zero crossing then we just need to use the entire array.
!      zitmpslice=zitmp(:,lx2,:)
!      if (any(zitmp(:,lx2,:) < 0)) then
!        ixs13 = minloc(zitmpslice, mask=zitmpslice < 0)
!        ix1=ixs13(1); ix3=ixs13(2);
!
!        ! northernmost edge is defined by the zero crossing (if any)
!      else
!        ix1=size(yitmp,1)     ! default in this case to last grid point
!      end if
!      ynrange(2)=yitmp(ix1,lx2,ix3)
!    end if
!
!    deallocate(xitmp,yitmp,zitmp,zitmpslice)
!  end procedure slabrange


  ! Attempt at a more general approach to try to detect overlap in grids
  ! ZZZ - also this is specific to dipole grids right now...
  module procedure slabrange
    real(wp), dimension(:,:,:), allocatable :: xitmp,yitmp,zitmp
    integer :: lx1tmp
    logical :: flagSH
    integer :: ix1,ix2,ix3
    integer, dimension(3) :: ixs
    integer :: lx1,lx2,lx3
    real(wp), dimension(4) :: xrange0,xrangeref,yrange0,yrangeref

    ! compute sizes from input arrays
    lx1=size(zimat,1); lx2=size(zimat,2); lx3=size(zimat,3);

    !in what hemisphere is our source?
    if (sourcemlat<=0) then
      flagSH=.true.
    else
      flagSH=.false.
    end if

    !peel the grid in half (source hemisphere if closed dipole)
    if (gridflag==0) then    !closed dipole grid
      ixs = maxloc(zimat(:,:,:))    !apex is by definition the highest altitude along a given field line
      ix1=ixs(1); ix2=ixs(2); ix3=ixs(3);

      if (flagSH) then
        lx1tmp=ix1                  !first piece of arrays
      else
        lx1tmp=lx1-ix1    !second (end) piece of arrays
      end if
      allocate(xitmp(lx1tmp,lx2,lx3), &
               yitmp(lx1tmp,lx2,lx3), &
               zitmp(lx1tmp,lx2,lx3))
      !! could this be done more less wastefully with pointers???

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
      allocate(xitmp(lx1tmp,lx2,lx3), &
               yitmp(lx1tmp,lx2,lx3), &
               zitmp(lx1tmp,lx2,lx3))
      !! could this be done more less wastefully with pointers?
      xitmp=ximat(1:lx1,1:lx2,1:lx3)
      yitmp=yimat(1:lx1,1:lx2,1:lx3)
      zitmp=zimat(1:lx1,1:lx2,1:lx3)
    !  flagSH=.true.    !treat is as southern, doesn't really matter in this case...
    end if


    !! Now perform the search for grid overlap in the appropriate hemisphere
    ! Find the 8 corners of the dipole mesh which fully encapsulate the neutral mesh and then
    ! define a "convex hull" of sorts for neutrals ranges that we need.
    call find_corners(xitmp,yitmp,zitmp,maxzn,flagSH,xrange0,xrangeref,yrange0,yrangeref)
    xnrange(1)=min(minval(xrange0),minval(xrangeref))
    xnrange(2)=max(maxval(xrange0),maxval(xrangeref))
    ynrange(1)=min(minval(yrange0),minval(yrangeref))
    ynrange(2)=max(maxval(yrange0),maxval(yrangeref))

    deallocate(xitmp,yitmp,zitmp)
  end procedure slabrange


  !> Find the y/lat. coordinates of eight "corners" of a dipole mesh intersecting a reference altitude and the ground
  subroutine find_corners(xi,yi,zi,zref,flagSH,xrange0,xrangeref,yrange0,yrangeref)
    real(wp), dimension(:,:,:), intent(in) :: xi,yi,zi
    real(wp), intent(in) :: zref         ! altitude of projection to indentify
    logical, intent(in) :: flagSH        ! source in southern hemisphere?
    real(wp), dimension(4), intent(inout) :: xrange0,xrangeref,yrange0,yrangeref
    integer :: ix1,ix2,ix3
    integer :: lx1,lx2,lx3
    integer, dimension(4) :: ix2s
    integer, dimension(4) :: ix3s
    integer :: icorner
    real(wp), dimension(1:size(zi,1)) :: zislice

    ! extents based on inputs
    lx1=size(xi,1); lx2=size(xi,2); lx3=size(xi,3);
    ix2s=[1,lx2,1,lx2]
    ix3s=[1,1,lx3,lx3]

    !! First check intersection with reference altitude specified in input arguments
    do icorner=1,4
      ix2=ix2s(icorner)
      ix3=ix3s(icorner)
      zislice=zi(:,ix2,ix3)
      if (any(zislice - zref > 0)) then
        ix1 = minloc(zislice-zref, dim=1, mask=zislice-zref > 0)
        ! find the min distance from maxzn subject to constraint that it is > 0 (so it *fully encloses* the neutral domain)
      else           ! default to highest possible altitude
        if (flagSH) then
          ix1=lx1
        else
          ix1=1
        end if
      end if
      xrangeref(icorner)=xi(ix1,ix2,ix3)
      yrangeref(icorner)=yi(ix1,ix2,ix3)
    end do

    !! Now we must check intersection with ground...
    !! an issue here is that the behavior in the case that the mask condition it not met is not well-defined so
    !!    we really need to check this separately and have the code do something sensible in this case.  I.e. if there is no
    !!    zero crossing then we just need to use the entire array.
    do icorner=1,4
      ix2=ix2s(icorner)
      ix3=ix3s(icorner)
      zislice=zi(:,ix2,ix3)
      if (any(zislice < 0)) then
        ix1 = maxloc(zislice, dim=1, mask=zislice < 0)
        ! northernmost edge is defined by the zero crossing (if any)
      else
        ! default to lowest possible altitude
        if (flagSH) then
          ix1=1
        else 
          ix1=lx1
        end if
      end if
      xrange0(icorner)=xi(ix1,ix2,ix3) 
      yrange0(icorner)=yi(ix1,ix2,ix3) 
    end do
  end subroutine find_corners


!    ! ZZZ - also this is specific to dipole grids right now...
!  subroutine slabrange_orig
!    real(wp), dimension(:,:,:), allocatable :: xitmp,yitmp,zitmp
!    integer :: lx1tmp
!    logical :: flagSH
!    integer :: ix1
!    integer :: lx1,lx2,lx3
!
!
!    ! compute sizes from input arrays
!    lx1=size(zimat,1); lx2=size(zimat,2); lx3=size(zimat,3);
!
!    !in what hemisphere is our source?
!    if (sourcemlat<=0) then
!      flagSH=.true.
!    else
!      flagSH=.false.
!    end if
!
!    !peel the grid in half (source hemisphere if closed dipole)
!    if (gridflag==0) then    !closed dipole grid
!
!      ix1 = maxloc(pack(zimat(:,1,1),.true.), dim=1)    !apex is by definition the highest altitude along a given field line
!      if (flagSH) then
!        lx1tmp=ix1                  !first piece of arrays
!      else
!        lx1tmp=lx1-ix1    !second (end) piece of arrays
!      end if
!      allocate(xitmp(lx1tmp,lx2,lx3), &
!               yitmp(lx1tmp,lx2,lx3), &
!               zitmp(lx1tmp,lx2,lx3))
!      !! could this be done more less wastefully with pointers???
!
!      if(flagSH) then    !southern hemisphere
!        xitmp=ximat(1:ix1,1:lx2,1:lx3)          !select beginning of the array - the southern half
!        yitmp=yimat(1:ix1,1:lx2,1:lx3)
!        zitmp=zimat(1:ix1,1:lx2,1:lx3)
!      else               !northern hemisphere
!          xitmp=ximat(ix1+1:lx1,1:lx2,1:lx3)    !select end half of the array
!          yitmp=yimat(ix1+1:lx1,1:lx2,1:lx3)
!          zitmp=zimat(ix1+1:lx1,1:lx2,1:lx3)
!      end if
!    else     !this is not an interhemispheric grid so our approach is to just use all of the data
!      lx1tmp=lx1
!      allocate(xitmp(lx1tmp,lx2,lx3), &
!               yitmp(lx1tmp,lx2,lx3), &
!               zitmp(lx1tmp,lx2,lx3))
!      !! could this be done more less wastefully with pointers?
!      xitmp=ximat(1:lx1,1:lx2,1:lx3)
!      yitmp=yimat(1:lx1,1:lx2,1:lx3)
!      zitmp=zimat(1:lx1,1:lx2,1:lx3)
!    !  flagSH=.true.    !treat is as southern, doesn't really matter in this case...
!    end if
!
!    !the min and max x are simply determined by longitude...
!    xnrange(1) = minval(xitmp)
!    xnrange(2) = maxval(xitmp)
!
!
!    !situation is more complicated for latitude due to dipole grid, need to determine by L-shell
!    if (flagSH) then
!      if (any(zitmp(:,1,1) - maxzn > 0)) then
!        ix1 = minloc(zitmp(:,1,1)-maxzn, dim=1, mask=zitmp(:,1,1) - maxzn > 0)
!      !! find the min distance from maxzn subject to constraint that it is > 0,
!      !! just use the first longitude slice since they will all have the same L-shell-field line relations
!      else
!        ix1 = lx1
!      end if
!      ynrange(2) = yitmp(ix1,1,1)
!      if (any(zitmp(:,lx2,1) < 0)) then
!        ix1 = minloc(zitmp(:,lx2,1), dim=1, mask=zitmp(:,lx2,1) < 0)
!      else
!        ix1 = 1
!      end if
!      !ix1=max(ix1,1)
!      ynrange(1)=yitmp(ix1,lx2,1)
!    else    !things are swapped around in NH
!      if (any(zitmp(:,1,1) - maxzn > 0)) then
!        ix1 = minloc(zitmp(:,1,1)-maxzn, dim=1, mask=zitmp(:,1,1) - maxzn > 0)
!        ! find the min distance from maxzn subject to constraint that it is > 0; this is the southernmost edge of the neutral slab we need
!      else
!        ix1=1    ! default to first grid point
!      end if
!      ynrange(1)=yitmp(ix1,1,1)
!      !! an issue here is that the behavior in the case that the mask condition it not met is not well-defined so
!      !!    we really need to check this separately and have the code do something sensible in this case.  I.e. if there is no
!      !!    zero crossing then we just need to use the entire array.
!      if (any(zitmp(:,lx2,1) < 0)) then
!        ix1 = minloc(zitmp(:,lx2,1), dim=1, mask=zitmp(:,lx2,1) < 0)
!        ! northernmost edge is defined by the zero crossing (if any)
!      else
!        ix1=size(yitmp,1)     ! default in this case to last grid point
!      end if
!      ynrange(2)=yitmp(ix1,lx2,1)
!    end if
!
!    deallocate(xitmp,yitmp,zitmp)
!  end subroutine slabrange_orig


  !> determine where the slab described by ranges falls within the global neutral grid
  module procedure range2inds
    real(wp) :: minzn,maxzn,minxn,maxxn,minyn,maxyn
    integer :: ixn,iyn
    integer :: lzn,lxnall,lynall

    ! pick off sizes for later use
    lzn=size(zn,1); lxnall=size(xnall,1); lynall=size(ynall,1);

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
    print*, mpi_cfg%myid
    print*, ranges
    print*, indices
    print*, lxnall,lynall
    print*, xnall(indices(3)),xnall(indices(4))
    print*, ynall(indices(5)),ynall(indices(6))
    print*, '!!!!!!!!!!!!!!!!!'

    !! corner cases - range is not at all within the neutral grid...
    !! Manifests as both indices being either 1 or lxi, interpolation should zero these out...
  end procedure range2inds


  !> transfer single state parameter data from root to workers (viz. "broadcast")
  module procedure dneu_root2workers
    integer :: iid
    real(wp), dimension(:,:,:), allocatable :: parmtmp
    integer :: lzn

    lzn=size(paramall,1)

    do iid=1,mpi_cfg%lid-1
      allocate(parmtmp(lzn,slabsizes(iid,1),slabsizes(iid,2)))    !get space for the parameters for this worker

      parmtmp=paramall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
      call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag,MPI_COMM_WORLD)

      deallocate(parmtmp)
    end do
    param=paramall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  end procedure dneu_root2workers


  !> get a chunk of neutral data from root
  module procedure dneu_workers_from_root
    integer :: lzn,lxn,lyn

    lzn=size(param,1); lxn=size(param,2); lyn=size(param,3);
    call mpi_recv(param,lzn*lxn*lyn,mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
  end procedure dneu_workers_from_root
end submodule neuslab_mpi
