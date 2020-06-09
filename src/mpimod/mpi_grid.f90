submodule (mpimod) mpi_gridsub

implicit none (type, external)

contains

module procedure mpi_manualgrid

integer, dimension(2) :: inds

if (lx2all/lid2in*lid2in /= lx2all) error stop 'user input grid split in x2 will not work'
if (lx3all/lid3in*lid3in /= lx3all) error stop 'user input grid split in x3 will not work'
if (lid2in*lid3in /= lid) error stop 'total number of processes not commensurate with x2 and x3 split'

lid2=lid2in
lid3=lid3in

!THIS PROCESS' LOCATION ON THE GRID
inds=ID2grid(myid)
myid2=inds(1)
myid3=inds(2)

end procedure mpi_manualgrid


module procedure mpigrid
!! Automatically determine the PROCESS GRID
!! sets value of lid2,lid3 globally
!! FIXME: should use derived type for these
!! FIXME: improve algorithm to use more CPU cores (be more effective in finding factors for x2 and x3)

integer, dimension(2) :: inds

if (lx3all==1) then
  !! 2D simulation, NOT swapped, divide in x3
  lid3 = min(lid, lx2all)
  lid2 = 1
elseif (lx2all==1) then
  !! 2D simulation, SWAP x2 to x3, divide in x3
  lid3 = min(lid, lx3all)
  lid2 = 1
else
  !! 3D simulation
  !if (modulo(lid,lx3all) /= 0) then
  !   write (stderr, *) 'mpigrid: 3D grid setup: lx2all, lx3all, lid:',lx2all,lx3all,lid
  !   error stop 'Grid is not divisible by number of processes (lx3all/lid*lid /= lx3all).'
  ! end if


  lid = min(lid, lx3all)
  !! more CPUs than lx3all, reduce used MPI images

  do while(modulo(lx3all, lid) /= 0)
    lid = lid-1
  end do
  !! make number of MPI images a factor of lx3all

  lid2=1
  lid3=lid
  do while( ((lid3/2)*2==lid3) .and. (lid3-lid2>lid3 .or. lid3-lid2>lid2) .and. &
            lx3all/(lid3/2)*(lid3/2)==lx3all .and. lx2all/(lid2*2)*(lid2*2)==lx2all .and. &
            lid3/2>1)
  !! ensure that lx3 is divisible by lid3 and lx2 by lid2 and lid3 must be > 1

    lid3=lid3/2
    lid2=lid2*2
  end do
end if


!FORCE THE CODE TO USE 1D PROCESS GRID
!lid2=1; lid3=lid;


!> THIS PROCESS' LOCATION ON THE GRID
inds = ID2grid(myid)
myid2 = inds(1)
myid3 = inds(2)

end procedure mpigrid


module procedure grid2id
!! COMPUTES A PROCESS ID FROM A LOCATION ON THE PROCESS GRID

grid2ID = i3 * lid2 + i2
!! this formula assumes that the first element is (i2,i3)=(0,0)

end procedure grid2id


module procedure id2grid
!! COMPUTES GRID LOCATION FROM A PROCESS ID

ID2grid(2) = ID / lid2
!! x3 index into process grid
ID2grid(1) = ID - ID2grid(2) * lid2
!! x2 index into process grid

end procedure ID2grid

end submodule mpi_gridsub
