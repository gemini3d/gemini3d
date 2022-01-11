submodule(io:io_aurora) io_aurora_raw

implicit none (type, external)

contains

module procedure output_aur_root_raw
  !! COLLECT COMPLETE DATA FROM WORKERS AND PROCESS FOR OUTPUT.
  !! NO GHOST CELLS (I HOPE)

  !real(wp), dimension(1:lx2,1:lwave,1:lx3) :: ivertmp
  !real(wp), dimension(1:lx2all,1:lwave,1:lx3all) :: iverall
  real(wp), dimension(1:lwave,1:lx2,1:lx3) :: ivertmp
  real(wp), dimension(1:lwave,1:lx2all,1:lx3all) :: iverall

  real(wp), dimension(1:lx2,1:lx3) :: emistmp                !single emission subgrid
  real(wp), dimension(1:lx2all,1:lx3all) :: emisall          !single emission total grid
  real(wp), dimension(1:lx2all,1:lx3all,1:lwave) :: iverout  !output array in the order scripts expect

  integer :: iwave

  !!ivertmp=reshape(iver,[lx2,lwave,lx3],order=[1,3,2])
  !ivertmp=reshape(iver,[lwave,lx2,lx3],order=[3,1,2])
  !call gather_recv(ivertmp,tag%Aur,iverall)
  do iwave=1,lwave
    emistmp=iver(:,:,iwave)
    call gather_recv(emistmp,tag%Aur,emisall)
    iverout(:,:,iwave)=emisall
  end do


  !FORM THE INPUT FILE NAME

  print *, '  Output file name (auroral maps):  ',filename
  block
    integer :: u
    open(newunit=u,file=filename,status='replace',form='unformatted',access='stream',action='write')
    !  write(u) reshape(iverall,[lx2all,lwave,lx3all],order=[2,1,3])
    !  write(u) reshape(iverall,[lx2all,lx3all,lwave],order=[2,3,1])
    write(u) iverout
    close(u)
  end block
end procedure output_aur_root_raw

end submodule io_aurora_raw
