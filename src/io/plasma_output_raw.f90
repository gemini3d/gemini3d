submodule (io:plasma_output) plasma_output_raw

use timeutils, only : date_filename

implicit none (type, external)

contains

module procedure output_root_stream_mpi_raw

!! COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE USING STREAM I/O.
!! STATE VARS ARE EXPECTED INCLUDE GHOST CELLS

integer :: isp
character(:), allocatable :: filenamefull
integer(8) :: recordlength   !can be 8 byte with compiler flag -frecord-marker=8

!FIGURE OUT THE FILENAME
filenamefull=date_filename(outdir,ymd,UTsec) // '.dat'
print *, 'Output file name:  ',filenamefull
! call logger(filenamefull,'filename.log')
! call logger(UTsec, 'UTsec.log')


!SOME DEBUG OUTPUT ON FILE SIZE
recordlength=int(8,8)+int(8,8)*int(3,8)*int(lx1,8)*int(lx2all,8)*int(lx3all,8)*int(lsp,8)+ &
             int(8,8)*int(5,8)*int(lx1,8)*int(lx2all,8)*int(lx3all,8)+ &
             int(8,8)*int(lx2,8)*int(lx3all,8)
print *, 'Output bit length:  ',recordlength,lx1,lx2all,lx3all,lsp

!WRITE THE DATA
block
integer :: u
open(newunit=u,file=filenamefull,status='replace',form='unformatted',access='stream',action='write')    !has no problem with > 2GB output files
write(u) real(ymd,wp),UTsec/3600._wp    !no matter what we must output date and time

select case (flagoutput)
  case (2)    !output ISR-like average parameters
    write(u) &
      neall(1:lx1,1:lx2all,1:lx3all), &
      v1avgall(1:lx1,1:lx2all,1:lx3all), &    !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
      Tavgall(1:lx1,1:lx2all,1:lx3all),&
      Teall(1:lx1,1:lx2all,1:lx3all),&
      J1all(1:lx1,1:lx2all,1:lx3all), &
      J2all(1:lx1,1:lx2all,1:lx3all), &
      J3all(1:lx1,1:lx2all,1:lx3all),&
      v2avgall(1:lx1,1:lx2all,1:lx3all),&
      v3avgall(1:lx1,1:lx2all,1:lx3all)
  case (3)     !just electron density
    print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'
    write(u) neall(1:lx1,1:lx2all,1:lx3all)
  case default    !output everything
    print *, '!!!NOTE:  Input file has selected full output, large files may result!'
    write(u) &
      nsall(1:lx1,1:lx2all,1:lx3all,:),&
      vs1all(1:lx1,1:lx2all,1:lx3all,:), &    !this is full output of all parameters in 3D
      Tsall(1:lx1,1:lx2all,1:lx3all,:),&
      J1all(1:lx1,1:lx2all,1:lx3all),&
      J2all(1:lx1,1:lx2all,1:lx3all), &
      J3all(1:lx1,1:lx2all,1:lx3all),&
      v2avgall(1:lx1,1:lx2all,1:lx3all),&
      v3avgall(1:lx1,1:lx2all,1:lx3all)
end select

if (gridflag==1) then
  print *, 'Writing topside boundary conditions for inverted-type grid...'
  write(u)  Phiall(1,1:lx2all,1:lx3all)
else
  print *, 'Writing topside boundary conditions for non-inverted-type grid...'
  write(u)  Phiall(lx1,1:lx2all,1:lx3all)
end if

close(u)
end block

end procedure output_root_stream_mpi_raw


end submodule plasma_output_raw
