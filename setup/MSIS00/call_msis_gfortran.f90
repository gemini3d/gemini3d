program neuatmos

implicit none
integer :: iyd,sec,mass=48,iz,lz
real :: f107a,f107,ap(7),stl,apday,ap3
real :: d(9),t(2)
real, allocatable :: glat(:),glon(:),alt(:)
character(256) :: infile,outfile

!read in msis inputs
call getarg(1, infile)
!open(unit=42,file=infile)

!time date etc.
!read(42,*) iyd,sec
!read(42,*) f107a,f107
!read(42,*) apday,ap3
!read(42,*) lz

open(unit=42,file=infile,status='old',form='unformatted',access='stream')    !use binary to reduce file size and read times
read(42) iyd
read(42) sec
read(42) f107a
read(42) f107
read(42) apday
read(42) ap3
read(42) lz

write(*,*) iyd,sec,f107a,f107,apday,ap3,lz


!allocate correct amount of memory
allocate(glat(lz),glon(lz),alt(lz))

!lat,lon,alt data
!do iz=1,lz
!  read(42,*) glat(iz),glon(iz),alt(iz)
!end do

read(42) glat,glon,alt

close(42)

 ap(1:7)=apday
 ap(2)=ap3

!switch to mksa units
 call meters(.true.)

!output file
call getarg(2,outfile)
open(unit=43,file=outfile,status='replace',form='unformatted',access='stream')    !use binary to reduce file size and read times

!call to msis routine
do iz=1,lz
  stl=sec/3600.+glon(iz)/15.
  call gtd7(iyd,sec,alt(iz),glat(iz),glon(iz),stl,f107a,f107,ap,mass,d,t)
  write(43) alt(iz),d(1:9),t(2)
!  write(*,'(1f20.4, 9e20.8, 1f20.4)') alt(iz),d(1:9),t(2)
! write(*,1000),alt(iz),d(1:9),t(2)
! 1000 format(1F,9E,1F)
end do

close(43)
deallocate(glat,glon,alt)

end program neuatmos
