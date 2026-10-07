program igrf_reference
implicit none
real(8) :: year,radius,colat,lon,x,y,z,f
integer :: ios
external :: igrf14syn
do
  read(*,*,iostat=ios) year,radius,colat,lon
  if(ios/=0)exit
  call igrf14syn(0,year,2,radius,colat,lon,x,y,z,f)
  write(*,'(4es26.17)') x,y,z,f
enddo
end program
