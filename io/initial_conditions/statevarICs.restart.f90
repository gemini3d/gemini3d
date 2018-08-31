  subroutine statevarICs(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    implicit none
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp
    integer :: lx1in,lx2in,lx3in
    integer :: ix1in,ix2in,ix3in


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lx1in=size(nsin,1)
    lx2in=size(nsin,2)
    lx3in=size(nsin,3)
    lsp=size(nsall,4)


    do isp=1,lsp
      do ix3=1,lx3all
        do ix2=1,lx2
          do ix1=1,lx1
            !DETERMINE INDICES IN INPUT DATA TO BE USED FOR EACH GRID POINT
            if (ix3 <= lx3in) then
              ix3in=ix3
            else
              ix3in=lx3in   !we are outside domain of input file, just copy last point
            end if

            if (ix2 <= lx2in) then
              ix2in=ix2
            else
              ix2in=lx2in
            end if

            if (ix1 <= lx1in) then
              ix1in=ix1
            else
              ix1in=lx1in
            end if

            nsall(ix1,ix2,ix3,isp)=nsin(ix1in,ix2in,ix3in,isp)
            vs1all(ix1,ix2,ix3,isp)=vs1in(ix1in,ix2in,ix3in,isp)
            Tsall(ix1,ix2,ix3,isp)=Tsin(ix1in,ix2in,ix3in,isp)
          end do
        end do
      end do
    end do

  end subroutine statevarICs

