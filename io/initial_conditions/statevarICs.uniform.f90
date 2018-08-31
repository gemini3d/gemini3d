  subroutine statevarICs(x1,x2,x3all,nsin,vs1in,Tsin,nsall,vs1all,Tsall)

    implicit none
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:), intent(in) :: x3all
    real(8), dimension(:,:,:,:), intent(in) :: nsin,vs1in,Tsin    !for loading equilibrium results
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: nsall,vs1all,Tsall

    real(8) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve
    integer :: ix1,ix2,ix3,isp,lx1,lx2,lx3all,lsp


    !SYSTEM SIZES
    lx1=size(nsall,1)-4
    lx2=size(nsall,2)-4
    lx3all=size(nsall,3)-4
    lsp=size(nsall,4)

      do isp=1,lsp
        nsall(1:lx1,1:lx2,1:lx3all,isp)=spread(spread(nsin(:,1,1,isp),2,lx2),3,lx3all)
      end do
      do isp=1,lsp
        vs1all(1:lx1,1:lx2,1:lx3all,isp)=spread(spread(vs1in(:,1,1,isp),2,lx2),3,lx3all)
      end do
      do isp=1,lsp
        Tsall(1:lx1,1:lx2,1:lx3all,isp)=spread(spread(Tsin(:,1,1,isp),2,lx2),3,lx3all)
      end do

  end subroutine statevarICs

