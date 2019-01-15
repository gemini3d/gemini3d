!  function curl3D1(f1,f2,f3,dx1,dx2,dx3)
!
!    !------------------------------------------------------------
!    !-------COMPUTE A 3D CURL, X1 COMPONENT.  IT IS EXPECTED THAT 
!    !-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!    !-------THEY ARE PASSED INTO THIS ROUTINE.  DX(I) IS PRESUMED
!    !-------TO BE THE *BACKWARD* DIFFERENCE AT POINT I.
!    !-------
!    !-------THE F1 COMPONENT AND DIFF COULD BE OMITTED.
!    !------------------------------------------------------------
!
!    real(wp), dimension(:,:,:), intent(in) :: f1,f2,f3
!    real(wp), dimension(1:size(f1,1)), intent(in) :: dx1
!    real(wp), dimension(1:size(f1,2)), intent(in) :: dx2
!    real(wp), dimension(1:size(f1,3)), intent(in) :: dx3
!
!    integer :: ix1,ix2,ix3,lx1,lx2,lx3
!
!    real(wp), dimension(1:size(f1,1),1:size(f1,2),1:size(f1,3)) :: curl3D1
!
!    lx1=size(f1,1)
!    lx2=size(f1,2)
!    lx3=size(f1,3)
!
!    if (lx2>1) then
!      do ix3=1,lx3
!        do ix1=1,lx1
!          curl3D1(ix1,1,ix3)=(f3(ix1,2,ix3)-f3(ix1,1,ix3))/dx2(2)   !fwd diff
!          curl3D1(ix1,2:lx2-1,ix3)=(f3(ix1,3:lx2,ix3)-f3(ix1,1:lx2-2,ix3)) &
!                                       /(dx2(3:lx2)+dx2(2:lx2-1))
!          curl3D1(ix1,lx2,ix3)=(f3(ix1,lx2,ix3)-f3(ix1,lx2-1,ix3))/dx2(lx2)    !bwd diff
!        end do
!      end do
!    else
!      curl3D1=0d0
!    end if
!
!    do ix2=1,lx2
!      do ix1=1,lx1
!        curl3D1(ix1,ix2,1)=curl3D1(ix1,ix2,1)-(f2(ix1,ix2,2)-f2(ix1,ix2,1))/dx3(2)
!        curl3D1(ix1,ix2,2:lx3-1)=curl3D1(ix1,ix2,2:lx3-1)- &
!                                 (f2(ix1,ix2,3:lx3)-f2(ix1,ix2,1:lx3-2)) &
!                                     /(dx3(3:lx3)+dx3(2:lx3-1))
!        curl3D1(ix1,ix2,lx3)=curl3D1(ix1,ix2,lx3)-(f2(ix1,ix2,lx3)-f2(ix1,ix2,lx3-1))/dx3(lx3)
!      end do
!    end do
!
!  end function curl3D1
!
!
!  function curl3D2(f1,f2,f3,dx1,dx2,dx3)
!
!    !------------------------------------------------------------
!    !-------COMPUTE A 3D CURL, X2 COMPONENT.  IT IS EXPECTED THAT 
!    !-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!    !-------THEY ARE PASSED INTO THIS ROUTINE.  DX(I) IS PRESUMED
!    !-------TO BE THE *BACKWARD* DIFFERENCE AT POINT I
!    !-------
!    !-------THE F2 COMPONENT AND DIFF COULD BE OMITTED.
!    !-------ARE THE SIGNS RIGHT?
!    !------------------------------------------------------------
!
!    real(wp), dimension(:,:,:), intent(in) :: f1,f2,f3
!    real(wp), dimension(1:size(f1,1)), intent(in) :: dx1
!    real(wp), dimension(1:size(f1,2)), intent(in) :: dx2
!    real(wp), dimension(1:size(f1,3)), intent(in) :: dx3
!
!    integer :: ix1,ix2,ix3,lx1,lx2,lx3
!
!    real(wp), dimension(1:size(f1,1),1:size(f1,2),1:size(f1,3)) :: curl3D2
!
!    lx1=size(f1,1)
!    lx2=size(f1,2)
!    lx3=size(f1,3)
!
!    do ix2=1,lx2
!      do ix1=1,lx1
!        curl3D2(1,ix2,ix3)=(f3(2,ix2,ix3)-f3(1,ix2,ix3))/dx1(2)   !fwd diff
!        curl3D2(2:lx1-1,ix2,ix3)=(f3(3:lx1,ix2,ix3)-f3(1:lx1-2,ix2,ix3)) &
!                                     /(dx1(3:lx1)+dx1(2:lx1-1))
!        curl3D2(lx1,ix2,ix3)=(f3(lx1,ix2,ix3)-f3(lx1-1,ix2,ix3))/dx1(lx1)    !bwd diff
!      end do
!    end do
!
!    do ix2=1,lx2
!      do ix1=1,lx1
!        curl3D2(ix1,ix2,1)=curl3D2(ix1,ix2,1)-(f1(ix1,ix2,2)-f1(ix1,ix2,1))/dx3(2)
!        curl3D2(ix1,ix2,2:lx3-1)=curl3D2(ix1,ix2,2:lx3-1)- &
!                                 (f1(ix1,ix2,3:lx3)-f1(ix1,ix2,1:lx3-2)) &
!                                     /(dx3(3:lx3)+dx3(2:lx3-1))
!        curl3D2(ix1,ix2,lx3)=curl3D2(ix1,ix2,lx3)-(f1(ix1,ix2,lx3)-f1(ix1,ix2,lx3-1))/dx3(lx3)
!      end do
!    end do
!
!  end function curl3D2
!
!
!  function curl3D3(f1,f2,f3,dx1,dx2,dx3)
!
!    !------------------------------------------------------------
!    !-------COMPUTE A 3D CURL, X3 COMPONENT.  IT IS EXPECTED THAT 
!    !-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!    !-------THEY ARE PASSED INTO THIS ROUTINE.  DX(I) IS PRESUMED
!    !-------TO BE THE *BACKWARD* DIFFERENCE AT POINT I
!    !-------
!    !-------THE F3 COMPONENT AND DIFF COULD BE OMITTED.
!    !------------------------------------------------------------
!
!    real(wp), dimension(:,:,:), intent(in) :: f1,f2,f3
!    real(wp), dimension(1:size(f1,1)), intent(in) :: dx1
!    real(wp), dimension(1:size(f1,2)), intent(in) :: dx2
!    real(wp), dimension(1:size(f1,3)), intent(in) :: dx3
!
!    integer :: ix1,ix2,ix3,lx1,lx2,lx3
!
!    real(wp), dimension(1:size(f1,1),1:size(f1,2),1:size(f1,3)) :: curl3D3
!
!    lx1=size(f1,1)
!    lx2=size(f1,2)
!    lx3=size(f1,3)
!
!    do ix2=1,lx2
!      do ix1=1,lx1
!        curl3D3(1,ix2,ix3)=(f2(2,ix2,ix3)-f2(1,ix2,ix3))/dx1(2)   !fwd diff
!        curl3D3(2:lx1-1,ix2,ix3)=(f2(3:lx1,ix2,ix3)-f2(1:lx1-2,ix2,ix3)) &
!                                     /(dx1(3:lx1)+dx1(2:lx1-1))
!        curl3D3(lx1,ix2,ix3)=(f2(lx1,ix2,ix3)-f2(lx1-1,ix2,ix3))/dx1(lx1)    !bwd diff
!      end do
!    end do
!
!    if (lx2>1) then
!      do ix3=1,lx3
!        do ix1=1,lx1
!          curl3D3(ix1,1,ix3)=curl3D3(ix1,1,ix3)-(f1(ix1,2,ix3)-f1(ix1,1,ix3))/dx2(2)   !fwd diff
!          curl3D3(ix1,2:lx2-1,ix3)=curl3D3(ix1,2:lx2-1,ix3)- &
!                                    (f1(ix1,3:lx2,ix3)-f1(ix1,1:lx2-2,ix3)) &
!                                       /(dx2(3:lx2)+dx2(2:lx2-1))
!          curl3D3(ix1,lx2,ix3)=curl3D3(ix1,lx2,ix3)-(f1(ix1,lx2,ix3)-f1(ix1,lx2-1,ix3))/dx2(lx2)    !bwd diff
!        end do
!      end do
!    end if
!
!  end function curl3D3
