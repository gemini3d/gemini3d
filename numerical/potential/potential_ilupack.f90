module potential_ilupack

use calculus
implicit none

contains


  function elliptic3D(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dx1,dx1i,dx2,dx2i,dx3,dx3i,Phi0,flagdirich)

    !------------------------------------------------------------
    !-------SOLVE POTENTIAL EQUATION IN 3D USING ILUPACK.  THIS PARTICULAR
    !-------IMPLEMENTATION IS A PURELY STATIC SOLUTION.
    !------------------------------------------------------------

    real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP,sigH,Phi0
    real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
    real(wp), dimension(0:), intent(in) :: dx1
    real(wp), dimension(:), intent(in) :: dx1i
    real(wp), dimension(0:), intent(in) :: dx2
    real(wp), dimension(:), intent(in) :: dx2i
    real(wp), dimension(0:), intent(in) :: dx3
    real(wp), dimension(:), intent(in) :: dx3i
    integer(8), intent(in) :: flagdirich

    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigP2,gradsigP3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigH2,gradsigH3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsig01
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: Ac,Bc,Cc,Dc,Ec,Fc

    !EXTERNAL ILUPACK ROUTINES
    integer, external :: dgnlamgfactor,dgnlamgsolver,dgnlamgnnz
    external dgnlamginit,dgnlamgsol,dgnlamgdelete,dgnlamginfo

    !ILUPACK external parameters
    integer :: mem
    integer :: matching,maxit,lfil,lfilS,nrestart,ierr,mixedprecision
    integer, dimension(:), allocatable :: ind
    character(20) :: ordering
    real(wp) :: droptol,droptolS,condest,restol,elbow

    !Variables that cover the and pass the C-pointers
    integer(8) :: param,PREC

    !LOCAL STORAGE
    integer :: ix1,ix2,ix3,lx1,lx2,lx3
    integer :: ik,ient,lk,lent
    integer, dimension(:), allocatable :: ir_ptr
    integer, dimension(:), allocatable :: ic
    real(wp), dimension(:), allocatable :: M
    real(wp), dimension(:), allocatable :: b
    integer, dimension(:), allocatable :: ir_ptr_copy
    integer, dimension(:), allocatable :: ic_copy
    real(wp), dimension(:), allocatable :: M_copy
    real(wp), dimension(:), allocatable :: b_copy
    real(wp), dimension(:), allocatable :: sol

    real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D


    !SYSTEM SIZES AND ALLOCATION FOR COMPRESSED ROW STORAGE
    write(*,*) 'Allocating space for STATIC elliptic equation solution...'
    lx1=size(srcterm,1)
    lx2=size(srcterm,2)
    lx3=size(srcterm,3)
    lk=lx1*lx2*lx3
!    write(*,*) 'System size:  ',lx1,lx2,lx3

    lent=7*(lx1-2)*(lx2-2)*(lx3-2)                                                 !interior entries
    lent=lent+2*(lx1-2)*(lx2-2)+2*(lx2-2)*(lx3-2)+2*(lx1-2)*(lx3-2)                !6 faces of cube
    lent=lent+4*(lx1-2)+4*(lx2-2)+4*(lx3-2)                                        !12 edges
    lent=lent+8                                                                    !8 corners, now we have total nonzero entries
    lent=lent+lx2*lx3                                                              !entries to deal with Neumann conditions on bottom
    if (flagdirich==0) then
      lent=lent+lx2*lx3
    end if

    write(*,*) 'Total unknowns and nonzero entries:  ',lk,lent
!    allocate(ir_ptr(lk+1),ic(lent),M(lent),b(lent))   !should be b(lk)?
!    allocate(ir_ptr_copy(lk+1),ic_copy(lent),M_copy(lent),b_copy(lent))
    allocate(ir_ptr(lk+1),ic(lent),M(lent),b(lk))
    allocate(ir_ptr_copy(lk+1),ic_copy(lent),M_copy(lent),b_copy(lk))
    allocate(ind(lk))      !not sure what this is used for
!    allocate(sol(lent))    !should be sol(lk) but extra entries okay???  Or are they needed for ilupack temp storage
    allocate(sol(lk))


    !INITIALIZATION
    M(1:lent)=0d0
    b(1:lk)=pack(srcterm,.true.)


    !COMPUTE AUXILIARY COEFFICIENTS
    write(*,*) 'Prepping coefficients for elliptic equation...'    
    gradsig01=grad3D1(sig0,dx1(1:lx1))
    gradsigP2=grad3D2(sigP,dx2(1:lx2)); gradsigP3=grad3D3(sigP,dx3(1:lx3));
    gradsigH2=grad3D2(sigH,dx2(1:lx2)); gradsigH3=grad3D3(sigH,dx3(1:lx3));

    Ac=sigP; Bc=sigP; Cc=sig0; 
    Dc=gradsigP2+gradsigH3
    Ec=gradsigP3-gradsigH2
    Fc=gradsig01


    !LOAD UP MATRIX ELEMENTS
    write(*,*) 'Building CRS sparse matrix...'
    ient=1
    ir_ptr(1)=1
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          ik=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix2,ix3) as a column vector.  Also row no. of big matrix

          if (ix1==1) then    !Check whether we are a bondary point and, if so, enforce boundary conditions
            ic(ient)=ik
            M(ient)=-1d0
            b(ik)=0d0
            ient=ient+1
            ic(ient)=ik+1
            M(ient)=1d0
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+2    !Neumann conditions on bottom (zero current, parallel field, normal derivative)
          elseif (ix1==lx1) then
            if (flagdirich==1) then
              ic(ient)=ik
              M(ient)=1d0
              b(ik)=Vmaxx1(ix2,ix3)    !Dirichlet on top BCs
              ient=ient+1
              ir_ptr(ik+1)=ir_ptr(ik)+1
            else
              ic(ient)=ik-1
              M(ient)=-1d0/dx1(lx1)
              b(ik)=Vmaxx1(ix2,ix3)
              ient=ient+1
              ic(ient)=ik
              M(ient)=1d0/dx1(lx1)
              ient=ient+1
              ir_ptr(ik+1)=ir_ptr(ik)+2
            end if
          elseif (ix2==1) then
            ic(ient)=ik  
            M(ient)=1d0
            b(ik)=Vminx2(ix1,ix3)
!            b(ik)=Vmaxx1(1,ix3)   !match the top BCs
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+1
          elseif (ix2==lx2) then
            ic(ient)=ik
            M(ient)=1d0
            b(ik)=Vmaxx2(ix1,ix3)
!            b(ik)=Vmaxx1(lx2,ix3)
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+1
          elseif (ix3==1) then
            ic(ient)=ik  
            M(ient)=1d0
            b(ik)=Vminx3(ix1,ix2)
!            b(ik)=Vmaxx1(ix2,1)
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+1
          elseif (ix3==lx3) then
            ic(ient)=ik
            M(ient)=1d0
            b(ik)=Vmaxx3(ix1,ix2)
!            b(ik)=Vmaxx1(ix2,lx3)
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+1
          else                       !INTERIOR
            !ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
            ic(ient)=ik-lx1*lx2
            M(ient)=Bc(ix1,ix2,ix3)/dx3(ix3)/dx3i(ix3)-Ec(ix1,ix2,ix3)/(dx3(ix3+1)+dx3(ix3))
            ient=ient+1

            !ix1,ix2-1,ix3
            ic(ient)=ik-lx1
            M(ient)=Ac(ix1,ix2,ix3)/dx2(ix2)/dx2i(ix2)-Dc(ix1,ix2,ix3)/(dx2(ix2+1)+dx2(ix2))
            ient=ient+1

            !ix1-1,ix2,ix3
            ic(ient)=ik-1
            M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1)/dx1i(ix1)-Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))
            ient=ient+1

            !ix1,ix2,ix3
            ic(ient)=ik
            M(ient)=-1d0*Ac(ix1,ix2,ix3)*(1d0/dx2(ix2+1)/dx2i(ix2)+1d0/dx2(ix2)/dx2i(ix2))- &
                         Bc(ix1,ix2,ix3)*(1d0/dx3(ix3+1)/dx3i(ix3)+1d0/dx3(ix3)/dx3i(ix3))- &
                         Cc(ix1,ix2,ix3)*(1d0/dx1(ix1+1)/dx1i(ix1)+1d0/dx1(ix1)/dx1i(ix1))
            ient=ient+1

            !ix1+1,ix2,ix3
            ic(ient)=ik+1
            M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1+1)/dx1i(ix1)+Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))
            ient=ient+1

            !ix1,ix2+1,ix3
            ic(ient)=ik+lx1
            M(ient)=Ac(ix1,ix2,ix3)/dx2(ix2+1)/dx2i(ix2)+Dc(ix1,ix2,ix3)/(dx2(ix2+1)+dx2(ix2))
            ient=ient+1

            !ix1,ix2,ix3+1
            ic(ient)=ik+lx1*lx2
            M(ient)=Bc(ix1,ix2,ix3)/dx3(ix3+1)/dx3i(ix3)+Ec(ix1,ix2,ix3)/(dx3(ix3+1)+dx3(ix3))
            ient=ient+1

            ir_ptr(ik+1)=ir_ptr(ik)+7
          end if
        end do
      end do
    end do
!    write(*,*) 'nnz+1=',ir_ptr(lk+1),ient,lent,size(M),size(b),size(ic),size(ir_ptr),maxval(ic),maxval(ir_ptr)


    !MAKE A COPY OF THE MATRIX
    ir_ptr_copy=ir_ptr; ic_copy=ic; M_copy=M; b_copy=b;


!     For `_' in {S, D, C, Z}: _GNLAMGinit              
!     init default parameter
    write(*,*) 'Initiating call to ILUPACK...'
    call dgnlamginit(lk,ir_ptr,ic,M,matching,ordering,droptol,droptolS,condest,restol, &
                       maxit,elbow,lfil,lfilS,nrestart,mixedprecision,ind)
    write(*,*) 'ILUPACK init successful...'
!    write(*,*) ordering,droptol,droptolS,condest,restol,maxit,elbow,lfil,lfilS,nrestart,mixedprecision,size(ind)

!     now the use may vary the parameters to gain optimal performance

!     maximum weight matching
!     default value is different from zero, matching turned on
!     matching=1

!     multilevel orderings
!     'amd' (default) Approximate Minimum Degree (1.25x factorization time of rcm)
!     'mmd'           Minimum Degree  (doesn't work)          
!     'rcm'           Reverse Cuthill-McKee
!     'metisn'        Metis multilevel nested dissection by nodes (3x rcm)
!     'metise'        Metis multilevel nested dissection by edges (seg. faults)
!     'pq'            ddPQ strategy by Saad (doesn't work)
    ordering='rcm'

!     threshold for ILU, default: 1e-2
!    droptol=0.1    !original program
!    droptol=0.02    !some of the Neumann problems benefit from this
!    droptol=0.01
    droptol=1d-3

!     threshold for the approximate Schur complements, default: 0.1*droptol
    droptolS=0.1*droptol

!     norm bound for the inverse factors L^{-1}, U^{-1}, default: 1e2
!    condest=10    !orig
    condest=100d0   !slight performance boost (~10 pct. faster)

!     relative error for the backward error (SPD case: relative energy
!     norm) used during the iterative solver, default: sqrt(eps)
!      restol=1e-12    !original code and works fairly well on ionospheric problem
    if (flagdirich==1) then
!      restol=1d-14    !dirichlet problem will converge to very tight tolerance in only a few interations
      restol=1d-16    !the highest that will converge in a reasonable time
    else
!      restol=1d-10    !neumann problem requires a much looser tolerance to converge at all...
      restol=5d-8
    end if
!     maximum number of iterative steps, default: 500
    maxit=5000

!     elbow space factor for the fill computed during the ILU, default: 10
!    elbow=100    !causes an incredible amount of memory to be block off...
    elbow=15

!     maximum number of nonzeros per column in L/ row in U, default: n+1
!     lfil=10

!     maximum number of nonzeros per row in the approximate Schur complement,
!     default: n+1
!     lfilS=10

!     restart length for GMRES, default: 30
!      nrestart=50

!     do you want to use a single precision preconditioner
!      mixedprecision=1   !example default
    mixedprecision=0

!     underlying block structure, only partially supported
!     default: right now turn it off!
    ind(1:lk)=0

!     compute multilevel ILU
!     cccccccccccccccccccccc
!     Note that the initial input matrix A will be rescaled by rows and
!     by columns (powers of 2.0) and that the order in the array might have
!     been altered
!     if you do need the original matrix (ia,ja,a) in for differen purposes,
!     you should use a copy (ib,jb,b) instead

!     For `_' in {S, D, C, Z}: _GNLAMGfactor 
!     compute multilevel ILU `PREC'
    write(*,*) 'Beginning incomplete LU factorization...'
    ierr=dgnlamgfactor(param,PREC,lk,ir_ptr_copy,ic_copy,M_copy,matching,ordering,droptol,droptolS,condest,restol, &
                         maxit,elbow,lfil,lfilS,nrestart,mixedprecision,ind)

    if (ierr.eq.-1) then
       write (6,'(A)') 'Error. input matrix may be wrong.'
    elseif (ierr.eq.-2) then
       write (6,'(A)') 'matrix L overflow, increase elbow and retry'
    elseif (ierr.eq.-3) then
       write (6,'(A)') 'matrix U overflow, increase elbow and retry'
    elseif (ierr.eq.-4) then
       write (6,'(A)') 'Illegal value for lfil'
    elseif (ierr.eq.-5) then
       write (6,'(A)') 'zero row encountered'
    elseif (ierr.eq.-6) then
       write (6,'(A)') 'zero column encountered'
    elseif (ierr.eq.-7) then
       write (6,'(A)') 'buffers are too small'
    elseif (ierr.ne.0) then
       write (6,'(A,I3)') 'zero pivot encountered at step number',ierr
    end if
      if (ierr.ne.0) stop


!     Just for fun: display multilevel information on the screen
!     For `_' in {S, D, C, Z}: _GNLAMGinfo
    write (6,'(A,F8.2)') '   final elbow space factor=',elbow+0.005
    write (6,'(A,F8.2)') '   final condest on level 1=',condest+0.005
    write (6,'(A)')      'ILUPACK,   multilevel structure'
    call dgnlamginfo(param,PREC,lk,ir_ptr_copy,ic_copy,M_copy)

!     Just for fun: if you want to know the logical number of nonzeros only
    mem=dgnlamgnnz(param,PREC)
    write (6,'(A,1P,E8.1)') 'fill-in factor nnz(LU)/nnz(A)',dble(mem)/dble(ir_ptr(lk+1)-1)


!     solve a single system with `PREC'
!     ccccccccccccccccccccccccccccccccc
!     This might be of interest if you want to apply ILUPACK inside your
!     own iterative method without referring to the convenient driver
!     artificial right hand side b=A*1
!      write (6,'(A)') 'right hand side'
!      do ik=1,lk
!         write(6,'(1P,E12.4)') b(ik)
!      enddo


!     For `_' in {S, D, C, Z}: _GNLAMGsol
!     solve a single linear system with `PREC'
!      call dgnlamgsol(param,PREC,b,sol,lk)


!      write (6,'(A)') 'ilu approximate solution'
!      do ik=1,lk
!         write(6,'(1P,E12.4)') sol(ik)
!      enddo


!     ok, the preconditioner is usually not exact.
!     as convenient ALTERNATIVE, ILUPACK offers an all-in-one iterative solver


!     solve Ax=b  until the desired accuracy is achieved
!     cccccccccccccccccccccccccccccccccccccccccccccccccc
!     provide an initial solution, e.g. 0
    sol(1:lk)=pack(Phi0,.true.)
!     For `_' in {S, D, C, Z}: _GNLAMGsolver
!     solve Ax=b iteratively
    ierr=dgnlamgsolver(param,PREC,b,sol,lk,ir_ptr_copy,ic_copy,M_copy,matching,ordering,droptol,droptolS,condest,restol, &
                         maxit,elbow,lfil,lfilS,nrestart,mixedprecision,ind)

    if (ierr.eq.-1) then
       write (6,'(A)') 'too many iterations'
    elseif (ierr.eq.-2) then
       write (6,'(A)') 'not enough work space'
    elseif (ierr.eq.-3) then
       write (6,'(A)') 'algorithm breaks down'
    elseif (ierr.ne.0) then
       write (6,'(A,I4)') 'solver exited with error code',ierr
    end if

    write (6,'(I4,A)') maxit,' iterations needed'


    !SAVE SOLUTION AS OUTPUT VARIABLE
    elliptic3D=reshape(sol(1:lk),[lx1,lx2,lx3])


    !DEALLOCATE MEMORY
    call dgnlamgdelete(param,PREC)
    deallocate(ir_ptr,ic,M,b)
    deallocate(ir_ptr_copy,ic_copy,M_copy,b_copy)
    deallocate(ind)
    deallocate(sol)

  end function elliptic3D


  function elliptic3D_pol_conv(srcterm,sig0,sigP,sigH,cm,v2,v3, &
                               Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                               dt,dx1,dx1i,dx2,dx2i,dx3,dx3i,Phi0,flagdirich)

    !------------------------------------------------------------
    !-------SOLVE POTENTIAL EQUATION IN 3D USING ILUPACK.  THIS
    !-------IMPLEMENTATION INCLUDE POLARIZATION CURRENT IN THE MATRIX
    !-------SOLUTION (BOTH TIME DERIV. AND CONVECTIVE TERMS).  VELOCITIES
    !-------SHOULD BE TRIMMED OF GHOST CELLS.  THIS CODE SHOULD PROBABLY
    !-------ALWAYS USE DIRICHLET CONDITIONS.
    !------------------------------------------------------------

    real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP,sigH,Phi0,cm,v2,v3
    real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
    real(wp), intent(in) :: dt
    real(wp), dimension(0:), intent(in) :: dx1
    real(wp), dimension(:), intent(in) :: dx1i
    real(wp), dimension(0:), intent(in) :: dx2
    real(wp), dimension(:), intent(in) :: dx2i
    real(wp), dimension(0:), intent(in) :: dx3
    real(wp), dimension(:), intent(in) :: dx3i
    integer(8), intent(in) :: flagdirich

    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigP2,gradsigP3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigH2,gradsigH3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsig01
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: cmh2,cmh3
    real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: Ac,Bc,Cc,Dc,Ec,Fc

    real(wp) :: coeff

    !EXTERNAL ILUPACK ROUTINES
    integer, external :: dgnlamgfactor,dgnlamgsolver,dgnlamgnnz
    external dgnlamginit,dgnlamgsol,dgnlamgdelete,dgnlamginfo

    !ILUPACK external parameters
    integer :: mem
    integer :: matching,maxit,lfil,lfilS,nrestart,ierr,mixedprecision
    integer, dimension(:), allocatable :: ind
    character(20) :: ordering
    real(wp) :: droptol,droptolS,condest,restol,elbow

    !Variables that cover the and pass the C-pointers
    integer(8) :: param,PREC

    !LOCAL STORAGE
    integer :: ix1,ix2,ix3,lx1,lx2,lx3
    integer :: ik,ient,lk,lent
    integer, dimension(:), allocatable :: ir_ptr
    integer, dimension(:), allocatable :: ic
    real(wp), dimension(:), allocatable :: M
    real(wp), dimension(:), allocatable :: b
    integer, dimension(:), allocatable :: ir_ptr_copy
    integer, dimension(:), allocatable :: ic_copy
    real(wp), dimension(:), allocatable :: M_copy
    real(wp), dimension(:), allocatable :: b_copy
    real(wp), dimension(:), allocatable :: sol

    real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D_pol_conv


    !SYSTEM SIZES AND ALLOCATION FOR COMPRESSED ROW STORAGE
    write(*,*) 'Allocating space for DYNAMIC elliptic equation solution...'
    lx1=size(srcterm,1)
    lx2=size(srcterm,2)
    lx3=size(srcterm,3)
    lk=lx1*lx2*lx3

!    lent=7*(lx1-2)*(lx2-2)*(lx3-2)                                                 !interior entries for static problem
    lent=19*(lx1-2)*(lx2-2)*(lx3-2)                                                !interior entries for polarization current problem
    lent=lent+2*(lx1-2)*(lx2-2)+2*(lx2-2)*(lx3-2)+2*(lx1-2)*(lx3-2)                !6 faces of cube
    lent=lent+4*(lx1-2)+4*(lx2-2)+4*(lx3-2)                                        !12 edges
    lent=lent+8                                                                    !8 corners
    lent=lent+lx2*lx3                                                              !entries to deal with Neumann conditions on bottom
    lent=lent-3*2*(lx2-2)*(lx1-2)-3*2*(lx3-2)*(lx1-2)                              !-x2_adj-x3_adj (3 entries for each)
    if (flagdirich==0) then                                                        !if we want Neumann on top boundary conditions
      lent=lent+lx2*lx3
    end if

    write(*,*) 'Total unknowns and nonzero entries:  ',lk,lent
    allocate(ir_ptr(lk+1),ic(lent),M(lent),b(lk))
    allocate(ir_ptr_copy(lk+1),ic_copy(lent),M_copy(lent),b_copy(lk))
    allocate(ind(lk))      !not sure what this is used for...
    allocate(sol(lk))


    !INITIALIZATION
    M(:)=0d0
    b(:)=pack(srcterm,.true.)


    !COMPUTE AUXILIARY COEFFICIENTS
    write(*,*) 'Prepping coefficients for elliptic equation...'    
    gradsig01=grad3D1(sig0,dx1(1:lx1))
    gradsigP2=grad3D2(sigP,dx2(1:lx2)); gradsigP3=grad3D3(sigP,dx3(1:lx3));
    gradsigH2=grad3D2(sigH,dx2(1:lx2)); gradsigH3=grad3D3(sigH,dx3(1:lx3));

    cmh2(:,1,:)=0d0
    cmh2(:,2:lx2,:)=0.5d0*(cm(:,1:lx2-1,:)+cm(:,2:lx2,:))
    cmh3(:,:,1)=0d0
    cmh3(:,:,2:lx3)=0.5d0*(cm(:,:,1:lx3-1)+cm(:,:,2:lx3))

    Ac=sigP; Bc=sigP; Cc=sig0; 
    Dc=gradsigP2+gradsigH3
    Ec=gradsigP3-gradsigH2
    Fc=gradsig01


    !LOAD UP MATRIX ELEMENTS
    write(*,*) 'Building CRS sparse matrix...'
    ient=1
    ir_ptr(1)=1
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          ik=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix2,ix3) as a column vector.  Also row no. of big matrix

          if (ix1==1) then    !Check whether we are a bondary point and, if so, enforce boundary conditions
            ic(ient)=ik
            M(ient)=-1d0
            b(ik)=0d0
            ient=ient+1
            ic(ient)=ik+1
            M(ient)=1d0
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+2    !Neumann conditions on bottom (zero current, parallel field, normal derivative)
          elseif (ix1==lx1) then
            if (flagdirich==1) then
              ic(ient)=ik
              M(ient)=1d0
              b(ik)=Vmaxx1(ix2,ix3)    !Dirichlet on top BCs
              ient=ient+1
              ir_ptr(ik+1)=ir_ptr(ik)+1
            else
              ic(ient)=ik-1
              M(ient)=-1d0/dx1(lx1)
              b(ik)=Vmaxx1(ix2,ix3)
              ient=ient+1
              ic(ient)=ik
              M(ient)=1d0/dx1(lx1)
              ient=ient+1
              ir_ptr(ik+1)=ir_ptr(ik)+2
            end if
          elseif (ix2==1) then
            ic(ient)=ik  
            M(ient)=1d0
            b(ik)=Vminx2(ix1,ix3)
!            b(ik)=Vmaxx1(1,ix3)   !match the top BCs
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+1
          elseif (ix2==lx2) then
            ic(ient)=ik
            M(ient)=1d0
            b(ik)=Vmaxx2(ix1,ix3)
!            b(ik)=Vmaxx1(lx2,ix3)
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+1
          elseif (ix3==1) then
            ic(ient)=ik  
            M(ient)=1d0
            b(ik)=Vminx3(ix1,ix2)
!            b(ik)=Vmaxx1(ix2,1)
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+1
          elseif (ix3==lx3) then
            ic(ient)=ik
            M(ient)=1d0
            b(ik)=Vmaxx3(ix1,ix2)
!            b(ik)=Vmaxx1(ix2,lx3)
            ient=ient+1
            ir_ptr(ik+1)=ir_ptr(ik)+1
          else                       !INTERIOR
            !ix1,ix2-1,ix3-2 grid point
            coeff=-1d0*cm(ix1,ix2,ix3-1)*v2(ix1,ix2,ix3-1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3-1)+dx3(ix3))*(dx2(ix2)+dx2(ix2+1)) )
            if (ix3==2) then    !out of bounds, use nearest BC, and add to known vector
              b(ik)=b(ik)-coeff*Vminx3(ix1,ix2-1)
            else    !in bounds, add to matrix
              ic(ient)=ik-lx1-2*lx2*lx1
              M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
              ient=ient+1
            end if


            !ix1,ix2,ix3-2 grid point
            coeff=-1d0*cm(ix1,ix2,ix3-1)*v3(ix1,ix2,ix3-1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3-1)*dx3i(ix3-1)) )
            if (ix3==2) then    !bit of intentional code duplication here and in the following sections to keep things organized in a way I can debug
              b(ik)=b(ik)-coeff*Vminx3(ix1,ix2)
            else
              ic(ient)=ik-2*lx2*lx1
              M(ient)=coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term
              ient=ient+1
            end if


            !ix1,ix2+1,ix3-2 grid point
            coeff=cm(ix1,ix2,ix3-1)*v2(ix1,ix2,ix3-1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3-1)+dx3(ix3))*(dx2(ix2)+dx2(ix2+1)) )
            if (ix3==2) then
              b(ik)=b(ik)-coeff*Vminx3(ix1,ix2+1)
            else
              ic(ient)=ik+lx1-2*lx2*lx1
              M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
              ient=ient+1
            end if


            !ix1,ix2-2,ix3-1
            coeff=-1d0*cm(ix1,ix2-1,ix3)*v3(ix1,ix2-1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2-1)+dx2(ix2))*(dx3(ix3)+dx3(ix3+1)) )
            if (ix2==2) then
              b(ik)=b(ik)-coeff*Vminx2(ix1,ix3-1)              
            else
              ic(ient)=ik-2*lx1-lx1*lx2
              M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )
              ient=ient+1
            end if


            !ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
            ic(ient)=ik-lx1*lx2
            M(ient)=Bc(ix1,ix2,ix3)/dx3(ix3)/dx3i(ix3)-Ec(ix1,ix2,ix3)/(dx3(ix3+1)+dx3(ix3))    !static terms

            coeff=cmh3(ix1,ix2,ix3)/(dt*dx3i(ix3)*dx3(ix3))
            M(ient)=M(ient)+coeff   !polarization time derivative terms
            b(ik)=b(ik)+coeff*Phi0(ix1,ix2,ix3-1)   !add in polarziation terms that include previous time step potential at this grid point

            coeff=cm(ix1,ix2,ix3-1)*v3(ix1,ix2,ix3-1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3)*dx3i(ix3-1)) )+ &
                  cm(ix1,ix2,ix3-1)*v3(ix1,ix2,ix3-1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3-1)*dx3i(ix3-1)) )
            M(ient)=M(ient)+coeff   !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

            coeff=cm(ix1,ix2+1,ix3)*v3(ix1,ix2+1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2+1)+dx2(ix2+2))*(dx3(ix3)+dx3(ix3+1)) )+ &
                  cm(ix1,ix2-1,ix3)*v3(ix1,ix2-1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2-1)+dx2(ix2))*(dx3(ix3)+dx3(ix3+1)) )
            M(ient)=M(ient)+coeff   !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

            ient=ient+1


            !ix1,ix2+2,ix3-1 grid point
            coeff=-1d0*cm(ix1,ix2+1,ix3)*v3(ix1,ix2+1,ix3)/( (dx2(ix2+1)+dx2(ix2+2))*(dx2(ix2)+dx2(ix2+1))*(dx3(ix3)+dx3(ix3+1)) )
            if (ix2==lx2-1) then
              b(ik)=b(ik)-coeff*Vmaxx2(ix1,ix3-1)              
            else
              ic(ient)=ik+2*lx1-lx1*lx2
              M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )
              ient=ient+1
            end if


            !ix2-2,ix3 grid point
            coeff=-1d0*cm(ix1,ix2-1,ix3)*v2(ix1,ix2-1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2-1)*dx2i(ix2-1)) )
            if (ix2==2) then
              b(ik)=b(ik)-coeff*Vminx2(ix1,ix3)              
            else
              ic(ient)=ik-2*lx1
              M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term
              ient=ient+1
            end if


            !ix1,ix2-1,ix3
            ic(ient)=ik-lx1
            M(ient)=Ac(ix1,ix2,ix3)/dx2(ix2)/dx2i(ix2)-Dc(ix1,ix2,ix3)/(dx2(ix2+1)+dx2(ix2))    !static terms

            coeff=cmh2(ix1,ix2,ix3)/(dt*dx2i(ix2)*dx2(ix2))
            M(ient)=M(ient)+coeff    !pol. time deriv.
            b(ik)=b(ik)+coeff*Phi0(ix1,ix2-1,ix3)    !BC's and pol. time deriv.

            coeff=cm(ix1,ix2-1,ix3)*v2(ix1,ix2-1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2)*dx2i(ix2-1)) )+ &
                  cm(ix1,ix2-1,ix3)*v2(ix1,ix2-1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2-1)*dx2i(ix2-1)) )
            M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

            coeff=cm(ix1,ix2,ix3+1)*v2(ix1,ix2,ix3+1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3+1)+dx3(ix3+2))*(dx2(ix2)+dx2(ix2+1)) )+ &
                  cm(ix1,ix2,ix3-1)*v2(ix1,ix2,ix3-1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3-1)+dx3(ix3))*(dx2(ix2)+dx2(ix2+1)) )
            M(ient)=M(ient)+coeff     !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

            ient=ient+1


            !ix1-1,ix2,ix3
            ic(ient)=ik-1
            M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1)/dx1i(ix1)-Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))    !only static terms here (pol. current assumed perp. to geomag. field
            ient=ient+1


            !ix1,ix2,ix3 (main diagonal)
            ic(ient)=ik
            M(ient)=-1d0*Ac(ix1,ix2,ix3)*(1d0/dx2(ix2+1)/dx2i(ix2)+1d0/dx2(ix2)/dx2i(ix2))- &
                         Bc(ix1,ix2,ix3)*(1d0/dx3(ix3+1)/dx3i(ix3)+1d0/dx3(ix3)/dx3i(ix3))- &
                         Cc(ix1,ix2,ix3)*(1d0/dx1(ix1+1)/dx1i(ix1)+1d0/dx1(ix1)/dx1i(ix1))

            coeff=-1d0*cmh2(ix1,ix2+1,ix3)/(dt*dx2i(ix2)*dx2(ix2+1)) &
                  -1d0*cmh2(ix1,ix2,ix3)/(dt*dx2i(ix2)*dx2(ix2)) &
                  -1d0*cmh3(ix1,ix2,ix3+1)/(dt*dx3i(ix3)*dx3(ix3+1)) &
                  -1d0*cmh3(ix1,ix2,ix3)/(dt*dx3i(ix3)*dx3(ix3))
            M(ient)=M(ient)+coeff    !pol. time deriv.
            b(ik)=b(ik)+coeff*Phi0(ix1,ix2,ix3)    !BC's and pol. time deriv.

            coeff=cm(ix1,ix2+1,ix3)*v2(ix1,ix2+1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2+1)*dx2i(ix2+1)) )+ &
                  (-1d0)*cm(ix1,ix2-1,ix3)*v2(ix1,ix2-1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2)*dx2i(ix2-1)) )
            M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

            coeff=cm(ix1,ix2,ix3+1)*v3(ix1,ix2,ix3+1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3+1)*dx3(ix3+1)) )+ &
                  (-1d0)*Cm(ix1,ix2,ix3-1)*v3(ix1,ix2,ix3-1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3)*dx3i(ix3-1)) )
            M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

            ient=ient+1


            !ix1+1,ix2,ix3
            ic(ient)=ik+1
            M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1+1)/dx1i(ix1)+Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))    !static terms (all that are used here)
            ient=ient+1


            !ix1,ix2+1,ix3
            ic(ient)=ik+lx1
            M(ient)=Ac(ix1,ix2,ix3)/dx2(ix2+1)/dx2i(ix2)+Dc(ix1,ix2,ix3)/(dx2(ix2+1)+dx2(ix2))    !static terms

            coeff=cmh2(ix1,ix2+1,ix3)/(dt*dx2i(ix2)*dx2(ix2+1))
            M(ient)=M(ient)+coeff    !pol. time deriv. terms
            b(ik)=b(ik)+coeff*Phi0(ix1,ix2+1,ix3)    !BC's and pol. time deriv.  

            coeff=-1d0*cm(ix1,ix2+1,ix3)*v2(ix1,ix2+1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2+2)*dx2i(ix2+1)) )+ &
                  (-1d0)*cm(ix1,ix2+1,ix3)*v2(ix1,ix2+1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2+1)*dx2i(ix2+1)) )
            M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term    
  
            coeff=-1d0*cm(ix1,ix2,ix3+1)*v2(ix1,ix2,ix3+1)/ &
                         ( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3+1)+dx3(ix3+2))*(dx2(ix2)+dx2(ix2+1)) )+ &
                  (-1d0)*cm(ix1,ix2,ix3-1)*v2(ix1,ix2,ix3-1)/ &
                         ( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3-1)+dx3(ix3))*(dx2(ix2)+dx2(ix2+1)) )
            M(ient)=M(ient)+coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

            ient=ient+1


            !ix1,ix2+2,ix3 grid point
            coeff=cm(ix1,ix2+1,ix3)*v2(ix1,ix2+1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2+2)*dx2i(ix2+1)) )
            if (ix2==lx2-1) then
              b(ik)=b(ik)-coeff*Vmaxx2(ix1,ix3)              
            else
              ic(ient)=ik+2*lx1
              M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term
              ient=ient+1
            end if


            !ix1,ix2-2,ix3+1 grid point
            coeff=cm(ix1,ix2-1,ix3)*v3(ix1,ix2-1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2-1)+dx2(ix2))*(dx3(ix3)+dx3(ix3+1)) )
            if (ix2==2) then
              b(ik)=b(ik)-coeff*Vminx2(ix1,ix3+1)              
            else
              ic(ient)=ik-2*lx1+lx1*lx2
              M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )
              ient=ient+1
            end if


            !ix1,ix2,ix3+1
            ic(ient)=ik+lx1*lx2
            M(ient)=Bc(ix1,ix2,ix3)/dx3(ix3+1)/dx3i(ix3)+Ec(ix1,ix2,ix3)/(dx3(ix3+1)+dx3(ix3))

            coeff=cmh3(ix1,ix2,ix3+1)/(dt*dx3i(ix3)*dx3(ix3+1))
            M(ient)=M(ient)+coeff    !pol. time deriv.
            b(ik)=b(ik)+coeff*Phi0(ix1,ix2,ix3+1)    !BC's and pol. time deriv.  

            coeff=-1d0*cm(ix1,ix2,ix3+1)*v3(ix1,ix2,ix3+1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3+2)*dx3i(ix3+1)) )+ &
                  (-1d0)*cm(ix1,ix2,ix3+1)*v3(ix1,ix2,ix3+1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3+1)*dx3i(ix3+1)) )
            M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

            coeff=-1d0*cm(ix1,ix2+1,ix3)*v3(ix1,ix2+1,ix3)/ &
                         ( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2+1)+dx2(ix2+2))*(dx3(ix3)+dx3(ix3+1)) )+ &
                  (-1d0)*cm(ix1,ix2-1,ix3)*v3(ix1,ix2-1,ix3)/ &
                         ( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2-1)+dx2(ix2))*(dx3(ix3)+dx3(ix3+1)) )
            M(ient)=M(ient)+coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

            ient=ient+1


            !ix1,ix2+2,ix3+1 grid point
            coeff=cm(ix1,ix2+1,ix3)*v3(ix1,ix2+1,ix3)/( (dx2(ix2)+dx2(ix2+1))*(dx2(ix2+1)+dx2(ix2+2))*(dx3(ix3)+dx3(ix3+1)) )
            if (ix2==lx2-1) then
              b(ik)=b(ik)-coeff*Vmaxx2(ix1,ix3+1)              
            else
              ic(ient)=ik+2*lx1+lx1*lx2
              M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )
              ient=ient+1
            end if


            !ix1,ix2-1,ix3+2 grid point
            coeff=-1d0*cm(ix1,ix2,ix3+1)*v2(ix1,ix2,ix3+1)/( (dx3(ix3+1)+dx3(ix3+2))*(dx3(ix3)+dx3(ix3+1))*(dx2(ix2)+dx2(ix2+1)) )
            if (ix3==lx3-1) then
              b(ik)=b(ik)-coeff*Vmaxx3(ix1,ix2-1)              
            else
              ic(ient)=ik-lx1+2*lx1*lx2
              M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
              ient=ient+1
            end if


            !ix1,ix2,ix3+2 grid point
            coeff=cm(ix1,ix2,ix3+1)*v3(ix1,ix2,ix3+1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3+2)*dx3i(ix3+1)) )
            if (ix3==lx3-1) then
              b(ik)=b(ik)-coeff*Vmaxx3(ix1,ix2)              
            else
              ic(ient)=ik+2*lx1*lx2
              M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )
              ient=ient+1
            end if


            !ix1,ix2+1,ix3+2 grid point
            coeff=cm(ix1,ix2,ix3+1)*v2(ix1,ix2,ix3+1)/( (dx3(ix3)+dx3(ix3+1))*(dx3(ix3+1)+dx3(ix3+2))*(dx2(ix2)+dx2(ix2+1)) )
            if (ix3==lx3-1) then
              b(ik)=b(ik)-coeff*Vmaxx3(ix1,ix2+1)              
            else
              ic(ient)=ik+lx1+2*lx1*lx2
              M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
              ient=ient+1
            end if     


            ir_ptr(ik+1)=ir_ptr(ik)+19
            if (ix2==2 .or. ix2==lx2-1) then    !must remove entries not used b/c we are adjacent to boundary
              ir_ptr(ik+1)=ir_ptr(ik+1)-3
            end if
            if (ix3==2 .or. ix3==lx3-1) then
              ir_ptr(ik+1)=ir_ptr(ik+1)-3
            end if
          end if
        end do
      end do
    end do
!    write(*,*) 'nnz+1=',ir_ptr(lk+1),ient,lent,size(M),size(b),size(ic),size(ir_ptr),maxval(ic),maxval(ir_ptr)


    !MAKE A COPY OF THE MATRIX
    ir_ptr_copy=ir_ptr; ic_copy=ic; M_copy=M; b_copy=b;


!     For `_' in {S, D, C, Z}: _GNLAMGinit              
!     init default parameter
    write(*,*) 'Initiating call to ILUPACK...'
    call dgnlamginit(lk,ir_ptr,ic,M,matching,ordering,droptol,droptolS,condest,restol, &
                       maxit,elbow,lfil,lfilS,nrestart,mixedprecision,ind)
    write(*,*) 'ILUPACK init successful...'
!    write(*,*) ordering,droptol,droptolS,condest,restol,maxit,elbow,lfil,lfilS,nrestart,mixedprecision,size(ind)

!     now the use may vary the parameters to gain optimal performance

!     maximum weight matching
!     default value is different from zero, matching turned on
!     matching=1

!     multilevel orderings
!     'amd' (default) Approximate Minimum Degree (1.25x factorization time of rcm)
!     'mmd'           Minimum Degree  (doesn't work)          
!     'rcm'           Reverse Cuthill-McKee
!     'metisn'        Metis multilevel nested dissection by nodes (3x rcm)
!     'metise'        Metis multilevel nested dissection by edges (seg. faults)
!     'pq'            ddPQ strategy by Saad (doesn't work)
    ordering='rcm'

!     threshold for ILU, default: 1e-2
    droptol=0.1    !original program
!    droptol=0.02
!    droptol=0.005

!     threshold for the approximate Schur complements, default: 0.1*droptol
    droptolS=0.1*droptol

!     norm bound for the inverse factors L^{-1}, U^{-1}, default: 1e2
    condest=10    !orig
!    condest=100d0   !slight performance boost (~10 pct. faster)

!     relative error for the backward error (SPD case: relative energy
!     norm) used during the iterative solver, default: sqrt(eps)
!      restol=1e-12    !original code and works fairly well on ionospheric problem
    if (flagdirich==1) then
!      restol=1d-14    !dirichlet problem will converge to very tight tolerance in only a few interations
      restol=1d-16    !the highest that will converge in a reasonable time
    else
!      restol=1d-10    !neumann problem requires a much looser tolerance to converge at all...
      restol=5d-8    !polar cap arc problem may need looser tolerance to converge
    end if
!     maximum number of iterative steps, default: 500
    maxit=10000

!     elbow space factor for the fill computed during the ILU, default: 10
!    elbow=100    !causes an incredible amount of memory to be block off...
    elbow=15

!     maximum number of nonzeros per column in L/ row in U, default: n+1
!     lfil=10

!     maximum number of nonzeros per row in the approximate Schur complement,
!     default: n+1
!     lfilS=10

!     restart length for GMRES, default: 30
!      nrestart=50

!     do you want to use a single precision preconditioner
!      mixedprecision=1   !example default
    mixedprecision=0

!     underlying block structure, only partially supported
!     default: right now turn it off!
    ind(1:lk)=0

!     compute multilevel ILU
!     cccccccccccccccccccccc
!     Note that the initial input matrix A will be rescaled by rows and
!     by columns (powers of 2.0) and that the order in the array might have
!     been altered
!     if you do need the original matrix (ia,ja,a) in for differen purposes,
!     you should use a copy (ib,jb,b) instead

!     For `_' in {S, D, C, Z}: _GNLAMGfactor 
!     compute multilevel ILU `PREC'
    write(*,*) 'Beginning incomplete LU factorization...'
    ierr=dgnlamgfactor(param,PREC,lk,ir_ptr_copy,ic_copy,M_copy,matching,ordering,droptol,droptolS,condest,restol, &
                         maxit,elbow,lfil,lfilS,nrestart,mixedprecision,ind)

    if (ierr.eq.-1) then
       write (6,'(A)') 'Error. input matrix may be wrong.'
    elseif (ierr.eq.-2) then
       write (6,'(A)') 'matrix L overflow, increase elbow and retry'
    elseif (ierr.eq.-3) then
       write (6,'(A)') 'matrix U overflow, increase elbow and retry'
    elseif (ierr.eq.-4) then
       write (6,'(A)') 'Illegal value for lfil'
    elseif (ierr.eq.-5) then
       write (6,'(A)') 'zero row encountered'
    elseif (ierr.eq.-6) then
       write (6,'(A)') 'zero column encountered'
    elseif (ierr.eq.-7) then
       write (6,'(A)') 'buffers are too small'
    elseif (ierr.ne.0) then
       write (6,'(A,I3)') 'zero pivot encountered at step number',ierr
    end if
      if (ierr.ne.0) stop


!     Just for fun: display multilevel information on the screen
!     For `_' in {S, D, C, Z}: _GNLAMGinfo
    write (6,'(A,F8.2)') '   final elbow space factor=',elbow+0.005
    write (6,'(A,F8.2)') '   final condest on level 1=',condest+0.005
    write (6,'(A)')      'ILUPACK,   multilevel structure'
    call dgnlamginfo(param,PREC,lk,ir_ptr_copy,ic_copy,M_copy)

!     Just for fun: if you want to know the logical number of nonzeros only
    mem=dgnlamgnnz(param,PREC)
    write (6,'(A,1P,E8.1)') 'fill-in factor nnz(LU)/nnz(A)',dble(mem)/dble(ir_ptr(lk+1)-1)


!     solve a single system with `PREC'
!     ccccccccccccccccccccccccccccccccc
!     This might be of interest if you want to apply ILUPACK inside your
!     own iterative method without referring to the convenient driver
!     artificial right hand side b=A*1
!      write (6,'(A)') 'right hand side'
!      do ik=1,lk
!         write(6,'(1P,E12.4)') b(ik)
!      enddo


!     For `_' in {S, D, C, Z}: _GNLAMGsol
!     solve a single linear system with `PREC'
!      call dgnlamgsol(param,PREC,b,sol,lk)


!      write (6,'(A)') 'ilu approximate solution'
!      do ik=1,lk
!         write(6,'(1P,E12.4)') sol(ik)
!      enddo


!     ok, the preconditioner is usually not exact.
!     as convenient ALTERNATIVE, ILUPACK offers an all-in-one iterative solver


!     solve Ax=b  until the desired accuracy is achieved
!     cccccccccccccccccccccccccccccccccccccccccccccccccc
!     provide an initial solution, e.g. 0
    sol(1:lk)=pack(Phi0,.true.)
!     For `_' in {S, D, C, Z}: _GNLAMGsolver
!     solve Ax=b iteratively
    ierr=dgnlamgsolver(param,PREC,b,sol,lk,ir_ptr_copy,ic_copy,M_copy,matching,ordering,droptol,droptolS,condest,restol, &
                         maxit,elbow,lfil,lfilS,nrestart,mixedprecision,ind)

    if (ierr.eq.-1) then
       write (6,'(A)') 'too many iterations'
    elseif (ierr.eq.-2) then
       write (6,'(A)') 'not enough work space'
    elseif (ierr.eq.-3) then
       write (6,'(A)') 'algorithm breaks down'
    elseif (ierr.ne.0) then
       write (6,'(A,I4)') 'solver exited with error code',ierr
    end if

    write (6,'(I4,A)') maxit,' iterations needed'


    !SAVE SOLUTION AS OUTPUT VARIABLE
    elliptic3D_pol_conv=reshape(sol(1:lk),[lx1,lx2,lx3])


    !DEALLOCATE MEMORY
    call dgnlamgdelete(param,PREC)
    deallocate(ir_ptr,ic,M,b)
    deallocate(ir_ptr_copy,ic_copy,M_copy,b_copy)
    deallocate(ind)
    deallocate(sol)

  end function elliptic3D_pol_conv

end module potential_ilupack
