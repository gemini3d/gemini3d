!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Code below here possibly deprecated...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function poisson2D(rho,Vminx1,Vmaxx1,Vminx2,Vmaxx2,dx1,perflag)

!------------------------------------------------------------
!-------SOLVE POISSONS'S EQUATION IN 2D USING MUMPS.  THIS MAY
!-------BE USED BY SOME OF THE UNIT TEST PROGRAMS, BUT I FORGET...
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: rho
real(wp), dimension(:), intent(in) :: Vminx1,Vmaxx1
real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
real(wp), intent(in) :: dx1
!    real(wp), dimension(:), intent(in) :: dx2 !! I guess this is assumed to be teh same as dx1???
logical, intent(in) :: perflag

integer :: ix1,ix2,lx1,lx2
integer :: lPhi, lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
real(wp) :: tstart,tfin

#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#else
type (DMUMPS_STRUC) mumps_par
#endif
real(wp), dimension(size(rho,1),size(rho,2)) :: poisson2D


if (myid==0) then
!------------------------------------------------------------
!-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
!-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
!-------GUIDE).
!------------------------------------------------------------
lx1=size(rho,1)
lx2=size(rho,2)
lPhi=lx1*lx2
lent=5*(lx1-2)*(lx2-2)+2*lx1+2*(lx2-2)
allocate(ir(lent),ic(lent),M(lent),b(lPhi))


!LOAD UP MATRIX ELEMENTS
M(:)=0.0
b=pack(rho,.true.)           !boundaries overwritten later
ient=1
do ix2=1,lx2
 do ix1=1,lx1
   iPhi=lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix2) as a column vector.  Also row of big matrix

   if (ix1==1) then          !BOTTOM GRID POINTS + CORNER
     ir(ient)=iPhi
     ic(ient)=iPhi
     M(ient)=1.0
     b(iPhi)=Vminx1(ix2)
     ient=ient+1
   elseif (ix1==lx1) then    !TOP GRID POINTS + CORNER
     ir(ient)=iPhi
     ic(ient)=iPhi
     M(ient)=1.0
     b(iPhi)=Vmaxx1(ix2)
     ient=ient+1
   elseif (ix2==1) then      !LEFT BOUNDARY
     ir(ient)=iPhi
     ic(ient)=iPhi
     M(ient)=1.0
     b(iPhi)=Vminx2(ix1)
     ient=ient+1
   elseif (ix2==lx2) then    !RIGHT BOUNDARY
     ir(ient)=iPhi
     ic(ient)=iPhi
     M(ient)=1.0
     b(iPhi)=Vmaxx2(ix1)
     ient=ient+1
   else                      !INTERIOR
     !ix1,ix2-1 grid point in ix1,ix2 equation
     ir(ient)=iPhi
     ic(ient)=iPhi-lx1
     M(ient)=1.0
     ient=ient+1

     !ix1-1,ix2 grid point
     ir(ient)=iPhi
     ic(ient)=iPhi-1
     M(ient)=1.0
     ient=ient+1

     !ix1,ix2 grid point
     ir(ient)=iPhi
     ic(ient)=iPhi
     M(ient)=-4.0
     ient=ient+1

     !ix1+1,ix2 grid point
     ir(ient)=iPhi
     ic(ient)=iPhi+1
     M(ient)=1.0
     ient=ient+1

     !ix1,ix2+1 grid point
     ir(ient)=iPhi
     ic(ient)=iPhi+lx1
     M(ient)=1.0
     ient=ient+1
   end if
 end do
end do


!CORRECT FOR DX /= 1
b=b*dx1**2
end if

!FIRE UP MUMPS
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
#if REALBITS==32
call SMUMPS(mumps_par)
#else
call DMUMPS(mumps_par)
#endif

call quiet_mumps(mumps_par)

!LOAD OUR PROBLEM (ROOT ONLY)
if ( myid == 0 ) then
 mumps_par%N=lPhi
 mumps_par%NZ=lent
 allocate( mumps_par%IRN ( mumps_par%NZ ) )
 allocate( mumps_par%JCN ( mumps_par%NZ ) )
 allocate( mumps_par%A( mumps_par%NZ ) )
 allocate( mumps_par%RHS ( mumps_par%N  ) )
 mumps_par%IRN=ir
 mumps_par%JCN=ic
 mumps_par%A=M
 mumps_par%RHS=b
 deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

 if (perflag) then       !used cached permutation
   allocate(mumps_par%PERM_IN(mumps_par%N))
   mumps_par%PERM_IN=mumps_perm
   mumps_par%ICNTL(7)=1
 end if
end if


!SOLVE OUR PROBLEM (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
call cpu_time(tstart)
#if REALBITS==32
call SMUMPS(mumps_par)
#else
call DMUMPS(mumps_par)
#endif

!> check if Mumps error occurred
if (mumps_par%INFO(1) < 0 .or. mumps_par%INFOG(1) < 0) then
 write(stderr, *) 'Gemini:potential_mumps:poisson2D  MUMPS ERROR: details:'
 write(stderr, *) 'Mumps Error: INFO(1)',mumps_par%INFO(1),'INFO(2)',mumps_par%INFO(2)
 write(stderr, *) 'Mumps Error: INFOG(1)',mumps_par%INFOG(1),'INFOG(2)',mumps_par%INFOG(2)
 error stop
endif


call cpu_time(tfin)
if (debug) print *, 'Solve took ',tfin-tstart,' seconds...'


!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
if ( myid == 0 ) then
 mumps_perm=mumps_par%SYM_PERM
 poisson2D=reshape(mumps_par%RHS/dx1**2,[lx1,lx2])

 deallocate( mumps_par%IRN )
 deallocate( mumps_par%JCN )
 deallocate( mumps_par%A   )
 deallocate( mumps_par%RHS )
end if

mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#else
call DMUMPS(mumps_par)
#endif

end function poisson2D


function elliptic3D_curv(srcterm,sig0,sigP,sigH,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                 x,flagdirich,perflag,it)

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
!-------ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
!-------LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
!-------GRID POINTS.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: srcterm,sig0,sigP,sigH
real(wp), dimension(:,:), intent(in) :: Vminx1,Vmaxx1
real(wp), dimension(:,:), intent(in) :: Vminx2,Vmaxx2
real(wp), dimension(:,:), intent(in) :: Vminx3,Vmaxx3
type(curvmesh), intent(in) :: x
integer, intent(in) :: flagdirich
logical, intent(in) :: perflag
integer, intent(in) :: it

real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigP2,gradsigP3
real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsigH2,gradsigH3
real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: gradsig01
real(wp), dimension(1:size(srcterm,1),1:size(srcterm,2),1:size(srcterm,3)) :: Ac,Bc,Cc,Dc,Ec,Fc

integer :: ix1,ix2,ix3,lx1,lx2,lx3
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
real(wp) :: tstart,tfin

#if REALBITS==32
type (SMUMPS_STRUC) mumps_par
#else
type (DMUMPS_STRUC) mumps_par
#endif

real(wp), dimension(size(srcterm,1),size(srcterm,2),size(srcterm,3)) :: elliptic3D_curv

!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
if (myid==0) then
 lx1=size(sig0,1)
 lx2=size(sig0,2)
 lx3=size(sig0,3)
 lPhi=lx1*lx2*lx3

 lent=7*(lx1-2)*(lx2-2)*(lx3-2)                                                 !interior entries
 lent=lent+2*(lx1-2)*(lx2-2)+2*(lx2-2)*(lx3-2)+2*(lx1-2)*(lx3-2)                !6 faces of cube
 lent=lent+4*(lx1-2)+4*(lx2-2)+4*(lx3-2)                                        !12 edges
 lent=lent+8                                                                    !8 corners, now we have total nonzero entries
 lent=lent+lx2*lx3                                                              !entries to deal with Neumann conditions on bottom
 if (flagdirich==0) then                                                        !more entries if Neumann on top
   lent=lent+lx2*lx3
 end if

 allocate(ir(lent),ic(lent),M(lent),b(lPhi))

 if (debug) print *, 'MUMPS will attempt a solve of size:  ',lx1,lx2,lx3
 if (debug) print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


 !COMPUTE AUXILIARY COEFFICIENTS
 if (debug) print *, 'Prepping coefficients for elliptic equation...'
 gradsig01=grad3D1(sig0,x,1,lx1,1,lx2,1,lx3)
 gradsigP2=grad3D2(sigP,x,1,lx1,1,lx2,1,lx3)
 gradsigP3=grad3D3(sigP,x,1,lx1,1,lx2,1,lx3)
 gradsigH2=grad3D2(sigH,x,1,lx1,1,lx2,1,lx3)
 gradsigH3=grad3D3(sigH,x,1,lx1,1,lx2,1,lx3)

 Ac=sigP; Bc=sigP; Cc=sig0;
 Dc=gradsigP2+gradsigH3
 Ec=gradsigP3-gradsigH2
 Fc=gradsig01


 !------------------------------------------------------------
 !-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
 !-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
 !-------GUIDE).
 !------------------------------------------------------------
 !LOAD UP MATRIX ELEMENTS
 M(:)=0d0
 b=pack(srcterm,.true.)           !boundaries overwritten later
 ient=1
 do ix3=1,lx3
   do ix2=1,lx2
   do ix1=1,lx1
     iPhi=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix3) as a column vector.  Also row of big matrix

     if (ix1==1) then          !BOTTOM GRID POINTS + CORNER, USE NEUMANN HERE, PRESUMABLY ZERO
       ir(ient)=iPhi
       ic(ient)=iPhi
       M(ient)=-1d0
       ient=ient+1
       ir(ient)=iPhi
       ic(ient)=iPhi+1
       M(ient)=1d0
!            b(iPhi)=Vminx1(ix3)
       b(iPhi)=0d0    !force bottom current to zero
       ient=ient+1
     elseif (ix1==lx1) then    !TOP GRID POINTS + CORNER
       if (flagdirich/=0) then
         ir(ient)=iPhi
         ic(ient)=iPhi
         M(ient)=1d0
         b(iPhi)=Vmaxx1(ix2,ix3)
         ient=ient+1
       else
         ir(ient)=iPhi
         ic(ient)=iPhi-1
         M(ient)=-1d0/x%dx1(lx1)
         ient=ient+1
         ir(ient)=iPhi
         ic(ient)=iPhi
         M(ient)=1d0/x%dx1(lx1)
         b(iPhi)=Vmaxx1(ix2,ix3)
         ient=ient+1
       end if
     elseif (ix2==1) then      !LEFT BOUNDARY
       ir(ient)=iPhi
       ic(ient)=iPhi
       M(ient)=1.0
       b(iPhi)=Vminx2(ix1,ix3)
       ient=ient+1
     elseif (ix2==lx2) then    !RIGHT BOUNDARY
       ir(ient)=iPhi
       ic(ient)=iPhi
       M(ient)=1.0
       b(iPhi)=Vmaxx2(ix1,ix3)
       ient=ient+1
     elseif (ix3==1) then
       ir(ient)=iPhi
       ic(ient)=iPhi
       M(ient)=1.0
       b(iPhi)=Vminx3(ix1,ix2)
       ient=ient+1
     elseif (ix3==lx3) then
       ir(ient)=iPhi
       ic(ient)=iPhi
       M(ient)=1.0
       b(iPhi)=Vmaxx3(ix1,ix2)
       ient=ient+1
     else                      !INTERIOR
       !ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
       ir(ient)=iPhi
       ic(ient)=iPhi-lx1*lx2
       M(ient)=Bc(ix1,ix2,ix3)/x%dx3all(ix3)/x%dx3iall(ix3)-Ec(ix1,ix2,ix3)/(x%dx3all(ix3+1)+x%dx3all(ix3))
       ient=ient+1

       !ix1,ix2-1,ix3
       ir(ient)=iPhi
       ic(ient)=iPhi-lx1
       M(ient)=Ac(ix1,ix2,ix3)/x%dx2(ix2)/x%dx2i(ix2)-Dc(ix1,ix2,ix3)/(x%dx2(ix2+1)+x%dx2(ix2))
       ient=ient+1

       !ix1-1,ix2,ix3
       ir(ient)=iPhi
       ic(ient)=iPhi-1
       M(ient)=Cc(ix1,ix2,ix3)/x%dx1(ix1)/x%dx1i(ix1)-Fc(ix1,ix2,ix3)/(x%dx1(ix1+1)+x%dx1(ix1))
       ient=ient+1

       !ix1,ix2,ix3
       ir(ient)=iPhi
       ic(ient)=iPhi
       M(ient)=-1d0*Ac(ix1,ix2,ix3)*(1d0/x%dx2(ix2+1)/x%dx2i(ix2)+1d0/x%dx2(ix2)/x%dx2i(ix2))- &
                    Bc(ix1,ix2,ix3)*(1d0/x%dx3all(ix3+1)/x%dx3iall(ix3)+1d0/x%dx3all(ix3)/x%dx3iall(ix3))- &
                    Cc(ix1,ix2,ix3)*(1d0/x%dx1(ix1+1)/x%dx1i(ix1)+1d0/x%dx1(ix1)/x%dx1i(ix1))
       ient=ient+1

       !ix1+1,ix2,ix3
       ir(ient)=iPhi
       ic(ient)=iPhi+1
       M(ient)=Cc(ix1,ix2,ix3)/x%dx1(ix1+1)/x%dx1i(ix1)+Fc(ix1,ix2,ix3)/(x%dx1(ix1+1)+x%dx1(ix1))
       ient=ient+1

       !ix1,ix2+1,ix3
       ir(ient)=iPhi
       ic(ient)=iPhi+lx1
       M(ient)=Ac(ix1,ix2,ix3)/x%dx2(ix2+1)/x%dx2i(ix2)+Dc(ix1,ix2,ix3)/(x%dx2(ix2+1)+x%dx2(ix2))
       ient=ient+1

       !ix1,ix2,ix3+1
       ir(ient)=iPhi
       ic(ient)=iPhi+lx1*lx2
       M(ient)=Bc(ix1,ix2,ix3)/x%dx3all(ix3+1)/x%dx3iall(ix3)+Ec(ix1,ix2,ix3)/(x%dx3all(ix3+1)+x%dx3all(ix3))
       ient=ient+1
     end if
   end do
 end do
 end do
end if
if (debug) print *, 'Number of entries used:  ',ient-1


! INIT MUMPS

mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

#if REALBITS==32
call SMUMPS(mumps_par)
#else
call DMUMPS(mumps_par)
#endif

call quiet_mumps(mumps_par)


!LOAD OUR PROBLEM
if ( myid==0 ) then
 mumps_par%N=lPhi
 mumps_par%NZ=lent
 allocate( mumps_par%IRN ( mumps_par%NZ ) )
 allocate( mumps_par%JCN ( mumps_par%NZ ) )
 allocate( mumps_par%A( mumps_par%NZ ) )
 allocate( mumps_par%RHS ( mumps_par%N  ) )
 mumps_par%IRN=ir
 mumps_par%JCN=ic
 mumps_par%A=M
 mumps_par%RHS=b
 deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

 if (perflag .and. it/=1) then       !used cached permutation
   allocate(mumps_par%PERM_IN(mumps_par%N))
   mumps_par%PERM_IN=mumps_perm
   mumps_par%ICNTL(7)=1
 end if
end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6
#if REALBITS==32
call SMUMPS(mumps_par)
#else
call DMUMPS(mumps_par)
#endif

!> check if Mumps error occurred
if (mumps_par%INFO(1) < 0 .or. mumps_par%INFOG(1) < 0) then
 write(stderr, *) 'Gemini:potential_mumps:elliptic3D_curv  MUMPS ERROR: details:'
 write(stderr, *) 'Mumps Error: INFO(1)',mumps_par%INFO(1),'INFO(2)',mumps_par%INFO(2)
 write(stderr, *) 'Mumps Error: INFOG(1)',mumps_par%INFOG(1),'INFOG(2)',mumps_par%INFOG(2)
 error stop
endif


!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
if ( myid==0 ) then
 if (debug) print *, 'Now organizing results...'

 if (perflag .and. it==1) then
   allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
   mumps_perm=mumps_par%SYM_PERM
 end if

 elliptic3D_curv=reshape(mumps_par%RHS,[lx1,lx2,lx3])

 if (debug) print *, 'Now attempting deallocations...'

 deallocate( mumps_par%IRN )
 deallocate( mumps_par%JCN )
 deallocate( mumps_par%A   )
 deallocate( mumps_par%RHS )
end if

mumps_par%JOB = -2
#if REALBITS==32
call SMUMPS(mumps_par)
#else
call DMUMPS(mumps_par)
#endif

end function elliptic3D_curv
