module gradshafranov

  implicit none
  
  vectype, allocatable :: psi(:)
  real, allocatable :: fun1(:), fun2(:), fun3(:), fun4(:)
  real :: psimin, dpsii
  real :: gamma2, gamma3, gamma4  
  integer :: itnum

  integer, parameter :: numvargs = 1

  integer, parameter :: NSTX_coils = 158
  real, dimension(NSTX_coils) :: NSTX_r, NSTX_z, NSTX_I

  data NSTX_r &
       / 0.1750, 0.1750, 0.1750, 0.1750, 0.1750, & ! PF1aU
         0.1750, 0.1750, 0.1750, 0.1750, 0.1750, & ! PF1aL
         0.2902, 0.2902, 0.2902, 0.2902,         & ! PF1b
         0.3386, 0.3386, 0.3386, 0.3386,         &
         0.7255, 0.7739, 0.8223, 0.8707,         & ! PF2U
         0.7255, 0.7739, 0.8223, 0.8707,         &
         0.7255, 0.7739, 0.8223, 0.8707,         &
         0.7255, 0.7739, 0.8223, 0.8707,         & ! PF2L
         0.7255, 0.7739, 0.8223, 0.8707,         &
         0.7255, 0.7739, 0.8223, 0.8707,         &
         1.4511, 1.5478, 1.4995,                 & ! PF3U
         1.4511, 1.5478, 1.4995,                 &
         1.4511, 1.5478, 1.4995, 1.4027, 1.4027, 1.4027, &
         1.4511, 1.5478, 1.4995,                 & ! PF3L
         1.4511, 1.5478, 1.4995,                 &
         1.4511, 1.5478, 1.4995, 1.4027, 1.4027, 1.4027, &
         1.9832, 2.0315, 1.9832, 2.0315, 1.9832, 2.0315, & ! PF5U
         1.9832, 2.0315, 1.9832, 2.0315, 1.9832, 2.0315, & ! PF5L /
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &   ! OH
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, &
         0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240, 0.1240 /
         
  data NSTX_z &
       / 1.5000, 1.5500, 1.6000, 1.6500, 1.7000, & ! PF1aU
        -1.5000,-1.5500,-1.6000,-1.6500,-1.7000, & ! PF1aL
        -1.9000,-1.8500,-1.8000,-1.7500,         & ! PF1b
        -1.9000,-1.8500,-1.8000,-1.7500,         &
         1.8500, 1.8500, 1.8500, 1.8500,         & ! PF2U
         1.9000, 1.9000, 1.9000, 1.9000,         &
         1.9500, 1.9500, 1.9500, 1.9500,         &
        -1.8500,-1.8500,-1.8500,-1.8500,         & ! PF2L
        -1.9000,-1.9000,-1.9000,-1.9000,         &
        -1.9500,-1.9500,-1.9500,-1.9500,         &
         1.6500, 1.6500, 1.5500,         & ! PF3U
         1.5500, 1.5500, 1.6500,         &
         1.6000, 1.6000, 1.6000, 1.5500, 1.6000, 1.6500, &
        -1.6500,-1.6500,-1.5500,         & ! PF3L
        -1.5500,-1.5500,-1.6500,         &
        -1.6000,-1.6000,-1.6000,-1.5500,-1.6000,-1.6500, &
         0.6500, 0.6500, 0.6000, 0.6000, 0.5500, 0.5500, & ! PF5U
        -0.6500,-0.6500,-0.6000,-0.6000,-0.5500,-0.5500, & ! PF5L
         0.0000, 0.0500, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.3500, & ! OH
         0.4000, 0.4500, 0.5000, 0.5500, 0.6000, 0.6500, 0.7000, 0.7500, &
         0.8000, 0.8500, 0.9000, 0.9500, 1.0000, 1.0500, 1.1000, 1.1500, &
         1.2000, 1.2500, 1.3000, 1.3500, 1.4000, 1.4500, 1.5000, 1.5500, &
         1.6000, 1.6500, 1.7000, 1.7500, 1.8000, 1.8500, 1.9000, 1.9500, &
         2.0000,-0.0500,-0.1000,-0.1500,-0.2000,-0.2500,-0.3000,-0.3500, &
        -0.4000,-0.4500,-0.5000,-0.5500,-0.6000,-0.6500,-0.7000,-0.7500, &
        -0.8000,-0.8500,-0.9000,-0.9500,-1.0000,-1.0500,-1.1000,-1.1500, &
        -1.2000,-1.2500,-1.3000,-1.3500,-1.4000,-1.4500,-1.5000,-1.5500, &
        -1.6000,-1.6500,-1.7000,-1.7500,-1.8000,-1.8500,-1.9000,-1.9500 /

!!$  data NSTX_I &
!!$       / 0.0288, 0.0288, 0.0288, 0.0288, 0.0288, & ! PF1aU
!!$         0.0493, 0.0493, 0.0493, 0.0493, 0.0493, & ! PF1aL
!!$         0.0064, 0.0064, 0.0064, 0.0064,         & ! PF1b
!!$         0.0064, 0.0064, 0.0064, 0.0064,         &
!!$         3.9e-4, 3.9e-4, 3.9e-4, 3.9e-4,         & ! PF2U
!!$         3.9e-4, 3.9e-4, 3.9e-4, 3.9e-4,         &
!!$         3.9e-4, 3.9e-4, 3.9e-4, 3.9e-4,         &
!!$         0.0000, 0.0000, 0.0000, 0.0000,         & ! PF2L
!!$         0.0000, 0.0000, 0.0000, 0.0000,         &
!!$         0.0000, 0.0000, 0.0000, 0.0000,         &
!!$        -0.0147,-0.0147,-0.0147,                 & ! PF3U
!!$        -0.0147,-0.0147,-0.0147,                 & 
!!$        -0.0147,-0.0147,-0.0147,-0.0147,-0.0147,-0.0147, &         
!!$        -0.0205,-0.0205,-0.0205,                 & ! PF3L
!!$        -0.0205,-0.0205,-0.0205,                 &
!!$        -0.0205,-0.0205,-0.0205,-0.0205,-0.0205,-0.0205, &
!!$        -0.0520,-0.0520,-0.0520,-0.0520,-0.0520,-0.0520, & ! PF5U 
!!$        -0.0520,-0.0520,-0.0520,-0.0520,-0.0520,-0.0520 /  ! PF5L
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, & ! OH
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, &
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, &
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, &
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, &
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, &
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, &
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, &
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647, &
!!$!        -0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647,-0.2647 /

  data NSTX_I &
       / 0.0700, 0.0700, 0.0700, 0.0700, 0.0700, & ! PF1aU
         0.0700, 0.0700, 0.0700, 0.0700, 0.0700, & ! PF1aL
         0.0100, 0.0100, 0.0100, 0.0100,         & ! PF1b
         0.0100, 0.0100, 0.0100, 0.0100,         &
         0.0200, 0.0200, 0.0200, 0.0200,         & ! PF2U
         0.0200, 0.0200, 0.0200, 0.0200,         &
         0.0200, 0.0200, 0.0200, 0.0200,         &
         0.0200, 0.0200, 0.0200, 0.0200,         & ! PF2L
         0.0200, 0.0200, 0.0200, 0.0200,         &
         0.0200, 0.0200, 0.0200, 0.0200,         &
        -0.0150,-0.0150,-0.0150,                 & ! PF3U
        -0.0150,-0.0150,-0.0150,                 & 
        -0.0150,-0.0150,-0.0150,-0.0150,-0.0150,-0.0150, &         
        -0.0150,-0.0150,-0.0150,                 & ! PF3L
        -0.0150,-0.0150,-0.0150,                 &
        -0.0150,-0.0150,-0.0150,-0.0150,-0.0150,-0.0150, &
        -0.0350,-0.0350,-0.0350,-0.0350,-0.0350,-0.0350, & ! PF5U 
        -0.0350,-0.0350,-0.0350,-0.0350,-0.0350,-0.0350, & ! PF5L
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, & ! OH
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082, &
        -0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082,-0.0082 /


contains

subroutine gradshafranov_init()

  use basic
  use arrays
  use diagnostics

  implicit none

  integer :: l, numnodes, ibegin, iendplusone, ibegin1, iendplusone1
  real :: tstart, tend, alx, alz, xmin, zmin, x, z
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call gradshafranov_solve()
  if(myrank.eq.0 .and. itimer.eq.1) then 
     call second(tend)
     t_gs = tend - tstart
  endif

  call numnod(numnodes)
  do l=1, numnodes

     call entdofs(vecsize, l, 0, ibegin, iendplusone)
     call xyznod(l, coords)

     x = coords(1) - xmin
     z = coords(2) - zmin

     call assign_local_pointers(l)

     u0_l = 0.
     vz0_l = 0.
     chi0_l = 0.

     call random_per(x,z,23)
  enddo
  
end subroutine gradshafranov_init

!============================================================
subroutine gradshafranov_solve

  use t_data
  use basic
  use arrays
  use sparse
  use diagnostics
  use nintegrate_mod

  implicit none

  include 'mpif.h'
  
  real :: gsint1,gsint4,gsint2,gsint3,lhs,cfac(18)
  real, allocatable :: temp(:)
  vectype, allocatable :: b1vecini(:)

  integer, parameter :: maxcoils = NSTX_coils
  integer :: numcoils

  integer :: itri,i,i1,j,j1,jone, k
  integer :: numelms, numnodes
  integer :: ibegin, iendplusone
  integer :: ibegin1, iendplusone1
  integer :: ineg, ier
  real :: dterm(18,18), sterm(18,18)
  real :: fac, aminor, bv, fintl(-6:maxi,-6:maxi)
  real, dimension(6,maxcoils) :: g
  real, dimension(maxcoils) :: xp, zp, xc, zc
  real :: x, z, xmin, zmin, xrel, zrel, xguess, zguess, error
  real :: sum, rhs, ajlim, curr, norm, rnorm, g0
  real, dimension(5) :: temp1, temp2
  real :: alx, alz
  
  double precision :: coords(3)
   
  real :: tstart, tend


  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, "Calculating Grad-Shafranov Equilibrium"

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart) ! t_gs_init

  t_gs_magaxis = 0.
  t_gs_solve = 0.
  t_gs_fundef = 0.

  call getmincoord(xmin, zmin)
  call numnod(numnodes)
  call numfac(numelms)

  ! allocate memory for arrays
  call createrealvec(temp, numvargs)
  call createvec(b1vecini, numvargs)
  call createvec(psi, numvargs)
  call createrealvec(fun1, numvargs)
  call createrealvec(fun4, numvargs)
  call createrealvec(fun2, numvargs)
  call createrealvec(fun3, numvargs)  

  ! form the grad-sharfranov matrix
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " forming the GS matrix..."

  call zerosuperlumatrix(gsmatrix_sm, icomplex, numvar1_numbering)

  ! populate the matrix
  do itri=1,numelms

     ! calculate the local sampling points and weights 
     ! for numerical integration
     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcpos(itri, si_79, eta_79, 79, x_79, z_79)
     r_79 = x_79
     ri_79 = 1./r_79

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, &
             ttri(itri), ri_79, 79, g79(:,:,i))
     end do
    
     do j=1,18
        j1 = isval1(itri,j)
        do i=1,18
           i1 = isval1(itri,i)

           temp79a = -ri_79* &
                (g79(:,OP_DR,i)*g79(:,OP_DR,j) &
                +g79(:,OP_DZ,i)*g79(:,OP_DZ,j))
           sum = int1(temp79a,weight_79,79)

           call insertval2(gsmatrix_sm, sum, 0, i1,j1,1)
        enddo
     enddo
  enddo

  ! insert boundary conditions
  call boundary_gs(gsmatrix_sm, b1vecini)
  call finalizematrix(gsmatrix_sm)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_gs_init = tend - tstart
  endif


  ! Define initial values of psi
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(myrank.eq.0 .and. iprint.gt.0) &
       print *, " calculating boundary conditions..."

  ! based on filiment with current tcuro
  ! and vertical field of strength bv given by shafranov formula
  ! NOTE:  This formula assumes (li/2 + beta_P) = 1.2
  fac = tcuro/(2.*pi)
  ! minor radius
  aminor = abs(xmag-xlim)
  bv = (1./(4.*pi*xmag))*(alog(8.*xmag/aminor) - 1.5 + 1.2)
  call getboundingboxsize(alx,alz)
  rnorm = xzero + alx/2.
  psi = 0.

  do i=1,numnodes

     call xyznod(i,coords)
     xp = coords(1) - xmin + xzero
     zp = coords(2) - zmin + zzero

     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     call assign_local_pointers(i)

     ! Field due to plasma current
     xc(1) = xmag
     zc(1) = zmag
     call gvect(xp,zp,xc,zc,1,g,0,ineg)
     psi(ibegin:ibegin+5) =   g(:,1)*fac

     ! Field due to external coils
     select case(idevice)
     case(1) ! CDX-U
        numcoils = 4
        xc(1) = 0.846
        zc(1) = 0.360
        xc(2) = 0.846
        zc(2) =-0.360
        xc(3) = 0.381
        zc(3) = 0.802
        xc(4) = 0.381
        zc(4) =-0.802
        call gvect(xp,zp,xc,zc,numcoils,g,0,ineg)     
        g = -g*.2*fac
     case(2) ! NSTX
        numcoils = NSTX_coils
        xc(1:NSTX_coils) = NSTX_r
        zc(1:NSTX_coils) = NSTX_z
        call gvect(xp,zp,xc,zc,numcoils,g,0,ineg)     
        do k=1,NSTX_coils 
           g(:,k) = g(:,k)*fac*NSTX_I(k)
        enddo
     case default ! Generic
        numcoils = 1
        xc(1) = 102.
        zc(1) = rnorm
        call gvect(xp,zp,xc,zc,numcoils,g,1,ineg)     
        g = g*bv*fac
     end select

     
     do k=1,numcoils 
        psi(ibegin:ibegin+5) = psi(ibegin:ibegin+5) + g(:,k)
     end do

     ! Add fields from divertor coils
     if(divertors.ge.1) then
        xc = xdiv
        zc(1) = zdiv
        if(divertors.eq.2) zc(2) = -zdiv
        call gvect(xp,zp,xc,zc,divertors,g,0,ineg)
        do k=1,divertors
           psi(ibegin:ibegin+5) = psi(ibegin:ibegin+5) + fac*divcur*g(:,k)
        end do
     endif

     ! store boundary conditions on psi
     psis_l = psi(ibegin:ibegin+5)
  enddo

  ! define initial b1vecini associated with delta-function source
  !     corresponding to current tcuro at location (xmag,zmag)
  xrel = xmag-xzero
  zrel = zmag-zzero


  if(myrank.eq.0 .and. iprint.gt.0) print *, " initializing current..."

  b1vecini = 0
  call deltafun(xrel,zrel,b1vecini,tcuro)
 
  !-------------------------------------------------------------------
  ! start of iteration loop on plasma current
  do itnum=1, igs

     if(myrank.eq.0 .and. iprint.eq.1) print *, "GS: iteration = ", itnum
     
     ! apply boundary conditions
     call boundary_gs(0, b1vecini)

     ! perform LU backsubstitution to get psi solution
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call solve(gsmatrix_sm,b1vecini,ier)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_gs_solve = t_gs_solve + tend - tstart
     endif

     if(itnum.eq.1) then
        psi = b1vecini
     else
        psi = th_gs*b1vecini + (1.-th_gs)*psi
     endif

    
     ! Find new magnetic axis (extremum of psi)
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     xguess = xmag - xzero
     zguess = zmag - zzero    
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call magaxis(xguess,zguess,psi,numvargs,psimin)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_gs_magaxis = t_gs_magaxis + tend - tstart
     endif
     xmag = xguess + xzero
     zmag = zguess + zzero
     
     ! calculate psi at the limiter
     xrel = xlim - xzero
     zrel = zlim - zzero
     itri = 0.

     call evaluate(xrel,zrel,psilim,ajlim,psi,1,numvargs,itri)

     ! define the pressure and toroidal field functions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call fundef
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_gs_fundef = t_gs_fundef + tend - tstart
     endif


     ! Calculate error in new solution
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! (nodes shared between processes are currenly overcounted)
     sum = 0.
     norm = 0.
     do i=1,numnodes
        
        call xyznod(i,coords)
        x = coords(1) - xmin + xzero
        z = coords(2) - zmin + zzero
        
        call entdofs(numvargs, i, 0, ibegin, iendplusone)

        if(real(psi(ibegin)).le.psilim) then
           lhs = (psi(ibegin+3)-psi(ibegin+1)/x+psi(ibegin+5))/x
           rhs =  -(fun1(ibegin)+gamma2*fun2(ibegin)+                        &
                gamma3*fun3(ibegin)+gamma4*fun4(ibegin))
           sum = sum + (lhs-rhs)**2
           norm = norm + lhs**2
        endif
     enddo

     if(maxrank.gt.1) then
        temp1(1) = sum
        temp1(2) = norm
        call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, &
             MPI_SUM, MPI_COMM_WORLD, ier)
        error = sqrt(temp2(1)/temp2(2))
     else
        error = sqrt(sum/norm)
     endif
     
     if(myrank.eq.0 .and. iprint.gt.0) print *, "GS: error = ", error
     
     ! start of loop over triangles to compute integrals needed to keep
     !     total current and q_0 constant using gamma4, gamma2, gamma3
     gsint1 = 0.
     gsint4 = 0.
     gsint2 = 0.
     gsint3 = 0.
     curr = 0.
     do itri=1,numelms

        call calcfint(fintl,maxi,atri(itri),btri(itri),ctri(itri))
        call calcsterm(itri, sterm, fintl)

        do j=1,18
           cfac(j) = 0.
           do k=1,20
              cfac(j) = cfac(j) + gtri(k,j,itri)*fintl(mi(k),ni(k))
           enddo
        enddo

        do i=1,18
           i1 = isval1(itri,i)
           
           gsint1 = gsint1 + cfac(i)*fun1(i1)
           gsint4 = gsint4 + cfac(i)*fun4(i1)
           gsint2 = gsint2 + cfac(i)*fun2(i1)
           gsint3 = gsint3 + cfac(i)*fun3(i1)

           do j=1,18
              jone = isval1(itri,j)
              curr = curr + sterm(i,j)*psi(i1)*rinv(jone)
           enddo
          
        enddo
     enddo

     if(maxrank.gt.1) then
        temp1(1) = curr
        temp1(2) = gsint1
        temp1(3) = gsint2
        temp1(4) = gsint3
        temp1(5) = gsint4
        call mpi_allreduce(temp1, temp2, 5, MPI_DOUBLE_PRECISION, &
             MPI_SUM, MPI_COMM_WORLD, ier)
        curr   = temp2(1)
        gsint1 = temp2(2)
        gsint2 = temp2(3)
        gsint3 = temp2(4)
        gsint4 = temp2(5)
     end if

     if(myrank.eq.0 .and. iprint.ge.1) then 
        print *, "GS: curr = ", curr
     endif


     ! choose gamma2 to fix q0/qstar.  Note that there is an additional
     ! degree of freedom in gamma3.  Could be used to fix qprime(0)
     g0 = bzero*xzero

     if(numvar.eq.1) then
        gamma2 = 0.
        gamma3 = 0.
        gamma4 = 0.
     else
        gamma2 = -2.*xmag*(xmag*p0*p1 + (2.*g0/(xmag**2*q0*dpsii)))
        gamma3 = -(xmag*djdpsi/dpsii + 2.*xmag**2*p0*p2)
        gamma4 = -(tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4
     endif
         
     ! start loop over elements to define RHS vector
     b1vecini = 0.

     do itri=1,numelms

        call calcfint(fintl,maxi,atri(itri),btri(itri),ctri(itri))
        call calcdterm(itri, dterm, fintl)
        
        do i=1,18
           i1 = isval1(itri,i)
           sum = 0.
           
           do j=1,18
              j1 = isval1(itri,j)
              
              sum = sum - dterm(i,j)* &
                   (       fun1(j1) + gamma4*fun4(j1)               &
                   +gamma2*fun2(j1) + gamma3*fun3(j1))
           enddo
           
           b1vecini(i1) =  b1vecini(i1) + sum
        enddo
     enddo
     call sumsharedppplvecvals(b1vecini)

  end do ! on itnum


  ! populate phi0 array
  ! ~~~~~~~~~~~~~~~~~~~
  do i=1,numnodes
     !.....defines the source functions for the GS equation:
     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     call assign_local_pointers(i)

     psi0_l = psi(ibegin:ibegin+5)

     ! temp = Psi (normalized flux)
     dpsii = 1./(psilim - psimin)
     temp(ibegin) = (psi(ibegin) - psimin)*dpsii
     temp(ibegin+1:ibegin+5) = psi(ibegin+1:ibegin+5)*dpsii

     ! pp = p(Psi) in pp
     call calc_pressure(temp(ibegin:ibegin+5),p0_l)
     call calc_toroidal_field(temp(ibegin:ibegin+5), bz0_l)

     pe0_l = (1. - ipres*pi0/p0)*p0_l

     den0_l(1) = p0_l(1)**expn
     den0_l(2) = p0_l(1)**(expn-1.)*p0_l(2)*expn
     den0_l(3) = p0_l(1)**(expn-1.)*p0_l(3)*expn
     den0_l(4) = p0_l(1)**(expn-1.)*p0_l(4)*expn &
          + p0_l(1)**(expn-2.)*p0_l(2)**2.*expn*(expn-1.)
     den0_l(5) = p0_l(1)**(expn-1.)*p0_l(5)*expn & 
          + p0_l(1)**(expn-2.)*p0_l(2)*p0_l(3) &
          * expn*(expn-1.)
     den0_l(6) = p0_l(1)**(expn-1.)*p0_l(6)*expn &
          + p0_l(1)**(expn-2.)*p0_l(3)**2.*expn*(expn-1.)
     den0_l = den0_l/p0**expn
  end do

  ! correct for left-handedness
  psimin = -psimin
  psilim = -psilim

  ! free memory
  call deleterealvec(temp)
  call deletevec(b1vecini)
  call deletevec(psi)
  call deleterealvec(fun1)
  call deleterealvec(fun4)
  call deleterealvec(fun2)
  call deleterealvec(fun3)

  return

end subroutine gradshafranov_solve

!============================================================
subroutine gvect(r,z,xi,zi,n,g,nmult,ineg)
  ! calculates derivatives wrt first argument

  implicit none

  integer, intent(in) :: n, nmult
  integer, intent(out) :: ineg
  real, dimension(n), intent(in) :: r, z, xi, zi
  real, dimension(6,n), intent(out) :: g
  
  real :: a0,a1,a2,a3,a4
  real :: b0,b1,b2,b3,b4
  real :: c1,c2,c3,c4
  real :: d1,d2,d3,d4
  real :: pi,tpi

  real :: rpxi, rxi, zmzi, rk, ce, ck, term1, term2, rz, co
  real :: rksq, sqrxi, x

  integer :: i, imult, imultp

  data a0,a1,a2,a3,a4/1.38629436112,9.666344259e-2,                 &
       3.590092383e-2,3.742563713e-2,1.451196212e-2/
  data b0,b1,b2,b3,b4/.5,.12498593597,6.880248576e-2,               &
       3.328355346e-2,4.41787012e-3/
  data c1,c2,c3,c4/.44325141463,6.260601220e-2,                     &
       4.757383546e-2,1.736506451e-2/
  data d1,d2,d3,d4/.24998368310,9.200180037e-2,                     &
       4.069697526e-2,5.26449639e-3/
  data pi,tpi/3.1415926535,6.283185308/
  
  if(nmult.gt.0) go to 101
  do i=1,n
     rpxi=r(i)+xi(i)
     rxi=r(i)*xi(i)
     zmzi=z(i)-zi(i)
     rksq=4.*rxi/(rpxi**2+zmzi**2)
     rk=sqrt(rksq)
     sqrxi=sqrt(rxi)
     x=1.-rksq
     ce=1.+x*(c1+x*(c2+x*(c3+x*c4)))+                                  &
          x*(d1+x*(d2+x*(d3+x*d4)))*(-alog(x))
     ck=a0+x*(a1+x*(a2+x*(a3+x*a4)))+                                  &
          (b0+x*(b1+x*(b2+x*(b3+x*b4))))*(-alog(x))
     
     term1=2.*ck-2.*ce-ce*rksq/x
     term2=2.*xi(i)-rksq*rpxi
     
     g(1,i) =- sqrxi*(2.*ck-2.*ce-ck*rksq)/rk
     g(2,i)=-rk*0.25/sqrxi*(rpxi*term1                                  &
          +2.*xi(i)*(ce/x-ck))
     g(3,i)=-rk*0.25*zmzi/sqrxi*term1
     g(5,i)=0.0625*zmzi*(rk/sqrxi)**3*(rpxi*term1                      &
          +(ce-ck+2.*ce*rksq/x)*                                            &
          (term2)/x)
     g(6,i)=-rk*0.25/sqrxi*(term1*                                     &
          (1.-rksq*zmzi**2/(4.*rxi))+zmzi**2*rksq**2/(4.*rxi*x)             &
          *(ce-ck+2.*ce*rksq/x))
     g(4,i)=-rk*0.25/sqrxi*(-rksq*rpxi/(4.*rxi)*                       &
          (rpxi*term1+2.*xi(i)*(ce/x-ck))+term1-                            &
          rksq*rpxi/(4.*rxi*x)*(ce-ck+2.*ce*rksq/x)*                        &
          (term2)+rksq/(2.*r(i)*x)*(2.*ce/x-ck)*term2)
100 end do
 
  return

101 continue
     
     ! check for multipolar coils
  do i=1,n
     if(xi(i) .lt. 100.) go to 200
     rz = zi(i)
     imult = ifix(xi(i) - 100.)
     if(imult .lt. 0 .or. imult.gt.10) go to 250
     imultp = imult + 1
     go to(10,11,12,13,14,15,16,17,18,19,20),imultp
10   continue
        
     ! even nullapole
     g(1,i) = tpi*rz**2
     g(2,i) = 0.
     g(3,i) = 0.
     g(4,i) = 0.
     g(5,i) = 0.
     g(6,i) = 0.
     go to 200
11   continue

     ! odd nullapole
     g(1,i) = 0.
     g(2,i) = 0.
     g(3,i) = 0.
     g(4,i) = 0.
     g(5,i) = 0.
     g(6,i) = 0.
     go to 200
12   continue

     ! even dipole
     g(1,i) = tpi*(r(i)**2 - rz**2)/2.
     g(2,i) = tpi*r(i)
     g(3,i) = 0.
     g(5,i) = 0.
     g(6,i) = 0.
     g(4,i) = tpi
     go to 200
13   continue
     
     ! odd dipole
     co=tpi/rz
     g(1,i) = co*(r(i)**2*z(i))
     g(2,i) = co*(2.*r(i)*z(i))
     g(3,i) = co*(r(i)**2)
     g(5,i) = co*2.*r(i)
     g(6,i) = 0.
     g(4,i) = co*2*z(i)
     go to 200
14   continue
        
     ! even quadrapole
     co=pi/(4.*rz**2)
     g(1,i) = co*(r(i)**4-4.*r(i)**2*z(i)**2 - 2.*r(i)**2*rz**2+rz**4)
     g(2,i) = co*(4.*r(i)**3-8.*r(i)*z(i)**2-4.*r(i)*rz**2)
     g(3,i) = co*(-8.*r(i)**2*z(i))
     g(5,i) = co*(-16.*r(i)*z(i))
     g(6,i) = co*(-8.*r(i)**2)
     g(4,i) = co*(12.*r(i)**2 - 8.*z(i)**2 - 4.*rz**2)
     go to 200
15   continue

     ! odd quadrapole
     co=pi/(3.*rz**3)
     g(1,i) = co*r(i)**2*z(i)*(3.*r(i)**2-4.*z(i)**2-3.*rz**2)
     g(2,i) = co*(12.*r(i)**3*z(i)-8.*r(i)*z(i)**3-6.*r(i)*z(i)*rz**2)
     g(3,i) = co*(3.*r(i)**4 - 12.*r(i)**2*z(i)**2 - 3.*r(i)**2*rz**2)
     g(5,i) = co*(12.*r(i)**3-24.*r(i)*z(i)**2 - 6.*r(i)*rz**2)
     g(6,i) = co*(-24.*r(i)**2*z(i))
     g(4,i) = co*(36.*r(i)**2*z(i)-8.*z(i)**3-6.*z(i)*rz**2)
     go to 200
16   continue

     ! even hexapole
     co=pi/(12.*rz**4)
     g(1,i) = co*(r(i)**6 - 12.*r(i)**4*z(i)**2 - 3.*r(i)**4*rz**2       &
          + 8.*r(i)**2*z(i)**4 + 12.*r(i)**2*z(i)**2*rz**2           &
          + 3.*r(i)**2*rz**4 - rz**6 )
     g(2,i)= co*(6.*r(i)**5 - 48.*r(i)**3*z(i)**2                       &
          - 12.*r(i)**3*rz**2 + 16.*r(i)*z(i)**4                     &
          + 24.*r(i)*z(i)**2*rz**2 + 6.*r(i)*rz**4 )
     g(3,i)= co*(-24.*r(i)**4*z(i) + 32.*r(i)**2*z(i)**3                &
          + 24.*r(i)**2*z(i)*rz**2)
     g(5,i)=co*(-96.*r(i)**3*z(i)+64.*r(i)*z(i)**3                     &
          + 48.*r(i)*z(i)*rz**2 )
     g(4,i)=co*(30.*r(i)**4-144.*r(i)**2*z(i)**2-36.*r(i)**2*rz**2     &
          + 16.*z(i)**4 + 24.*z(i)**2*rz**2 + 6.*rz**4)
     g(6,i)=co*(-24.*r(i)**4 + 96.*r(i)**2*z(i)**2 + 24.*r(i)**2*rz**2)
     go to 200
17   continue
        
     ! odd hexapole
     co=pi/(30.*rz**5)
     g(1,i) = co*(15.*r(i)**6*z(i) - 60.*r(i)**4*z(i)**3                 &
          - 30.*r(i)**4*z(i)*rz**2 + 24.*r(i)**2*z(i)**5             &
          + 40.*r(i)**2*z(i)**3*rz**2 + 15.*r(i)**2*z(i)*rz**4)
     g(2,i)= co*(90.*r(i)**5*z(i) - 240.*r(i)**3*z(i)**3                &
          - 120.*r(i)**3*z(i)*rz**2 + 48.*r(i)*z(i)**5               &
          + 80.*r(i)*z(i)**3*rz**2 + 30.*r(i)*z(i)*rz**4)
     g(3,i)= co*(15.*r(i)**6 - 180.*r(i)**4*z(i)**2                     &
          - 30.*r(i)**4*rz**2 + 120.*r(i)**2*z(i)**4                 &
          +120.*r(i)**2*z(i)**2*rz**2 + 15.*r(i)**2*rz**4)
     g(5,i)=co*(90.*r(i)**5 - 720.*r(i)**3*z(i)**2                     &
          - 120.*r(i)**3*rz**2 + 240.*r(i)*z(i)**4                   &
          + 240.*r(i)*z(i)**2*rz**2 + 30.*r(i)*rz**4)
     g(6,i)=co*(-360.*r(i)**4*z(i) + 480.*r(i)**2*z(i)**3              &
          + 240.*r(i)**2*z(i)*rz**2)
     g(4,i)=co*(450.*r(i)**4*z(i) - 720.*r(i)**2*z(i)**3               &
          - 360.*r(i)**2*z(i)*rz**2 + 48.*z(i)**5                    &
          + 80.*z(i)**3*rz**2 + 30.*z(i)*rz**4)
     go to 200
18   continue
        
     ! even octapole
     co=pi/(160.*rz**6)
     g(1,i) = co*(5.*r(i)**8 - 120.*r(i)**6*z(i)**2 - 20.*r(i)**6*rz**2  &
          + 240.*r(i)**4*z(i)**4 + 240.*r(i)**4*z(i)**2*rz**2        &
          + 30.*r(i)**4*rz**4 - 64.*r(i)**2*z(i)**6                  &
          - 160.*r(i)**2*z(i)**4*rz**2 - 120.*r(i)**2*z(i)**2*rz**4  &
          - 20.*r(i)**2*rz**6 + 5.*rz**8)
     g(2,i)= co*(40.*r(i)**7 - 720.*r(i)**5*z(i)**2 - 120.*r(i)**5*rz**2 &
          + 960.*r(i)**3*z(i)**4 + 960.*r(i)**3*z(i)**2*rz**2        &
          + 120.*r(i)**3*rz**4 - 128.*r(i)*z(i)**6                   &
          - 320.*r(i)*z(i)**4*rz**2 - 240.*r(i)*z(i)**2*rz**4        &
          - 40.*r(i)*rz**6)
     g(3,i)= co*(-240.*r(i)**6*z(i) + 960.*r(i)**4*z(i)**3              &
          + 480.*r(i)**4*z(i)*rz**2 - 384.*r(i)**2*z(i)**5           &
          - 640.*r(i)**2*z(i)**3*rz**2 - 240.*r(i)**2*z(i)*rz**4)
     g(5,i)=co*(-1440.*r(i)**5*z(i) + 3840.*r(i)**3*z(i)**3            &
          + 1920.*r(i)**3*z(i)*rz**2 - 768.*r(i)*z(i)**5             &
          - 1280.*r(i)*z(i)**3*rz**2 - 480.*r(i)*z(i)*rz**4)
     g(6,i)=co*(-240.*r(i)**6 + 2880.*r(i)**4*z(i)**2                  &
          + 480.*r(i)**4*rz**2 - 1920.*r(i)**2*z(i)**4               &
          - 1920.*r(i)**2*z(i)**2*rz**2 - 240.*r(i)**2*rz**4)
     g(4,i)=co*(280.*r(i)**6 - 3600.*r(i)**4*z(i)**2                   &
          - 600.*r(i)**4*rz**2 + 2880.*r(i)**2*z(i)**4               &
          + 2880.*r(i)**2*z(i)**2*rz**2                              &
          + 360.*r(i)**2*rz**4 - 128.*z(i)**6                        &
          - 320.*z(i)**4*rz**2 - 240.*z(i)**2*rz**4 - 40.*rz**6)
     go to 200
19   continue

     ! odd octapole
     co=pi/(140.*rz**7)
     g(1,i) = co*r(i)**2*z(i)*(35.*r(i)**6 - 280.*r(i)**4*z(i)**2        &
          - 105.*r(i)**4*rz**2 + 336.*r(i)**2*z(i)**4                &
          + 420.*r(i)**2*z(i)**2*rz**2 + 105.*r(i)**2*rz**4          &
          - 64.*z(i)**6 - 168.*z(i)**4*rz**2                         &
          - 140.*z(i)**2*rz**4 - 35.*rz**6)
     g(2,i)= co*(280.*r(i)**7*z(i) - 1680.*r(i)**5*z(i)**3              &
          - 630.*r(i)**5*z(i)*rz**2 + 1344.*r(i)**3*z(i)**5          &
          + 1680.*r(i)**3*z(i)**3*rz**2 + 420.*r(i)**3*z(i)*rz**4    &
          - 128.*r(i)*z(i)**7 - 336.*r(i)*z(i)**5*rz**2              &
          - 280.*r(i)*z(i)**3*rz**4 - 70.*r(i)*z(i)*rz**6)
     g(3,i)= co*(35.*r(i)**8-840.*r(i)**6*z(i)**2-105.*r(i)**6*rz**2    &
          + 1680.*r(i)**4*z(i)**4 + 1260.*r(i)**4*z(i)**2*rz**2      &
          + 105.*r(i)**4*rz**4 - 448.*r(i)**2*z(i)**6                &
          - 840.*r(i)**2*z(i)**4*rz**2 - 420.*r(i)**2*z(i)**2*rz**4  &
          - 35.*r(i)**2*rz**6)
     g(5,i)=co*(280.*r(i)**7 - 5040.*r(i)**5*z(i)**2                   &
          - 630.*r(i)**5*rz**2 + 6720.*r(i)**3*z(i)**4               &
          + 5040.*r(i)**3*z(i)**2*rz**2 + 420.*r(i)**3*rz**4         &
          - 896.*r(i)*z(i)**6 - 1680.*r(i)*z(i)**4*rz**2             &
          - 840.*r(i)*z(i)**2*rz**4 - 70.*r(i)*rz**6)
     g(6,i)=co*(-1680.*r(i)**6*z(i) + 6720.*r(i)**4*z(i)**3            &
          + 2520.*r(i)**4*z(i)*rz**2 - 2688.*r(i)**2*z(i)**5         &
          - 3360.*r(i)**2*z(i)**3*rz**2 - 840.*r(i)**2*z(i)*rz**4)
     g(4,i)=co*(1960.*r(i)**6*z(i) - 8400.*r(i)**4*z(i)**3             &
          - 3150.*r(i)**4*z(i)*rz**2 + 4032.*r(i)**2*z(i)**5         &
          + 5040.*r(i)**2*z(i)**3*rz**2 + 1260.*r(i)**2*z(i)*rz**4   &
          - 128.*z(i)**7 - 336.*z(i)**5*rz**2                        &
          - 280.*z(i)**3*rz**4 - 70.*z(i)*rz**6)
     go to 200
20   continue
        
     ! even decapole
     co=pi/(560.*rz**8)
     g(1,i) = co*(7.*r(i)**10 - 280.*r(i)**8*z(i)**2 - 35.*r(i)**8*rz**2 &
          + 1120.*r(i)**6*z(i)**4 + 840.*r(i)**6*z(i)**2*rz**2       &
          + 70.*r(i)**6*rz**4 - 896.*r(i)**4*z(i)**6                 &
          - 1680.*r(i)**4*z(i)**4*rz**2 - 840.*r(i)**4*z(i)**2*rz**4 &
          - 70.*r(i)**4*rz**6 + 128.*r(i)**2*z(i)**8                 &
          + 448.*r(i)**2*z(i)**6*rz**2 + 560.*r(i)**2*z(i)**4*rz**4  &
          + 280.*r(i)**2*z(i)**2*rz**6 + 35.*r(i)**2*rz**8           &
          - 7.*rz**10)
     g(2,i)= co*(70.*r(i)**9 - 2240.*r(i)**7*z(i)**2                    &
          - 280.*r(i)**7*rz**2 + 6720.*r(i)**5*z(i)**4               &
          + 5040.*r(i)**5*z(i)**2*rz**2 + 420.*r(i)**5*rz**4         &
          - 3584.*r(i)**3*z(i)**6 - 6720.*r(i)**3*z(i)**4*rz**2      &
          - 3360.*r(i)**3*z(i)**2*rz**4 - 280.*r(i)**3*rz**6         &
          + 256.*r(i)*z(i)**8 + 896.*r(i)*z(i)**6*rz**2              &
          + 1120.*r(i)*z(i)**4*rz**4 + 560.*r(i)*z(i)**2*rz**6       &
          + 70.*r(i)*rz**8)
     g(3,i)= co*(-560.*r(i)**8*z(i) + 4480.*r(i)**6*z(i)**3             &
          + 1680.*r(i)**6*z(i)*rz**2 - 5376.*r(i)**4*z(i)**5         &
          - 6720.*r(i)**4*z(i)**3*rz**2 - 1680.*r(i)**4*z(i)*rz**4   &
          + 1024.*r(i)**2*z(i)**7 + 2688.*r(i)**2*z(i)**5*rz**2      &
          + 2240.*r(i)**2*z(i)**3*rz**4 + 560.*r(i)**2*z(i)*rz**6)
     g(5,i)=co*(-4480.*r(i)**7*z(i) + 26880.*r(i)**5*z(i)**3           &
          + 10080.*r(i)**5*z(i)*rz**2 - 21504.*r(i)**3*z(i)**5       &
          - 26880.*r(i)**3*z(i)**3*rz**2 - 6720.*r(i)**3*z(i)*rz**4  &
          + 2048.*r(i)*z(i)**7 + 5376.*r(i)*z(i)**5*rz**2            &
          + 4480.*r(i)*z(i)**3*rz**4 + 1120.*r(i)*z(i)*rz**6)
     g(6,i)=co*(-560.*r(i)**8 + 13440.*r(i)**6*z(i)**2                 &
          + 1680.*r(i)**6*rz**2 - 26880.*r(i)**4*z(i)**4             &
          - 20160.*r(i)**4*z(i)**2*rz**2 - 1680.*r(i)**4*rz**4       &
          + 7168.*r(i)**2*z(i)**6 + 13440.*r(i)**2*z(i)**4*rz**2     &
          + 6720.*r(i)**2*z(i)**2*rz**4 + 560.*r(i)**2*rz**6)
     g(4,i)=co*(630.*r(i)**8 - 15680.*r(i)**6*z(i)**2                  &
          - 1960.*r(i)**6*rz**2 + 33600*r(i)**4*z(i)**4              &
          + 25200.*r(i)**4*z(i)**2*rz**2 + 2100.*r(i)**4*rz**4       &
          - 10752.*r(i)**2*z(i)**6 - 20160.*r(i)**2*z(i)**4*rz**2    &
          - 10080.*r(i)**2*z(i)**2*rz**4 - 840.*r(i)**2*rz**6        &
          + 256.*z(i)**8 + 896.*z(i)**6*rz**2                        &
          + 1120.*z(i)**4*rz**4 + 560.*z(i)**2*rz**6                 &
          + 70.*rz**8)
     go to 200
200 end do

  return
250 continue
   
  ! error
  ineg=39
  
  return
end subroutine gvect

! ===========================================================
subroutine deltafun(x,z,dum,val)

  use t_data
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z, val
  vectype, intent(out) :: dum(*)

  integer :: itri, i, k, index
  real :: x1, z1, b, theta, si, eta, sum
  
  call whattri(x,z,itri,x1,z1)

  if(itri.gt.0) then

     ! calculate local coordinates
     theta = ttri(itri)
     b = btri(itri)
     si  = (x-x1)*cos(theta) + (z-z1)*sin(theta) - b
     eta =-(x-x1)*sin(theta) + (z-z1)*cos(theta)

     ! calculate the contribution to b1vecini
     do i=1,18
        index = isval1(itri,i)

        sum = 0.
        do k=1,20
           sum = sum + gtri(k,i,itri)*si**mi(k)*eta**ni(k)
        enddo
        dum(index) = dum(index) + sum*val
     enddo
  end if

  call sumsharedppplvecvals(dum)

end subroutine deltafun
!============================================================
subroutine fundef

!.....defines the source functions for the GS equation:
  ! fun1 = r*p'
  ! fun4 = G1'/2r
  ! fun2 = G2'/2r
  ! fun3 = G3'/2r

  use basic
  use diagnostics
  
  implicit none 
  integer :: l, numnodes, i, ibegin, iendplusone
  real :: x, z, xmin, zmin, pso, psox, psoy, psoxx, psoxy, psoyy, fbig
  real :: fbigp, fbigpp, g4big, g4bigp, g4bigpp, g2big, g2bigp
  real :: g2bigpp, g3big, g3bigp, g3bigpp
  double precision :: coords(3)

  call getmincoord(xmin, zmin)
  dpsii = 1./(psilim - psimin)

  call numnod(numnodes)
  do l=1,numnodes

     call xyznod(l,coords)
     x = coords(1) - xmin + xzero
     z = coords(2) - zmin + zzero

     call entdofs(numvargs, l, 0, ibegin, iendplusone)
     pso =  (psi(ibegin)-psimin)*dpsii

!     if(pso .lt. 0. .or. pso .gt. 1.) then
     if(pso .gt. 1.) then
        do i=0,5
           fun1(ibegin+i) = 0.
           fun4(ibegin+i) = 0.
           fun2(ibegin+i) = 0.
           fun3(ibegin+i) = 0.
        enddo
     else
        psox = psi(ibegin+1)*dpsii
        psoy = psi(ibegin+2)*dpsii
        psoxx= psi(ibegin+3)*dpsii
        psoxy= psi(ibegin+4)*dpsii
        psoyy= psi(ibegin+5)*dpsii
     
        fbig = p0*dpsii*(p1 + 2.*p2*pso - 3.*(20. + 10.*p1+4.*p2)*pso**2     &
             + 4.*(45.+20.*p1+6.*p2)*pso**3 - 5.*(36.+15.*p1+4.*p2)*pso**4   &
             + 6.*(10.+4.*p1+p2)*pso**5)
        fbigp = p0*dpsii*(2.*p2 - 6.*(20. + 10.*p1+4.*p2)*pso                &
             + 12.*(45.+20.*p1+6*p2)*pso**2 - 20.*(36.+15*p1+4*p2)*pso**3    &
             + 30.*(10.+4.*p1+p2)*pso**4)
        fbigpp= p0*dpsii*(- 6.*(20. + 10.*p1+4.*p2)                          &
             + 24.*(45.+20.*p1+6*p2)*pso - 60.*(36.+15*p1+4*p2)*pso**2       &
             + 120.*(10.+4.*p1+p2)*pso**3)

        fun1(ibegin)   = x*fbig
        fun1(ibegin+1) = fbig + x*fbigp*psox
        fun1(ibegin+2) =        x*fbigp*psoy
        fun1(ibegin+3) = 2.*fbigp*psox + x*(fbigpp*psox**2+fbigp*psoxx)
        fun1(ibegin+4) = fbigp*psoy + x*(fbigpp*psox*psoy +fbigp*psoxy)
        fun1(ibegin+5) = x*(fbigpp*psoy**2 + fbigp*psoyy)

!!$     g4big =   dpsii*(-pso**2+3.*pso**3-3.*pso**4+pso**5)
!!$     g4bigp =  dpsii*(-2*pso+9.*pso**2-12.*pso**3+5.*pso**4)
!!$     g4bigpp = dpsii*(-2 + 18.*pso-36.*pso**2+20*pso**3)
        g4big = dpsii*(-60*pso**2+180*pso**3-180*pso**4+60*pso**5)
        g4bigp= dpsii*(-120*pso+540*pso**2-720*pso**3+300*pso**4)
        g4bigpp=dpsii*(-120   +1080*pso  -2160*pso**2+1200*pso**3)

        fun4(ibegin)  = g4big/x
        fun4(ibegin+1)= g4bigp*psox/x - g4big/x**2
        fun4(ibegin+2)= g4bigp*psoy/x
        fun4(ibegin+3)= (g4bigpp*psox**2 + g4bigp*psoxx)/x                  &
             - 2.*g4bigp*psox/x**2 + 2.*g4big/x**3
        fun4(ibegin+4)= (g4bigpp*psox*psoy+g4bigp*psoxy)/x                  &
             - g4bigp*psoy/x**2
        fun4(ibegin+5)=  (g4bigpp*psoy**2 + g4bigp*psoyy)/x
        
        g2big =  dpsii*(1. - 30.*pso**2 + 80.*pso**3                     &
             - 75.*pso**4 + 24.*pso**5)
        g2bigp =  dpsii*(-60.*pso + 240.*pso**2                         &
             - 300.*pso**3 + 120.*pso**4)
        g2bigpp =  dpsii*(-60. + 480.*pso                               &
             - 900.*pso**2 + 480.*pso**3)
        fun2(ibegin)  =  g2big/x
        fun2(ibegin+1)=  g2bigp*psox/x - g2big/x**2
        fun2(ibegin+2)=  g2bigp*psoy/x
        fun2(ibegin+3)=  (g2bigpp*psox**2 + g2bigp*psoxx)/x                 &
             - 2.*g2bigp*psox/x**2 + 2.*g2big/x**3
        fun2(ibegin+4)=(g2bigpp*psox*psoy+g2bigp*psoxy)/x                   &
             - g2bigp*psoy/x**2
        fun2(ibegin+5)= (g2bigpp*psoy**2 + g2bigp*psoyy)/x
        
        g3big =  dpsii*(2.*pso - 12.*pso**2 + 24.*pso**3                &
             - 20.*pso**4 + 6.*pso**5)
        g3bigp =  dpsii*(2. - 24.*pso + 72.*pso**2                      &
             - 80.*pso**3 + 30.*pso**4)
        g3bigpp =  dpsii*(- 24. + 144.*pso                              &
             - 240.*pso**2 + 120.*pso**3)
        fun3(ibegin)= g3big/x
        fun3(ibegin+1)= g3bigp*psox/x - g3big/x**2
        fun3(ibegin+2)= g3bigp*psoy/x
        fun3(ibegin+3)= (g3bigpp*psox**2 + g3bigp*psoxx)/x                  &
             - 2.*g3bigp*psox/x**2 + 2.*g3big/x**3
        fun3(ibegin+4)= (g3bigpp*psox*psoy+g3bigp*psoxy)/x                  &
             - g3bigp*psoy/x**2
        fun3(ibegin+5)=  (g3bigpp*psoy**2 + g3bigp*psoyy)/x
     endif
  enddo

  fun2 = fun2 / 2.
  fun3 = fun3 / 2.
  fun4 = fun4 / 2.
  
  return
end subroutine fundef


!======================================================================
! calc_toroidal_field
! ~~~~~~~~~~~~~~~~~~~
!
! calculates the toroidal field (I) as a function of the 
! normalized flux
!======================================================================
subroutine calc_toroidal_field(psii,tf)
  use basic

  real, intent(in), dimension(6)  :: psii     ! normalized flux
  vectype, intent(out), dimension(6) :: tf    ! toroidal field (I)

  real :: g1, g1x, g1z, g1xx, g1xz, g1zz
  real :: g2, g2x, g2z, g2xx, g2xz, g2zz
  real :: g3, g3x, g3z, g3xx, g3xz, g3zz
  
!  if(psii(1) .lt. 0. .or. psii(1) .gt. 1.) then
  if(psii(1) .gt. 1.) then
     call constant_field(tf, bzero*xzero)
  else
     g1  = psii(1) - 10.*psii(1)**3 + 20.*psii(1)**4 &
          - 15.*psii(1)**5 + 4.*psii(1)**6
     g1x = psii(2)*(1. - 30.*psii(1)**2 + 80.*psii(1)**3 &
          - 75.*psii(1)**4 + 24.*psii(1)**5)
     g1z = psii(3)*(1. - 30.*psii(1)**2 + 80.*psii(1)**3 &
          - 75.*psii(1)**4 + 24.*psii(1)**5)
     g1xx= psii(4)  *(1. - 30.*psii(1)**2 + 80.*psii(1)**3 &
          - 75.*psii(1)**4 + 24.*psii(1)**5) + &
          psii(1+1)**2*(-60.*psii(1)+240.*psii(1)**2 &
          -300.*psii(1)**3 +120.*psii(1)**4)
     g1xz= psii(5)  *(1. - 30.*psii(1)**2 + 80.*psii(1)**3 &
          - 75.*psii(1)**4 + 24.*psii(1)**5) + &
          psii(2)*psii(3)* &
          (-60.*psii(1)   +240.*psii(1)**2 &
          -300.*psii(1)**3 +120.*psii(1)**4)
     g1zz= psii(6)  *(1. - 30.*psii(1)**2 + 80.*psii(1)**3 &
          - 75.*psii(1)**4 + 24.*psii(1)**5) + &
          psii(3)**2*(-60.*psii(1)+240.*psii(1)**2 &
          -300.*psii(1)**3 +120.*psii(1)**4)
     
     g2  = psii(1)**2 - 4.*psii(1)**3 + 6.*psii(1)**4 &
          - 4.*psii(1)**5 + psii(1)**6
     g2x = psii(2)*(2.*psii(1) - 12.*psii(1)**2 &
          + 24.*psii(1)**3 - 20.*psii(1)**4 +  6.*psii(1)**5)
     g2z = psii(3)*(2.*psii(1) - 12.*psii(1)**2 &
          + 24.*psii(1)**3 - 20.*psii(1)**4 +  6.*psii(1)**5)
     g2xx= psii(4)*(2.*psii(1) - 12.*psii(1)**2 &
          + 24.*psii(1)**3 - 20.*psii(1)**4 +  6.*psii(1)**5) + &
          psii(2)**2*(2. - 24.*psii(1) &
          + 72.*psii(1)**2 - 80.*psii(1)**3 + 30.*psii(1)**4)
     g2xz= psii(5)*(2.*psii(1) - 12.*psii(1)**2 &
          + 24.*psii(1)**3 - 20.*psii(1)**4 +  6.*psii(1)**5) + &
          psii(2)*psii(3)*(2. - 24.*psii(1) &
          + 72.*psii(1)**2 - 80.*psii(1)**3 + 30.*psii(1)**4)
     g2zz= psii(6)*(2.*psii(1) - 12.*psii(1)**2 &
          + 24.*psii(1)**3 - 20.*psii(1)**4 +  6.*psii(1)**5) + &
          psii(3)**2*(2. - 24.*psii(1) &
          + 72.*psii(1)**2 - 80.*psii(1)**3 + 30.*psii(1)**4)
     
     g3 = 1. - 20.*psii(1)**3 + 45.*psii(1)**4 &
                - 36.*psii(1)**5 + 10.*psii(1)**6
     g3x = psii(2)*(-60.*psii(1)**2 +180.*psii(1)**3 &
          -180.*psii(1)**4 + 60.*psii(1)**5)
     g3z = psii(3)*(-60.*psii(1)**2 +180.*psii(1)**3 &
          -180.*psii(1)**4 + 60.*psii(1)**5)
     g3xx= psii(4)*(-60.*psii(1)**2 +180.*psii(1)**3 &
          -180.*psii(1)**4 + 60.*psii(1)**5) + &
          psii(2)**2*(-120.*psii(1) +540.*psii(1)**2 &
          -720.*psii(1)**3 +300.*psii(1)**4)
     g3xz= psii(5)*(-60.*psii(1)**2 +180.*psii(1)**3 &
          -180.*psii(1)**4 + 60.*psii(1)**5) + &
          psii(2)*psii(3)*(-120.*psii(1) +540.*psii(1)**2 &
          -720.*psii(1)**3 +300.*psii(1)**4)
     g3zz= psii(6)*(-60.*psii(1)**2 +180.*psii(1)**3 &
          -180.*psii(1)**4 + 60.*psii(1)**5) + &
          psii(3)**2*(-120.*psii(1) +540.*psii(1)**2 &
          -720.*psii(1)**3 +300.*psii(1)**4)

     
     tf(1) = sqrt((bzero*xzero)**2 + gamma2*g1 + gamma3*g2 + gamma4*g3)
     tf(2) = 0.5*(gamma2*g1x + gamma3*g2x + gamma4*g3x) / tf(1)
     tf(3) = 0.5*(gamma2*g1z + gamma3*g2z + gamma4*g3z) / tf(1)
     tf(4) = 0.5*(gamma2*g1xx + gamma3*g2xx + gamma4*g3xx) / tf(1) &
          - (0.5*(gamma2*g1x + gamma3*g2x + gamma4*g3x))**2 / tf(1)**3
     tf(5) = 0.5*(gamma2*g1xz + gamma3*g2xz + gamma4*g3xz) / tf(1) &
          -  0.5*(gamma2*g1x + gamma3*g2x + gamma4*g3x) &
          *  0.5*(gamma2*g1z + gamma3*g2z + gamma4*g3z) / tf(1)**3
     tf(6) = 0.5*(gamma2*g1zz + gamma3*g2zz + gamma4*g3zz) / tf(1) &
          - (0.5*(gamma2*g1z + gamma3*g2z + gamma4*g3z))**2 / tf(1)**3
  endif
end subroutine calc_toroidal_field




!======================================================================
! calc_pressure
! ~~~~~~~~~~~~~
!
! calculates the pressure as a function of the normalized flux
!======================================================================
subroutine calc_pressure(psii,pres)
  
  use basic

  real, intent(in), dimension(6)  :: psii     ! normalized flux
  vectype, intent(out), dimension(6) :: pres     ! pressure

!  if(psii(1) .lt. 0 .or. psii(1) .gt. 1) then
  if(psii(1) .gt. 1) then
     pres = 0.
  else
     pres(1) = 1.+p1*psii(1)+p2*psii(1)**2 &
          -(20. + 10.*p1 + 4.*p2)*psii(1)**3 &
          +(45. + 20.*p1 + 6.*p2)*psii(1)**4 &
          -(36. + 15.*p1 + 4.*p2)*psii(1)**5 &
          +(10. +  4.*p1 +    p2)*psii(1)**6
     pres(2) = psii(2)* &
          (p1+2.*p2*psii(1) &
          -3.*(20. + 10.*p1 + 4.*p2)*psii(1)**2 &
          +4.*(45. + 20.*p1 + 6.*p2)*psii(1)**3 &
          -5.*(36. + 15.*p1 + 4.*p2)*psii(1)**4 &
          +6.*(10. +  4.*p1 +    p2)*psii(1)**5)
     pres(3) = psii(3)* &
          (p1+2.*p2*psii(1) &
          -3.*(20. + 10.*p1 + 4.*p2)*psii(1)**2 &
          +4.*(45. + 20.*p1 + 6.*p2)*psii(1)**3 &
          -5.*(36. + 15.*p1 + 4.*p2)*psii(1)**4 &
          +6.*(10. +  4.*p1 +    p2)*psii(1)**5)
     pres(4) = psii(4)* &
          (p1+2.*p2*psii(1) &
          -3.*(20. + 10.*p1 + 4.*p2)*psii(1)**2 &
          +4.*(45. + 20.*p1 + 6.*p2)*psii(1)**3 &
          -5.*(36. + 15.*p1 + 4.*p2)*psii(1)**4 &
          +6.*(10. +  4.*p1 +    p2)*psii(1)**5) + &
          psii(2)**2* &
          (2.*p2 &
          - 6.*(20. + 10.*p1 + 4.*p2)*psii(1) &
          +12.*(45. + 20.*p1 + 6.*p2)*psii(1)**2 &
          -20.*(36. + 15.*p1 + 4.*p2)*psii(1)**3 &
          +30.*(10. +  4.*p1 +    p2)*psii(1)**4)
     pres(5) = psii(5)* &
          (p1+2.*p2*psii(1) &
          -3.*(20. + 10.*p1 + 4.*p2)*psii(1)**2 &
          +4.*(45. + 20.*p1 + 6.*p2)*psii(1)**3 &
          -5.*(36. + 15.*p1 + 4.*p2)*psii(1)**4 &
          +6.*(10. +  4.*p1 +    p2)*psii(1)**5) + &
          psii(2)*psii(3)* &
          (2.*p2 &
          - 6.*(20. + 10.*p1 + 4.*p2)*psii(1) &
          +12.*(45. + 20.*p1 + 6.*p2)*psii(1)**2 &
          -20.*(36. + 15.*p1 + 4.*p2)*psii(1)**3 &
          +30.*(10. +  4.*p1 +    p2)*psii(1)**4)
     pres(6) = psii(6)* &
          (p1+2.*p2*psii(1) &
          -3.*(20. + 10.*p1 + 4.*p2)*psii(1)**2 &
          +4.*(45. + 20.*p1 + 6.*p2)*psii(1)**3 &
          -5.*(36. + 15.*p1 + 4.*p2)*psii(1)**4 &
          +6.*(10. +  4.*p1 +    p2)*psii(1)**5) + &
          psii(3)**2* &
          (2.*p2 &
          - 6.*(20. + 10.*p1 + 4.*p2)*psii(1) &
          +12.*(45. + 20.*p1 + 6.*p2)*psii(1)**2 &
          -20.*(36. + 15.*p1 + 4.*p2)*psii(1)**3 &
          +30.*(10. +  4.*p1 +    p2)*psii(1)**4)
  endif

  pres = p0*pres
  pres(1) = pres(1) + pedge

end subroutine calc_pressure

end module gradshafranov
