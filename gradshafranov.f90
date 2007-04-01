module gradshafranov

  implicit none
  
  real, allocatable :: psi(:), fun1(:), fun2(:), fun3(:), fun4(:)
  real :: psimin, psilim, dpsii
  real :: gamma2, gamma3, gamma4  

  integer, parameter :: numvargs = 1

contains

subroutine gradshafranov_init()

  use basic
  use arrays
  use diagnostics

  implicit none

  integer :: l, numnodes, ibegin, iendplusone
  real :: tstart, tend

  call numnod(numnodes)
  do l=1, numnodes
     call entdofs(numvar, l, 0, ibegin, iendplusone)

     call static_equ(ibegin)
     
     if(idens.eq.1) then
        call entdofs(1, l, 0, ibegin, iendplusone)
        call constant_field(den0(ibegin:ibegin+5), 1.)
     endif
  enddo

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call gradshafranov_solve()
  if(myrank.eq.0 .and. itimer.eq.1) then 
     call second(tend)
     t_gs = tend - tstart
  endif
  
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
  
  integer, parameter :: iterations = 80

  real   gsint1,gsint4,gsint2,gsint3,lhs,cfac(18)
  real, allocatable :: temp(:), b1vecini(:)
  integer, allocatable :: iboundgs(:)

  integer :: itri,i,i1,j,j1,jone, k
  integer :: numelms, numnodes
  integer :: ibegin, iendplusone, ibeginn, iendplusonen
  integer :: ineg, itnum, ier, nbcgs
  real :: dterm(18,18), sterm(18,18)
  real :: fac, aminor, bv, fintl(-6:maxi,-6:maxi)
  real :: g, gx, gz, gxx, gxz, gzz, g0
  real :: gv, gvx, gvz, gvxx, gvxz, gvzz
  real :: x,z, xmin, zmin, xrel, zrel, xguess, zguess, error
  real :: sum, rhs, ajlim, curr, q0, qstar, norm, rnorm
  real :: g1, g1x, g1z, g1xx, g1xz, g1zz
  real :: g2, g2x, g2z, g2xx, g2xz, g2zz
  real :: g3, g3x, g3z, g3xx, g3xz, g3zz
  real, dimension(5) :: temp1, temp2
  real :: alx, alz

  double precision :: coords(3)
   
  real :: tstart, tend, tsol, tmagaxis, tfundef, tplot


  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "gradshafranov called"

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart) ! t_gs_init

  t_gs_magaxis = 0.
  t_gs_solve = 0.
  t_gs_fundef = 0.

  call getmincoord(xmin, zmin)
  call numnod(numnodes)
  call numfac(numelms)

  ! allocate memory for arrays
  call createvec(temp, numvargs)
  call createvec(b1vecini, numvargs)
  call createvec(psi, numvargs)
  call createvec(fun1, numvargs)
  call createvec(fun4, numvargs)
  call createvec(fun2, numvargs)
  call createvec(fun3, numvargs)
  allocate(iboundgs(6*numnodes*numvargs))
  

  ! form the grad-sharfranov matrix
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call zeroarray4solve(gsmatrix_sm,numvar1_numbering)

  ! populate the matrix
  do itri=1,numelms

     ! calculate the local sampling points and weights for numerical integration
     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcr(itri, si_79, eta_79, 79, r_79)
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

           call insertval(gsmatrix_sm, sum, i1,j1,1)
        enddo
     enddo
  enddo

  ! insert boundary conditions
  iboundgs = 0.
  call boundarygs(iboundgs,nbcgs)
  do i=1,nbcgs
     call setdiribc(gsmatrix_sm, iboundgs(i))
  enddo
  call finalizearray4solve(gsmatrix_sm)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_gs_init = tend - tstart
  endif


  ! Define initial values of psi
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
     x = coords(1) - xmin + xzero
     z = coords(2) - zmin + zzero

     call entdofs(numvargs, i, 0, ibegin, iendplusone)

     call gvect1(x,z,xmag,zmag,g,gx,gz,gxx,gxz,gzz,0,ineg)
     call gvect1(x,z,102.,rnorm,gv,gvx,gvz,gvxx,gvxz,gvzz,1,ineg)

     psi(ibegin  ) = (g  +  gv*bv)*fac
     psi(ibegin+1) = (gx + gvx*bv)*fac
     psi(ibegin+2) = (gz + gvz*bv)*fac
     psi(ibegin+3) = (gxx+gvxx*bv)*fac
     psi(ibegin+4) = (gxz+gvxz*bv)*fac
     psi(ibegin+5) = (gzz+gvzz*bv)*fac
  enddo

  ! store boundary conditions on psi
  psibounds = 0.
  do i=1,nbcgs
     psibounds(i) = psi(iboundgs(i))
  enddo

  ! define initial b1vecini associated with delta-function source
  !     corresponding to current tcuro at location (xmag,zmag)
  xrel = xmag-xzero
  zrel = zmag-zzero

  b1vecini = 0
  call deltafun(xrel,zrel,b1vecini,tcuro)
 
  !-------------------------------------------------------------------
  ! start of iteration loop on plasma current
  do itnum=1, iterations

     if(myrank.eq.0 .and. iprint.eq.1) print *, "GS: iteration = ", itnum
     
     if(myrank.eq.0 .and. maxrank .eq. 1) call oneplot(b1vecini,1,numvargs,"b1vecini",0)

     ! apply boundary conditions
     do i=1,nbcgs
        b1vecini(iboundgs(i)) = psibounds(i)
     enddo

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
        psi = 0.5*(b1vecini + psi)
     endif

     if(myrank.eq.0 .and. maxrank .eq. 1) call oneplot(psi,1,numvargs,"psi ",0)
    
     ! Find new magnetic axis (extremum of psi)
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     xguess = xmag - xzero
     zguess = zmag - zzero    
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call magaxis(xguess,zguess)
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

        if(psi(ibegin).le.psilim) then
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
              cfac(j) = cfac(j) + gtri(k,j,itri)*fint(mi(k),ni(k))
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
     q0 = 1.
     qstar = 4.
     g0 = bzero*xzero
!!$     gamma2 = -2.*xmag*(xmag*p0*p1 + (2.*g0/(xmag**2*q0*dpsii)))
!!$     gamma3 = -(xmag*djdpsi/dpsii + 2*xmag**2*p0*p2)

!!$     gamma2 = gamma2 / 2.
!!$     gamma3 = gamma3 / 2.
!!$     
!!$     ! Nate:
     gamma2 =  - xmag*(xmag*p0*p1 + 2.*g0*(psilim-psimin)/(xmag**2*q0))
     gamma3 =  - (djdpsi*(psilim-psimin)/2. + p0*p2)

     gamma4 = -(tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4
     
!!$     if(myrank.eq.0 .and. iprint.gt.0) then
!!$        write(*,*) "gsint1, gsint2, gsint3, gsint4", gsint1, gsint2, gsint3, gsint4
!!$        write(*,*) "gamma2, gamma3, gamma4", gamma2, gamma3, gamma4
!!$     endif
     
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
     call sumshareddofs(b1vecini)

  end do ! on itnum



  ! populate phi0 array
  ! ~~~~~~~~~~~~~~~~~~~
  do i=1,numnodes
     !.....defines the source functions for the GS equation:
     call entdofs(numvargs, i, 0, ibegin, iendplusone)
     call entdofs(numvar, i, 0, ibeginn, iendplusonen)

     phi0(ibeginn:ibeginn+5) = psi(ibegin:ibegin+5)

!    I = sqrt(g0**2 + gamma_i*G_i)
     if(numvar.ge.2) then
        call xyznod(i,coords)
        x = coords(1) - xmin + xzero

        dpsii = 1./(psilim - psimin)
        temp(ibegin) = (psi(ibegin) - psimin)*dpsii
        temp(ibegin+1:ibegin+5) = psi(ibegin+1:ibegin+5)*dpsii
        if(temp(ibegin) .lt. 0. .or. temp(ibegin) .gt. 1.) then
           call constant_field(phi0(ibeginn+6 :ibeginn+11), g0)
        else
           g1  = temp(ibegin) - 10.*temp(ibegin)**3 + 20.*temp(ibegin)**4 &
                - 15.*temp(ibegin)**5 + 4.*temp(ibegin)**6
           g1x = temp(ibegin+1)*(1. - 30.*temp(ibegin)**2 + 80.*temp(ibegin)**3 &
                - 75.*temp(ibegin)**4 + 24.*temp(ibegin)**5)
           g1z = temp(ibegin+2)*(1. - 30.*temp(ibegin)**2 + 80.*temp(ibegin)**3 &
                - 75.*temp(ibegin)**4 + 24.*temp(ibegin)**5)
           g1xx= temp(ibegin+3)  *(1. - 30.*temp(ibegin)**2 + 80.*temp(ibegin)**3 &
                - 75.*temp(ibegin)**4 + 24.*temp(ibegin)**5) + &
                temp(ibegin+1)**2*(-60.*temp(ibegin)+240.*temp(ibegin)**2 &
                -300.*temp(ibegin)**3 +120.*temp(ibegin)**4)
           g1xz= temp(ibegin+4)  *(1. - 30.*temp(ibegin)**2 + 80.*temp(ibegin)**3 &
                - 75.*temp(ibegin)**4 + 24.*temp(ibegin)**5) + &
                temp(ibegin+1)*temp(ibegin+2)* &
                (-60.*temp(ibegin)   +240.*temp(ibegin)**2 &
                -300.*temp(ibegin)**3 +120.*temp(ibegin)**4)
           g1zz= temp(ibegin+5)  *(1. - 30.*temp(ibegin)**2 + 80.*temp(ibegin)**3 &
                - 75.*temp(ibegin)**4 + 24.*temp(ibegin)**5) + &
                temp(ibegin+2)**2*(-60.*temp(ibegin)+240.*temp(ibegin)**2 &
                -300.*temp(ibegin)**3 +120.*temp(ibegin)**4)

           g2  = temp(ibegin)**2 - 4.*temp(ibegin)**3 + 6.*temp(ibegin)**4 &
                - 4.*temp(ibegin)**5 + temp(ibegin)**6
           g2x = temp(ibegin+1)*(2.*temp(ibegin) - 12.*temp(ibegin)**2 &
                + 24.*temp(ibegin)**3 - 20.*temp(ibegin)**4 +  6.*temp(ibegin)**5)
           g2z = temp(ibegin+2)*(2.*temp(ibegin) - 12.*temp(ibegin)**2 &
                + 24.*temp(ibegin)**3 - 20.*temp(ibegin)**4 +  6.*temp(ibegin)**5)
           g2xx= temp(ibegin+3)*(2.*temp(ibegin) - 12.*temp(ibegin)**2 &
                + 24.*temp(ibegin)**3 - 20.*temp(ibegin)**4 +  6.*temp(ibegin)**5) + &
                temp(ibegin+1)**2*(2. - 24.*temp(ibegin) &
                + 72.*temp(ibegin)**2 - 80.*temp(ibegin)**3 + 30.*temp(ibegin)**4)
           g2xz= temp(ibegin+4)*(2.*temp(ibegin) - 12.*temp(ibegin)**2 &
                + 24.*temp(ibegin)**3 - 20.*temp(ibegin)**4 +  6.*temp(ibegin)**5) + &
                temp(ibegin+1)*temp(ibegin+2)*(2. - 24.*temp(ibegin) &
                + 72.*temp(ibegin)**2 - 80.*temp(ibegin)**3 + 30.*temp(ibegin)**4)
           g2zz= temp(ibegin+5)*(2.*temp(ibegin) - 12.*temp(ibegin)**2 &
                + 24.*temp(ibegin)**3 - 20.*temp(ibegin)**4 +  6.*temp(ibegin)**5) + &
                temp(ibegin+2)**2*(2. - 24.*temp(ibegin) &
                + 72.*temp(ibegin)**2 - 80.*temp(ibegin)**3 + 30.*temp(ibegin)**4)


           g3 = 1. - 20.*temp(ibegin)**3 + 45.*temp(ibegin)**4 &
                - 36.*temp(ibegin)**5 + 10.*temp(ibegin)**6
           g3x = temp(ibegin+1)*(-60.*temp(ibegin)**2 +180.*temp(ibegin)**3 &
                -180.*temp(ibegin)**4 + 60.*temp(ibegin)**5)
           g3z = temp(ibegin+2)*(-60.*temp(ibegin)**2 +180.*temp(ibegin)**3 &
                -180.*temp(ibegin)**4 + 60.*temp(ibegin)**5)
           g3xx= temp(ibegin+3)*(-60.*temp(ibegin)**2 +180.*temp(ibegin)**3 &
                -180.*temp(ibegin)**4 + 60.*temp(ibegin)**5) + &
                 temp(ibegin+1)**2*(-120.*temp(ibegin) +540.*temp(ibegin)**2 &
                -720.*temp(ibegin)**3 +300.*temp(ibegin)**4)
           g3xz= temp(ibegin+4)*(-60.*temp(ibegin)**2 +180.*temp(ibegin)**3 &
                -180.*temp(ibegin)**4 + 60.*temp(ibegin)**5) + &
                 temp(ibegin+1)*temp(ibegin+2)*(-120.*temp(ibegin) +540.*temp(ibegin)**2 &
                -720.*temp(ibegin)**3 +300.*temp(ibegin)**4)
           g3zz= temp(ibegin+5)*(-60.*temp(ibegin)**2 +180.*temp(ibegin)**3 &
                -180.*temp(ibegin)**4 + 60.*temp(ibegin)**5) + &
                 temp(ibegin+2)**2*(-120.*temp(ibegin) +540.*temp(ibegin)**2 &
                -720.*temp(ibegin)**3 +300.*temp(ibegin)**4)


           phi0(ibeginn+6) = sqrt(g0**2 + gamma2*g1 + gamma3*g2 + gamma4*g3)
           phi0(ibeginn+7) = 0.5*(gamma2*g1x + gamma3*g2x + gamma4*g3x) &
                / phi0(ibeginn+6)
           phi0(ibeginn+8) = 0.5*(gamma2*g1z + gamma3*g2z + gamma4*g3z) &
                / phi0(ibeginn+6)
           phi0(ibeginn+9) = 0.5*(gamma2*g1xx + gamma3*g2xx + gamma4*g3xx) &
                / phi0(ibeginn+6) - &
                (0.5*(gamma2*g1x + gamma3*g2x + gamma4*g3x))**2 &
                / phi0(ibeginn+6)**3
           phi0(ibeginn+10) = 0.5*(gamma2*g1xz + gamma3*g2xz + gamma4*g3xz) &
                / phi0(ibeginn+6) - &
                0.5*(gamma2*g1x + gamma3*g2x + gamma4*g3x)* &
                0.5*(gamma2*g1z + gamma3*g2z + gamma4*g3z)  &
                / phi0(ibeginn+6)**3
           phi0(ibeginn+11) = 0.5*(gamma2*g1zz + gamma3*g2zz + gamma4*g3zz) &
                / phi0(ibeginn+6) - &
                (0.5*(gamma2*g1z + gamma3*g2z + gamma4*g3z))**2 &
                / phi0(ibeginn+6)**3
        endif
     end if

     if(numvar.ge.3) then
        sum = p0 - pi0*ipres

        if(temp(ibegin) .lt. 0 .or. temp(ibegin) .gt. 1) then
           call constant_field(phi0(ibeginn+12:ibeginn+17), pedge * sum/p0)
        else
           phi0(ibeginn+12) = pedge * sum/p0 + sum* &
                (1.+p1*temp(ibegin)+p2*temp(ibegin)**2 &
                -(20. + 10.*p1 + 4.*p2)*temp(ibegin)**3 &
                +(45. + 20.*p1 + 6.*p2)*temp(ibegin)**4 &
                -(36. + 15.*p1 + 4.*p2)*temp(ibegin)**5 &
                +(10. +  4.*p1 +    p2)*temp(ibegin)**6)
           phi0(ibeginn+13) = sum*temp(ibegin+1)* &
                (p1+2.*p2*temp(ibegin) &
                -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                +6.*(10. +  4.*p1 +    p2)*temp(ibegin)**5)
           phi0(ibeginn+14) = sum*temp(ibegin+2)* &
                (p1+2.*p2*temp(ibegin) &
                -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                +6.*(10. +  4.*p1 +    p2)*temp(ibegin)**5)
           phi0(ibeginn+15) = sum*temp(ibegin+3)* &
                (p1+2.*p2*temp(ibegin) &
                -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                +6.*(10. +  4.*p1 +    p2)*temp(ibegin)**5) + &
                sum*temp(ibegin+1)**2* &
                (2.*p2 &
                - 6.*(20. + 10.*p1 + 4.*p2)*temp(ibegin) &
                +12.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**2 &
                -20.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**3 &
                +30.*(10. +  4.*p1 +    p2)*temp(ibegin)**4)
           phi0(ibeginn+16) = sum*temp(ibegin+4)* &
                (p1+2.*p2*temp(ibegin) &
                -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                +6.*(10. +  4.*p1 +    p2)*temp(ibegin)**5) + &
                sum*temp(ibegin+1)*temp(ibegin+2)* &
                (2.*p2 &
                - 6.*(20. + 10.*p1 + 4.*p2)*temp(ibegin) &
                +12.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**2 &
                -20.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**3 &
                +30.*(10. +  4.*p1 +    p2)*temp(ibegin)**4)
           phi0(ibeginn+17) = sum*temp(ibegin+5)* &
                (p1+2.*p2*temp(ibegin) &
                -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                +6.*(10. +  4.*p1 +    p2)*temp(ibegin)**5) + &
                sum*temp(ibegin+2)**2* &
                (2.*p2 &
                - 6.*(20. + 10.*p1 + 4.*p2)*temp(ibegin) &
                +12.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**2 &
                -20.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**3 &
                +30.*(10. +  4.*p1 +    p2)*temp(ibegin)**4)
        endif

        if(ipres.eq.1) then
           call entdofs(1, i, 0, ibeginn, iendplusonen)

           sum = p0

           if(temp(ibegin) .lt. 0 .or. temp(ibegin) .gt. 1) then
              call constant_field(pres0(ibeginn:ibeginn+5), pedge)
           else
              pres0(ibeginn) = pedge + sum* &
                   (1.+p1*temp(ibegin)+p2*temp(ibegin)**2 &
                   -(20. + 10.*p1 + 4.*p2)*temp(ibegin)**3 &
                   +(45. + 20.*p1 + 6.*p2)*temp(ibegin)**4 &
                   -(36. + 15.*p1 + 4.*p2)*temp(ibegin)**5 &
                   +(10. + 4.*p1 + p2)*temp(ibegin)**6)
              pres0(ibeginn+1) = sum*temp(ibegin+1)* &
                   (p1+2.*p2*temp(ibegin) &
                   -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                   +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                   -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                   +6.*(10. + 4.*p1 + p2)*temp(ibegin)**5)
              pres0(ibeginn+2) = sum*temp(ibegin+2)* &
                   (p1+2.*p2*temp(ibegin) &
                   -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                   +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                   -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                   +6.*(10. + 4.*p1 + p2)*temp(ibegin)**5)
              pres0(ibeginn+3) = sum*temp(ibegin+3)* &
                   (p1+2.*p2*temp(ibegin) &
                   -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                   +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                   -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                   +6.*(10. + 4.*p1 + p2)*temp(ibegin)**5) + &
                   sum*temp(ibegin+1)**2* &
                   (2.*p2 &
                   - 6.*(20. + 10.*p1 + 4.*p2)*temp(ibegin) &
                   +12.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**2 &
                   -20.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**3 &
                   +30.*(10. + 4.*p1 + p2)*temp(ibegin)**4)
              pres0(ibeginn+4) = sum*temp(ibegin+4)* &
                   (p1+2.*p2*temp(ibegin) &
                   -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                   +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                   -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                   +6.*(10. + 4.*p1 + p2)*temp(ibegin)**5) + &
                   sum*temp(ibegin+1)*temp(ibegin+2)* &
                   (2.*p2 &
                   - 6.*(20. + 10.*p1 + 4.*p2)*temp(ibegin) &
                   +12.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**2 &
                   -20.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**3 &
                   +30.*(10. + 4.*p1 + p2)*temp(ibegin)**4)
              pres0(ibeginn+5) = sum*temp(ibegin+5)* &
                   (p1+2.*p2*temp(ibegin) &
                   -3.*(20. + 10.*p1 + 4.*p2)*temp(ibegin)**2 &
                   +4.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**3 &
                   -5.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**4 &
                   +6.*(10. + 4.*p1 + p2)*temp(ibegin)**5) + &
                   sum*temp(ibegin+2)**2* &
                   (2.*p2 &
                   - 6.*(20. + 10.*p1 + 4.*p2)*temp(ibegin) &
                   +12.*(45. + 20.*p1 + 6.*p2)*temp(ibegin)**2 &
                   -20.*(36. + 15.*p1 + 4.*p2)*temp(ibegin)**3 &
                   +30.*(10. + 4.*p1 + p2)*temp(ibegin)**4)
           end if
        endif
     end if

  end do


  ! free memory
  call deletevec(temp)
  call deletevec(b1vecini)
  call deletevec(psi)
  call deletevec(fun1)
  call deletevec(fun4)
  call deletevec(fun2)
  call deletevec(fun3)
  deallocate(iboundgs)

  if(myrank.eq.0 .and. itimer.eq.1) then
     write(*,*) "gradshafranov: Time spent in magaxis: ", tmagaxis
     write(*,*) "gradshafranov: Time spent in fundef: ", tfundef
     write(*,*) "gradshafranov: Time spent in solve: ", tsol
     write(*,*) "gradshafranov: Time spent in plot: ", tplot
  endif

  return

end subroutine gradshafranov_solve


!==================================
subroutine magaxis(xguess,zguess)
  use basic
  use t_data
  use nintegrate_mod

  implicit none

  include 'mpif.h'

  real, intent(inout) :: xguess, zguess

  integer, parameter :: iterations = 5

  integer :: itri, itrit, itrinew, inews
  integer :: i, ier
  real :: x1, z1, x, z, theta, b, co, sn, si, eta
  real :: sum, sum1, sum2, sum3, sum4, sum5
  real :: term1, term2, term3, term4, term5
  real :: pt, pt1, pt2, p11, p22, p12
  real :: xnew, znew, denom, sinew, etanew
  real :: alx, alz
  real, dimension(20) :: avector
  real, dimension(3) :: temp1, temp2

  !     locates the magnetic axis and the value of psi there

  if(myrank.eq.0 .and. iprint.gt.0) &
       print *,  "magaxis: guess=", xguess, zguess

  call getboundingboxsize(alx, alz)

  itrit = 0
  inews = 0

  call whattri(xguess,zguess,itrit,x1,z1)

  itri = itrit
  x = xguess
  z = zguess
  
  do inews=1, iterations

     ! calculate position of minimum
     if(itri.gt.0) then
        call calcavector(itri, psi, 1, numvargs, avector)
         
        ! calculate local coordinates
        theta = ttri(itri)
        b = btri(itri)
        co = cos(theta)
        sn = sin(theta)
        si  = (x-x1)*co + (z-z1)*sn - b
        eta =-(x-x1)*sn + (z-z1)*co
  
        ! evaluate the polynomial and second derivative
        sum = 0.
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        sum5 = 0.
        do i=1,20
           sum = sum + avector(i)*si**mi(i)*eta**ni(i)
           term1 = 0.
           if(mi(i).ge.1) term1 = mi(i)*si**(mi(i)-1)*eta**ni(i)
           term2 = 0.
           if(ni(i).ge.1) term2 = ni(i)*si**mi(i)*eta**(ni(i)-1)
           term3 = 0.
           if(mi(i).ge.2) term3 = mi(i)*(mi(i)-1)*si**(mi(i)-2)*eta**ni(i)
           term4 = 0.
           if(ni(i).ge.2) term4 = ni(i)*(ni(i)-1)*si**mi(i)*eta**(ni(i)-2)
           term5 = 0.
           if(ni(i)*mi(i) .ge. 1)                                          &
                term5 = mi(i)*ni(i)*si**(mi(i)-1)*eta**(ni(i)-1)
           
           sum1 = sum1 + avector(i)*term1
           sum2 = sum2 + avector(i)*term2
           sum3 = sum3 + avector(i)*term3
           sum4 = sum4 + avector(i)*term4
           sum5 = sum5 + avector(i)*term5
        enddo
        pt  = sum
        pt1 = sum1
        pt2 = sum2
        p11 = sum3
        p22 = sum4
        p12 = sum5

        denom = p22*p11 - p12**2
        sinew = si -  ( p22*pt1 - p12*pt2)/denom
        etanew= eta - (-p12*pt1 + p11*pt2)/denom

        xnew = x1 + co*(b+sinew) - sn*etanew
        znew = z1 + sn*(b+sinew) + co*etanew
     else
        xnew = 0.
        znew = 0.
     endif  ! on itri.gt.0
     
     ! communicate new minimum to all processors
     if(maxrank.gt.1) then
        temp1(1) = xnew
        temp1(2) = znew
        call mpi_allreduce(temp1, temp2, 2, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
        xnew = temp2(1)
        znew = temp2(2)
     endif
     
     ! check to see whether the new minimum is outside the simulation domain
     if(xnew .lt. 0 .or. xnew.gt.alx .or. &
          znew .lt. 0 .or. znew.gt.alz .or. &
          xnew.ne.xnew .or. znew.ne.znew) then
        ! if not within the domain, safestop.

        write(*,3333) inews,x,z,xnew,znew
3333    format("magaxis: new minimum outside domain. ",i3,1p4e12.4)
        call safestop(27)       
     else
        x = xnew
        z = znew
     endif

     call whattri(xnew,znew,itrinew,x1,z1) 
     
     if(itrinew.le.0) then
        x = 0
        z = 0
        pt = 0
     endif
     itri = itrinew
  end do

  xguess = x
  zguess = z
  psimin = pt

  if(maxrank.gt.1) then
     temp1(1) = xguess
     temp1(2) = zguess
     temp1(3) = psimin
     call mpi_allreduce(temp1, temp2, 3, &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
     xguess = temp2(1)
     zguess = temp2(2)
     psimin = temp2(3)
  endif

  if(myrank.eq.0 .and. iprint.gt.0) print *, "magaxis: minimum at ", xguess, zguess
  
end subroutine magaxis

!============================================================
subroutine gvect1(r,z,xi,zi,g,gr,gz,grr,grz,gzz,nmult,ineg)
  integer, intent(in) :: nmult
  integer, intent(out) :: ineg
  real, intent(in) :: r, z, xi, zi
  real, intent(out) :: g,gr,gz,grr,grz,gzz

  real, dimension(1) :: rv, zv, xiv, ziv
  real, dimension(1) :: gv,grv,gzv,grrv,grzv,gzzv

  rv(1) = r
  zv(1) = z
  xiv(1) = xi
  ziv(1) = zi

  call gvect(rv,zv,xiv,ziv,1,gv,grv,gzv,grrv,grzv,gzzv,nmult,ineg)

  g = gv(1)
  gr = grv(1)
  gz = gzv(1)
  grr = grrv(1)
  grz = grzv(1)
  gzz = gzzv(1)

end subroutine gvect1


!============================================================
subroutine gvect(r,z,xi,zi,n,g,gr,gz,grr,grz,gzz,nmult,ineg)
  ! calculates derivatives wrt first argument

  implicit none

  integer, intent(in) :: n, nmult
  integer, intent(out) :: ineg
  real, dimension(n), intent(in) :: r, z, xi, zi
  real, dimension(n), intent(out) :: g,gr,gz,grr,grz,gzz

!!$  real, intent(in)  :: r, z, xi, zi
!!$  real, intent(out) :: g,gr,gz,grr,grz,gzz
  
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
     
     g(i) =- sqrxi*(2.*ck-2.*ce-ck*rksq)/rk
     gr(i)=-rk*0.25/sqrxi*(rpxi*term1                                  &
          +2.*xi(i)*(ce/x-ck))
     gz(i)=-rk*0.25*zmzi/sqrxi*term1
     grz(i)=0.0625*zmzi*(rk/sqrxi)**3*(rpxi*term1                      &
          +(ce-ck+2.*ce*rksq/x)*                                            &
          (term2)/x)
     gzz(i)=-rk*0.25/sqrxi*(term1*                                     &
          (1.-rksq*zmzi**2/(4.*rxi))+zmzi**2*rksq**2/(4.*rxi*x)             &
          *(ce-ck+2.*ce*rksq/x))
     grr(i)=-rk*0.25/sqrxi*(-rksq*rpxi/(4.*rxi)*                       &
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
     g(i) = tpi*rz**2
     gr(i) = 0.
     gz(i) = 0.
     grz(i) = 0.
     gzz(i) = 0.
     grr(i) = 0.
     go to 200
11   continue

     ! odd nullapole
     g(i) = 0.
     gr(i) = 0.
     gz(i) = 0.
     grz(i) = 0.
     gzz(i) = 0.
     grr(i) = 0.
     go to 200
12   continue

     ! even dipole
     g(i) = tpi*(r(i)**2 - rz**2)/2.
     gr(i) = tpi*r(i)
     gz(i) = 0.
     grz(i) = 0.
     gzz(i) = 0.
     grr(i) = tpi
     go to 200
13   continue
     
     ! odd dipole
     co=tpi/rz
     g(i) = co*(r(i)**2*z(i))
     gr(i) = co*(2.*r(i)*z(i))
     gz(i) = co*(r(i)**2)
     grz(i) = co*2.*r(i)
     gzz(i) = 0.
     grr(i) = co*2*z(i)
     go to 200
14   continue
        
     ! even quadrapole
     co=pi/(4.*rz**2)
     g(i) = co*(r(i)**4-4.*r(i)**2*z(i)**2 - 2.*r(i)**2*rz**2+rz**4)
     gr(i) = co*(4.*r(i)**3-8.*r(i)*z(i)**2-4.*r(i)*rz**2)
     gz(i) = co*(-8.*r(i)**2*z(i))
     grz(i) = co*(-16.*r(i)*z(i))
     gzz(i) = co*(-8.*r(i)**2)
     grr(i) = co*(12.*r(i)**2 - 8.*z(i)**2 - 4.*rz**2)
     go to 200
15   continue

     ! odd quadrapole
     co=pi/(3.*rz**3)
     g(i) = co*r(i)**2*z(i)*(3.*r(i)**2-4.*z(i)**2-3.*rz**2)
     gr(i) = co*(12.*r(i)**3*z(i)-8.*r(i)*z(i)**3-6.*r(i)*z(i)*rz**2)
     gz(i) = co*(3.*r(i)**4 - 12.*r(i)**2*z(i)**2 - 3.*r(i)**2*rz**2)
     grz(i) = co*(12.*r(i)**3-24.*r(i)*z(i)**2 - 6.*r(i)*rz**2)
     gzz(i) = co*(-24.*r(i)**2*z(i))
     grr(i) = co*(36.*r(i)**2*z(i)-8.*z(i)**3-6.*z(i)*rz**2)
     go to 200
16   continue

     ! even hexapole
     co=pi/(12.*rz**4)
     g(i) = co*(r(i)**6 - 12.*r(i)**4*z(i)**2 - 3.*r(i)**4*rz**2       &
          + 8.*r(i)**2*z(i)**4 + 12.*r(i)**2*z(i)**2*rz**2           &
          + 3.*r(i)**2*rz**4 - rz**6 )
     gr(i)= co*(6.*r(i)**5 - 48.*r(i)**3*z(i)**2                       &
          - 12.*r(i)**3*rz**2 + 16.*r(i)*z(i)**4                     &
          + 24.*r(i)*z(i)**2*rz**2 + 6.*r(i)*rz**4 )
     gz(i)= co*(-24.*r(i)**4*z(i) + 32.*r(i)**2*z(i)**3                &
          + 24.*r(i)**2*z(i)*rz**2)
     grz(i)=co*(-96.*r(i)**3*z(i)+64.*r(i)*z(i)**3                     &
          + 48.*r(i)*z(i)*rz**2 )
     grr(i)=co*(30.*r(i)**4-144.*r(i)**2*z(i)**2-36.*r(i)**2*rz**2     &
          + 16.*z(i)**4 + 24.*z(i)**2*rz**2 + 6.*rz**4)
     gzz(i)=co*(-24.*r(i)**4 + 96.*r(i)**2*z(i)**2 + 24.*r(i)**2*rz**2)
     go to 200
17   continue
        
     ! odd hexapole
     co=pi/(30.*rz**5)
     g(i) = co*(15.*r(i)**6*z(i) - 60.*r(i)**4*z(i)**3                 &
          - 30.*r(i)**4*z(i)*rz**2 + 24.*r(i)**2*z(i)**5             &
          + 40.*r(i)**2*z(i)**3*rz**2 + 15.*r(i)**2*z(i)*rz**4)
     gr(i)= co*(90.*r(i)**5*z(i) - 240.*r(i)**3*z(i)**3                &
          - 120.*r(i)**3*z(i)*rz**2 + 48.*r(i)*z(i)**5               &
          + 80.*r(i)*z(i)**3*rz**2 + 30.*r(i)*z(i)*rz**4)
     gz(i)= co*(15.*r(i)**6 - 180.*r(i)**4*z(i)**2                     &
          - 30.*r(i)**4*rz**2 + 120.*r(i)**2*z(i)**4                 &
          +120.*r(i)**2*z(i)**2*rz**2 + 15.*r(i)**2*rz**4)
     grz(i)=co*(90.*r(i)**5 - 720.*r(i)**3*z(i)**2                     &
          - 120.*r(i)**3*rz**2 + 240.*r(i)*z(i)**4                   &
          + 240.*r(i)*z(i)**2*rz**2 + 30.*r(i)*rz**4)
     gzz(i)=co*(-360.*r(i)**4*z(i) + 480.*r(i)**2*z(i)**3              &
          + 240.*r(i)**2*z(i)*rz**2)
     grr(i)=co*(450.*r(i)**4*z(i) - 720.*r(i)**2*z(i)**3               &
          - 360.*r(i)**2*z(i)*rz**2 + 48.*z(i)**5                    &
          + 80.*z(i)**3*rz**2 + 30.*z(i)*rz**4)
     go to 200
18   continue
        
     ! even octapole
     co=pi/(160.*rz**6)
     g(i) = co*(5.*r(i)**8 - 120.*r(i)**6*z(i)**2 - 20.*r(i)**6*rz**2  &
          + 240.*r(i)**4*z(i)**4 + 240.*r(i)**4*z(i)**2*rz**2        &
          + 30.*r(i)**4*rz**4 - 64.*r(i)**2*z(i)**6                  &
          - 160.*r(i)**2*z(i)**4*rz**2 - 120.*r(i)**2*z(i)**2*rz**4  &
          - 20.*r(i)**2*rz**6 + 5.*rz**8)
     gr(i)= co*(40.*r(i)**7 - 720.*r(i)**5*z(i)**2 - 120.*r(i)**5*rz**2 &
          + 960.*r(i)**3*z(i)**4 + 960.*r(i)**3*z(i)**2*rz**2        &
          + 120.*r(i)**3*rz**4 - 128.*r(i)*z(i)**6                   &
          - 320.*r(i)*z(i)**4*rz**2 - 240.*r(i)*z(i)**2*rz**4        &
          - 40.*r(i)*rz**6)
     gz(i)= co*(-240.*r(i)**6*z(i) + 960.*r(i)**4*z(i)**3              &
          + 480.*r(i)**4*z(i)*rz**2 - 384.*r(i)**2*z(i)**5           &
          - 640.*r(i)**2*z(i)**3*rz**2 - 240.*r(i)**2*z(i)*rz**4)
     grz(i)=co*(-1440.*r(i)**5*z(i) + 3840.*r(i)**3*z(i)**3            &
          + 1920.*r(i)**3*z(i)*rz**2 - 768.*r(i)*z(i)**5             &
          - 1280.*r(i)*z(i)**3*rz**2 - 480.*r(i)*z(i)*rz**4)
     gzz(i)=co*(-240.*r(i)**6 + 2880.*r(i)**4*z(i)**2                  &
          + 480.*r(i)**4*rz**2 - 1920.*r(i)**2*z(i)**4               &
          - 1920.*r(i)**2*z(i)**2*rz**2 - 240.*r(i)**2*rz**4)
     grr(i)=co*(280.*r(i)**6 - 3600.*r(i)**4*z(i)**2                   &
          - 600.*r(i)**4*rz**2 + 2880.*r(i)**2*z(i)**4               &
          + 2880.*r(i)**2*z(i)**2*rz**2                              &
          + 360.*r(i)**2*rz**4 - 128.*z(i)**6                        &
          - 320.*z(i)**4*rz**2 - 240.*z(i)**2*rz**4 - 40.*rz**6)
     go to 200
19   continue

     ! odd octapole
     co=pi/(140.*rz**7)
     g(i) = co*r(i)**2*z(i)*(35.*r(i)**6 - 280.*r(i)**4*z(i)**2        &
          - 105.*r(i)**4*rz**2 + 336.*r(i)**2*z(i)**4                &
          + 420.*r(i)**2*z(i)**2*rz**2 + 105.*r(i)**2*rz**4          &
          - 64.*z(i)**6 - 168.*z(i)**4*rz**2                         &
          - 140.*z(i)**2*rz**4 - 35.*rz**6)
     gr(i)= co*(280.*r(i)**7*z(i) - 1680.*r(i)**5*z(i)**3              &
          - 630.*r(i)**5*z(i)*rz**2 + 1344.*r(i)**3*z(i)**5          &
          + 1680.*r(i)**3*z(i)**3*rz**2 + 420.*r(i)**3*z(i)*rz**4    &
          - 128.*r(i)*z(i)**7 - 336.*r(i)*z(i)**5*rz**2              &
          - 280.*r(i)*z(i)**3*rz**4 - 70.*r(i)*z(i)*rz**6)
     gz(i)= co*(35.*r(i)**8-840.*r(i)**6*z(i)**2-105.*r(i)**6*rz**2    &
          + 1680.*r(i)**4*z(i)**4 + 1260.*r(i)**4*z(i)**2*rz**2      &
          + 105.*r(i)**4*rz**4 - 448.*r(i)**2*z(i)**6                &
          - 840.*r(i)**2*z(i)**4*rz**2 - 420.*r(i)**2*z(i)**2*rz**4  &
          - 35.*r(i)**2*rz**6)
     grz(i)=co*(280.*r(i)**7 - 5040.*r(i)**5*z(i)**2                   &
          - 630.*r(i)**5*rz**2 + 6720.*r(i)**3*z(i)**4               &
          + 5040.*r(i)**3*z(i)**2*rz**2 + 420.*r(i)**3*rz**4         &
          - 896.*r(i)*z(i)**6 - 1680.*r(i)*z(i)**4*rz**2             &
          - 840.*r(i)*z(i)**2*rz**4 - 70.*r(i)*rz**6)
     gzz(i)=co*(-1680.*r(i)**6*z(i) + 6720.*r(i)**4*z(i)**3            &
          + 2520.*r(i)**4*z(i)*rz**2 - 2688.*r(i)**2*z(i)**5         &
          - 3360.*r(i)**2*z(i)**3*rz**2 - 840.*r(i)**2*z(i)*rz**4)
     grr(i)=co*(1960.*r(i)**6*z(i) - 8400.*r(i)**4*z(i)**3             &
          - 3150.*r(i)**4*z(i)*rz**2 + 4032.*r(i)**2*z(i)**5         &
          + 5040.*r(i)**2*z(i)**3*rz**2 + 1260.*r(i)**2*z(i)*rz**4   &
          - 128.*z(i)**7 - 336.*z(i)**5*rz**2                        &
          - 280.*z(i)**3*rz**4 - 70.*z(i)*rz**6)
     go to 200
20   continue
        
     ! even decapole
     co=pi/(560.*rz**8)
     g(i) = co*(7.*r(i)**10 - 280.*r(i)**8*z(i)**2 - 35.*r(i)**8*rz**2 &
          + 1120.*r(i)**6*z(i)**4 + 840.*r(i)**6*z(i)**2*rz**2       &
          + 70.*r(i)**6*rz**4 - 896.*r(i)**4*z(i)**6                 &
          - 1680.*r(i)**4*z(i)**4*rz**2 - 840.*r(i)**4*z(i)**2*rz**4 &
          - 70.*r(i)**4*rz**6 + 128.*r(i)**2*z(i)**8                 &
          + 448.*r(i)**2*z(i)**6*rz**2 + 560.*r(i)**2*z(i)**4*rz**4  &
          + 280.*r(i)**2*z(i)**2*rz**6 + 35.*r(i)**2*rz**8           &
          - 7.*rz**10)
     gr(i)= co*(70.*r(i)**9 - 2240.*r(i)**7*z(i)**2                    &
          - 280.*r(i)**7*rz**2 + 6720.*r(i)**5*z(i)**4               &
          + 5040.*r(i)**5*z(i)**2*rz**2 + 420.*r(i)**5*rz**4         &
          - 3584.*r(i)**3*z(i)**6 - 6720.*r(i)**3*z(i)**4*rz**2      &
          - 3360.*r(i)**3*z(i)**2*rz**4 - 280.*r(i)**3*rz**6         &
          + 256.*r(i)*z(i)**8 + 896.*r(i)*z(i)**6*rz**2              &
          + 1120.*r(i)*z(i)**4*rz**4 + 560.*r(i)*z(i)**2*rz**6       &
          + 70.*r(i)*rz**8)
     gz(i)= co*(-560.*r(i)**8*z(i) + 4480.*r(i)**6*z(i)**3             &
          + 1680.*r(i)**6*z(i)*rz**2 - 5376.*r(i)**4*z(i)**5         &
          - 6720.*r(i)**4*z(i)**3*rz**2 - 1680.*r(i)**4*z(i)*rz**4   &
          + 1024.*r(i)**2*z(i)**7 + 2688.*r(i)**2*z(i)**5*rz**2      &
          + 2240.*r(i)**2*z(i)**3*rz**4 + 560.*r(i)**2*z(i)*rz**6)
     grz(i)=co*(-4480.*r(i)**7*z(i) + 26880.*r(i)**5*z(i)**3           &
          + 10080.*r(i)**5*z(i)*rz**2 - 21504.*r(i)**3*z(i)**5       &
          - 26880.*r(i)**3*z(i)**3*rz**2 - 6720.*r(i)**3*z(i)*rz**4  &
          + 2048.*r(i)*z(i)**7 + 5376.*r(i)*z(i)**5*rz**2            &
          + 4480.*r(i)*z(i)**3*rz**4 + 1120.*r(i)*z(i)*rz**6)
     gzz(i)=co*(-560.*r(i)**8 + 13440.*r(i)**6*z(i)**2                 &
          + 1680.*r(i)**6*rz**2 - 26880.*r(i)**4*z(i)**4             &
          - 20160.*r(i)**4*z(i)**2*rz**2 - 1680.*r(i)**4*rz**4       &
          + 7168.*r(i)**2*z(i)**6 + 13440.*r(i)**2*z(i)**4*rz**2     &
          + 6720.*r(i)**2*z(i)**2*rz**4 + 560.*r(i)**2*rz**6)
     grr(i)=co*(630.*r(i)**8 - 15680.*r(i)**6*z(i)**2                  &
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

!============================================================
subroutine boundarygs(ibound,nbc)

  use basic

  implicit none
  
  integer, intent(out) :: ibound(*), nbc

  integer :: numnodes, i, izone, izonedim
  integer :: ibottom, iright, itop, ileft, ibegin, iendplusone
  
  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)
  nbc = 0

  do i=1,numnodes

     call zonenod(i,izone, izonedim)
     call entdofs(numvargs, i, 0, ibegin, iendplusone)

     if(izonedim .eq. 1) then ! an edge...
        if(izone .eq. itop .or. izone .eq. ibottom) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+1
           ibound(nbc+3) = ibegin+3
           nbc =  nbc+3
        else
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+2
           ibound(nbc+3) = ibegin+5
           nbc =  nbc+3
        endif
     else if(izonedim .eq. 0) then
        ibound(nbc+1) = ibegin
        ibound(nbc+2) = ibegin+1
        ibound(nbc+3) = ibegin+2
        ibound(nbc+4) = ibegin+3
        ibound(nbc+5) = ibegin+4
        ibound(nbc+6) = ibegin+5
        nbc =  nbc+6
     endif
  enddo

end subroutine boundarygs

! ===========================================================
subroutine deltafun(x,z,dum,val)

  use t_data
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z, val
  real, intent(out) :: dum(*)

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

  call sumshareddofs(dum)

end subroutine deltafun
!============================================================
subroutine fundef

!.....defines the source functions for the GS equation:
  ! fun1 = r*p'
  ! fun4 = G1'/2r
  ! fun2 = G2'/2r
  ! fun3 = G3'/2r

  use basic
  
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

     if(pso .lt. 0. .or. pso .gt. 1.) then
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
             + 4.*(45.+20.*p1+6.*p2)*pso**3 - 5*(36.+15.*p1+4.*p2)*pso**4    &
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


end module gradshafranov
