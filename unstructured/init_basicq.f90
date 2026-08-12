module basicq
  use LZeqbm

  implicit none

  real, private :: q0_qp, rzero_qp, p0_qp, bz_qp, r0_qp, r1_qp, q2_qp, q4_qp, pedge_qp
  real, private :: q6_qp, q8_qp, q10_qp, q12_qp, q14_qp
  real, private :: kappa_qp, kappae_qp, coolrate_qp, v0_qp, v1_qp, beta_qp
  integer, private :: myrank_qp, iprint_qp, itaylor_qp

  integer, private, parameter :: nint=10000  !  (number of intervals)
  real, private, dimension (0:nint) :: bpsi, btor, bpolor, psi, jphi, jthor, gradpor, equor, pary, vary


contains

  subroutine init_qp
    use basic

    implicit none
    !input variables:
    bz_qp = bzero
    r0_qp = alpha0
    r1_qp = th_gs
    q0_qp = q0
    q2_qp = alpha1
    q4_qp = alpha2
    rzero_qp = rzero
    p0_qp = p0
    pedge_qp = pedge
    kappa_qp = kappa0
    kappae_qp = alpha3
    iprint_qp = iprint
    myrank_qp = myrank
    itaylor_qp = itaylor
    coolrate_qp = coolrate
    q6_qp = libetap
    q8_qp = p1
    q10_qp = p2
    q12_qp = djdpsi
    q14_qp = divcur
    v0_qp = v0_cyl
    v1_qp = v1_cyl
    beta_qp = beta
  end subroutine init_qp

  subroutine fixed_q_profiles()

    use basic
    use math
    use mesh_mod
    use sparse
    use arrays
    use m3dc1_nint
    use newvar_mod
    use boundary_conditions
    use model
    use gradshafranov
    use init_common
    implicit none

    vectype, dimension (dofs_per_element) :: dofsps, dofsbz, dofspr, dofsvz, dofsden
    real , dimension(MAX_PTS) :: rtemp79a, rtemp79b, rtemp79c, rtemp79d
    real :: r, dum1, dum2, dum3
    integer :: nelms, itri, j
    type (field_type) :: psi_vec, bz_vec, p_vec, vz_vec, den_vec

    call create_field(psi_vec)
    call create_field(bz_vec)
    call create_field(p_vec)
    call create_field(vz_vec)
    call create_field(den_vec)

    psi_vec = 0.
    bz_vec = 0.
    p_vec = 0.
    vz_vec = 0.
    den_vec = 0.
    
    if(myrank.eq.0 .and. iprint.ge.1) then 
       write (*,'(A,1pe12.4)') 'bz_qp =', bz_qp 
       write (*,'(A,1pe12.4)') 'r0_qp =', r0_qp 
       write (*,'(A,1pe12.4)') 'q0_qp =', q0_qp 
       write (*,'(A,1pe12.4)') 'q2_qp =', q2_qp 
       write (*,'(A,1pe12.4)') 'q4_qp =', q4_qp 
       write (*,'(A,1pe12.4)') 'rzero_qp =', rzero_qp 
       write (*,'(A,1pe12.4)') 'p0_qp =', p0_qp 
       write (*,'(A,1pe12.4)') 'pedge_qp =', pedge_qp 
       write (*,'(A,1pe12.4)') 'kappa_qp =', kappa_qp 
       write (*,'(A,1pe12.4)') 'kappae_qp =', kappae_qp 
       write (*,'(A,i5)') 'iprint_qp =', iprint_qp 
       write (*,'(A,i5)') 'myrank_qp =', myrank_qp 
       write (*,'(A,i5)') 'itaylor_qp =', itaylor_qp 
       write (*,'(A,1pe12.4)') 'coolrate_qp =', coolrate_qp 
       write (*,'(A,1pe12.4)') 'q6_qp =', q6_qp 
       write (*,'(A,1pe12.4)') 'q8_qp =', q8_qp 
       write (*,'(A,1pe12.4)') 'q10_qp =', q10_qp 
       write (*,'(A,1pe12.4)') 'q12_qp =', q12_qp 
       write (*,'(A,1pe12.4)') 'q14_qp =', q14_qp 
       write (*,'(A,1pe12.4)') 'rq_qp =', r1_qp
       write (*,'(A,1pe12.4)') 'v0_qp =', v0_qp 
       write (*,'(A,1pe12.4)') 'v1_qp =', v1_qp
       write (*,'(A,1pe12.4)') 'beta_qp =', beta_qp
  end if

    if(itaylor.eq.22) call setupLZeqbm

    if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "before loop over elements"

    call setup_qsolver

    nelms = local_elements()
!!$OMP PARALLEL DO &
!!$OMP& PRIVATE(dofsps,dofsbz,dofspr,dofsvz,dofsden,r,rtemp79a,rtemp79b,rtemp79c,rtemp79d,dum1,dum2,dum3,j)
    do itri=1,nelms
       
       call define_element_quadrature(itri,int_pts_diag, int_pts_tor)
       call define_fields(itri,0,1,0) ! defines x_79,z_79,mu,nu
       
       !  assemble matrix
       do j=1,npoints
          r = (sqrt((x_79(j)-xmag)**2 + (z_79(j)-zmag)**2))/r0_qp  ! normalized radius
          call getvals_qsolver(r,rtemp79a(j),rtemp79b(j),rtemp79c(j),rtemp79d(j))
       enddo
       
#ifdef USECOMPLEX
       temp79a = cmplx(rtemp79a)
       temp79b = cmplx(rtemp79b)
       temp79c = cmplx(rtemp79c)
       temp79d = cmplx(rtemp79d)
#else
       temp79a = rtemp79a
       temp79b = rtemp79b
       temp79c = rtemp79c
       temp79d = rtemp79d
#endif

       dofsps = intx2(mu79(:,:,OP_1),temp79a)
       dofsbz = intx2(mu79(:,:,OP_1),temp79b)
       dofspr = intx2(mu79(:,:,OP_1),temp79c)
       dofsvz = intx2(mu79(:,:,OP_1),temp79d)
       ! dofsden = den0*intx1(mu79(:,:,OP_1))
      temp79e = den0
      if(linear.eq.0) temp79e = den0*(temp79c/p0_qp)**expn + den_edge
      dofsden = intx2(mu79(:,:,OP_1),temp79e)

!!$OMP CRITICAL
       call vector_insert_block(psi_vec%vec,itri,1,dofsps,VEC_ADD)
       call vector_insert_block(bz_vec%vec ,itri,1,dofsbz,VEC_ADD)
       call vector_insert_block(p_vec%vec  ,itri,1,dofspr,VEC_ADD)
       call vector_insert_block(vz_vec%vec ,itri,1,dofsvz,VEC_ADD)
       call vector_insert_block(den_vec%vec,itri,1,dofsden,VEC_ADD)
!!$OMP END CRITICAL
    enddo
!!$OMP END PARALLEL DO
    
    ! solve for psi
    if(myrank.eq.0 .and. iprint.ge.1) print *, "solving psi"
    call newvar_solve(psi_vec%vec,mass_mat_lhs)
    if(myrank.eq.0 .and. iprint.ge.1) print *, "solving bz"
    call newvar_solve(bz_vec%vec ,mass_mat_lhs)
    if(myrank.eq.0 .and. iprint.ge.1) print *, "solving p"
    call newvar_solve(p_vec%vec  ,mass_mat_lhs)
    if(myrank.eq.0 .and. iprint.ge.1) print *, "solving vz"
    call newvar_solve(vz_vec%vec  ,mass_mat_lhs)
    if(myrank.eq.0 .and. iprint.ge.1) print *, "solving den"
    call newvar_solve(den_vec%vec  ,mass_mat_lhs)
    
    psi_field(0) = psi_vec
    bz_field(0)  = bz_vec
    p_field(0)   = p_vec
    vz_field(0)  = vz_vec
    den_field(0) = den_vec
    pe_field(0) = p_field(0)
    call mult(pe_field(0),pefac)
    
    call destroy_field(psi_vec)
    call destroy_field(bz_vec)
    call destroy_field(p_vec)
    call destroy_field(vz_vec)
    call destroy_field(den_vec)
    
    call init_perturbations
    
    if(itaylor.eq.27) then
       call getvals_qsolver(0.,psimin,dum1,dum2,dum3)
       call getvals_qsolver(q4_qp,psibound,dum1,dum2,dum3)
       if(myrank.eq.0) write(*,3001) psimin,psibound
3001   format("psimin,psibound = ", 1p2e12.4)
    endif
    
    if(myrank.eq.0 .and. iprint.ge.1) print *, "end fixed_q_profiles"
  end subroutine fixed_q_profiles

  subroutine setup_qsolver()
    implicit none
    
    real :: dpsi
    real :: psimid, qmid, qpmid, ppmid, denom,fterm,gterm,aquad,bquad,cquad,disc,A_qp
    integer :: j
    
    A_qp = rzero_qp/r0_qp
    
    !           small psi is normalized r**2
    dpsi = 1./nint
    do j=0,nint
       psi(j) = j*dpsi
       !  DEBUG
       !   if(iprint_qp .ge.1 .and. myrank_qp .eq.0) write(*,4000) j, psi(j), qfunc(psi(j))
4000   format('j   psi   qfunc(psi)', i5, 1p2e12.4)
    enddo
    
    !  boundary condition at edge
    btor(nint) = bz_qp
    bpsi(nint) = 0.
    bpolor(nint) = btor(nint)/(2.*A_qp*qfunc(psi(nint)))
    if(myrank_qp.eq.0 .and. iprint_qp.ge.1) write(*,3000) btor(nint),bpsi(nint),bpolor(nint)
3000 format( 'btor(nint), bpsi(nint), bpolor(nint) =',1p3e12.4)
    if(myrank_qp.eq.0 .and. iprint_qp.ge.1) write(*,*) "diagnostics to follow"
    !  integrate first order ode from boundary in
    do j=nint,1,-1
       psimid = (j-.5)*dpsi
       qmid = A_qp*qfunc(psimid)
       qpmid= A_qp*qpfunc(psimid)
       ppmid= ppfunc(psimid)
       denom = psimid + qmid**2
       fterm = -(1 -psimid*qpmid/qmid)/denom
       gterm = -ppmid*qmid**2/denom
       
       aquad = 1. + .5*dpsi*fterm
       bquad = dpsi*fterm*btor(j)
       cquad = -btor(j)**2*(1.-.5*dpsi*fterm) + 2.*dpsi*gterm
       
       disc = bquad**2 - 4.*aquad*cquad
       btor(j-1) = (-bquad + sqrt(disc))/(2.*aquad)
       bpolor(j-1) =btor(j-1)*r0_qp**2/(2.*rzero_qp*qfunc(psi(j-1)))
       bpsi(j-1)   = bpsi(j) - .5*dpsi*(bpolor(j)+bpolor(j-1)) 
!!$          if(myrank_qp.eq.0 .and. iprint_qp.ge.1) &
!!$               write(*,1002) qmid,denom,fterm,gterm,aquad,bquad,cquad,disc,btor(j-1)
    enddo
1002 format(1p9e9.1)
    
    !  calculate poloidal and toroidal fields in cell centers
    do j=1,nint-1
       jphi(j) = 4.*((j+.5)*(bpsi(j+1)-bpsi(j))-(j-.5)*(bpsi(j)-bpsi(j-1)))/(dpsi*r0_qp**2) 
       jthor(j) =  2*((bpsi(j+1)-bpsi(j))*rzero_qp*qfunc((j+0.5)*dpsi) &
            - (bpsi(j)-bpsi(j-1))*rzero_qp*qfunc((j-0.5)*dpsi))/(dpsi**2*r0_qp**2)
       gradpor(j) =  ppfunc(j*dpsi)
       pary(j) = pfunc(psi(j))
       vary(j) = vfunc(psi(j))
       !    error in equilibrium equation
       equor(j) = (jphi(j)*bpolor(j)+jthor(j)*btor(j)+gradpor(j))*sqrt(j*dpsi)
    enddo
    pary(0) =  pfunc(psi(0))
    pary(nint) =  pfunc(psi(nint))
    vary(0) =  vfunc(psi(0))
    vary(nint) =  vfunc(psi(nint))
    !
!!$       if(myrank_qp .eq. 0 .and. iprint_qp .ge. 2) then
!!$          write(6,1001)
!!$1001      format(" j       r**2       bpsi        btor         p            v         equil")
!!$          do j=0,nint
!!$             write(6,1000) j,psi(j),bpsi(j),btor(j),pary(j),vary(j), equor(j)
!!$          enddo
!!$1000      format(i4,1p8e12.4)
!!$       endif
 
  end subroutine setup_qsolver
  
  subroutine getvals_qsolver(rval,bpsival,ival,pval,vval)
    implicit none
    
    real, intent(in) :: rval
    real, intent(out) :: bpsival, ival, pval, vval

    real :: psival
    
    psival = rval*rval
    bpsival = cubicinterp(psival,psi,bpsi,nint)
    ival = cubicinterp(psival,psi,btor,nint)
    pval = cubicinterp(psival,psi,pary,nint)
    vval = cubicinterp(psival,psi,vary,nint)
  end subroutine getvals_qsolver

  real function cubicinterp(x,xary,yary,N)
    implicit none
    real, intent(in) :: x
    integer, intent(in) :: N
    real, intent(in), dimension(0:N) :: xary, yary
    integer :: i
    real :: xt,del,m1,m2,a,b,c,d

    xt = 0
    a = 0
    b = 0
    c = 0
    d = 0
    if      (x .le. xary(0))   then
       a = yary(0)
    else if (x .ge. xary(N))   then
       a = yary(N)
    else if (x .le. xary(1))   then
       xt = x - xary(0)
       del = xary(1) - xary(0)
       m2 =  (yary(2)-yary(0))/(xary(2)-xary(0))
       a = yary(0)
       b = 2.*(yary(1) - yary(0))/del    - m2
       c =   -(yary(1) - yary(0))/del**2 + m2/del
    else if (x .ge. xary(N-1)) then
       xt = x - xary(N-1)
       del = xary(N) - xary(N-1)
       m1 =  (yary(N)-yary(N-2))/(xary(N)-xary(N-2))
       a = yary(N-1)
       b = m1
       c = (yary(N) - yary(N-1) - m1*del)/del**2
    else
       
       do i=1,N-2
          if(x.ge.xary(i) .and. x.le.xary(i+1)) then
             xt = x - xary(i)
             del = xary(i+1) - xary(i)
             m1 = (yary(i+1)-yary(i-1))/(xary(i+1)-xary(i-1))
             m2 = (yary(i+2)-yary(i  ))/(xary(i+2)-xary(i  ))
             a = yary(i)
             b = m1
             c =  3.*(yary(i+1) - yary(i) - m1*del)/del**2 - (m2 - m1)/del
             d = -2.*(yary(i+1) - yary(i) - m1*del)/del**3 + (m2 - m1)/del**2
             exit
          endif
       enddo
       
    endif
    cubicinterp = a + b*xt + c*xt**2 + d*xt**3
  end function cubicinterp

  real function qfunc(psi)    !   q  (safety factor)
    implicit none

    real, intent(in) :: psi !  note:  psi = r**2
    real :: c0,c1,c2,c3,c4 
    real :: asq, bigA, bigB
    real :: ra0, x, xbtheta
    
    select case(itaylor_qp)
       
    case(21)
       c0 = 4.179343
       c1 = -0.080417
       c2=-8.659146
       c3 = 10.668674
       c4 = -4.108323
       qfunc = (q0_qp) + psi**2*(c0+c1*psi+c2*psi**2+c3*psi**3+c4*psi**4)
       
    case(22)
       qfunc = q_LZ(psi)
       
    case(25)
       qfunc = (q0_qp) + psi*(q2_qp + q4_qp*psi)
       
    case(26)
       qfunc = q0_qp*(1. + (psi/q2_qp)**q4_qp )**(1./q4_qp)
       
    case(27)
       !new coding
       asq = q4_qp**2
       bigA = (-2. + 3.*q0_qp/q2_qp)/asq 
       bigB = (1. -  2.*q0_qp/q2_qp)/asq**2
       if(psi .le. asq) then
          qfunc = q0_qp/(1. + bigA*psi + bigB*psi**2)
       else
          qfunc = q2_qp*psi/asq
       endif
       
    case(28)
       ra0 = q4_qp*abs(((q8_qp/q10_qp/q0_qp)**(q12_qp+q14_qp*q4_qp**2)-(1.,0.))**(-0.5/(q12_qp+q14_qp*q4_qp**2)))
       qfunc = (1+(psi/ra0**2)**(q12_qp+psi*q14_qp))**(1/(q12_qp+psi*q14_qp))*q0_qp*(1+q6_qp/exp((sqrt(psi)-r1_qp)**2/q2_qp**2))
       
    case(30)

       x = sqrt(psi)
       xbtheta=(atan(psi)/4.+(x/2.+(x-1)**2/4.-1./4.)/(4.*x+6.*(x-1)**2+4.*(x-1)**3 &
            +(x-1)**4-2.)-5.*(1-x)**(1./2.)+(70.*(1-x)**(3./2.))/3.-91.*(1-x)**(5./2.) &
            +260.*(1-x)**(7./2.)-(5005.*(1-x)**(9./2.))/9.+910.*(1-x)**(11./2.) &
            -1155.*(1-x)**(13./2.)+1144.*(1-x)**(15./2.)-(15015.*(1-x)**(17./2.))/17. &
            +(10010.*(1-x)**(19./2.))/19.-(715.*(1-x)**(21./2.))/3.+(1820.*(1-x)**(23./2.))/23. &
            -(91.*(1-x)**(25./2.))/5.+(70.*(1-x)**(27./2.))/27.-(5.*(1-x)**(29./2.))/29.)
       
       qfunc=-psi/xbtheta

    case(32,34)
    if(psi .lt. q2_qp) then
       qfunc = q0_qp
    else
       qfunc = q0_qp + q4_qp*(psi - q2_qp)**2
    endif
       
    end select
  end function qfunc

  real function qpfunc(psi)   !   derivative of q wrt psi
    implicit none

    real, intent(in) :: psi !  note:  psi=r^2
    real :: c0,c1,c2,c3,c4   
    real :: asq, bigA, bigB  
    real :: ra0, psis, x, xbtheta
    
    select case (itaylor_qp)
       
    case(21)
       c0 = 4.179343
       c1 = -0.080417
       c2=-8.659146
       c3 = 10.668674
       c4 = -4.108323
       qpfunc =  psi*(2.*c0+3.*c1*psi+4.*c2*psi**2+5.*c3*psi**3+6.*c4*psi**4)
       
    case(22)
       qpfunc = qprime_LZ(psi)
       
    case(25)
       qpfunc = (q2_qp + 2.*q4_qp*psi)
       
    case(26)
       qpfunc = q0_qp*(1. + (psi/q2_qp)**q4_qp )**((1.-q4_qp)/q4_qp)       &
            *(1./q2_qp)**q4_qp*q4_qp*psi**(q4_qp-1)
    case(27)
       !new coding
       asq = q4_qp**2
       bigA = (-2. + 3.*q0_qp/q2_qp)/asq 
       bigB = (1. -  2.*q0_qp/q2_qp)/asq**2
       if(psi .le. asq) then
          qpfunc = -q0_qp*(bigA + 2.*bigB*psi)/(1. + bigA*psi + bigB*psi**2)**2
       else
          qpfunc = q2_qp/asq
       endif
       
    case(28)
       psis = max(1.e-5,psi)
       ra0 = q4_qp*abs(((q8_qp/q10_qp/q0_qp)**(q12_qp+q14_qp*q4_qp**2)-(1.,0.))**(-0.5/(q12_qp+q14_qp*q4_qp**2)))
       qpfunc = (1+(psis/ra0**2)**(q12_qp+psis*q14_qp))**(1/(q12_qp+psis*q14_qp)) &
            *q0_qp*(-((log(1+(psis/ra0**2)**(q12_qp+psis*q14_qp))*q14_qp)/(q12_qp+psis*q14_qp)**2) &
            +((psis/ra0**2)**(q12_qp+psis*q14_qp)*(log(psis/ra0**2)*q14_qp &
            +(q12_qp+psis*q14_qp)/psis))/((1+(psis/ra0**2)**(q12_qp+psis*q14_qp))*(q12_qp+psis*q14_qp))) &
            *(1+q6_qp/exp((sqrt(psis)-r1_qp)**2/q2_qp**2)) &
            -((1+(psis/ra0**2)**(q12_qp+psis*q14_qp))**(1/(q12_qp+psis*q14_qp)) &
            *q0_qp*q6_qp*(sqrt(psis)-r1_qp))/(exp((sqrt(psis)-r1_qp)**2/q2_qp**2)*sqrt(psis)*q2_qp**2)
       
    case(30)
       x = sqrt(psi)
       xbtheta=(atan(psi)/4.+(x/2.+(x-1)**2/4.-1./4.)/(4.*x+6.*(x-1)**2+4.*(x-1)**3+(x-1)**4-2.) &
            -5.*(1-x)**(1./2.)+(70.*(1-x)**(3./2.))/3.-91.*(1-x)**(5./2.)+260.*(1-x)**(7./2.) &
            -(5005.*(1-x)**(9./2.))/9.+910.*(1-x)**(11./2.)-1155.*(1-x)**(13./2.) &
            +1144.*(1-x)**(15./2.)-(15015.*(1-x)**(17./2.))/17.+(10010.*(1-x)**(19./2.))/19. &
            -(715.*(1-x)**(21./2.))/3.+(1820.*(1-x)**(23./2.))/23.-(91.*(1-x)**(25./2.))/5. &
            +(70.*(1-x)**(27./2.))/27.-(5.*(1-x)**(29./2.))/29.)

       qpfunc= (2.*x*xbtheta-psi*(2.*x/(4.*(1+psi))+(-2.+6.*x-4.*x**2+6.*(x-1)**2+(x-1)**3 &
            -3.*(x-1)**4-(x-1)**5-6.*x*(x-1)**2+(x*(x-1)**4)/2.)/((4.*x+6.*(x-1)**2 &
            +4.*(x-1)**3+(x-1)**4-2.)**2)+(5.*(1-x)**(-1./2.))/2.-(70.*(1-x)**(1./2.))/2. &
            +(5.*91.*(1-x)**(3./2.))/2.-(7.*260.*(1-x)**(5./2.))/2.+(5005.*(1-x)**(7./2.))/2. &
            -(910.*11.*(1-x)**(9./2.))/2.+(1155.*13.*(1-x)**(11./2.))/2. &
            -(1144.*15.*(1-x)**(13./2.))/2.+(15015.*(1-x)**(15./2.))/2.-(10010.*(1-x)**(17./2.))/2. &
            +(715.*7.*(1-x)**(19./2.))/2.-(1820.*(1-x)**(21./2.))/2.+(91.*5.*(1-x)**(23./2.))/2. &
            -(70.*(1-x)**(25./2.))/2.+(5.*(1-x)**(27./2.))/2.))/(2.*x*xbtheta**2) 

    case(32,34)
       if(psi .lt. q2_qp) then
          qpfunc = 0
       else
          qpfunc = 2.*q4_qp*(psi - q2_qp)
       endif
    end select
  end function qpfunc

  real function pfunc(psi)    !   p  (pressure)
    implicit none

    real, intent(in) :: psi !  note:  psi=r^2
    real :: asq, bigA, bigB 
    select case(itaylor_qp)
       
    case(21,25,26)
       pfunc = p0_qp * (1. - 3.2*psi + 4.16*psi**2 - 2.56*psi**3 + 0.64*psi**4)
       
    case(22)
       pfunc = p_LZ(psi)
       
    case(27)
       !new coding
       asq = q4_qp**2
       bigA = (-4. + 6.*q0_qp/q2_qp)/asq 
       bigB = (3. -  6.*q0_qp/q2_qp)/asq**2
       if(psi .le. asq) then
          pfunc = p0_qp*(1+ bigA*psi + bigB*psi**2)**(2./3.) + pedge_qp
       else
          pfunc = pedge_qp
       endif

    case(32)
       pfunc = (1.-kappae_qp)*p0_qp*(1. - psi   )    &
             +     kappae_qp *p0_qp*(1. - psi**2)    & 
             + pedge_qp
    case(34)
       pfunc = p0_qp*(1. - psi )**2 + pedge_qp
    end select
    return
  end function pfunc
  
  real function ppfunc(psi)    !  derivative of p wrt psi
    implicit none
    
    real, intent(in) :: psi !  note:  psi=r^2
    real :: asq, bigA, bigB

    select case(itaylor_qp)
       
    case(21,25,26)
       ppfunc = p0_qp * (-3.2 + 8.32*psi - 7.68*psi**2 + 2.56*psi**3)
    case(22)
       ppfunc = pprime_LZ(psi)
    case(27)
       asq = q4_qp**2
       bigA = (-4. + 6.*q0_qp/q2_qp)/asq 
       bigB = ( 3. - 6.*q0_qp/q2_qp)/asq**2
       if(psi .lt. asq) then
          ppfunc = p0_qp*(2./3.)*(1+ bigA*psi + bigB*psi**2)**(-1./3.)    &
               *(bigA + 2.*bigB*psi)
       else
          ppfunc = 0.
       endif

    case(32)
       ppfunc = - (1.- kappae_qp)*p0_qp    &
                -    kappae_qp   *p0_qp*2.*psi
    case(34)
       ppfunc = -2.*p0_qp*(1.-psi)
       
    end select
  end function ppfunc
  real function vfunc(psi)
    implicit none

    real, intent(in) :: psi !  note:  psi=r^2
    
    vfunc = v0_qp + v1_qp*psi**beta_qp

    return
   end function vfunc

  real function get_kappa(psi)  ! thermal conductivity for itaylor=27, ikappafunc=12   
    implicit none
    
    real, intent(in) :: psi   !  note:  psi=r^2
    real :: asq, bigA, bigB, num1, num2, denom, jedge, psin
    
    psin = psi / r0_qp**2
    
    asq = q4_qp**2
    bigA = (1. - 1.5*q0_qp/q2_qp)/asq 
    bigB = (1. -  2.*q0_qp/q2_qp)/asq**2
    jedge = (pedge_qp/p0_qp)**1.5
    !   if(psi .gt. asq) psi = asq    !     temporary fix
    if(psin .le. asq) then
       num1 = 1. - 2*psin*bigA + psin**2*bigB + jedge
       num2 = (1. - 4.*psin*bigA + 3.*psin**2*bigB + jedge )**(1./3.)
       denom =  asq*(bigA - 1.5*psin*bigB)
       get_kappa = kappa_qp*num1*num2/denom
    else
       get_kappa = kappae_qp
    endif
  end function get_kappa
  
  function hsink_qp(psi)
    implicit none
    real :: psi, hsink_qp
    hsink_qp = coolrate_qp*pfunc(psi)

  end function hsink_qp


end module basicq
