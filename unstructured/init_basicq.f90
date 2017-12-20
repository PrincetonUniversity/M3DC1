module basicq
  implicit none

  real, private :: q0_qp, rzero_qp, p0_qp, bz_qp, r0_qp, r1_qp, q2_qp, q4_qp, pedge_qp
  real, private :: q6_qp, q8_qp, q10_qp, q12_qp, q14_qp
  real, private :: kappa_qp, kappae_qp, coolrate_qp, v0_qp, v1_qp
  integer, private :: myrank_qp, iprint_qp, itaylor_qp

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
    real , dimension(npoints) :: rtemp79a, rtemp79b, rtemp79c, rtemp79d
    real :: r, dum1, dum2, dum3
    integer :: nelms, itri, i, j
    type (field_type) :: psi_vec, bz_vec, p_vec, vz_vec, den_vec

    call create_field(psi_vec)
    call create_field(bz_vec)
    call create_field(p_vec)
    call create_field(vz_vec)
    call create_field(den_vec)
    
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2000) bz_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2001) r0_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2002) q0_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2003) q2_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2004) q4_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2005) rzero_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2006) p0_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2007) pedge_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2008) kappa_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2009) kappae_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2010) iprint_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2011) myrank_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2012) itaylor_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2013) coolrate_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2014) q6_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2015) q8_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2016) q10_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2017) q12_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2018) q14_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2019) r1_qp
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2020) v0_qp 
    if(myrank.eq.0 .and. iprint.ge.1) write (*,2021) v1_qp

2000 format( 'bz_qp =', 1pe12.4)
2001 format( 'r0_qp =', 1pe12.4)
2002 format( 'q0_qp =', 1pe12.4)
2003 format( 'q2_qp =', 1pe12.4)
2004 format( 'q4_qp =', 1pe12.4)
2005 format( 'rzero_qp =', 1pe12.4)
2006 format( 'p0_qp =', 1pe12.4)
2007 format( 'pedge_qp =', 1pe12.4)
2008 format( 'kappa_qp =', 1pe12.4)
2009 format( 'kappae_qp =', 1pe12.4)
2010 format( 'iprint_qp =', i5)
2011 format( 'myrank_qp =', i5)
2012 format( 'itaylor_qp =', i5)
2013 format( 'coolrate_qp =', 1pe12.4)
2014 format( 'q6_qp =', 1pe12.4)
2015 format( 'q8_qp =', 1pe12.4)
2016 format( 'q10_qp =', 1pe12.4)
2017 format( 'q12_qp =', 1pe12.4)
2018 format( 'q14_qp =', 1pe12.4)
2019 format( 'r1_qp  =', 1pe12.4)
2020 format( 'v0_qp =', 1pe12.4)
2021 format( 'v1_qp  =', 1pe12.4)

    if(itaylor.eq.22) call setupLZeqbm

    if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "before loop over elements"

    nelms = local_elements()
    do itri=1,nelms
       
       call define_element_quadrature(itri,int_pts_diag, int_pts_tor)
       call define_fields(itri,0,1,0) ! defines x_79,z_79,mu,nu
       
       !  assemble matrix
       
       do i=1,dofs_per_element
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
          
          dofsps(i) = int2(mu79(:,OP_1,i),temp79a)
          dofsbz(i) = int2(mu79(:,OP_1,i),temp79b)
          dofspr(i) = int2(mu79(:,OP_1,i),temp79c)
          dofsvz(i) = int2(mu79(:,OP_1,i),temp79d)
          dofsden(i) = den0*int1(mu79(:,OP_1,i))
       enddo
       call vector_insert_block(psi_vec%vec,itri,1,dofsps,VEC_ADD)
       call vector_insert_block(bz_vec%vec ,itri,1,dofsbz,VEC_ADD)
       call vector_insert_block(p_vec%vec  ,itri,1,dofspr,VEC_ADD)
       call vector_insert_block(vz_vec%vec  ,itri,1,dofsvz,VEC_ADD)
       call vector_insert_block(den_vec%vec,itri,1,dofsden,VEC_ADD)
    enddo
    
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
  
  subroutine getvals_qsolver(rval,bpsival,ival,pval,vval)
    implicit none
    
    integer, parameter :: N=1000  !  (number of intervals)
    integer :: ifirstq = 1
    real :: dpsi,rval,psival,bpsival,ival,pval,vval
    real :: psimid, qmid, qpmid, ppmid, denom,fterm,gterm,aquad,bquad,cquad,disc,A_qp
    integer :: j
    real, dimension (0:N) :: bpsi, btor, bpolor, psi, jphi, jthor, gradpor, equor, pary, vary
    
    if(ifirstq.eq.1) then
       ifirstq = 0
       A_qp = rzero_qp/r0_qp
       
       !           small psi is normalized r**2
       dpsi = 1./N
       do j=0,N
          psi(j) = j*dpsi
          !  DEBUG
          !   if(iprint_qp .ge.1 .and. myrank_qp .eq.0) write(*,4000) j, psi(j), qfunc(psi(j))
4000      format('j   psi   qfunc(psi)', i5, 1p2e12.4)
       enddo
       
       !  boundary condition at edge
       btor(N) = bz_qp
       bpsi(N) = 0.
       bpolor(N) = btor(N)/(2.*A_qp*qfunc(psi(N)))
       if(myrank_qp.eq.0 .and. iprint_qp.ge.1) write(*,3000) btor(N),bpsi(N),bpolor(N)
3000   format( 'btor(N), bpsi(N), bpolor(N) =',1p3e12.4)
       if(myrank_qp.eq.0 .and. iprint_qp.ge.1) write(*,*) "diagnostics to follow"
       !  integrate first order ode from boundary in
       do j=N,1,-1
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
          if(myrank_qp.eq.0 .and. iprint_qp.ge.1) &
               write(*,1002) qmid,denom,fterm,gterm,aquad,bquad,cquad,disc,btor(j-1)
       enddo
1002   format(1p9e9.1)
       
       !  calculate poloidal and toroidal fields in cell centers
       do j=1,N-1
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
       pary(N) =  pfunc(psi(N))
       vary(0) =  vfunc(psi(0))
       vary(N) =  vfunc(psi(N))
       !
       if(myrank_qp .eq. 0 .and. iprint_qp .ge. 2) then
          write(6,1001)
1001      format(" j       r**2       bpsi        btor         p            v         equil")
          do j=0,N
             write(6,1000) j,psi(j),bpsi(j),btor(j),pary(j),vary(j), equor(j)
          enddo
1000      format(i4,1p8e12.4)
       endif
       
    endif !   end of initialization
    
    psival = rval*rval
    bpsival = cubicinterp(psival,psi,bpsi,N)
    ival = cubicinterp(psival,psi,btor,N)
    pval = cubicinterp(psival,psi,pary,N)
    vval = cubicinterp(psival,psi,vary,N)
  end subroutine getvals_qsolver

  real function cubicinterp(x,xary,yary,N)
    implicit none
    real :: x,xt,del,m1,m2,a,b,c,d
    real, dimension(0:N) :: xary, yary
    integer :: i,N
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
    real :: q_LZ  
    real :: c0,c1,c2,c3,c4 
    real :: asq, bigA, bigB
    real :: ra0
    
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
       
    end select
  end function qfunc

  real function qpfunc(psi)   !   derivative of q wrt psi
    implicit none

    real, intent(in) :: psi !  note:  psi=r^2
    real :: qprime_LZ   
    real :: c0,c1,c2,c3,c4   
    real :: asq, bigA, bigB  
    real :: ra0, psis
    
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
       qpfunc = (1+(psis/ra0**2)**(q12_qp+psis*q14_qp))**(1/(q12_qp+psis*q14_qp))*q0_qp*(-((log(1+(psis/ra0**2)**(q12_qp+psis*q14_qp))*q14_qp)/(q12_qp+psis*q14_qp)**2)+((psis/ra0**2)**(q12_qp+psis*q14_qp)*(log(psis/ra0**2)*q14_qp+(q12_qp+psis*q14_qp)/psis))/((1+(psis/ra0**2)**(q12_qp+psis*q14_qp))*(q12_qp+psis*q14_qp)))*(1+q6_qp/exp((sqrt(psis)-r1_qp)**2/q2_qp**2))-((1+(psis/ra0**2)**(q12_qp+psis*q14_qp))**(1/(q12_qp+psis*q14_qp))*q0_qp*q6_qp*(sqrt(psis)-r1_qp))/(exp((sqrt(psis)-r1_qp)**2/q2_qp**2)*sqrt(psis)*q2_qp**2)
       
    end select
  end function qpfunc

  real function pfunc(psi)    !   p  (pressure)
    implicit none

    real, intent(in) :: psi !  note:  psi=r^2
    real :: p_LZ   
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
       
    end select
    return
  end function pfunc
  
  real function ppfunc(psi)    !  derivative of p wrt psi
    implicit none
    
    real, intent(in) :: psi !  note:  psi=r^2
    real :: asq, bigA, bigB, pprime_LZ

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
       
    end select
  end function ppfunc
  real function vfunc(psi)
    implicit none

    real, intent(in) :: psi !  note:  psi=r^2
    
    vfunc = v0_qp + psi*v1_qp

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
    real :: psi, hsink_qp
    hsink_qp = coolrate_qp*pfunc(psi)

  end function hsink_qp


end module basicq
