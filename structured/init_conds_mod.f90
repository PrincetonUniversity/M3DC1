#define REAL64 real

module init_conds_mod

implicit none

! subroutine velequ
! subroutine phiequ
! subroutine denequ
! subroutine presequ
! subroutine velinit
! subroutine phiinit
! subroutine deninit
! subroutine presinit

contains

real function sech(x)
  real, intent(in) :: x
  real :: cosh
  sech = 1./cosh(x)
  return
end function sech

!============================================================
subroutine velequ(dum)
  use basic
  
  implicit none

#ifndef BIT64
  real, intent(out), dimension(*) :: dum
#else
  REAL64, intent(out), dimension(*) :: dum
#endif

  real :: x, z
  integer :: lx, lz, l, i0, i
  
  ! initialize equilibrium velocity variables to zero
  do lx=1,n
     x = (lx-1)*deex - alx/2.*(1-isym)
     do lz=1,m+jper
        z = (lz-1)*deez - alz/2.*(1-jsym)
        l = lx + (lz-1)*n
        i0 = (l-1)*numvar*6
!
        if(itaylor.eq.6) go to 375
        do i=i0+1,i0+6*numvar
           dum(i) = 0.
        enddo

        go to 502
!
!    itest=6   blob propagation
375     continue
        do i=i0+1,i0+6*numvar
           dum(i) = 0.
        enddo
!       go to 502
502     continue
     enddo
  enddo
  
  return
end subroutine velequ

!============================================================
subroutine phiequ(dum)
  use basic

  implicit none

#ifndef BIT64
  real, intent(out), dimension(*) :: dum
#else
  REAL64, intent(out), dimension(*) :: dum
#endif

  integer :: lx, lz, i0, l, i, ifail1, ifail2, ifail3
  real :: x, z, ff, fp, fpp, r, alam
  real :: k, j0, j1, kb, rt, ri, arg, befo, akz
  real :: n0, pn0, ppn0, pezero
  real :: s17aef, s17aff

  real :: loc, locp, r0

  k = 3.8317059702

  ! start a loop over the number of vertices to define initial conditions

  ! note coordinates are centered at 0,0 for isym=jsym=0
  do lx=1,n
     x = (lx-1)*deex - alx/2.*(1-isym)
     ! NATE: used to be 1,m+1
     do lz=1,m+jper
        z = (lz-1)*deez - alz/2.*(1-jsym)
        l = lx + (lz-1)*n
        i0 = (l-1)*numvar*6
        
        if(itaylor.eq.1) go to 500
        if(itaylor.eq.2) go to 450
        if(itaylor.eq.3) go to 400
        if(itaylor.eq.4) go to 350
        if(itaylor.eq.5) go to 300
        if(itaylor.eq.6) go to 350
        if(itaylor.eq.7) go to 250

        ! for itaylor = 0, use the equilibrium given in R. Richard, et al,
        !     Phys. Fluids B, 2 (3) 1990 p 489
        ! (also Strauss and Longcope)
        !
        r = sqrt(x**2 + z**2)
        ri = 0.
        if(r.ne.0) ri = 1./r
        arg = k*r
        if(r.le.1) then
           befo = 2./s17aef(k,ifail1)
           j0 = s17aef(arg,ifail2)
           j1 = s17aff(arg,ifail3)
           ff = .5*befo
           fp = 0
           fpp = -.125*k**2*befo
           if(arg.ne.0)  then
              ff   = befo*j1/arg
              fp  = befo*(j0 - 2.*j1/arg)/r
              fpp = befo*(-3*j0 + 6.*j1/arg - arg*j1)/r**2
           endif
        else
           ff   = (1.-1./r**2)
           fp  = 2./r**3
           fpp = -6/r**4
        endif

        dum(i0+1) = ff* z
        dum(i0+3) = fp*z**2*ri + ff
        dum(i0+2) = fp*z*x *ri
        dum(i0+6) = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)
        dum(i0+5) = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
        dum(i0+4) = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)

        if(numvar.lt.2) go to 501
        
        if(r.gt.1) go to 501
        !.....define toroidal field
        kb = k**2*(1.-beta)
        dum(i0+7) = sqrt(kb*dum(i0+1)**2+bzero**2)
        dum(i0+8) = kb/dum(i0+7)*dum(i0+1)*dum(i0+2)
        dum(i0+9) = kb/dum(i0+7)*dum(i0+1)*dum(i0+3)
        dum(i0+10) = kb/dum(i0+7)                                   &
             *(-kb/dum(i0+7)**2*(dum(i0+1)*dum(i0+2))**2            &
             + dum(i0+2)**2+dum(i0+1)*dum(i0+4))
        dum(i0+11) = kb/dum(i0+7)                                   &
             *(-kb/dum(i0+7)**2*dum(i0+1)**2*dum(i0+2)*dum(i0+3)    &
             + dum(i0+2)*dum(i0+3)+dum(i0+1)*dum(i0+5))
        dum(i0+12) = kb/dum(i0+7)                                   &
             *(-kb/dum(i0+7)**2*(dum(i0+1)*dum(i0+3))**2            &
             + dum(i0+3)**2+dum(i0+1)*dum(i0+6))
        !.......define pressure
        if(numvar.ge.3) then
           kb = k**2*beta
           if(ipres.eq.0) then
              dum(i0+13) = 0.5*kb*dum(i0+1)**2 + p0
           else
              dum(i0+13) = 0.5*kb*dum(i0+1)**2 + p0 - pi0
           endif
           dum(i0+14) = kb*dum(i0+1)*dum(i0+2)
           dum(i0+15) = kb*dum(i0+1)*dum(i0+3)
           dum(i0+16) = kb*(dum(i0+2)**2+dum(i0+1)*dum(i0+4))
           dum(i0+17) = kb*(dum(i0+2)*dum(i0+3)+dum(i0+1)*dum(i0+5))
           dum(i0+18) = kb*(dum(i0+3)**2+dum(i0+1)*dum(i0+6))
        endif
        go to 502
500     continue

        ! Taylor reconnection (itaylor=1)

        dum(i0+1) = -0.5*z**2
        dum(i0+2) = 0.
        dum(i0+3) = -z
        dum(i0+4) = 0
        dum(i0+5) = 0
        dum(i0+6) = -1
        go to 501

450     continue
!
        ! Force-free taylor state (itaylor=2)
!        x = (lx-1)*deex
!        z = (lz-1)*deez 
        alam = sqrt((pi/alx)**2+(pi/alz)**2)
        dum(i0+1) = sin(pi*x/alx)*sin(pi*z/alz)
        dum(i0+2) = (pi/alx)*cos(pi*x/alx)*sin(pi*z/alz)
        dum(i0+3) = (pi/alz)*sin(pi*x/alx)*cos(pi*z/alz)
        dum(i0+4) = -(pi/alx)**2*sin(pi*x/alx)*sin(pi*z/alz)
        dum(i0+5) = (pi/alx)*(pi/alz)*cos(pi*x/alx)*cos(pi*z/alz)
        dum(i0+6) = -(pi/alz)**2*sin(pi*x/alx)*sin(pi*z/alz)

        if(numvar.ge.2) then
           dum(i0+7) = alam*dum(i0+1)
           dum(i0+8) = alam*dum(i0+2)
           dum(i0+9) = alam*dum(i0+3)
           dum(i0+10) = alam*dum(i0+4)
           dum(i0+11) = alam*dum(i0+5)
           dum(i0+12) = alam*dum(i0+6)
        endif
        go to 502

400     continue

        ! GEM reconnection (itaylor=3)
        
        ! periodic option added 6/17/05
        if(jper.eq.0) then
           dum(i0+1) = 0.5*alog(cosh(2.*z))
           dum(i0+2) = 0.0
           dum(i0+3) = tanh(2.*z)  
           dum(i0+4) = 0.0
           dum(i0+5) = 0.0
           dum(i0+6) = 2.*sech(2.*z)**2

           ! if numvar = 2, then use Bz to satisfy force balance
           if(numvar.eq.2) then
              dum(i0+7) =  sqrt(bzero**2 + sech(2.*z)**2)
              dum(i0+8) =  0.
              dum(i0+9) =  -2.*tanh(2.*z)*sech(2.*z)**2/dum(i0+7)
              dum(i0+10) = 0.
              dum(i0+11) = 0.
              dum(i0+12) = (2.*sech(2.*z)**2/dum(i0+7))                      &
                   *(dum(i0+9)*tanh(2.*z)/dum(i0+7)                          &
                   - 2.*sech(2.*z)**2 + 4.*tanh(2.*z)**2)
           endif

           ! if numvar >= 3, then use pressure to satisfy force balance
           if(numvar.ge.3) then
              
              dum(i0+7) = bzero
              do i=i0+8,i0+12
                 dum(i) = 0.
              enddo

              if(ipres.eq.1) then
                 pezero = p0-pi0
              else
                 pezero = p0
              endif
              dum(i0+13) = pezero*(sech(2.*z)**2 + 0.2)
              dum(i0+14) = 0.
              dum(i0+15) = pezero*(-4.*sech(2.*z)**2*tanh(2.*z))
              dum(i0+16) = 0.
              dum(i0+17) = 0.
              dum(i0+18) = pezero                                         &
                   *(-8.*sech(2.*z)**2*(sech(2.*z)**2-2.*tanh(2.*z)**2))
           endif                     ! if on numvar.ge.3

        else ! jper.ne.0
           akz = 2.*pi/alz
           dum(i0+1) = cos(akz*z)
           dum(i0+2) = 0.0
           dum(i0+3) = -akz*sin(akz*z)
           dum(i0+4) = 0.0
           dum(i0+5) = 0.0
           dum(i0+6) = -akz**2*cos(akz*z)
           if(numvar.ge.2) then
              dum(i0+7) = bzero
              do i=i0+8,i0+12
                 dum(i) = 0.
              enddo
           endif
           if(numvar.ge.3) then
              dum(i0+13) = p0 - pi0
              do i=i0+14,i0+18
                 dum(i) = 0.
              enddo
           endif
        endif !  end on jper.eq.0  

        go to 502
350     continue

        !  Wave propagation (itaylor=4 and itaylor=6)
        
        dum(i0+1) = z*del
        dum(i0+2) = 0.0
        dum(i0+3) = del
        dum(i0+4) = 0.0
        dum(i0+5) = 0.0
        dum(i0+6) = 0.0
 
        go to 501
        
300     continue
        ! Gravitational Instability (itaylor = 5)

        n0 = exp(z/ln)
        pn0 = exp(z/ln)/ln
        ppn0 = exp(z/ln)/ln**2

        ! psi
        dum(i0+1) = 0.
        dum(i0+2) = 0.
        dum(i0+3) = 0.
        dum(i0+4) = 0.
        dum(i0+5) = 0.
        dum(i0+6) = 0.

        !.....bz
        if(numvar.ge.2) then
           if(itaylorw.eq.2) then   ! simple equilibrium
              rt = 0.
           else
              rt = grav*ln+gam*p0   ! physical equilibrium
           endif

           dum(i0+7) = sqrt(bzero**2 - 2.*rt*(n0-1.))
           dum(i0+8) = 0.
           dum(i0+9) = -pn0*rt/dum(i0+7)
           dum(i0+10) = 0.
           dum(i0+11) = 0.
           dum(i0+12) = (rt/dum(i0+7))*(pn0*dum(i0+9)/dum(i0+7) - ppn0)
        endif
        ! pe
        if(numvar.ge.3) then

           if(ipres.eq.1) then
              pezero = p0-pi0
           else
              pezero = p0
           endif

           if(itaylorw.eq.2 .and. p0.gt.0) then  ! simple equilibrium
              dum(i0+13) = (p0 - grav*ln*n0)*pezero/p0
              dum(i0+14) = 0.
              dum(i0+15) = -grav*ln*pn0*pezero/p0
              dum(i0+16) = 0.
              dum(i0+17) = 0.
              dum(i0+18) = -grav*ln*ppn0*pezero/p0
           else                                  ! physical equilibrium
              dum(i0+13) = pezero+gam*pezero*(n0-1.)
              dum(i0+14) = 0.
              dum(i0+15) = gam*pezero*pn0
              dum(i0+16) = 0.
              dum(i0+17) = 0.
              dum(i0+18) = gam*pezero*ppn0
           endif
        endif
        go to 502

250     continue
        ! MTI (itaylor=7)

        if(ipres.eq.1) then
           pezero = p0-pi0
        else
           pezero = p0
        endif

        select case(itaylorw)
        case(2)
           r = sqrt((x/alx)**2+(z/alz)**2+eps)
           dum(i0+1) = del*r
           dum(i0+2) = del*x/(r*alx**2)
           dum(i0+3) = del*z/(r*alz**2)
           dum(i0+4) = del*(1.-(x/(r*alx**2))**2)/(r*alz**2)
           dum(i0+5) =-del*(x/alx**2)*(z/alz**2)/r**3
           dum(i0+6) = del*(1.-(z/(r*alz**2))**2)/(r*alz**2)
   
        case default 
           dum(i0+1) = del*z
           dum(i0+2) = 0.
           dum(i0+3) = del
           dum(i0+4) = 0.
           dum(i0+5) = 0.
           dum(i0+6) = 0.
        end select

        ! bz
        if(numvar.ge.2) then
           dum(i0+7) = bzero
           dum(i0+8) = 0.
           dum(i0+9) = 0.
           dum(i0+10) = 0.
           dum(i0+11) = 0.
           dum(i0+12) = 0.
        endif

        if(numvar.ge.3) then
           if(ipres.eq.1) then
              pezero = p0-pi0
           else
              pezero = p0
           endif

           select case(itaylorw)
           case(0)
              dum(i0+13) = pezero*(1.-z/ln)**3
              dum(i0+14) = 0.
              dum(i0+15) =-3.*pezero*(1.-z/ln)**2/ln
              dum(i0+16) = 0.
              dum(i0+17) = 0.
              dum(i0+18) = 6.*pezero*(1.-z/ln)/ln**2

           case(3)

              r = sqrt(x**2+z**2+eps)
              r0 = 0.0

              loc = alx/ln
              locp = alz/ln

              dum(i0+13) = eps*exp(-((x-r0)/loc)**2-(z/locp)**2)
              dum(i0+14) = -2.*(x-r0)*dum(i0+13)/loc**2
              dum(i0+15) = -2.*z*dum(i0+13)/locp**2
              dum(i0+16) = 2.*(2.*((x-r0)/loc)**2 - 1.)*dum(i0+13)/loc**2
              dum(i0+17) = 4.*(x-r0)*z*dum(i0+13)/(loc*locp)**2
              dum(i0+18) = 2.*(2.*(z/locp)**2 - 1.)*dum(i0+13)/locp**2

              dum(i0+13) = dum(i0+13)+pezero

           case default
              dum(i0+13) = pezero
              dum(i0+14) = 0.
              dum(i0+15) = 0.
              dum(i0+16) = 0.
              dum(i0+17) = 0.
              dum(i0+18) = 0.
           end select
        endif

        goto 502
        
501     continue

        if(numvar.gt.1) then
           dum(i0+7) = bzero
           do i=i0+8,i0+12
              dum(i) = 0.
           enddo
        endif
        if(numvar.gt.2) then
           dum(i0+13) = p0 - pi0
           do i=i0+14,i0+18
              dum(i) = 0.
           enddo
        endif

 502  continue

     enddo
  enddo

  return
end subroutine phiequ
!============================================================
subroutine denequ(dum)
  use basic

  implicit none

#ifndef BIT64
  real, intent(out), dimension(*) :: dum
#else
  REAL64, intent(out), dimension(*) :: dum
#endif

  integer :: lx, lz, l, i0, ii
  real :: x, z, rsq,  n0, pn0, ppn0
  

  ! initialize equilibrium to be the Strauss and Longcope solution
  
  ! start a loop over the number of vertices to define initial conditions

  ! note coordinates are centered at 0,0 for isym=jsym=0
  do lx=1,n
  x = (lx-1)*deex - alx/2.*(1-isym)
  do lz=1,m
     z = (lz-1)*deez - alz/2.*(1-jsym)
     l = lx + (lz-1)*n
     i0 = (l-1)*6
     
     if(itaylor.eq.1) go to 500
     if(itaylor.eq.2) go to 500
     if(itaylor.eq.3) go to 400
     if(itaylor.eq.4) go to 500
     if(itaylor.eq.5) go to 300
     if(itaylor.eq.6) go to 375
     if(itaylor.eq.7) go to 250

     ! for itaylor = 0, use the equilibrium given in R. Richard, et al,
     !     Phys. Fluids B, 2 (3) 1990 p 489
     ! ..... uniform density

     dum(i0+1) = 1.0
     dum(i0+2) = 0.
     dum(i0+3) = 0.
     dum(i0+4) = 0.
     dum(i0+5) = 0.
     dum(i0+6) = 0.
     
     go to 502
500  continue
     dum(i0+1) = 1.
     do ii=2,6
        dum(i0+ii) = 0.
     enddo
     
     go to 502
400  continue
     
     ! GEM reconnection (itaylor=3)
!            (NOTE: GEM density for facden=1, constant density for facden=0.)
     dum(i0+1) = facden*(sech(2.*z)**2 + 0.2) + (1.-facden)*0.2
     dum(i0+2) = 0.
     dum(i0+3) = -facden*4.*sech(2.*z)**2*tanh(2.*z)
     dum(i0+4) = 0.
     dum(i0+5) = 0.
     dum(i0+6) = -facden*8.*sech(2.*z)**2*(sech(2.*z)**2-2.*tanh(2.*z)**2)
     go to 502     
     
300  continue
     
     ! Gravitational Instability (itaylor=5)
     n0 = exp(z/ln)
     pn0 = exp(z/ln)/ln
     ppn0 = exp(z/ln)/ln**2
     dum(i0+1) = n0
     dum(i0+2) = 0.
     dum(i0+3) = pn0
     dum(i0+4) = 0.
     dum(i0+5) = 0.
     dum(i0+6) = ppn0
     go to 502

375     continue
!
!    blob transport
     rsq = (x**2+z**2)*2.
     dum(i0+1) = p0*exp(-rsq)+ 0.01
     dum(i0+2) = -2.*x*p0*exp(-rsq)*2.
     dum(i0+3) = -2.*z*p0*exp(-rsq)*2.
     dum(i0+4) = (-2.+4.*x**2)*p0*exp(-rsq)*4.
     dum(i0+5) = 4*x*z*p0*exp(-rsq)*4.
     dum(i0+6) = (-2.+4.*z**2)*p0*exp(-rsq)*4.
     go to 502

250 continue
     ! MTI (itaylor=7)

     select case(itaylorw)
     case(0)
        dum(i0+1) = (1.-z/ln)**2
        dum(i0+2) = 0.
        dum(i0+3) = -2.*(1.-z/ln)/ln
        dum(i0+4) = 0.
        dum(i0+5) = 0.
        dum(i0+6) = 2./ln**2
     case default
        dum(i0+1) = 1.
        dum(i0+2) = 0.
        dum(i0+3) = 0.
        dum(i0+4) = 0.
        dum(i0+5) = 0.
        dum(i0+6) = 0.
     end select
     go to 502
     
    
502  continue
     
  enddo
  enddo
  
  return
end subroutine denequ
!============================================================
subroutine presequ(dum)
  use basic
  
  implicit none

#ifndef BIT64
  real, intent(out), dimension(*) :: dum
#else
  REAL64, intent(out), dimension(*) :: dum
#endif

  integer :: lx, lz, l, i0
  real :: x, z, rsq


  ! initialize equilibrium to be the Strauss and Longcope solution

  ! start a loop over the number of vertices to define initial conditions

  ! note coordinates are centered at 0,0
  do lx=1,n
     x = (lx-1)*deex - alx/2.*(1-isym)
     do lz=1,m
        z = (lz-1)*deez - alz/2.*(1-jsym)
        l = lx + (lz-1)*n
        i0 = (l-1)*6
        
        if(itaylor.eq.3) go to 300
        if(itaylor.eq.5) go to 500
        if(itaylor.eq.6) go to 375
        if(itaylor.eq.7) go to 250
      
        dum(i0+1) = p0
        dum(i0+2) = 0.
        dum(i0+3) = 0.
        dum(i0+4) = 0.
        dum(i0+5) = 0.
        dum(i0+6) = 0.
        go to 502

300     continue

        ! GEM reconnection (itaylor=3)
        dum(i0+1) = p0*(sech(2.*z)**2 + 0.2)
        dum(i0+2) = 0.
        dum(i0+3) = p0*(-4.*sech(2.*z)**2*tanh(2.*z))
        dum(i0+4) = 0.
        dum(i0+5) = 0.
        dum(i0+6) = p0*(-8.*sech(2.*z)**2*(sech(2.*z)**2                  &
             -2.*tanh(2.*z)**2))

        go to 502     
        
500     continue

        ! Gravitational instability (itaylor=5)
        if(itaylorw.eq.2) then                     ! simple equilibrium
           dum(i0+1) = p0-grav*ln*exp(z/ln)
           dum(i0+2) = 0.
           dum(i0+3) = -grav*exp(z/ln)
           dum(i0+4) = 0.
           dum(i0+5) = 0.
           dum(i0+6) = -grav*exp(z/ln)/ln
        else                                     ! physical equilibrium
           dum(i0+1) = p0+gam*p0*(exp(z/ln)-1.)
           dum(i0+2) = 0.
           dum(i0+3) = gam*p0*exp(z/ln)/ln
           dum(i0+4) = 0.
           dum(i0+5) = 0.
           dum(i0+6) = gam*p0*exp(z/ln)/ln**2
        endif
        
        go to 502     

375     continue
!
!       blob transport
        rsq = x**2+z**2
        dum(i0+1) = p0*exp(-rsq)+.01
        dum(i0+2) = -2.*x*p0*exp(-rsq)
        dum(i0+3) = -2.*z*p0*exp(-rsq)
        dum(i0+4) = (-2.+4.*x**2)*p0*exp(-rsq) 
        dum(i0+5) = 4*x*z*p0*exp(-rsq)
        dum(i0+6) = (-2.+4.*z**2)*p0*exp(-rsq)
        go to 502

250     continue
        ! p

        select case(itaylorw)
        case(0)
           dum(i0+1) = p0*(1.-z/ln)**3
           dum(i0+2) = 0.
           dum(i0+3) = -3.*p0*(1.-z/ln)**2/ln
           dum(i0+4) = 0.
           dum(i0+5) = 0.
           dum(i0+6) = 6.*p0*(1.-z/ln)/ln**2
        case default
           dum(i0+1) = p0
           dum(i0+2) = 0.
           dum(i0+3) = 0.
           dum(i0+4) = 0.
           dum(i0+5) = 0.
           dum(i0+6) = 0.
        end select

        go to 502

502     continue

        
     enddo
  enddo

  
  return
end subroutine presequ

!========================================================================
subroutine velinit(dum)
  use basic

  implicit none

#ifndef BIT64
  real, intent(out), dimension(*) :: dum
#else
  REAL64, intent(out), dimension(*) :: dum
#endif

  integer :: lx, lz, i0, i, l
  real :: phiper, vzper, chiper
  real :: f,fx,fxx,fxz,fz,fzz,g,gx,gz,gxz,gxx,gzz,uzero,czero,befo
  real :: x, z, akx, akz, loc, locp, locpp
  real rsq

  ! initialize perturbed velocity variables
  
  if(itaylor.eq.4) then
     call wave_perturbation(phiper, vzper, chiper, 2)
  endif
  
  do lx=1,n
  do lz=1,m
     x = (lx-1)*deex - alx/2.*(1-isym)
     z = (lz-1)*deez - alz/2.*(1-jsym)
     l = lx + (lz-1)*n
     i0 = (l-1)*numvar*6
     
     if(itaylor.eq.4) go to 350
     if(itaylor.eq.5) go to 300
     if(itaylor.eq.6) go to 375
     if(itaylor.eq.7) go to 250
!
     if(itaylor.eq.0) then
       befo = eps*exp(-(x**2+z**2))
       uzero = befo*cos(z)
       czero = befo*sin(z)
       dum(i0+1) = uzero
       dum(i0+2) = - 2.*x*uzero
       dum(i0+3) = - 2.*z*uzero - czero
       dum(i0+4) = (-2.+4.*x**2)*uzero
       dum(i0+5) = 4.*x*z*uzero + 2.*x*czero
       dum(i0+6) = (-3.+4.*z**2)*uzero + 4.*z*czero

       if(numvar.ge.2) then
         do i=7,12
           dum(i0+i) = 0
         enddo
       endif
       if(numvar.ge.3) then
         dum(i0+13) = czero
         dum(i0+14) = -2*x*czero
         dum(i0+15) = -2*z*czero + uzero
         dum(i0+16) = (-2.+4*x**2)*czero
         dum(i0+17) = 4.*x*z*czero - 2.*x*uzero
         dum(i0+18) = (-3.+4.*z**2)*czero - 4*z*uzero
       endif

       
go to 502
     endif
     
     
     do i=i0+1,i0+6*numvar
        dum(i) = 0.
     enddo
     go to 502

350  continue
     ! initialize initial perturbation for wave propagation (itaylor=4)
     akx = 2.*pi/alx

     ! phi
     dum(i0+1) =  phiper*cos(akx*x)
     dum(i0+2) = -phiper*akx*sin(akx*x)
     dum(i0+3) =  0.
     dum(i0+4) = -phiper*akx**2*cos(akx*x)
     dum(i0+5) =  0.
     dum(i0+6) =  0.
     
     ! vz
     if(numvar.ge.2) then
        dum(i0+7) =  vzper*cos(akx*x)
        dum(i0+8) = -vzper*akx*sin(akx*x)
        dum(i0+9) =  0.
        dum(i0+10)= -vzper*akx**2*cos(akx*x)
        dum(i0+11)=  0.
        dum(i0+12)=  0.
     endif
    
     ! chi
     if(numvar.ge.3) then
        dum(i0+13) = -chiper*sin(akx*x)
        dum(i0+14) = -chiper*akx*cos(akx*x)
        dum(i0+15) =  0.
        dum(i0+16) =  chiper*akx**2*sin(akx*x)
        dum(i0+17) =  0.
        dum(i0+18) =  0.
     endif
     go to 502

!
 375 continue
     ! initialize initial perturbation for blob transport (itaylor=7)
     akx = 2.*pi/alx
     f =           sin(akx*x)
     fx=    akx*   cos(akx*x)
     fz = 0
     fxx = -akx**2*sin(akx*x)
     fzz = 0.
     fxz = 0.
     rsq = (x**2+z**2)*2.
     g = eps*exp(-rsq)                
     gx= -2.*x*eps*exp(-rsq)*2.
     gz= -2.*z*eps*exp(-rsq)*2.
     gxx= (-2.+4.*x**2)*eps*exp(-rsq)*4.
     gxz= 4*x*z*eps*exp(-rsq)*4.
     gzz= (-2.+4.*z**2)*eps*exp(-rsq)*4.
!
     dum(i0+1) = f*g
     dum(i0+2) = fx*g + f*gx
     dum(i0+3) = fz*g + f*gz
     dum(i0+4) = fxx*g + 2*fx*gx + f*gxx
     dum(i0+5) = fxz*g + fx*gz + fz*gx + f*gxz
     dum(i0+6) = fzz*g + 2*fz*gz + f*gzz


     ! phi
!     dum(i0+1) =   eps*       cos(akx*x)*cos(akx*z)
!     dum(i0+2) =  -eps*akx*   sin(akx*x)*cos(akx*z)
!     dum(i0+3) =  -eps*akx*   cos(akx*x)*sin(akx*z)
!     dum(i0+4) =  -eps*akx**2*cos(akx*x)*cos(akx*z)
!     dum(i0+5) =  +eps*akx**2*sin(akx*x)*sin(akx*z)
!     dum(i0+6) =  -eps*akx**2*cos(akx*x)*cos(akx*z)
     
     ! vz
     if(numvar.ge.2) then
        dum(i0+7) =  0.
        dum(i0+8) =  0.
        dum(i0+9) =  0.
        dum(i0+10)=  0.
        dum(i0+11)=  0.
        dum(i0+12)=  0.
     endif
    
     ! chi
     if(numvar.ge.3) then
        dum(i0+13) = 0.
        dum(i0+14) = 0.
        dum(i0+15) =  0.
        dum(i0+16) =  0.
        dum(i0+17) =  0.
        dum(i0+18) =  0.
     endif
     go to 502

     ! initialize velocity perturbation for grav. inst. (itaylor=5)
300  continue

     akx = 2.*del*pi/alx
     akz = 2.*pi/alz
     
     loc = cos(akz*z/2.)
     locp = -sin(akz*z/2.)*akz/2.
     locpp = -cos(akz*z/2.)*akz**2/4.
     
     if(grav.eq.0.) then
        phiper = eps
     else
        phiper = 0.
     endif

     dum(i0+1) = loc  *( phiper*cos(akx*x))
     dum(i0+2) = loc  *(-phiper*sin(akx*x))*akx
     dum(i0+3) = locp *( phiper*cos(akx*x))
     dum(i0+4) = loc  *(-phiper*cos(akx*x))*akx**2
     dum(i0+5) = locp *(-phiper*sin(akx*x))*akx
     dum(i0+6) = locpp*( phiper*cos(akx*x))
     if(numvar.ge.2) then 
        do i=i0+7,i0+12
           dum(i) = 0.
        enddo
     endif
     if(numvar.ge.3) then
        if(itaylorw.eq.1) then
           chiper = eps
        else
           chiper = 0.
        endif
        dum(i0+13) = loc  *( chiper*cos(akx*x))
        dum(i0+14) = loc  *(-chiper*sin(akx*x))*akx
        dum(i0+15) = locp *( chiper*cos(akx*x))
        dum(i0+16) = loc  *(-chiper*cos(akx*x))*akx**2
        dum(i0+17) = locp *(-chiper*sin(akx*x))*akx
        dum(i0+18) = locpp*( chiper*cos(akx*x))
     endif
     go to 502
   
     ! initialize velocity perturbation for MTI (itaylor=7)

250  continue


     select case(itaylorw)
!     case(0)
!        akx = 2.*pi/alx
!        akz = 2.*pi/alz
!
!!        loc = 27.5453*eps*cos(akz*z/2.)
!!        locp = -27.5453*eps*sin(akz*z/2.)*akz/2.
!!        locpp = -27.5453*eps*cos(akz*z/2.)*akz**2/4.
!        loc = eps*cos(akz*z/2.)
!        locp = -eps*sin(akz*z/2.)*akz/2.
!        locpp = -eps*cos(akz*z/2.)*akz**2/4.
!
     case default
        loc = 0.
        locp = 0.
        locpp = 0.
     end select
    
     dum(i0+1) = loc  *( sin(akx*x))
     dum(i0+2) = loc  *( cos(akx*x))*akx
     dum(i0+3) = locp *( sin(akx*x))
     dum(i0+4) = loc  *(-sin(akx*x))*akx**2
     dum(i0+5) = locp *( cos(akx*x))*akx
     dum(i0+6) = locpp*( sin(akx*x))

     if(numvar.ge.2) then 
        do i=i0+7,i0+12
           dum(i) = 0.
        enddo
     endif

     if(numvar.ge.3) then
        select case(itaylorw)
!        case(0)
!           akx = 2.*pi/alx
!           akz = 2.*pi/alz
!
!           loc = 0.131509*eps*cos(akz*z/2.)
!           locp = -0.131509*eps*sin(akz*z/2.)*akz/2.
!           locpp = -0.131509*eps*cos(akz*z/2.)*akz**2/4.

        case default
           loc = 0.
           locp = 0.
           locpp = 0.
        end select
    
        dum(i0+13) = loc  *( cos(akx*x))
        dum(i0+14) = loc  *(-sin(akx*x))*akx
        dum(i0+15) = locp *( cos(akx*x))
        dum(i0+16) = loc  *(-cos(akx*x))*akx**2
        dum(i0+17) = locp *(-sin(akx*x))*akx
        dum(i0+18) = locpp*( cos(akx*x))
     endif
     go to 502

502  continue

  enddo
  enddo
  return
end subroutine velinit

!============================================================
subroutine phiinit(dum)
  use basic

  implicit none

#ifndef BIT64
  real, intent(out), dimension(*) :: dum
#else
  REAL64, intent(out), dimension(*) :: dum
#endif

  integer :: lx, lz, i0, l, i
  real :: x, z, r, r0
  real :: psiper, bzper, peper, akx, akz, loc, locp, locpp

  ! start a loop over the number of vertices to define initial conditions
  
  ! note coordinates are centered at 0,0 for isym=jsym=0
  
  if(itaylor.eq.4) then
     call wave_perturbation(psiper, bzper, peper, 1)
  endif
  
  do lx=1,n
  x = (lx-1)*deex - alx/2.*(1-isym)
  do lz=1,m
     z = (lz-1)*deez - alz/2.*(1-jsym)
     l = lx + (lz-1)*n
     i0 = (l-1)*numvar*6
     
     if(itaylor.eq.1) go to 500
     if(itaylor.eq.2) go to 450
     if(itaylor.eq.3) go to 400
     if(itaylor.eq.4) go to 350
     if(itaylor.eq.5) go to 300
     if(itaylor.eq.7) go to 250
     
     dum(i0+1) = 0.
     dum(i0+3) = 0.
     dum(i0+2) = 0.
     dum(i0+6) = 0.
     dum(i0+5) = 0.
     dum(i0+4) = 0.
     go to 501
     
500  continue

     ! Taylor reconnection (itaylor=1)

     dum(i0+1) = 0.
     dum(i0+2) = 0.
     dum(i0+3) = 0.
     dum(i0+4) = 0
     dum(i0+5) = 0
     dum(i0+6) = 0.
     go to 501
450  continue
     
     ! define a force-free taylor state for itaylor=2
     dum(i0+1) = 0.
     dum(i0+2) = 0.
     dum(i0+3) = 0.
     dum(i0+4) = 0.
     dum(i0+5) = 0.
     dum(i0+6) = 0.
     
     if(numvar.ge.2) then
        dum(i0+7) = 0.
        dum(i0+8) = 0.
        dum(i0+9) = 0.
        dum(i0+10) = 0.
        dum(i0+11) = 0.
        dum(i0+12) = 0.
     endif
     go to 502
     
400  continue

     ! GEM reconnection (itaylor=3)
     akx = 2.*pi/alx
     akz = pi/alz
     if(jper.eq.1) akz = 2.*pi/alz
     
     ! flux
     dum(i0+1) =  eps*cos(akx*x)*cos(akz*z)
     dum(i0+2) = -eps*akx*sin(akx*x)*cos(akz*z)
     dum(i0+3) = -eps*akz*cos(akx*x)*sin(akz*z)
     dum(i0+4) = -eps*akx**2*cos(akx*x)*cos(akz*z)
     dum(i0+5) =  eps*akx*akz*sin(akx*x)*sin(akz*z)
     dum(i0+6) = -eps*akz**2*cos(akx*x)*cos(akz*z)
     
     if(numvar.ge.2) then
        do i=i0+7,i0+12
           dum(i) = 0.
        enddo
     endif
     if(numvar.ge.3) then
        do i=i0+13,i0+18
           dum(i) = 0.
        enddo
     endif
     go to 502
350  continue
     ! Wave propagation (itaylor=4)
     akx = 2.*pi/alx
     
     ! psi
     dum(i0+1) =  psiper*cos(akx*x)
     dum(i0+2) = -psiper*akx*sin(akx*x)
     dum(i0+3) =  0.
     dum(i0+4) = -psiper*akx**2*cos(akx*x)
     dum(i0+5) =  0.
     dum(i0+6) =  0.
     
     ! Bz
     if(numvar.ge.2) then 
        dum(i0+7) =   bzper*cos(akx*x)
        dum(i0+8) =  -bzper*akx*sin(akx*x)
        dum(i0+9) =   0.
        dum(i0+10) = -bzper*akx**2*cos(akx*x)
        dum(i0+11) =  0.
        dum(i0+12) =  0.
     endif
     
     ! pe
     if(numvar.ge.3) then
        dum(i0+13) =  peper*cos(akx*x)
        dum(i0+14) = -peper*akx*sin(akx*x)
        dum(i0+15) =  0.
        dum(i0+16) = -peper*akx**2*cos(akx*x)
        dum(i0+17) =  0.
        dum(i0+18) =  0.
     endif
     go to 502
     
     ! Gravitational Instability (itaylor=5)
300   continue
     
     do i=i0+1,i0+6
        dum(i) = 0.
     enddo
     
     if(numvar.ge.2) then 
        do i=i0+7, i0+12
           dum(i) = 0.
        enddo
     endif
     
     if(numvar.ge.3) then
        akx = 2.*del*pi/alx
        akz = 2.*pi/alz
        
        loc = cos(akz*z/2.)
        locp = -sin(akz*z/2.)*akz/2.
        locpp = -cos(akz*z/2.)*akz**2/4.
        
        if(grav.eq.0. .or. itaylorw.eq.1) then 
           peper = 0. 
        else
           peper = eps*(p0-pi0)
        endif
        
        dum(i0+13) = peper*loc*cos(akx*x)
        dum(i0+14) =-peper*loc*sin(akx*x)*akx
        dum(i0+15) = peper*locp*cos(akx*x)
        dum(i0+16) =-peper*loc*cos(akx*x)*akx**2
        dum(i0+17) =-peper*locp*sin(akx*x)*akx
        dum(i0+18) = peper*locpp*cos(akx*x)
     endif
     go to 502

250 continue
     
     do i=i0+1,i0+6
        dum(i) = 0.
     enddo

     akx = 2.*pi/alx
     akz = 2.*pi/alz

     select case(itaylorw) 
!     case(0)
!        loc = eps*cos(akz*z/2.)
!        locp = -eps*sin(akz*z/2.)*akz/2.
!        locpp = -eps*cos(akz*z/2.)*akz**2/4.
     case default
        loc = 0.
        locp = 0.
        locpp = 0.
     end select
    
     dum(i0+1) = loc  *( cos(akx*x))
     dum(i0+2) = loc  *(-sin(akx*x))*akx
     dum(i0+3) = locp *( cos(akx*x))
     dum(i0+4) = loc  *(-cos(akx*x))*akx**2
     dum(i0+5) = locp *(-sin(akx*x))*akx
     dum(i0+6) = locpp*( cos(akx*x))
     
     if(numvar.ge.2) then 
        do i=i0+7, i0+12
           dum(i) = 0.
        enddo
     endif
     
     if(numvar.ge.3) then

        akx = 2.*pi/alx
        akz = 2.*pi/alz
        
        peper = eps*(p0-pi0)

!        peper = 0.

        select case(itaylorw)
        case(1)
           r = sqrt(x**2+z**2+eps)
           r0 = 0.0

           loc = alx/ln
           locp = alz/ln

           dum(i0+13) = peper*exp(-((x-r0)/loc)**2-(z/locp)**2)
           dum(i0+14) = -2.*(x-r0)*dum(i0+13)/loc**2
           dum(i0+15) = -2.*z*dum(i0+13)/locp**2
           dum(i0+16) = 2.*(2.*((x-r0)/loc)**2 - 1.)*dum(i0+13)/loc**2
           dum(i0+17) = 4.*(x-r0)*z*dum(i0+13)/(loc*locp)**2
           dum(i0+18) = 2.*(2.*(z/locp)**2 - 1.)*dum(i0+13)/locp**2

        case(2)
           r = sqrt(x**2+z**2+eps)
           r0 = alx/4.

           loc = alx/ln
           locp = alz/ln

           dum(i0+13) = peper*exp(-((x-r0)/loc)**2-(z/locp)**2)
           dum(i0+14) = -2.*(x-r0)*dum(i0+13)/loc**2
           dum(i0+15) = -2.*z*dum(i0+13)/locp**2
           dum(i0+16) = 2.*(2.*((x-r0)/loc)**2 - 1.)*dum(i0+13)/loc**2
           dum(i0+17) = 4.*(x-r0)*z*dum(i0+13)/(loc*locp)**2
           dum(i0+18) = 2.*(2.*(z/locp)**2 - 1.)*dum(i0+13)/locp**2

        case default
           dum(i0+13) = 0.
           dum(i0+14) = 0.
           dum(i0+15) = 0.
           dum(i0+16) = 0.
           dum(i0+17) = 0.
           dum(i0+18) = 0.
        end select



     endif
     go to 502

501  continue
     
     if(numvar.gt.1) then
        dum(i0+7) = 0.
        do i=i0+8,i0+6*numvar
           dum(i) = 0.
        enddo
     endif
502  continue
     
  enddo
  enddo

  return
end subroutine phiinit

!============================================================
subroutine deninit(dum)
  use basic

  implicit none

#ifndef BIT64
  real, intent(out), dimension(*) :: dum
#else
  REAL64, intent(out), dimension(*) :: dum
#endif

  integer :: lx, lz, l, i0
  real :: x, z, akx, akz
  real :: nper, loc, locp, locpp 

  ! start a loop over the number of vertices to define initial conditions

  ! note coordinates are centered at 0,0 for isym=jsym=0

  if(itaylor.eq.4) then
     call wave_perturbation(nper, nper, nper, 3)
  endif
  
  do lx=1,n
  x = (lx-1)*deex - alx/2.*(1-isym)
  do lz=1,m
     z = (lz-1)*deez - alz/2.*(1-jsym)
     l = lx + (lz-1)*n
     i0 = (l-1)* 6

     if(itaylor.eq.4) go to 400
     if(itaylor.eq.5) go to 500
     if(itaylor.eq.7) go to 600


     dum(i0+1) = 0.
     dum(i0+3) = 0.
     dum(i0+2) = 0.
     dum(i0+6) = 0.
     dum(i0+5) = 0.
     dum(i0+4) = 0.

     go to 1000
     
     ! initialize density for wave propagation (itaylor=4)
400  continue
     akx = 2.*pi/alx
     dum(i0+1) =  nper*cos(akx*x)
     dum(i0+2) = -nper*akx*sin(akx*x)
     dum(i0+3) =  0.
     dum(i0+4) = -nper*akx**2*cos(akx*x)
     dum(i0+5) =  0.
     dum(i0+6) =  0.
     
     go to 1000
     
     ! initialize density for gravitational instability (itaylor=5)
500  continue
     akx = 2.*del*pi/alx
     akz = 2.*pi/alz
     
     loc = cos(akz*z/2.)
     locp = -sin(akz*z/2.)*akz/2.
     locpp = -cos(akz*z/2.)*akz**2/4.
     
     if(grav.eq.0. .or. itaylorw.eq.1) then 
        nper = 0
     else
        nper = eps
     endif

     dum(i0+1) = nper*loc*cos(akx*x)
     dum(i0+2) =-nper*loc*sin(akx*x)*akx
     dum(i0+3) = nper*locp*cos(akx*x)
     dum(i0+4) =-nper*loc*cos(akx*x)*akx**2
     dum(i0+5) =-nper*locp*sin(akx*x)*akx
     dum(i0+6) = nper*locpp*cos(akx*x)
     
     go to 1000

600  continue

     select case(itaylorw)
     case(0)
        akx = 2.*pi/alx
        akz = 2.*pi/alz
     
        loc = cos(akz*z/2.)
        locp = -sin(akz*z/2.)*akz/2.
        locpp = -cos(akz*z/2.)*akz**2/4.

!        nper = 0.999766*eps
        nper = eps
!        nper = 0

        dum(i0+1) = loc  *( nper*cos(akx*x))
        dum(i0+2) = loc  *(-nper*sin(akx*x))*akx
        dum(i0+3) = locp *( nper*cos(akx*x))
        dum(i0+4) = loc  *(-nper*cos(akx*x))*akx**2
        dum(i0+5) = locp *(-nper*sin(akx*x))*akx
        dum(i0+6) = locpp*( nper*cos(akx*x))

     case default
        dum(i0+1) = 0.
        dum(i0+2) = 0.
        dum(i0+3) = 0.
        dum(i0+4) = 0.
        dum(i0+5) = 0.
        dum(i0+6) = 0.
     end select


     go to 1000


1000 continue
  enddo
  enddo

  return
end subroutine deninit
!============================================================
subroutine presinit(dum)
  use basic

  implicit none
  
#ifndef BIT64
  real, intent(out), dimension(*) :: dum
#else
  REAL64, intent(out), dimension(*) :: dum
#endif

  integer :: lx, lz, l, i0
  real :: x, z, akz
  real :: pper,  akx, loc, locp, locpp
  
  ! start a loop over the number of vertices to define initial conditions
  
  ! note coordinates are centered at 0,0
  
  if(itaylor.eq.4) then
     call wave_perturbation(pper, pper, pper, 4)
  endif

  do lx=1,n
  x = (lx-1)*deex - alx/2.
  do lz=1,m
     z = (lz-1)*deez - alz/2.
     l = lx + (lz-1)*n
     i0 = (l-1)* 6

     if(itaylor.eq.4) go to 400
     if(itaylor.eq.5) go to 500

     dum(i0+1) = 0.
     dum(i0+2) = 0.
     dum(i0+3) = 0.
     dum(i0+4) = 0.
     dum(i0+5) = 0.
     dum(i0+6) = 0.

     go to 1000
     
400  continue

     ! initialize pressure perturbation for wave propagation (itaylor=4)
     akx = 2.*pi/alx
     dum(i0+1) =  pper*cos(akx*x)
     dum(i0+2) = -pper*akx*sin(akx*x)
     dum(i0+3) =  0.
     dum(i0+4) = -pper*akx**2*cos(akx*x)
     dum(i0+5) =  0.
     dum(i0+6) =  0.

     go to 1000
     
     ! initialize pressure perturbation for gravitational inst. (itaylor=5)
500  continue
     akx = 2.*del*pi/alx
     akz = 2.*pi/alz
     
     loc = cos(akz*z/2.)
     locp = -sin(akz*z/2.)*akz/2.
     locpp = -cos(akz*z/2.)*akz**2/4.

     if(grav.eq.0. .or. itaylorw.eq.1) then
        pper = 0 
     else
        pper = eps*p0
     endif

     dum(i0+1) = loc  *( pper*cos(akx*x))
     dum(i0+2) = loc  *(-pper*sin(akx*x))*akx
     dum(i0+3) = locp *( pper*cos(akx*x))
     dum(i0+4) = loc  *(-pper*cos(akx*x))*akx**2
     dum(i0+5) = locp *(-pper*sin(akx*x))*akx
     dum(i0+6) = locpp*( pper*cos(akx*x))

     go to 1000

1000 continue
  enddo
  enddo
  
  return
end subroutine presinit

end module init_conds_mod
