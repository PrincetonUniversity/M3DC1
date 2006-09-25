!=====================================
subroutine velequ(dum, numberingid)
  
  use basic
  
  implicit none
  integer :: numnodes, l, ibegin, iendplusone, i, numberingid
  real :: dum(*)

  real :: alx, alz, akx, akz, x, z, xmin, zmin
  double precision :: coords(3)

  call getboundingboxsize(alx, alz)

  ! initialize equilibrium velocity variables to zero 
  call numnod(numnodes)
  do l=1, numnodes
     call entdofs(numberingid, l, 0, ibegin, iendplusone)

     call xyznod(l, coords)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     do i=ibegin, iendplusone-1
        dum(i) = 0.
     enddo

!     select case(itaylor)

!     end select
  enddo

  return
end subroutine velequ

!============================================================
subroutine phiequ(dum, numberingid)
  
  use basic
  implicit none
  integer :: l, i, ifail1, ifail2, ifail3, numnodes, ibegin, iendplusone, numberingid
  real :: dum(*), x, z, r, ri, arg, befo, s17aef, s17aff
  real :: ff, fp, fpp, alam, sech, alx, alz, xmin, zmin
  real :: k, j0, j1, kb, akx, akz
  double precision :: coords(3)

  k = 3.8317059702

  ! start a loop over the number of vertices to define initial conditions

  ! note coordinates are centered at 0,0

  call numnod(numnodes)
  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)
  
  do l=1,numnodes
     call entdofs(numberingid, l, 0, ibegin, iendplusone)
     call xyznod(l, coords)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5
     
     select case(itaylor)
     
     case(0)
!
!.....for itaylor = 0, use the equilibrium given in R. Richard, et al,
!     Phys. Fluids B, 2 (3) 1990 p 489
!.....(also Strauss and Longcope)
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
   
        dum(ibegin) = ff* z
        dum(ibegin+1) = fp*z*x *ri
        dum(ibegin+2) = fp*z**2*ri + ff
        dum(ibegin+3) = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)
        dum(ibegin+4) = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
        dum(ibegin+5) = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)


        if(numvar.ge.2) then

           if(r.gt.1) then

              if(numvar.ge.2) then
                 dum(ibegin+6) = bzero
                 do i=ibegin+7,ibegin+11
                    dum(i) = 0.
                 enddo
              endif
              if(numvar.ge.3) then
                 dum(ibegin+12) = p0 - pi0
                 do i=ibegin+13,ibegin+17
                    dum(i) = 0.
                 enddo
              endif

           else 

              kb = k**2*(1.-beta)
              dum(ibegin+6) = sqrt(kb*dum(ibegin)**2+bzero**2)
              dum(ibegin+7) = kb/dum(ibegin+6)*dum(ibegin)*dum(ibegin+1)
              dum(ibegin+8) = kb/dum(ibegin+6)*dum(ibegin)*dum(ibegin+2)
              dum(ibegin+9) = kb/dum(ibegin+6)                               &
                   *(-kb/dum(ibegin+6)**2*(dum(ibegin)*dum(ibegin+1))**2     &
                   + dum(ibegin+1)**2+dum(ibegin)*dum(ibegin+3))
              dum(ibegin+10) = kb/dum(ibegin+6)*(-kb/dum(ibegin+6)**2*       &
                   dum(ibegin)**2*dum(ibegin+1)*dum(ibegin+2)                &
                   + dum(ibegin+1)*dum(ibegin+2)+dum(ibegin)*dum(ibegin+4))
              dum(ibegin+11) = kb/dum(ibegin+6)                              &
                   *(-kb/dum(ibegin+6)**2*(dum(ibegin)*dum(ibegin+2))**2     &
                   + dum(ibegin+2)**2+dum(ibegin)*dum(ibegin+5))
              
              if(numvar.ge.3) then
                 kb = k**2*beta
                 dum(ibegin+12) = 0.5*kb*dum(ibegin)**2
                 dum(ibegin+13) = kb*dum(ibegin)*dum(ibegin+1)
                 dum(ibegin+14) = kb*dum(ibegin)*dum(ibegin+2)
                 dum(ibegin+15) = kb*(dum(ibegin+1)**2+                      &
                      dum(ibegin)*dum(ibegin+3))
                 dum(ibegin+16) = kb*(dum(ibegin+1)*dum(ibegin+2)+           &
                      dum(ibegin)*dum(ibegin+4))
                 dum(ibegin+17) = kb*(dum(ibegin+2)**2+                      &
                      dum(ibegin)*dum(ibegin+5))
              endif
           endif
        endif

      case(1) !  Taylor reconnection (itaylor=1)

         dum(ibegin) = -0.5*z**2
         dum(ibegin+1) = 0.
         dum(ibegin+2) = -z
         dum(ibegin+4) = 0
         dum(ibegin+4) = 0
         dum(ibegin+5) = -1

      case(2)  !  force-free taylor state (itaylor=2)
         call xyznod(l, coords)
         x = coords(1)
         z = coords(2)
         alam = sqrt((pi/alx)**2+(pi/alz)**2)

         dum(ibegin) = sin(pi*x/alx)*sin(pi*z/alz)
         dum(ibegin+1) = (pi/alx)*cos(pi*x/alx)*sin(pi*z/alz)
         dum(ibegin+2) = (pi/alz)*sin(pi*x/alx)*cos(pi*z/alz)
         dum(ibegin+3) = -(pi/alx)**2*sin(pi*x/alx)*sin(pi*z/alz)
         dum(ibegin+4) = (pi/alx)*(pi/alz)*cos(pi*x/alx)*cos(pi*z/alz)
         dum(ibegin+5) = -(pi/alz)**2*sin(pi*x/alx)*sin(pi*z/alz)

         dum(ibegin+6) = alam*dum(ibegin)
         dum(ibegin+7) = alam*dum(ibegin+1)
         dum(ibegin+8) = alam*dum(ibegin+2)
         dum(ibegin+9) = alam*dum(ibegin+3)
         dum(ibegin+10) = alam*dum(ibegin+4)
         dum(ibegin+11) = alam*dum(ibegin+5)

      case(3)  ! GEM reconnection (itaylor=3)
         
         dum(ibegin) = 0.5*alog(cosh(2.*z))
         dum(ibegin+1) = 0.0
         dum(ibegin+2) = tanh(2.*z)  
         dum(ibegin+3) = 0.0
         dum(ibegin+4) = 0.0
         dum(ibegin+5) = 2.*(1.-tanh(2.*z)**2)

         if(numvar.ge.2) then
            dum(ibegin+6) = bzero
            do i=ibegin+7,ibegin+11
               dum(i) = 0.
            enddo
         endif

         if(numvar.ge.3) then
            dum(ibegin+12) = 0.5*sech(2*z)**2 + 0.1
            dum(ibegin+13) = 0.
            dum(ibegin+14) = -2.*sech(2*z)**2*tanh(2*z)
            dum(ibegin+15) = 0.
            dum(ibegin+16) = 0.
            dum(ibegin+17) = -4.*sech(2*z)**2*             &
                 &        (sech(2*z)**2-2.*tanh(2*z)**2)
         endif ! if on numvar.ge.3
         
      case(4)  ! wave propagation (itaylor=4)
         
         dum(ibegin) = z
         dum(ibegin+1) = 0.0
         dum(ibegin+2) = 1.0
         dum(ibegin+3) = 0.0
         dum(ibegin+4) = 0.0
         dum(ibegin+5) = 0.0
         
         if(numvar.ge.2) then
            dum(ibegin+6) = bzero
            do i=ibegin+7,ibegin+11
               dum(i) = 0.
            enddo
         endif

      case(100) ! test case
         
         akx = 2.*pi/alx
         akz = 2.*2.*pi/alz

         dum(ibegin  ) = cos(akx*x)*cos(akz*z)
         dum(ibegin+1) =-sin(akx*x)*cos(akz*z)*akx
         dum(ibegin+2) =-cos(akx*x)*sin(akz*z)*akz
         dum(ibegin+3) =-cos(akx*x)*cos(akz*z)*akx**2
         dum(ibegin+4) = sin(akx*x)*sin(akz*z)*akx*akz
         dum(ibegin+5) =-cos(akx*x)*cos(akz*z)*akz**2
      end select
   end do

end subroutine phiequ
!============================================================
subroutine denequ(dum, numberingid)
  use basic
  implicit none

  real, intent(inout) :: dum(*)
  integer, intent(in) :: numberingid

  integer :: l, numnodes, ii, ifail1, ifail2
  integer :: ifail3, ibegin, iendplusone

  real :: x, z
  real :: j0, j1, k, kb, r, ri, arg, befo, ff, s17aef, s17aff
  real :: fp, fpp, d1, d3, d2, d6, d5, d4, sech
  double precision :: coords(3)
  real :: alx, alz, xmin, zmin

  k = 3.8317059702

  call getboundingboxsize(alx, alz)
  call getmincoord(xmin, zmin)

  ! initialize equilibrium to be the Strauss and Longcope solution
  
  ! start a loop over the number of vertices to define initial conditions

  ! note coordinates are centered at 0,0
  call numnod(numnodes)

  print *, 'reached dd', numnodes

  do l=1,numnodes
     call entdofs(numberingid, l, 0, ibegin, iendplusone)
     call xyznod(l,coords)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     if(idens.eq.0) then
        dum(ibegin) = 1.
        dum(ibegin+1) = 0.
        dum(ibegin+2) = 0.
        dum(ibegin+3) = 0.
        dum(ibegin+4) = 0.
        dum(ibegin+5) = 0.
        cycle
     endif

     select case(itaylor)

     case(0)         

        ! for itaylor = 0, use the equilibrium given in R. Richard, et al,
        !     Phys. Fluids B, 2 (3) 1990 p 489

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

        d1 = ff* z
        d3 = fp*z**2*ri + ff
        d2 = fp*z*x *ri
        d6 = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)
        d5 = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
        d4 = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)

        kb = k**2*beta
        dum(ibegin) = 0.5*kb*d1**2

        dum(ibegin+1) = kb*d1*d2
        dum(ibegin+2) = kb*d1*d3
        dum(ibegin+3) = kb*(d2**2+d1*d4)
        dum(ibegin+4) = kb*(d2*d3+d1*d5)
        dum(ibegin+5) = kb*(d3**2+d1*d6)

     case(3) ! GEM reconnection (itaylor=3)
        dum(ibegin) = sech(2*z)**2 + 0.2
        dum(ibegin+1) = 0.
        dum(ibegin+2) = -4.*sech(2*z)**2*tanh(2*z)
        dum(ibegin+3) = 0.
        dum(ibegin+4) = 0.
        dum(ibegin+5) = -8.*sech(2*z)**2*(sech(2*z)**2-2.*tanh(2*z)**2)

     case default
        dum(ibegin) = 1.
        do ii=1,5
           dum(ibegin+ii) = 0.
        enddo
        
     end select
     
     
  end do

end subroutine denequ
!========================
subroutine velinit(dum)

  use basic
  implicit none
  real dum(*), x, z, akx, akz, akx2, omega, befo, alx, alz, xmin, zmin
  integer numnodes, l, ibegin, iendplusone
  double precision coords(3) 

  call numnod(numnodes)
  call getboundingboxsize(alx, alz)
  call getmincoord(xmin, zmin)

  ! initialize perturbed velocity variables

  do l=1,numnodes
     call entdofs(numvar, l, 0, ibegin, iendplusone)
     call xyznod(l,coords)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5
   
     dum(ibegin:iendplusone-1) = 0.

     select case(itaylor)

     case(4) ! wave propagation (itaylor=4)
        akx = 2.*pi/alx
        akx2 = akx**2
        omega = sqrt(akx2 + .5*db**2*akx2**2                                 &
             + db*akx2*sqrt(akx2 + .25*db**2*akx2**2))
        befo = eps*akx/omega
        
        dum(ibegin) = befo*sin(akx*x)
        dum(ibegin+1) = befo*akx*cos(akx*x)
        dum(ibegin+2) =  0.
        dum(ibegin+3) =  -befo*akx**2*sin(akx*x)
        dum(ibegin+4) = 0.
        dum(ibegin+5) =  0.
        
        if(numvar.ge.2) then
           befo = eps*(akx2/omega**2-1)
           
           dum(ibegin+6) = befo*sin(akx*x)
           dum(ibegin+7) = befo*akx*cos(akx*x)
           dum(ibegin+8) =  0.
           dum(ibegin+9) =  -befo*akx**2*sin(akx*x)
           dum(ibegin+10) = 0.
           dum(ibegin+11) =  0.
        endif

     case(100) ! test case
         
        akx = 3.*2.*pi/alx
        akz = 4.*2.*pi/alz
         
        dum(ibegin  ) = cos(akx*x)*cos(akz*z)
        dum(ibegin+1) =-sin(akx*x)*cos(akz*z)*akx
        dum(ibegin+2) =-cos(akx*x)*sin(akz*z)*akz
        dum(ibegin+3) =-cos(akx*x)*cos(akz*z)*akx**2
        dum(ibegin+4) = sin(akx*x)*sin(akz*z)*akx*akz
        dum(ibegin+5) =-cos(akx*x)*cos(akz*z)*akz**2
        
     end select

  end do

end subroutine velinit
!===========================================================================
subroutine phiinit(dum)
  use basic
  implicit none
  
  integer :: l, numnodes, ibegin, iendplusone
  real x, z, akx, akz, akx2, omega, befo, alx, alz, xmin, zmin
  real dum(*)
  double precision :: coords(3)

! start a loop over the number of vertices to define initial conditions

! note coordinates are centered at 0,0
  call getmincoord(xmin, zmin)
  call getboundingboxsize(alx, alz)
  call numnod(numnodes)

 
  do l=1,numnodes
     call xyznod(l,coords)
     call entdofs(numvar, l, 0, ibegin, iendplusone)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     dum(ibegin:iendplusone-1) = 0.

     select case(itaylor)

     case(0,1,2)
     
     case(3) ! GEM reconnection (itaylor=3)
        akx = 2.*pi/alx
        akz = pi/alz

        dum(ibegin) = eps*cos(akx*x)*cos(akz*z)
        dum(ibegin+1) = -eps*akx*sin(akx*x)*cos(akz*z)
        dum(ibegin+2) =  - eps*akz*cos(akx*x)*sin(akz*z)
        dum(ibegin+3) = -eps*akx**2*cos(akx*x)*cos(akz*z)
        dum(ibegin+4) = eps*akx*akz*sin(akx*x)*sin(akz*z)
        dum(ibegin+5) =  - eps*akz**2*cos(akx*x)*cos(akz*z)

     case(4) ! wave propagation (itaylor=4)
        akx = 2.*pi/alx
        akx2 = akx**2
        omega = sqrt(akx2 + .5*db**2*akx2**2                                 &
             + db*akx2*sqrt(akx2 + .25*db**2*akx2**2))

        dum(ibegin) =  eps*sin(akx*x)
        dum(ibegin+1) =  eps*akx*cos(akx*x)
        dum(ibegin+2) =  0.
        dum(ibegin+3) = -eps*akx**2*sin(akx*x)
        dum(ibegin+4) = 0.
        dum(ibegin+5) =  0.
        
        if(numvar.ge.2) then
           befo = eps*(akx/omega-omega/akx)
           dum(ibegin+6) =  befo*sin(akx*x)
           dum(ibegin+7) =  befo*akx*cos(akx*x)
           dum(ibegin+8) =  0.
           dum(ibegin+9) = -befo*akx**2*sin(akx*x)
           dum(ibegin+10) = 0.
           dum(ibegin+11) =  0.
        endif ! on numvar.ge.2
   
     end select

  end do

end subroutine phiinit
!============================================================
subroutine deninit(dum)
  use basic

  implicit none

  real, intent(out) :: dum(*)
  
  integer :: l, numnodes, ibegin, iendplusone
  real :: x, z, alx, alz, xmin, zmin

  double precision :: coords(3)


  ! note coordinates are centered at 0,0
  call numnod(numnodes)
  call getboundingboxsize(alx, alz)
  call getmincoord(xmin, zmin)

  do l=1,numnodes
     call xyznod(l,coords)
     call entdofs(1, l, 0, ibegin, iendplusone)
     x = coords(1) - xmin - alx*.5
     z = coords(2) - zmin - alz*.5

     dum(ibegin) = 0.
     dum(ibegin+1) = 0.
     dum(ibegin+2) = 0.
     dum(ibegin+3) = 0.
     dum(ibegin+4) = 0.
     dum(ibegin+5) = 0.
  end do

end subroutine deninit
