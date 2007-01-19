!============================================================
subroutine gradshafranov
  use p_data
  use t_data
  use basic
  use arrays
  use sparse

  use nintegrate_mod

  implicit none
  
  real   gsint1,gsint4,gsint2,gsint3,lhs,cfac(18)
  real, allocatable::temp(:)

  integer :: itri,i,i1,j,j1,jone, k
  integer :: numelms, numnodes, ibegin, iendplusone
  integer :: ineg, itnum, ier, numvargs
  integer :: itop, iright, ibottom, ileft, izone, izonedim
  real :: dterm(18,18), sterm(18,18), error
  real :: fac, aminor, bv, fintl(-6:maxi,-6:maxi)
  real :: g, gx, gz, gxx, gxz, gzz, g0
  real :: gv, gvx, gvz, gvxx, gvxz, gvzz
  real :: x,z, xmin, zmin, xrel, zrel, xguess, zguess
  real :: th, sum, count, rhs, ajlim, curr, q0, qstar
  real :: dofs(20)
  double precision :: coords(3)
 
  real, dimension(20) :: avec
  
  real :: tstart, tend, tsol, tmagaxis, tfundef, tplot

  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "gradshafranov called"

  tmagaxis = 0.
  tsol = 0.
  tfundef = 0.
  tplot = 0.

  call getmodeltags(ibottom, iright, itop, ileft)
  call getmincoord(xmin, zmin)
  call numnod(numnodes)
  call numfac(numelms)
  
  write(*,*) "this function has not been updated to work in parallel"
  call safestop
  allocate (temp(maxdofs1))
  numvargs = 1
  
  ! compute LU decomposition only once

  ! form matrices
  call zeroarray4solve(gsmatrix_sm,numvar1_numbering)

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

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
        do i=1,18
           j1 = isval1(itri,j)
           i1 = isval1(itri,i)

           sum = int3(ri_79, g79(:,OP_DR,i), g79(:,OP_DR,j),weight_79,79) &
                +int3(ri_79, g79(:,OP_DZ,i), g79(:,OP_DZ,j),weight_79,79)

           call insertval(gsmatrix_sm, -sum, i1,j1,1)

        enddo
     enddo
  enddo

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     write(*,*) "Time defining GS matrix: ", tend - tstart
  endif
  
  ! modify the s-matrix, inserting the boundary conditions
  
  ! define indices for boundary arrays
  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "before call to boundarygs"
  call boundarygs(iboundgs,nbcgs)
  
  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "after call to boundarygs", nbcgs
  
  ! define initial values based on filiment with current tcuro
  ! and vertical field of strength bv given by shafranov formula
  ! NOTE:  This formula assumes (li/2 + beta_P) = 1.2
  fac = tcuro/(2.*pi)
  ! minor radius
  aminor = abs(xmag-xlim)
  bv = (1./(4.*pi*xmag))*(alog(8.*xmag/aminor) - 1.5 + 1.2)

  if(myrank.eq.0 .and. iprint.gt.0) then
     write(*,*) "xmin, xmag, xlim", xmin, xmag, xlim
     write(*,*) "zmin, zmag, zlim", zmin, zmag, zlim
  endif

  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "initializing phi"
  do i=1,numnodes

     call zonenod(i,izone, izonedim)

     if(izonedim .eq. 1 .or. izonedim.eq.0) then
        call xyznod(i,coords)
        x = coords(1) - xmin + xzero
        z = coords(2) - zmin + zzero

        call entdofs(numvargs, i, 0, ibegin, iendplusone)

        call gvect(x,z,xmag,zmag,1,g,gx,gz,gxx,gxz,gzz,0,ineg)
        call gvect(x,z,102.,xmag,1,gv,gvx,gvz,gvxx,gvxz,gvzz,1,ineg)

        phi(ibegin  ) = (g  +  gv*bv)*fac
        phi(ibegin+1) = (gx + gvx*bv)*fac
        phi(ibegin+2) = (gz + gvz*bv)*fac
        phi(ibegin+3) = (gxx+gvxx*bv)*fac
        phi(ibegin+4) = (gxz+gvxz*bv)*fac
        phi(ibegin+5) = (gzz+gvzz*bv)*fac

!        term1 = gxx - gx/x +gzz
!        term2 = gvxx-gvx/x +gvzz        
     endif
  enddo

  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "done initializing phi"

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(myrank.eq.0 .and. maxrank .eq. 1) call oneplot(phi,1,numvargs,"phi ",0)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tplot = tplot + tend - tstart
  endif
  
  psibounds = 0.
  do i=1,nbcgs
     psibounds(i) = phi(iboundgs(i))
     call setdiribc(gsmatrix_sm, iboundgs(i))
  enddo
  call finalizearray4solve(gsmatrix_sm)
 
  write(*,*) "before call to dsupralu", xzero
 
  b1vecini = 0
  
  ! define initial b1vecini associated with delta-function source
  !     corresponding to current tcuro at location (xmag,zmag)
  xrel = xmag-xzero
  zrel = zmag-zzero

  call deltafun(xrel,zrel,b1vecini,tcuro)
  
!!$  write(*,3004)
!!$3004 format("      curr        xmag        zmag",                      &
!!$          "      psimin      psilim      error")
  itnum = 0
  th = 1.
  error = 0.
  
  !-------------------------------------------------------------------
  ! start of iteration loop on plasma current
500 continue

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(myrank.eq.0 .and. maxrank .eq. 1) call oneplot(b1vecini,1,1,"b1vecini",0)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tplot = tplot + tend - tstart
  endif

  ! apply boundary conditions
  do i=1,nbcgs
     b1vecini(iboundgs(i)) = psibounds(i)
  enddo

  itnum = itnum + 1

  ! perform LU backsubstitution to get psi solution
  write(*,*) 'about to solve gradshaf',myrank
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call solve(gsmatrix_sm,b1vecini,ier)
  write(*,*) ' solved gradshaf',myrank
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tsol = tsol + tend - tstart
  endif

  phi = th*b1vecini + (1.-th)*phi

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(myrank.eq.0 .and. maxrank .eq. 1) call oneplot(phi,1,numvargs,"phi ",0)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tplot = tplot + tend - tstart
  endif

  th = 0.9

  ! calculate the error
  sum = 0.
  count = 0.

  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "calculating error:"
  do i=1,numnodes

     call xyznod(i,coords)
     x = coords(1) - xmin + xzero
     z = coords(2) - zmin + zzero
     
     call entdofs(numvargs, i, 0, ibegin, iendplusone)

     if(phi(ibegin).le.psilim) then
        lhs = (phi(ibegin+3)-phi(ibegin+1)/x+phi(ibegin+5))/x
        rhs =  -(fun1(ibegin)+gamma2*fun2(ibegin)+                        &
             gamma3*fun3(ibegin)+gamma4*fun4(ibegin))
        sum = sum + (lhs-rhs)**2
        count = count + 1
     endif
  enddo

  error = sqrt(sum/numnodes)

  if(myrank.eq.0 .and. iprint.gt.0) write(*,*) "gradshafranov: error = ", error
  
  ! calculate psi at the magnetic axis and the limiter
  xguess = xmag - xzero
  zguess = zmag - zzero

  write(*,*) xmag, xzero, zmag, zzero

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call magaxis(phi,xguess,zguess)
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tmagaxis = tmagaxis + tend - tstart
  endif


  xmag = xguess + xzero
  zmag = zguess + zzero
  
  xrel = xlim - xzero
  zrel = zlim - zzero
  call evaluate(xrel,zrel,psilim,ajlim,phi,1,numvargs,0)
  dpsii = 1./(psilim-psimin)
  
  ! define the pressure and toroidal field functions
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call fundef
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tfundef = tfundef + tend - tstart
  endif
  
  ! start of loop over triangles to compute integrals needed to keep
  !     total current and q_0 constant using gamma4, gamma2, gamma3
  gsint1 = 0.
  gsint4 = 0.
  gsint2 = 0.
  gsint3 = 0.
  curr = 0.
  itri = 0.
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
        i1 = isvaln(itri,i)
        
        gsint1 = gsint1 + cfac(i)*fun1(i1)
        gsint4 = gsint4 + cfac(i)*fun4(i1)
        gsint2 = gsint2 + cfac(i)*fun2(i1)
        gsint3 = gsint3 + cfac(i)*fun3(i1)
        do j=1,18
           jone = isval1(itri,j)
           curr = curr + sterm(i,j)*phi(i1)*rinv(jone)
        enddo
     enddo

  enddo

  ! end of loop to define gamma4, gamma2, gamma3
  write(*,1001) itnum, curr, xmag, zmag, psimin, psilim, error
1001 format(i5,1p3e12.4,1p2e20.12,1pe12.4) 

  ! choose gamma2 to fix q0/qstar.  Note that there is an additional
  ! degree of freedom in gamma3.  Could be used to fix qprime(0)
  q0 = 1.
  qstar = 4.
  g0 = 36.4
!  gamma2 =-xmag*(xmag*p0*p1 + (tcuro*qstar/(pi*aminor**2*q0*dpsii)))
!  gamma3 = -(.5*djdpsi/dpsii + p0*p2)
  gamma2 =- 2.*xmag*(xmag*p0*p1 + (2.*g0/(xmag**2*q0*dpsii)))
  gamma3 = -(xmag*djdpsi/dpsii + 2*xmag**2*p0*p2)
  gamma4 = -(tcuro + gamma2*gsint2 + gamma3*gsint3 + gsint1)/gsint4
  
  if(myrank.eq.0 .and. iprint.gt.0) then
     write(*,*) "gsint1, gsint2, gsint3, gsint4", gsint1, gsint2, gsint3, gsint4
     write(*,*) "gamma2, gamma3, gamma4", gamma2, gamma3, gamma4
  endif

  ! start loop over elements to define RHS vector
  b1vecini = 0.

  do itri=1,numelms

     call calcfint(fintl,maxi,atri(itri),btri(itri),ctri(itri))
     call calcdterm(itri, dterm, fintl)

     do i=1,18
        i1 = isvaln(itri,i)
        sum = 0.

        do j=1,18
           j1 = isvaln(itri,j)
           
           sum = sum - dterm(i,j)*(fun1(j1) + gamma4*fun4(j1)               &
                +  gamma2*fun2(j1) + gamma3*fun3(j1))
        enddo

        b1vecini(i1) =  b1vecini(i1) + sum
     enddo
  enddo

  ! end of loop to define RHS vector 
  !
  ! diagnostic plots
  if(itnum.le.10) go to 500
  ntime = itnum
  temp = -(fun1+gamma2*fun2+gamma3*fun3+gamma4*fun4)     

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(maxrank .eq. 1) then
     call plotit(phi,temp,0)
     call oneplot(temp,1,1,"temp",0)
     if(myrank.eq.0) call oneplot(fun1,1,1,"fun1",0)
     if(myrank.eq.0) call oneplot(fun2,1,1,"fun2",0)
     if(myrank.eq.0) call oneplot(fun3,1,1,"fun3",0)
     if(myrank.eq.0) call oneplot(fun4,1,1,"fun4",0)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        tplot = tplot + tend - tstart
     endif
  endif

  deallocate (temp)

  if(myrank.eq.0 .and. itimer.eq.1) then
     write(*,*) "gradshafranov: Time spent in magaxis: ", tmagaxis
     write(*,*) "gradshafranov: Time spent in fundef: ", tfundef
     write(*,*) "gradshafranov: Time spent in solve: ", tsol
     write(*,*) "gradshafranov: Time spent in plot: ", tplot
  endif

  return

end subroutine gradshafranov
!==================================
subroutine magaxis(phi,xguess,zguess)
  use basic
  use p_data
  use t_data
  use nintegrate_mod

  implicit none

  real, intent(in) :: phi(*)
  real, intent(inout) :: xguess, zguess

  integer :: itri, itrit, jrect, irect, inews, iii
  integer :: ii, i, index, k, inewt, itrinew, ibegin, iendplusone
  integer :: numvargs
  real :: x1, z1, x, z, theta, b, co, sn, si, eta
  real :: sum, sum1, sum2, sum3, sum4, sum5
  real :: term1, term2, term3, term4, term5
  real :: pt, pt1, pt2, p11, p22, p12, denom, sinew, etanew, xnew
  real :: znew, deex, alx, alz
  real, dimension(20) :: wlocal, avector

  !     locates the magnetic axis and the value of psi there

  numvargs = 1

  if(myrank.eq.0 .and. iprint.gt.0) then
     write(*,*) "magaxis called with xguess, yguess", xguess, zguess
  endif

  itrit = 0
  call safestop(9832)
  call whattri(xguess,zguess,itrit,x1,z1)

  call getboundingboxsize(alx, alz)
  call getdeex(itrit,deex)
  itri = itrit

  inews = 0
300 inews = inews+1

  call calcavector(itri, phi, 1, numvargs, avector)
  
  if(inews.le.1) then
     x = xguess
     z = zguess
  endif

  ! calculate local coordinates
  theta = ttri(itri)
  b = btri(itri)
  co = cos(theta)
  sn = sin(theta)
  si  = (x-x1)*co + (z-z1)*sn - b
  eta =-(x-x1)*sn + (z-z1)*co
  
  inewt = 0
301 continue
  inewt = inewt + 1

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
  pt1  = sum1
  pt2  = sum2
  p11 = sum3
  p22 = sum4
  p12 = sum5

  if(myrank.eq.0 .and. iprint.gt.0) then
     write(*,*) "pt,pt1,pt2,p11pp22,p12", pt,pt1,pt2,p11,p22,p12
  endif

!  denom = p22*p11 - p12**2
!  sinew = si -  ( p22*pt1 - p12*pt2)/denom
!  etanew= eta - (-p12*pt1 + p11*pt2)/denom
  denom = 4.*p22*p11 - p12**2
  sinew = si -  ( 2.*p22*pt1 - p12*pt2)/denom
  etanew= eta - (-p12*pt1 + 2.*p11*pt2)/denom
  xnew = x1 + co*(b+sinew) - sn*etanew
  znew = z1 + sn*(b+sinew) + co*etanew
  
  if(myrank.eq.0 .and. iprint.gt.0) then
     write(*,*), "Found maximum at: ", xnew, znew
  endif

  ! determine if this new point is nearby
  call safestop(482)
  call whattri(xnew,znew,itrinew,x1,z1)
  if( (xnew-x)**2 + (znew-z)**2 .le. 2.*deex**2                     &
       .and. xnew .gt. 0 .and. xnew.lt.alx                            &
       .and. znew .gt. 0 .and. znew.lt.alz) then
     x = xnew
     z = znew
     itri = itrinew
     irect = (x/deex) + 1
     jrect = (z/deex) + 1
     if(inews.le.4 )go to 300
     if(inews.eq.5) then
        xguess = x
        zguess = z
        psimin = pt
        return
     endif
3001 format(1p2e12.4)
  endif

  ! error exit
  write(*,3333) inews,x,z,xnew,znew
3333 format(" error exit from magaxis",i3,1p4e12.4)
  call safestop(27)

end subroutine magaxis

!============================================================
subroutine gvect(r,z,xi,zi,n,g,gr,gz,grr,grz,gzz,nmult,ineg)
  ! calculates derivatives wrt first argument

  implicit none

  integer, intent(in) :: n, nmult
  integer, intent(out) :: ineg
  real, dimension(n), intent(in) :: r, z, xi, zi
  real, dimension(n), intent(out) :: g,gr,gz,grr,grz,gzz
  
!!$  dimension r(n),z(n),xi(n),zi(n),g(n),gr(n),gz(n),                 &
!!$       grz(n),gzz(n),grr(n)

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
  
  integer, intent(out) :: ibound(*),nbc

  integer numnodes, i, izone, iplace, izonedim, j
  integer ibottom, iright, itop, ileft, ibegin, iendplusone
  
  call getmodeltags(ibottom, iright, itop, ileft)
  call numnod(numnodes)
  nbc = 0

  do i=1,numnodes

     call zonenod(i,izone, izonedim)
     call entdofs(1, i, 0, ibegin, iendplusone)

     if(izonedim .eq. 1) then ! an edge...
        if(izone .eq. itop .or. izone .eq. ibottom) then
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+1
           ibound(nbc+3) = ibegin+3
           nbc = nbc+3
        else
           ibound(nbc+1) = ibegin
           ibound(nbc+2) = ibegin+2
           ibound(nbc+3) = ibegin+5
           nbc = nbc+3
        endif
     else if(izonedim .eq. 0) then
        ibound(nbc+1) = ibegin
        ibound(nbc+2) = ibegin+1
        ibound(nbc+3) = ibegin+2
        ibound(nbc+4) = ibegin+3
        ibound(nbc+5) = ibegin+4
        ibound(nbc+6) = ibegin+5
        nbc = nbc + 6
     endif
  enddo

end subroutine boundarygs

subroutine deltafun(x,z,dum,val)

  use p_data
  use t_data
  use basic

  implicit none

  real, intent(in) :: x, z, val
  real, intent(out) :: dum(*)

  integer :: itri, i, ii, iii, k, ibegin, iendplusone, index
  real :: x1, z1, b, theta, si, eta, sum
  
  call safestop(4421)
  call whattri(x,z,itri,x1,z1)

  ! calculate local coordinates
  theta = ttri(itri)
  b = btri(itri)
  si  = (x-x1)*cos(theta) + (z-z1)*sin(theta) - b
  eta =-(x-x1)*sin(theta) + (z-z1)*cos(theta)

  ! calculate the contribution to b1vecini
  do iii=1,3     
     call entdofs(1, ist(itri, iii)+1, 0, ibegin, iendplusone)
     do ii=1,6
        i = (iii-1)*6 + ii
        index = ibegin+ii-1
        sum = 0.
        do k=1,20
           sum = sum + gtri(k,i,itri)*si**mi(k)*eta**ni(k)
        enddo
        dum(index) = dum(index) + sum*val
     enddo
  enddo

end subroutine deltafun
!============================================================
subroutine fundef

!.....defines the source functions for the GS equation:
!     fun1 = R*p' 
!     fun4 = G1/R
!     fun2 = G2/R
!     fun3 = G3/R

  use arrays
  use basic
  
  implicit none

  integer :: l, numnodes, k, ibegin, iendplusone
  real :: x, z, xmin, zmin, psi, psix, psiy, psixx, psixy, psiyy, fbig
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

     call entdofs(1, l, 0, ibegin, iendplusone)
     psi =  (phi(ibegin)-psimin)*dpsii
     if(psi .lt. 0 .or. psi .gt. 1) go to 500
     psix = phi(ibegin+1)*dpsii
     psiy = phi(ibegin+2)*dpsii
     psixx= phi(ibegin+3)*dpsii
     psixy= phi(ibegin+4)*dpsii
     psiyy= phi(ibegin+5)*dpsii
     
     fbig = p0*dpsii*(p1 + 2.*p2*psi - 3.*(20 + 10*p1+4.*p2)*psi**2       &
          + 4.*(45.+20.*p1+6*p2)*psi**3 - 5*(36.+15*p1+4*p2)*psi**4       &
          + 6.*(10.+4.*p1+p2)*psi**5)
     fbigp = p0*dpsii*(2.*p2 - 6.*(20 + 10*p1+4.*p2)*psi                  &
          + 12.*(45.+20.*p1+6*p2)*psi**2 - 20.*(36.+15*p1+4*p2)*psi**3    &
          + 30.*(10.+4.*p1+p2)*psi**4)
     fbigpp= p0*dpsii*(- 6.*(20 + 10*p1+4.*p2)                            &
          + 24.*(45.+20.*p1+6*p2)*psi - 60.*(36.+15*p1+4*p2)*psi**2       &
          + 120.*(10.+4.*p1+p2)*psi**3)
     fun1(ibegin) = x*fbig
     fun1(ibegin+1) = fbig + x*fbigp*psix
     fun1(ibegin+2) =        x*fbigp*psiy
     fun1(ibegin+3) = 2.*fbigp*psix + x*(fbigpp*psix**2+fbigp*psixx)
     fun1(ibegin+4) = fbigp*psiy + x*(fbigpp*psix*psiy +fbigp*psixy)
     fun1(ibegin+5) = x*(fbigpp*psiy**2 + fbigp*psiyy)

!!$     g4big =   dpsii*(-psi**2+3.*psi**3-3.*psi**4+psi**5)
!!$     g4bigp =  dpsii*(-2*psi+9.*psi**2-12.*psi**3+5.*psi**4)
!!$     g4bigpp = dpsii*(-2 + 18.*psi-36.*psi**2+20*psi**3)
     g4big = dpsii*(-60*psi**2+180*psi**3-180*psi**4+60*psi**5)
     g4bigp= dpsii*(-120*psi+540*psi**2-720*psi**3+300*psi**4)
     g4bigpp=dpsii*(-120   +1080*psi  -2160*psi**2+1200*psi**3)

     fun4(ibegin)  = g4big/x
     fun4(ibegin+1)= g4bigp*psix/x - g4big/x**2
     fun4(ibegin+2)= g4bigp*psiy/x
     fun4(ibegin+3)= (g4bigpp*psix**2 + g4bigp*psixx)/x                  &
          - 2*g4bigp*psix/x**2 + 2.*g4big/x**3
     fun4(ibegin+4)= (g4bigpp*psix*psiy+g4bigp*psixy)/x                  &
          - g4bigp*psiy/x**2
     fun4(ibegin+5)=  (g4bigpp*psiy**2 + g4bigp*psiyy)/x

     g2big =  dpsii*(1 - 30.*psi**2 + 80.*psi**3                     &
          - 75.*psi**4 + 24.*psi**5)
     g2bigp =  dpsii*(-60.*psi + 240.*psi**2                         &
          - 300.*psi**3 + 120.*psi**4)
     g2bigpp =  dpsii*(-60. + 480.*psi                               &
          - 900.*psi**2 + 480.*psi**3)
     fun2(ibegin)  =  g2big/x
     fun2(ibegin+1)=  g2bigp*psix/x - g2big/x**2
     fun2(ibegin+2)=  g2bigp*psiy/x
     fun2(ibegin+3)=  (g2bigpp*psix**2 + g2bigp*psixx)/x                 &
          - 2*g2bigp*psix/x**2 + 2.*g2big/x**3
     fun2(ibegin+4)=(g2bigpp*psix*psiy+g2bigp*psixy)/x                   &
          - g2bigp*psiy/x**2
     fun2(ibegin+5)= (g2bigpp*psiy**2 + g2bigp*psiyy)/x

     g3big =  dpsii*(2.*psi - 12.*psi**2 + 24.*psi**3                &
          - 20.*psi**4 + 6.*psi**5)
     g3bigp =  dpsii*(2. - 24.*psi + 72.*psi**2                      &
          - 80.*psi**3 + 30*psi**4)
     g3bigpp =  dpsii*(- 24. + 144.*psi                              &
          - 240.*psi**2 + 120*psi**3)
     fun3(ibegin)= g3big/x
     fun3(ibegin+1)= g3bigp*psix/x - g3big/x**2
     fun3(ibegin+2)= g3bigp*psiy/x
     fun3(ibegin+3)= (g3bigpp*psix**2 + g3bigp*psixx)/x                  &
          - 2*g3bigp*psix/x**2 + 2.*g3big/x**3
     fun3(ibegin+4)= (g3bigpp*psix*psiy+g3bigp*psixy)/x                  &
          - g3bigp*psiy/x**2
     fun3(ibegin+5)=  (g3bigpp*psiy**2 + g3bigp*psiyy)/x

     go to 501
500  continue
     do k=0,5
        fun1(ibegin+k) = 0.
        fun4(ibegin+k) = 0.
        fun2(ibegin+k) = 0.
        fun3(ibegin+k) = 0.
     enddo
501  continue
  enddo
  
  return
end subroutine fundef
