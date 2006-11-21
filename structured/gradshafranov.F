#define REAL64 real

!============================================================
      subroutine gradshafranov
      use p_data
      use t_data
      use basic
      use arrays
      use superlu
      use inout_mod

      implicit none

      real   int1,int4,int2,int3,lhs,cfac(18),g0
#ifndef BIT64
      integer, parameter ::  r8c = selected_real_kind(12,100)
      real(r8c), allocatable::temp(:),b1vecini(:)
#else
<<<<<<< .mine
      real, allocatable::temp(:),b1vecini(:)
=======
      REAL64, allocatable::temp(:)
>>>>>>> .r16
#endif
      integer :: mmnn6, nrags, numvargs, mbandsgs, msizegs
      integer :: m_lgs, nnz_lgs, ifs_rgs, ineg
      integer :: itri, iodd, ll, j, jj, jjj, i, ii, iii, lx, lz, j1, i1
      integer :: l, i0, jjp, jjjp, jp, k, ier, jer, itnum, jsp
      real :: terma, x, z, aminor, fac, bv, termd, term1, term2
      real :: g, gx, gz, gxx, gxz, gzz
      real :: gv, gvx, gvz, gvxx, gvxz, gvzz
      real :: th, sum, count, rhs, ajlim, curr, q0, qstar
      real :: xrel, zrel, error, xguess, zguess
      real :: f

#ifdef mpi
      type(sparseR8d_obj) :: gsmatrix_lu
#else
      type(sparseR8_obj) :: gsmatrix_lu
      integer*8 :: gshandle
#endif
!
      if(iprint.eq.1) write(*,*) "gradshafranov called"
      mmnn6 = m*n*6
      allocate (temp(mmnn6),b1vecini(mmnn6))
      numvargs = 1
      nrags = n*m*6
      call jbdecomp1(n*m, 6*numvargs, m_lgs, nnz_lgs, ifs_rgs)
      mbandsgs = (n+2)*6-1
      msizegs = 9*n*m*(6*numvargs)**2
      base = 0
!
!.....compute LU decomposition only once
!
!.....form matrices
      ss = 0
      itri = 0
      do iodd=1,ioddmx
      do ll=1,nreg
        itri = itri + 1
        do jjj=1,3
        do jj=1,6
          j = (jjj-1)*6 + jj
          do iii=1,3
          do ii=1,6
            i = (iii-1)*6 + ii
            terma = hterm(itri,i,j)
            j1 = 6*ist(itri,jjj)+jj
            i1 = 6*ist(itri,iii)+ii
            call increment(ss,m*n,i1,j1,terma,1,10)
!
          enddo
          enddo
        enddo
        enddo
      enddo
      enddo
!
!.....modify the s-matrix, inserting the boundary conditions
!
!
!......define indices for boundary arrays
      call boundarygs(iboundgs,m,n,nbcgs)
      if(iprint.eq.1) write(*,*) "after call to boundarygs", nbcgs
!
!.....define initial values based on filiment with current tcuro
!     and vertical field of strength bv given by shafranov formula
!.....NOTE:  This formula assumes (li/2 + beta_P) = 1.2
      fac = tcuro/(2.*pi)
!.....minor radius
      aminor = abs(xmag-xlim)
      bv = (1./(4.*pi*xmag))*(alog(8.*xmag/aminor) - 1.5 + 1.2)
      do lx=1,n
      x = xzero + (lx-1)*deex
        do lz=1,m,m-1
        z = zzero + (lz-1)*deez
          l = lx + (lz-1)*n
          i0 = (l-1)*6
          call gvect(x,z,xmag,zmag,1,g,gx,gz,gxx,gxz,gzz,0,ineg)
          call gvect(x,z,102.,xmag,1,gv,gvx,gvz,gvxx,gvxz,gvzz,1,ineg)
          phi(i0+1) = (g  +  gv*bv)*fac
          phi(i0+2) = (gx + gvx*bv)*fac
          phi(i0+3) = (gz + gvz*bv)*fac
          phi(i0+4) = (gxx+gvxx*bv)*fac
          phi(i0+5) = (gxz+gvxz*bv)*fac
          phi(i0+6) = (gzz+gvzz*bv)*fac
          term1 = gxx - gx/x +gzz
          term2 = gvxx-gvx/x +gvzz
      if(lx.eq.1 .or. lx.eq.n .or. lz .eq.1 .or. lz.eq.m) then
      write(88,1088) lx,lz,x,z,phi(i0+1)
      endif
 1088 format(2i3,1p3e12.4)
        enddo
      enddo
      do lx=1,n,n-1
      x = xzero + (lx-1)*deex
        do lz=2,m-1
        z = zzero + (lz-1)*deez
          l = lx + (lz-1)*n
          i0 = (l-1)*6
          call gvect(x,z,xmag,zmag,1,g,gx,gz,gxx,gxz,gzz,0,ineg)
          call gvect(x,z,102.,12.,1,gv,gvx,gvz,gvxx,gvxz,gvzz,1,ineg)
          phi(i0+1) = (g  +  gv*bv)*fac
          phi(i0+2) = (gx + gvx*bv)*fac
          phi(i0+3) = (gz + gvz*bv)*fac
          phi(i0+4) = (gxx+gvxx*bv)*fac
          phi(i0+5) = (gxz+gvxz*bv)*fac
          phi(i0+6) = (gzz+gvzz*bv)*fac
          term1 = gxx - gx/x +gzz
          term2 = gvxx-gvx/x +gvzz
      if(lx.eq.1 .or. lx.eq.n .or. lz .eq.1 .or. lz.eq.m) then
      write(88,1088) lx,lz,x,z,phi(i0+1)
      endif
        enddo
      enddo
      if(myrank.eq.0) call oneplot(0,0,phi,1,1,"phi  ")
!
      do l=1,nbcgs
        i = iboundgs(l)
        psibounds(l) = phi(iboundgs(l))
            if(i.le.0 .or. i.ge.m*n*6*numvargs+1) then
            if(iprint.eq.1) write(*,3344) i,l,iboundgs(l)
 3344       format(" error in iboundgs index",3i5)
            call safestop(3344)
          endif
        do j= 1,nrags
           call assign(ss,n*m,i,j,0.,numvargs)
        enddo
      enddo
      do l=1,nbcgs
          call assign(ss,n*m,iboundgs(l),iboundgs(l),1.,numvargs)
      enddo
! 
      if(iprint.eq.1) write(*,*) "before first call to sparseR8"
!
!....perform LU decomposition of a banded matrix "gsmatrix"
#ifdef mpi
      call sparseR8d_init(gsmatrix_lu,nrags,nnz_lgs,m_lgs,ifs_rgs,ier)
#else
      call sparseR8_init(gsmatrix_lu,nrags,msizegs,base,ier)
       if(ier.ne.0) call safestop(50)
!
!......experiment with different perm_spec
       gsmatrix_lu %permc_spec = 0
#endif
      call define(gsmatrix_lu,ss,m*n,1)
      jer = 0
      if(iprint.eq.1) write(*,*) "before call to dsupralu"
!
#ifdef mpi
       call sparseR8d_new(gsmatrix_lu, jer)
       if(jer.ne.0) call safestop(51)
#else
       call dsupralu_new(gshandle, gsmatrix_lu%values(1),               &
     &      gsmatrix_lu%irow(1), gsmatrix_lu%jcol_ptr(1),               &
     &      gsmatrix_lu%nnz, nrags , jer)
       if(jer.ne.0) call safestop(52)
       if(iprint.eq.1) write(*,*) "before call to dsupralu_colperm"
       call dsupralu_colperm(gshandle, gsmatrix_lu% permc_spec, jer)
       if(jer.ne.0) call safestop(53)
       if(iprint.eq.1) write(*,*) "before call to dsupralu_lu"
       call dsupralu_lu(gshandle, jer)
       if(jer.ne.0) call safestop(54)
#endif
!
      b1vecini = 0
!
!.....define initial b1vecini associated with delta-function source
!     corresponding to current tcuro at location (xmag,zmag)
      xrel = xmag-xzero
      zrel = zmag-zzero
      call deltafun(xrel,zrel,b1vecini,tcuro)
!
      if(iprint.eq.1) write(*,3004)
 3004 format("      curr        xmag        zmag",                      &
     &       "      psimin      psilim      error")
      itnum = 0
      th = 1.
      error = 0.
!
!-------------------------------------------------------------------
!.....start of iteration loop on plasma current
 500  continue
      itnum = itnum + 1
      if(myrank.eq.0) call oneplot(0,0,b1vecini,1,1,"b1vec")
!
!.....apply boundary conditions
      do l=1,nbcgs
         b1vecini(iboundgs(l)) = psibounds(l)
      enddo
!.....perform LU backsubstitution to get psi solution
      if(iprint.eq.1) write(*,*) "before call to solve"
#ifdef mpi
      call sparseR8d_solve(gsmatrix_lu,b1vecini,ier)
#else
      call dsupralu_solve(gshandle,b1vecini(1),ier)
#endif
      do i=1,nrags
      phi(i) = th*b1vecini(i) + (1.-th)*phi(i)
      enddo
      if(myrank.eq.0) call oneplot(0,0,phi,1,1,"phi  ")
      th = 0.9
!
!
!.....calculate the error
      sum = 0.
      count = 0.
      do lx=1,n
        x = xzero + (lx-1)*deex
        do lz=1,m
          l = lx + (lz-1)*n
          i0 = (l-1)*6
          if(phi(i0+1).gt.psilim) go to 202
          lhs = (phi(i0+4)-phi(i0+2)/x+phi(i0+6))/x
          rhs =  -(fun1(i0+1)+gamma2*fun2(i0+1)+                        &
     &             gamma3*fun3(i0+1)+gamma4*fun4(i0+1))
        sum = sum + (lhs-rhs)**2
        count = count + 1
 202    continue
        enddo
      enddo
      error = sqrt(sum/(n*m))
      write(*,*) "gradshafranov: error = ", error
!
!.....calculate psi at the magnetic axis and the limiter
      xguess = xmag - xzero
      zguess = zmag - zzero
      if(iprint.eq.1) write(*,*) "before magaxis",xguess,zguess,                 &
     &                 phi(6*(n*(m-1)/2 + (n-1)/2)+1)
      call magaxis(phi,xguess,zguess)
      xmag = xguess + xzero
      zmag = zguess + zzero
!
      xrel = xlim - xzero
      zrel = zlim - zzero
      call evaluate(xrel,zrel,psilim,ajlim,phi,1,numvargs)
      dpsii = 1./(psilim-psimin)
      if(iprint.eq.1) write(*,*) "xmag,zmag",xmag,zmag
      if(iprint.eq.1) write(*,*) "xlim,zlim",xlim,zlim
      if(iprint.eq.1) write(*,*) "psimin,psilim",psimin,psilim
!
!.....define the pressure and toroidal field functions
      call fundef
!
!.....start of loop over triangles to compute integrals needed to keep
!     total current and q_0 constant using gamma4, gamma2, gamma3
      int1 = 0.
      int4 = 0.
      int2 = 0.
      int3 = 0.
      curr = 0.
      itri = 0.
      do iodd=1,ioddmx
!
!.....calculate common factors
!     calculate matrix elements of the integration matrix fint
      do i=0,5
      do j=0,5
        fint(i,j) = f(i,j,atri(iodd),btri(iodd),ctri(iodd))
      enddo
      enddo
!
      do jjj=1,3
      do jj=1,6
        j = (jjj-1)*6 + jj
        cfac(j) = 0.
        do k=1,20
          cfac(j) = cfac(j) + gtri(k,j,iodd)*fint(mi(k),ni(k))
        enddo
      enddo
      enddo
!
      do ll=1,nreg
      itri = itri + 1
!
        do jjj=1,3
        do jj=1,6
          j = (jjj-1)*6 + jj
          j1 = 6*ist(itri,jjj)+jj
!
            int1 = int1 + cfac(j)*fun1(j1)
            int4 = int4 + cfac(j)*fun4(j1)
            int2 = int2 + cfac(j)*fun2(j1)
            int3 = int3 + cfac(j)*fun3(j1)
              do jjjp=1,3
              do jjp=1,6
                jp = (jjjp-1)*6 + jjp
                jsp= 6*ist(itri,jjjp)+jjp
                curr = curr + sterm(iodd,j,jp)*phi(j1)*rinv(jsp)
              enddo
              enddo
!
        enddo
        enddo
      enddo
      enddo
!.....end of loop to define gamma4, gamma2, gamma3
      write(*,1001) itnum, curr, xmag, zmag, psimin, psilim, error
 1001 format(i5,1p3e12.4,1p2e20.12,1pe12.4)
!
!
!.....choose gamma2 to fix q0/qstar.  Note that there is an additional
!     degree of freedom in gamma3.  Could be used to fix qprime(0)
      q0 = 1.
      qstar = 4.
      g0 = 36.4
      gamma2 =-2.*xmag*(xmag*p0*p1 + (2.*g0/(xmag**2*q0*dpsii)))
      gamma3 = -(xmag*djdpsi/dpsii + 2*xmag**2*p0*p2)
      gamma4 = -(tcuro + gamma2*int2 + gamma3*int3 + int1)/int4
!
!     start loop over elements to define RHS vector
      b1vecini = 0.
      itri = 0
      do iodd=1,ioddmx
      do ll=1,nreg
      itri = itri + 1
!
        do iii=1,3
        do ii=1,6
          i = (iii-1)*6 + ii
          i1 = 6*ist(itri,iii) + ii
          sum = 0.
          do jjj=1,3
          do jj=1,6
            j = (jjj-1)*6 + jj
            j1 = 6*ist(itri,jjj) + jj
            termd = dterm(iodd,i,j)
            sum = sum - termd*(fun1(j1) + gamma4*fun4(j1)               &
     &                + gamma2*fun2(j1) + gamma3*fun3(j1))
          enddo
          enddo
          b1vecini(i1) =  b1vecini(i1) + sum
        enddo
        enddo
      enddo
      enddo
!.....end of loop to define RHS vector 
!
!.....diagnostic plots
      if(itnum.le.30) go to 500
      ntime = itnum
      temp = -(fun1+gamma2*fun2+gamma3*fun3+gamma4*fun4)     
      if(myrank.eq.0) call plotit(phi,temp,0)
      if(myrank.eq.0) call oneplot(0,0,temp,1,1,"temp ")
      if(myrank.eq.0) call oneplot(0,4,fun1,1,1,"Rp-p ")
      if(myrank.eq.0) call oneplot(0,3,fun2,1,1,"G1-p ")
      if(myrank.eq.0) call oneplot(0,2,fun3,1,1,"G2-p ")
      if(myrank.eq.0) call oneplot(0,0,fun4,1,1,"G3-p ")
#ifdef mpi
      call sparseR8d_del(gsmatrix_lu,jer)
#else
      call dsupralu_del(gshandle,jer)
      call sparseR8_free(d2matrix_lu,ier)
#endif
      deallocate (temp)
      return
!
      end
      subroutine gvect(r,z,xi,zi,n,g,gr,gz,grr,grz,gzz,nmult,ineg)
!
!
!
!......calculates derivatives wrt first argument
      dimension r(n),z(n),xi(n),zi(n),g(n),gr(n),gz(n),                 &
     &          grz(n),gzz(n),grr(n)
!
      data a0,a1,a2,a3,a4/1.38629436112,9.666344259e-2,                 &
     &3.590092383e-2,3.742563713e-2,1.451196212e-2/
      data b0,b1,b2,b3,b4/.5,.12498593597,6.880248576e-2,               &
     &3.328355346e-2,4.41787012e-3/
      data c1,c2,c3,c4/.44325141463,6.260601220e-2,                     &
     &4.757383546e-2,1.736506451e-2/
      data d1,d2,d3,d4/.24998368310,9.200180037e-2,                     &
     &4.069697526e-2,5.26449639e-3/
      data pi,tpi/3.1415926535,6.283185308/
!
      if(nmult.gt.0) go to 101
      do 100 i=1,n
      rpxi=r(i)+xi(i)
      rxi=r(i)*xi(i)
      zmzi=z(i)-zi(i)
      rksq=4.*rxi/(rpxi**2+zmzi**2)
      rk=sqrt(rksq)
      sqrxi=sqrt(rxi)
      x=1.-rksq
      ce=1.+x*(c1+x*(c2+x*(c3+x*c4)))+                                  &
     &x*(d1+x*(d2+x*(d3+x*d4)))*(-alog(x))
      ck=a0+x*(a1+x*(a2+x*(a3+x*a4)))+                                  &
     &(b0+x*(b1+x*(b2+x*(b3+x*b4))))*(-alog(x))

      term1=2.*ck-2.*ce-ce*rksq/x
      term2=2.*xi(i)-rksq*rpxi

      g(i) =- sqrxi*(2.*ck-2.*ce-ck*rksq)/rk
      gr(i)=-rk*0.25/sqrxi*(rpxi*term1                                  &
     &+2.*xi(i)*(ce/x-ck))
      gz(i)=-rk*0.25*zmzi/sqrxi*term1
      grz(i)=0.0625*zmzi*(rk/sqrxi)**3*(rpxi*term1                      &
     &+(ce-ck+2.*ce*rksq/x)*                                            &
     &(term2)/x)
      gzz(i)=-rk*0.25/sqrxi*(term1*                                     &
     &(1.-rksq*zmzi**2/(4.*rxi))+zmzi**2*rksq**2/(4.*rxi*x)             &
     &*(ce-ck+2.*ce*rksq/x))
      grr(i)=-rk*0.25/sqrxi*(-rksq*rpxi/(4.*rxi)*                       &
     &(rpxi*term1+2.*xi(i)*(ce/x-ck))+term1-                            &
     &rksq*rpxi/(4.*rxi*x)*(ce-ck+2.*ce*rksq/x)*                        &
     &(term2)+rksq/(2.*r(i)*x)*(2.*ce/x-ck)*term2)
  100 continue
      return
 101  continue
!
!.....check for multipolar coils
      do 200 i=1,n
      if(xi(i) .lt. 100.) go to 200
      rz = zi(i)
      imult = ifix(xi(i) - 100.)
      if(imult .lt. 0 .or. imult.gt.10) go to 250
      imultp = imult + 1
      go to(10,11,12,13,14,15,16,17,18,19,20),imultp
   10 continue
!
!....even nullapole
      g(i) = tpi*rz**2
      gr(i) = 0.
      gz(i) = 0.
      grz(i) = 0.
      gzz(i) = 0.
      grr(i) = 0.
      go to 200
   11 continue
!
!....odd nullapole
      g(i) = 0.
      gr(i) = 0.
      gz(i) = 0.
      grz(i) = 0.
      gzz(i) = 0.
      grr(i) = 0.
      go to 200
   12 continue
!
!....even dipole
      g(i) = tpi*(r(i)**2 - rz**2)/2.
      gr(i) = tpi*r(i)
      gz(i) = 0.
      grz(i) = 0.
      gzz(i) = 0.
      grr(i) = tpi
      go to 200
   13 continue
!
!....odd dipole
      co=tpi/rz
      g(i) = co*(r(i)**2*z(i))
      gr(i) = co*(2.*r(i)*z(i))
      gz(i) = co*(r(i)**2)
      grz(i) = co*2.*r(i)
      gzz(i) = 0.
      grr(i) = co*2*z(i)
      go to 200
   14 continue
!
!....even quadrapole
      co=pi/(4.*rz**2)
      g(i) = co*(r(i)**4-4.*r(i)**2*z(i)**2 - 2.*r(i)**2*rz**2+rz**4)
      gr(i) = co*(4.*r(i)**3-8.*r(i)*z(i)**2-4.*r(i)*rz**2)
      gz(i) = co*(-8.*r(i)**2*z(i))
      grz(i) = co*(-16.*r(i)*z(i))
      gzz(i) = co*(-8.*r(i)**2)
      grr(i) = co*(12.*r(i)**2 - 8.*z(i)**2 - 4.*rz**2)
      go to 200
   15 continue
!
!....odd quadrapole
      co=pi/(3.*rz**3)
      g(i) = co*r(i)**2*z(i)*(3.*r(i)**2-4.*z(i)**2-3.*rz**2)
      gr(i) = co*(12.*r(i)**3*z(i)-8.*r(i)*z(i)**3-6.*r(i)*z(i)*rz**2)
      gz(i) = co*(3.*r(i)**4 - 12.*r(i)**2*z(i)**2 - 3.*r(i)**2*rz**2)
      grz(i) = co*(12.*r(i)**3-24.*r(i)*z(i)**2 - 6.*r(i)*rz**2)
      gzz(i) = co*(-24.*r(i)**2*z(i))
      grr(i) = co*(36.*r(i)**2*z(i)-8.*z(i)**3-6.*z(i)*rz**2)
      go to 200
   16 continue
!
!.....even hexapole
      co=pi/(12.*rz**4)
      g(i) = co*(r(i)**6 - 12.*r(i)**4*z(i)**2 - 3.*r(i)**4*rz**2       &
     &       + 8.*r(i)**2*z(i)**4 + 12.*r(i)**2*z(i)**2*rz**2           &
     &       + 3.*r(i)**2*rz**4 - rz**6 )
      gr(i)= co*(6.*r(i)**5 - 48.*r(i)**3*z(i)**2                       &
     &       - 12.*r(i)**3*rz**2 + 16.*r(i)*z(i)**4                     &
     &       + 24.*r(i)*z(i)**2*rz**2 + 6.*r(i)*rz**4 )
      gz(i)= co*(-24.*r(i)**4*z(i) + 32.*r(i)**2*z(i)**3                &
     &       + 24.*r(i)**2*z(i)*rz**2)
      grz(i)=co*(-96.*r(i)**3*z(i)+64.*r(i)*z(i)**3                     &
     &       + 48.*r(i)*z(i)*rz**2 )
      grr(i)=co*(30.*r(i)**4-144.*r(i)**2*z(i)**2-36.*r(i)**2*rz**2     &
     &       + 16.*z(i)**4 + 24.*z(i)**2*rz**2 + 6.*rz**4)
      gzz(i)=co*(-24.*r(i)**4 + 96.*r(i)**2*z(i)**2 + 24.*r(i)**2*rz**2)
      go to 200
   17 continue
!
!.....odd hexapole
      co=pi/(30.*rz**5)
      g(i) = co*(15.*r(i)**6*z(i) - 60.*r(i)**4*z(i)**3                 &
     &       - 30.*r(i)**4*z(i)*rz**2 + 24.*r(i)**2*z(i)**5             &
     &       + 40.*r(i)**2*z(i)**3*rz**2 + 15.*r(i)**2*z(i)*rz**4)
      gr(i)= co*(90.*r(i)**5*z(i) - 240.*r(i)**3*z(i)**3                &
     &       - 120.*r(i)**3*z(i)*rz**2 + 48.*r(i)*z(i)**5               &
     &       + 80.*r(i)*z(i)**3*rz**2 + 30.*r(i)*z(i)*rz**4)
      gz(i)= co*(15.*r(i)**6 - 180.*r(i)**4*z(i)**2                     &
     &       - 30.*r(i)**4*rz**2 + 120.*r(i)**2*z(i)**4                 &
     &       +120.*r(i)**2*z(i)**2*rz**2 + 15.*r(i)**2*rz**4)
      grz(i)=co*(90.*r(i)**5 - 720.*r(i)**3*z(i)**2                     &
     &       - 120.*r(i)**3*rz**2 + 240.*r(i)*z(i)**4                   &
     &       + 240.*r(i)*z(i)**2*rz**2 + 30.*r(i)*rz**4)
      gzz(i)=co*(-360.*r(i)**4*z(i) + 480.*r(i)**2*z(i)**3              &
     &       + 240.*r(i)**2*z(i)*rz**2)
      grr(i)=co*(450.*r(i)**4*z(i) - 720.*r(i)**2*z(i)**3               &
     &       - 360.*r(i)**2*z(i)*rz**2 + 48.*z(i)**5                    &
     &       + 80.*z(i)**3*rz**2 + 30.*z(i)*rz**4)
      go to 200
   18 continue
!
!....even octapole
      co=pi/(160.*rz**6)
      g(i) = co*(5.*r(i)**8 - 120.*r(i)**6*z(i)**2 - 20.*r(i)**6*rz**2  &
     &       + 240.*r(i)**4*z(i)**4 + 240.*r(i)**4*z(i)**2*rz**2        &
     &       + 30.*r(i)**4*rz**4 - 64.*r(i)**2*z(i)**6                  &
     &       - 160.*r(i)**2*z(i)**4*rz**2 - 120.*r(i)**2*z(i)**2*rz**4  &
     &       - 20.*r(i)**2*rz**6 + 5.*rz**8)
      gr(i)= co*(40.*r(i)**7 - 720.*r(i)**5*z(i)**2 - 120.*r(i)**5*rz**2 &
     &       + 960.*r(i)**3*z(i)**4 + 960.*r(i)**3*z(i)**2*rz**2        &
     &       + 120.*r(i)**3*rz**4 - 128.*r(i)*z(i)**6                   &
     &       - 320.*r(i)*z(i)**4*rz**2 - 240.*r(i)*z(i)**2*rz**4        &
     &       - 40.*r(i)*rz**6)
      gz(i)= co*(-240.*r(i)**6*z(i) + 960.*r(i)**4*z(i)**3              &
     &       + 480.*r(i)**4*z(i)*rz**2 - 384.*r(i)**2*z(i)**5           &
     &       - 640.*r(i)**2*z(i)**3*rz**2 - 240.*r(i)**2*z(i)*rz**4)
      grz(i)=co*(-1440.*r(i)**5*z(i) + 3840.*r(i)**3*z(i)**3            &
     &       + 1920.*r(i)**3*z(i)*rz**2 - 768.*r(i)*z(i)**5             &
     &       - 1280.*r(i)*z(i)**3*rz**2 - 480.*r(i)*z(i)*rz**4)
      gzz(i)=co*(-240.*r(i)**6 + 2880.*r(i)**4*z(i)**2                  &
     &       + 480.*r(i)**4*rz**2 - 1920.*r(i)**2*z(i)**4               &
     &       - 1920.*r(i)**2*z(i)**2*rz**2 - 240.*r(i)**2*rz**4)
      grr(i)=co*(280.*r(i)**6 - 3600.*r(i)**4*z(i)**2                   &
     &       - 600.*r(i)**4*rz**2 + 2880.*r(i)**2*z(i)**4               &
     &       + 2880.*r(i)**2*z(i)**2*rz**2                              &
     &       + 360.*r(i)**2*rz**4 - 128.*z(i)**6                        &
     &       - 320.*z(i)**4*rz**2 - 240.*z(i)**2*rz**4 - 40.*rz**6)
      go to 200
   19 continue
!
!....odd octapole
      co=pi/(140.*rz**7)
      g(i) = co*r(i)**2*z(i)*(35.*r(i)**6 - 280.*r(i)**4*z(i)**2        &
     &       - 105.*r(i)**4*rz**2 + 336.*r(i)**2*z(i)**4                &
     &       + 420.*r(i)**2*z(i)**2*rz**2 + 105.*r(i)**2*rz**4          &
     &       - 64.*z(i)**6 - 168.*z(i)**4*rz**2                         &
     &       - 140.*z(i)**2*rz**4 - 35.*rz**6)
      gr(i)= co*(280.*r(i)**7*z(i) - 1680.*r(i)**5*z(i)**3              &
     &       - 630.*r(i)**5*z(i)*rz**2 + 1344.*r(i)**3*z(i)**5          &
     &       + 1680.*r(i)**3*z(i)**3*rz**2 + 420.*r(i)**3*z(i)*rz**4    &
     &       - 128.*r(i)*z(i)**7 - 336.*r(i)*z(i)**5*rz**2              &
     &       - 280.*r(i)*z(i)**3*rz**4 - 70.*r(i)*z(i)*rz**6)
      gz(i)= co*(35.*r(i)**8-840.*r(i)**6*z(i)**2-105.*r(i)**6*rz**2    &
     &       + 1680.*r(i)**4*z(i)**4 + 1260.*r(i)**4*z(i)**2*rz**2      &
     &       + 105.*r(i)**4*rz**4 - 448.*r(i)**2*z(i)**6                &
     &       - 840.*r(i)**2*z(i)**4*rz**2 - 420.*r(i)**2*z(i)**2*rz**4  &
     &       - 35.*r(i)**2*rz**6)
      grz(i)=co*(280.*r(i)**7 - 5040.*r(i)**5*z(i)**2                   &
     &       - 630.*r(i)**5*rz**2 + 6720.*r(i)**3*z(i)**4               &
     &       + 5040.*r(i)**3*z(i)**2*rz**2 + 420.*r(i)**3*rz**4         &
     &       - 896.*r(i)*z(i)**6 - 1680.*r(i)*z(i)**4*rz**2             &
     &       - 840.*r(i)*z(i)**2*rz**4 - 70.*r(i)*rz**6)
      gzz(i)=co*(-1680.*r(i)**6*z(i) + 6720.*r(i)**4*z(i)**3            &
     &       + 2520.*r(i)**4*z(i)*rz**2 - 2688.*r(i)**2*z(i)**5         &
     &       - 3360.*r(i)**2*z(i)**3*rz**2 - 840.*r(i)**2*z(i)*rz**4)
      grr(i)=co*(1960.*r(i)**6*z(i) - 8400.*r(i)**4*z(i)**3             &
     &       - 3150.*r(i)**4*z(i)*rz**2 + 4032.*r(i)**2*z(i)**5         &
     &       + 5040.*r(i)**2*z(i)**3*rz**2 + 1260.*r(i)**2*z(i)*rz**4   &
     &       - 128.*z(i)**7 - 336.*z(i)**5*rz**2                        &
     &       - 280.*z(i)**3*rz**4 - 70.*z(i)*rz**6)
      go to 200
   20 continue
!
!....even decapole
      co=pi/(560.*rz**8)
      g(i) = co*(7.*r(i)**10 - 280.*r(i)**8*z(i)**2 - 35.*r(i)**8*rz**2 &
     &       + 1120.*r(i)**6*z(i)**4 + 840.*r(i)**6*z(i)**2*rz**2       &
     &       + 70.*r(i)**6*rz**4 - 896.*r(i)**4*z(i)**6                 &
     &       - 1680.*r(i)**4*z(i)**4*rz**2 - 840.*r(i)**4*z(i)**2*rz**4 &
     &       - 70.*r(i)**4*rz**6 + 128.*r(i)**2*z(i)**8                 &
     &       + 448.*r(i)**2*z(i)**6*rz**2 + 560.*r(i)**2*z(i)**4*rz**4  &
     &       + 280.*r(i)**2*z(i)**2*rz**6 + 35.*r(i)**2*rz**8           &
     &       - 7.*rz**10)
      gr(i)= co*(70.*r(i)**9 - 2240.*r(i)**7*z(i)**2                    &
     &       - 280.*r(i)**7*rz**2 + 6720.*r(i)**5*z(i)**4               &
     &       + 5040.*r(i)**5*z(i)**2*rz**2 + 420.*r(i)**5*rz**4         &
     &       - 3584.*r(i)**3*z(i)**6 - 6720.*r(i)**3*z(i)**4*rz**2      &
     &       - 3360.*r(i)**3*z(i)**2*rz**4 - 280.*r(i)**3*rz**6         &
     &       + 256.*r(i)*z(i)**8 + 896.*r(i)*z(i)**6*rz**2              &
     &       + 1120.*r(i)*z(i)**4*rz**4 + 560.*r(i)*z(i)**2*rz**6       &
     &       + 70.*r(i)*rz**8)
      gz(i)= co*(-560.*r(i)**8*z(i) + 4480.*r(i)**6*z(i)**3             &
     &       + 1680.*r(i)**6*z(i)*rz**2 - 5376.*r(i)**4*z(i)**5         &
     &       - 6720.*r(i)**4*z(i)**3*rz**2 - 1680.*r(i)**4*z(i)*rz**4   &
     &       + 1024.*r(i)**2*z(i)**7 + 2688.*r(i)**2*z(i)**5*rz**2      &
     &       + 2240.*r(i)**2*z(i)**3*rz**4 + 560.*r(i)**2*z(i)*rz**6)
      grz(i)=co*(-4480.*r(i)**7*z(i) + 26880.*r(i)**5*z(i)**3           &
     &       + 10080.*r(i)**5*z(i)*rz**2 - 21504.*r(i)**3*z(i)**5       &
     &       - 26880.*r(i)**3*z(i)**3*rz**2 - 6720.*r(i)**3*z(i)*rz**4  &
     &       + 2048.*r(i)*z(i)**7 + 5376.*r(i)*z(i)**5*rz**2            &
     &       + 4480.*r(i)*z(i)**3*rz**4 + 1120.*r(i)*z(i)*rz**6)
      gzz(i)=co*(-560.*r(i)**8 + 13440.*r(i)**6*z(i)**2                 &
     &       + 1680.*r(i)**6*rz**2 - 26880.*r(i)**4*z(i)**4             &
     &       - 20160.*r(i)**4*z(i)**2*rz**2 - 1680.*r(i)**4*rz**4       &
     &       + 7168.*r(i)**2*z(i)**6 + 13440.*r(i)**2*z(i)**4*rz**2     &
     &       + 6720.*r(i)**2*z(i)**2*rz**4 + 560.*r(i)**2*rz**6)
      grr(i)=co*(630.*r(i)**8 - 15680.*r(i)**6*z(i)**2                  &
     &       - 1960.*r(i)**6*rz**2 + 33600*r(i)**4*z(i)**4              &
     &       + 25200.*r(i)**4*z(i)**2*rz**2 + 2100.*r(i)**4*rz**4       &
     &       - 10752.*r(i)**2*z(i)**6 - 20160.*r(i)**2*z(i)**4*rz**2    &
     &       - 10080.*r(i)**2*z(i)**2*rz**4 - 840.*r(i)**2*rz**6        &
     &       + 256.*z(i)**8 + 896.*z(i)**6*rz**2                        &
     &       + 1120.*z(i)**4*rz**4 + 560.*z(i)**2*rz**6                 &
     &       + 70.*rz**8)
      go to 200
  200 continue
      return
  250 continue
!......error
      ineg=39

      return
      end
!============================================================
      subroutine boundarygs(ibound,m,n,nbc)
!
!.....defines the boundary conditions for a rectangular region
      dimension ibound(*)
!
!.....Define Boundary arrays for the grad shafranov equation
      index = 0
!.....start loop on top and bottom boundaries (x-direction)
      do lx=1,n
      index = index+2
      ibound(index-1) = 6*(          lx-1) + 1
      ibound(index  ) = 6*((m-1)*n + lx-1) + 1
!.....tangential derivative
      index = index+2
      ibound(index-1) = 6*(          lx-1) + 2
      ibound(index  ) = 6*((m-1)*n + lx-1) + 2
!.....normal derivative
!      index = index+2
!      ibound(index-1) = 6*(          lx-1) + 3
!      ibound(index  ) = 6*((m-1)*n + lx-1) + 3
!.....tangential second derivatives
      index = index+2
      ibound(index-1) = 6*(          lx-1) + 4
      ibound(index  ) = 6*((m-1)*n + lx-1) + 4
!.....cross derivative
!      index = index+2
!      ibound(index-1) = 6*(          lx-1) + 5
!      ibound(index  ) = 6*((m-1)*n + lx-1) + 5
!.....normal second derivatives
!      index = index+2
!      ibound(index-1) = 6*(          lx-1) + 6
!      ibound(index  ) = 6*((m-1)*n + lx-1) + 6
      enddo
!.....start loop on left and right boundaries (z-direction)
      if(m.le.2) go to 501
      do lz=1,m-2
      index = index+2
      ibound(index-1) = 6*( lz*n      ) + 1
      ibound(index  ) = 6*( lz*n + n-1) + 1
!.....normal derivative
!      index = index+2
!      ibound(index-1) = 6*( lz*n      ) + 2
!      ibound(index  ) = 6*( lz*n + n-1) + 2
!.....tangential derivative
      index = index+2
      ibound(index-1) = 6*( lz*n      ) + 3
      ibound(index  ) = 6*( lz*n + n-1) + 3
!.....normal second derivatives
!      index = index+2
!      ibound(index-1) = 6*( lz*n      ) + 4
!      ibound(index  ) = 6*( lz*n + n-1) + 4
!.....cross derivative
!      index = index+2
!      ibound(index-1) = 6*( lz*n      ) + 5
!      ibound(index  ) = 6*( lz*n + n-1) + 5
!.....tangential second derivatives
      index = index+2
      ibound(index-1) = 6*( lz*n      ) + 6
      ibound(index  ) = 6*( lz*n + n-1) + 6
      enddo
 501  continue
!.....corner points
      index = index+4
      ibound(index-3) =                3
      ibound(index-2) = 6*(n-1     ) + 3
      ibound(index-1) = 6*((m-1)*n ) + 3
      ibound(index  ) = 6*( m*n - 1) + 3
      index = index+4
      ibound(index-3) =                6
      ibound(index-2) = 6*(n-1     ) + 6
      ibound(index-1) = 6*((m-1)*n ) + 6
      ibound(index  ) = 6*( m*n - 1) + 6
!.....cross derivatives
      index = index+4
      ibound(index-3) =                5
      ibound(index-2) = 6*(n-1     ) + 5
      ibound(index-1) = 6*((m-1)*n ) + 5
      ibound(index  ) = 6*( m*n - 1) + 5
      nbc = index
      return
      end
! ================================================
      subroutine axis(dum,xsep,zsep,itime)
      use basic
      use p_data
      use t_data
      use mesh_mod
      dimension xsep(*),zsep(*)
      dimension dum(*),wlocal(20),avector(20)
!
!     itime=0 for first call
!          =1 for subsequent calls
!     
!     locates the magnetic axis and separatrix locations
!
!.....evaluate the solution to get the value [ans] at one point (x,z)

!     loop over all triangles
      nsep = 0
      itrit = 0
      nmax = n-1
      if(iper.eq.1) nmax = n
      do ioddt=1,2
      do jrectt=1,m-1
      do irectt=1,nmax
      iodd = ioddt
      jrect = jrectt
      irect = irectt
      itrit = itrit + 1
      itri = itrit
!
      inews = 0
 300  inews = inews+1
      z1 = (jrect-1)*deez
      x1 = (irect-1)*deex
!
!.....construct the 18 vector corresponding to this triangle
!.....calculate the index and local coordinates for this triangle
      sumx = 0.
      sumz = 0.
      do iii=1,3    
        do ii=1,6
          i = (iii-1)*6 + ii
          index = numvar*6*ist(itri,iii)+ii
          wlocal(i) = dum(index)
        enddo
        index = ist(itri,iii) + 1
        sumx = sumx + xcord(index)
        sumz = sumz + zcord(index)
      enddo
        if(inews.le.1) then
        x = sumx/3.
        z = sumz/3.
                       endif
!
!.....calculate local coordinates
      theta = ttri(iodd)
      b = btri(iodd)
      co = cos(theta)
      sn = sin(theta)
      si  = (x-x1)*co + (z-z1)*sn - b
      eta =-(x-x1)*sn + (z-z1)*co
!
!     calculate the function value corresponding to this point
      do i=1,20
        sum = 0.
        do k=1,18
          sum = sum + gtri(i,k,iodd)*wlocal(k)
        enddo
        avector(i) = sum
      enddo
      inewt = 0
 301  continue
      inewt = inewt + 1
!
!.....evaluate the polynomial and second derivative
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
     &                 term5 = mi(i)*ni(i)*si**(mi(i)-1)*eta**(ni(i)-1)
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
!
      denom =  p22*p11 - p12**2
      sinew = si -  ( p22*pt1 - p12*pt2)/denom
      etanew= eta - (-p12*pt1 + p11*pt2)/denom
      xnew = x1 + co*(b+sinew) - sn*etanew
      znew = z1 + sn*(b+sinew) + co*etanew
!
!.....determine if this new point is nearby
      call whattri(xnew,znew,itrinew,ioddnew,x1,z1)
      if( (xnew-x)**2 + (znew-z)**2 .le. deex**2+deez**2                &
     &   .and. xnew .gt. 0 .and. xnew.lt.alx                             &
     &   .and. znew .gt. 0 .and. znew.lt.alz) then
            x = xnew
            z = znew
            itri = itrinew
            irect = (x/deex) + 1
            jrect = (z/deez) + 1
            iodd = ioddnew
            if(inews.le.4 )go to 300
            if(inews.eq.5) then
!
              if(nsep.eq.0) then
                nsep = 1
                xsep(nsep) = x
                zsep(nsep) = z
              else
                do isep = 1,nsep
                if(sqrt((x-xsep(isep))**2+(z-zsep(isep))**2)            &
     &           .le. .001*deex) go to 400
                enddo
                nsep = nsep + 1
                xsep(nsep) = x
                zsep(nsep) = z
 400            continue
              endif
!
            endif
 3001       format(1p2e12.4)
      endif
!
      enddo
      enddo
      enddo
!
!.....end of loop over all triangles
!     return
      end
!
!============================================================
      subroutine magaxis(phi,xguess,zguess)
      use basic
      use p_data
      use t_data
      use mesh_mod
      dimension phi(*),wlocal(20),avector(20)

      if(myrank.eq.0 .and. iprint.gt.0) then
         write(*,*) "magaxis called with xguess, yguess", xguess, zguess
      endif
!
!     locates the magnetic axis and the value of psi there
!
      itrit = 0
      call whattri(xguess,zguess,itrit,ioddt,x1,z1)
      iodd = ioddt
      jrect = (zguess/deex) + 1
      irect = (xguess/deez) + 1
      itri = itrit
!
      inews = 0
 300  inews = inews+1
      z1 = (jrect-1)*deez
      x1 = (irect-1)*deex
!
!.....construct the 18 vector corresponding to this triangle
!.....calculate the index and local coordinates for this triangle
      sumx = 0.
      sumz = 0.
      do iii=1,3    
        do ii=1,6
          i = (iii-1)*6 + ii
          index = numvar*6*ist(itri,iii)+ii
          wlocal(i) = phi(index)
        enddo
      enddo
        if(inews.le.1) then
        x = xguess
        z = zguess
                       endif
!
!.....calculate local coordinates
      theta = ttri(iodd)
      b = btri(iodd)
      co = cos(theta)
      sn = sin(theta)
      si  = (x-x1)*co + (z-z1)*sn - b
      eta =-(x-x1)*sn + (z-z1)*co
!
!     calculate the function value corresponding to this point
      do i=1,20
        sum = 0.
        do k=1,18
          sum = sum + gtri(i,k,iodd)*wlocal(k)
        enddo
        avector(i) = sum
      enddo
      inewt = 0
 301  continue
      inewt = inewt + 1
!
!.....evaluate the polynomial and second derivative
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
     &                 term5 = mi(i)*ni(i)*si**(mi(i)-1)*eta**(ni(i)-1)
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
!
      if(myrank.eq.0 .and. iprint.gt.0) then
         write(*,*) "pt,pt1,pt2,p11pp22,p12", pt,pt1,pt2,p11,p22,p12
      endif
!
      denom = p22*p11 - p12**2
      sinew = si -  ( p22*pt1 - p12*pt2)/denom
      etanew= eta - (-p12*pt1 + p11*pt2)/denom
      xnew = x1 + co*(b+sinew) - sn*etanew
      znew = z1 + sn*(b+sinew) + co*etanew
      if(iprint.eq.1) write(*,*) "from magaxis, xnew,znew=",xnew,znew
!
!.....determine if this new point is nearby
      call whattri(xnew,znew,itrinew,ioddnew,x1,z1)
      if( (xnew-x)**2 + (znew-z)**2 .le. deex**2+deez**2                &
     &   .and. xnew .gt. 0 .and. xnew.lt.alx                             &
     &   .and. znew .gt. 0 .and. znew.lt.alz) then
            x = xnew
            z = znew
            itri = itrinew
            irect = (x/deex) + 1
            jrect = (z/deez) + 1
            iodd = ioddnew
            if(inews.le.4 )go to 300
            if(inews.eq.5) then
              xguess = x
              zguess = z
              psimin = pt
              return
            endif
 3001       format(1p2e12.4)
      endif
!
!.....error exit
      write(*,3333) inews,x,z,xnew,znew
 3333 format(" error exit from magaxis",i3,1p4e12.4)
      call safestop(27)
!
!
!.....end of loop over all triangles
!     return
      end
      subroutine fundef
      use arrays
      use basic
!
!.....defines the source functions for the GS equation:
!     fun1 = R*p' 
!     fun4 = G1/R
!     fun2 = G2/R
!     fun3 = G3/R
!
      dpsii = 1./(psilim - psimin)
!
      do lx=1,n
      x = xzero + (lx-1)*deex
        do lz=1,m
          l = lx + (lz-1)*n
          i0 = (l-1)*6
          psi =  (phi(i0+1)-psimin)*dpsii
          if(psi .lt. 0 .or. psi .gt. 1) go to 500
          psix = phi(i0+2)*dpsii
          psiy = phi(i0+3)*dpsii
          psixx= phi(i0+4)*dpsii
          psixy= phi(i0+5)*dpsii
          psiyy= phi(i0+6)*dpsii
!
        fbig = p0*dpsii*(p1 + 2.*p2*psi - 3.*(20 + 10*p1+4.*p2)*psi**2  &
     &  + 4.*(45.+20.*p1+6*p2)*psi**3 - 5*(36.+15*p1+4*p2)*psi**4       &
     &  + 6.*(10.+4.*p1+p2)*psi**5)
        fbigp = p0*dpsii*(2.*p2 - 6.*(20 + 10*p1+4.*p2)*psi             &
     &  + 12.*(45.+20.*p1+6*p2)*psi**2 - 20.*(36.+15*p1+4*p2)*psi**3    &
     &  + 30.*(10.+4.*p1+p2)*psi**4)
        fbigpp= p0*dpsii*(- 6.*(20 + 10*p1+4.*p2)                       &
     &  + 24.*(45.+20.*p1+6*p2)*psi - 60.*(36.+15*p1+4*p2)*psi**2       &
     &  + 120.*(10.+4.*p1+p2)*psi**3)
        fun1(i0+1) = x*fbig
        fun1(i0+2) = fbig + x*fbigp*psix
        fun1(i0+3) =        x*fbigp*psiy
        fun1(i0+4) = 2.*fbigp*psix + x*(fbigpp*psix**2+fbigp*psixx)
        fun1(i0+5) = fbigp*psiy + x*(fbigpp*psix*psiy +fbigp*psixy)
        fun1(i0+6) = x*(fbigpp*psiy**2 + fbigp*psiyy)
!
!       g4big =   dpsii*(-psi**2+3.*psi**3-3.*psi**4+psi**5)
!       g4bigp =  dpsii*(-2*psi+9.*psi**2-12.*psi**3+5.*psi**4)
!       g4bigpp = dpsii*(-2 + 18.*psi-36.*psi**2+20*psi**3)
      g4big = dpsii*(-60*psi**2+180*psi**3-180*psi**4+60*psi**5)
      g4bigp= dpsii*(-120*psi+540*psi**2-720*psi**3+300*psi**4)
      g4bigpp=dpsii*(-120   +1080*psi  -2160*psi**2+1200*psi**3)
        fun4(i0+1)= g4big/x
        fun4(i0+2)= g4bigp*psix/x - g4big/x**2
        fun4(i0+3)= g4bigp*psiy/x
        fun4(i0+4)= (g4bigpp*psix**2 + g4bigp*psixx)/x                  &
     &           - 2*g4bigp*psix/x**2 + 2.*g4big/x**3
        fun4(i0+5)= (g4bigpp*psix*psiy+g4bigp*psixy)/x                  &
     &           - g4bigp*psiy/x**2
        fun4(i0+6)=  (g4bigpp*psiy**2 + g4bigp*psiyy)/x
!
        g2big =  dpsii*(1 - 30.*psi**2 + 80.*psi**3                     &
     &        - 75.*psi**4 + 24.*psi**5)
        g2bigp =  dpsii*(-60.*psi + 240.*psi**2                         &
     &        - 300.*psi**3 + 120.*psi**4)
        g2bigpp =  dpsii*(-60. + 480.*psi                               &
     &        - 900.*psi**2 + 480.*psi**3)
        fun2(i0+1)=  g2big/x
        fun2(i0+2)=  g2bigp*psix/x - g2big/x**2
        fun2(i0+3)=  g2bigp*psiy/x
        fun2(i0+4)=  (g2bigpp*psix**2 + g2bigp*psixx)/x                 &
     &           - 2*g2bigp*psix/x**2 + 2.*g2big/x**3
        fun2(i0+5)=(g2bigpp*psix*psiy+g2bigp*psixy)/x                   &
     &           - g2bigp*psiy/x**2
        fun2(i0+6)= (g2bigpp*psiy**2 + g2bigp*psiyy)/x
!
        g3big =  dpsii*(2.*psi - 12.*psi**2 + 24.*psi**3                &
     &        - 20.*psi**4 + 6.*psi**5)
        g3bigp =  dpsii*(2. - 24.*psi + 72.*psi**2                      &
     &        - 80.*psi**3 + 30*psi**4)
        g3bigpp =  dpsii*(- 24. + 144.*psi                              &
     &        - 240.*psi**2 + 120*psi**3)
        fun3(i0+1)= g3big/x
        fun3(i0+2)= g3bigp*psix/x - g3big/x**2
        fun3(i0+3)= g3bigp*psiy/x
        fun3(i0+4)= (g3bigpp*psix**2 + g3bigp*psixx)/x                  &
     &           - 2*g3bigp*psix/x**2 + 2.*g3big/x**3
        fun3(i0+5)= (g3bigpp*psix*psiy+g3bigp*psixy)/x                  &
     &             - g3bigp*psiy/x**2
        fun3(i0+6)=  (g3bigpp*psiy**2 + g3bigp*psiyy)/x
!
        go to 501
 500    continue
        do k=1,6
          fun1(i0+k) = 0.
          fun4(i0+k) = 0.
          fun2(i0+k) = 0.
          fun3(i0+k) = 0.
        enddo
 501    continue
        enddo
      enddo
!
      return
      end
!============================================================
      subroutine deltafun(x,z,b1vecini,val)
      use p_data
      use t_data
      use basic
      use mesh_mod
      dimension b1vecini(*)
!
!.....calculate the source term associated with a delta
!     function source at x,z (relative coordinates)
!
!
!     first determine what rectangle x,z lies within
!      irect = ifix(x/(2*deex)) + 1
!      jrect = ifix(z/(2*deez)) + 1
!      if(irect.gt.(n-1)/2) irect = (n-1)/2
!      if(jrect.gt.(m-1)/2) jrect = (m-1)/2 
!      x1 = ((irect-1)*2+1)*deex
!      z1 = ((jrect-1)*2+1)*deez
!
!
!     next determine what region number
!      ll = irect + (jrect-1)*(n-1)/2
!
!      angle = atan2(z-z1,x-x1)
!      iodd = angle/(pi/4.) + 1
!      if(iodd.gt.8) iodd = 8
!      itri = (iodd-1)*(n-1)*(m-1)/4 + ll
!
      call whattri(x,z,itri,iodd,x1,z1)
!
!
!.....calculate local coordinates
      theta = ttri(iodd)
      b = btri(iodd)
      si  = (x-x1)*cos(theta) + (z-z1)*sin(theta) - b
      eta =-(x-x1)*sin(theta) + (z-z1)*cos(theta)
!
!.....calculate the contribution to b1vecini
      do iii=1,3     
        do ii=1,6
          i = (iii-1)*6 + ii
          index = numvar*6*ist(itri,iii)+ii
          sum = 0.
          do k=1,20
            sum = sum + gtri(k,i,iodd)*si**mi(k)*eta**ni(k)
          enddo
          b1vecini(index) = b1vecini(index) + sum*val
        enddo
      enddo
!
!
      return
      end
c$$$      subroutine deltafun(x0,z0,dum,val)
c$$$      use p_data
c$$$      use t_data
c$$$      use basic
c$$$      use mesh_mod
c$$$
c$$$      implicit none
c$$$
c$$$      real dum, val, x0, z0
c$$$      dimension dum(*)
c$$$
c$$$      real dr, x, z, a
c$$$      integer lx, lz, i0, l
c$$$
c$$$      dr = 0.01
c$$$      a = val/(2.*3.14159265*dr**2)
c$$$
c$$$      do lx=1,n
c$$$         x = (lx-1)*deex
c$$$
c$$$         do lz=1,m+jper
c$$$            z = (lz-1)*deez
c$$$            l = lx + (lz-1)*n
c$$$            i0 = (l-1)*numvar*6
c$$$
c$$$            dum(i0+1)   = a*exp(-0.5*((x-x0)**2 + (z-z0)**2)/dr**2)
c$$$            dum(i0+2) = -(x-x0)*dum(i0+1)/dr**2
c$$$            dum(i0+3) = -(z-z0)*dum(i0+1)/dr**2
c$$$            dum(i0+4) = (((x-x0)/dr)**2-1.)*dum(i0+1)/dr**2
c$$$            dum(i0+5) = (x-x0)*(z-z0)*dum(i0+1)/dr**4
c$$$            dum(i0+6) = (((z-z0)/dr)**2-1.)*dum(i0+1)/dr**2
c$$$         enddo
c$$$      enddo
c$$$
c$$$      end
!============================================================
      subroutine rinvdef(rinv)
      use basic
      dimension rinv(*)
      do lx=1,n
      x = xzero + (lx-1)*deex
        do lz=1,m
          l = lx + (lz-1)*n
          i0 = (l-1)*6
          rinv(i0+1) =  1./x
          rinv(i0+2) = -1./x**2
          rinv(i0+4) =  2./x**3
        enddo
      enddo
      return
      end
