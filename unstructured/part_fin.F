      subroutine mask(x,z,factor)
      use basic
      implicit none
      real factor, r1, r2, alpha, rmid, rnorm, r, x, z, alx, alz,
     &     xmin, zmin
      factor = 1.
! acbauer -- changed to use xmin and zmin since it apears that 
!            this routine may want to center the domain around (0,0)

      call getboundingboxsize(alx, alz)
      call getmincoord(xmin, zmin)
!
!     radius where the factor becomes 0
      r1 = 2.0
!
!     radius where the factor becomes 1
      r2 = 1.5
      alpha = 6./(r1 - r2)
      rmid = (r1+r2)/2.
!
      rnorm = sqrt(alx*alz)/4.
      r = sqrt( (x-xmin-alx*.5)**2 + (z-zmin-alz*.5)**2)/rnorm
!
      factor = (1 + tanh(-alpha*(r-rmid)))*.5
!
!     write(66,1001) x,z,r,factor
 1001 format(1p4e12.4)
      return
      end
!============================================================
      subroutine space
!
!     allocates space for big arrays
      use p_data
      use t_data
      use basic
      use arrays
      use nintegrate_mod

      implicit none
      integer numnodes, nbound, numelms, i, izone, izonedim 

      call numnod(numnodes)
      call numfac(numelms)
      ntri   = 2*numelms !should this be ntri = numelms???

      nbound = 0 ! number of boundary vertices
      do i=1,numnodes
         call zonenod(i, izone, izonedim)
         if(izonedim .ne. 2) nbound = nbound + 1
      enddo
      nbound = 12*nbound*numvar
!
!.....arrays defined at all vertices
c     created with ordering specified by numvar
      call createvec(b1vecini, numvar)
      call createvec(vel, numvar)
      call createvec(vels, numvar)
      call createvec(velold, numvar)
      call createvec(veln, numvar)
      call createvec(denold, numvar)
      call createvec(vel0, numvar)
      call createvec(vel1, numvar)
      call createvec(b2vecini, numvar)
      call createvec(phi, numvar)
      call createvec(phis, numvar)
      call createvec(phiold, numvar)
      call createvec(phi0, numvar)
      call createvec(phi1, numvar)
      call createvec(b1vector, numvar)
      call createvec(b2vector, numvar)
      call createvec(b3vector, numvar)
      call createvec(vtemp, numvar)
      call createvec(r4, numvar)
      call createvec(q4, numvar)
      call createvec(fun1, 1)
      call createvec(fun4, 1)
      call createvec(fun2, 1)
      call createvec(fun3, 1)
      call createvec(jphi, 1)
      jphi = 0.
      call createvec(jphi0, 1)
      call createvec(sb1, 1)
      call createvec(sb2, 1)
      call createvec(sb3, 1)
      call createvec(sp1, 1)
      call createvec(vor, 1)
      call createvec(vor0, 1)
      call createvec(com, 1)
      call createvec(com0, 1)
      call createvec(den, 1)
      call createvec(den0, 1)
      call createvec(deni, 1)
      call createvec(pres, 1)
      call createvec(pres0, 1)
      call createvec(qn4, 1)
!
!.....boundary arrays
      allocate(iboundgs(nbound),iboundv(nbound),iboundv2(nbound),       &
     &     iboundp(nbound),iboundn(nbound),                             &
     &     psibounds(nbound),velbounds(nbound))
!
!.....arrays associated with the triangles
      allocate(ist(ntri,3),isvaln(ntri,18),isval1(ntri,18),             &
     &         isval2(ntri,18))
!
!.....metric terms involved in toroidal geometry
      allocate(rinv(numnodes*6))
!
      return
      end
!=============================
      subroutine smoother3(inarray,outarray,numnodes,numvard,iplace)
!
!     This routine takes as input the Laplacian of the scalar field in
!     the array inarray, smooths it, and returns the smoothed scalar array
!     in location iplace in the stand-alone array outarray(mmnn6).  The LU
!     decomposition takes place the first-time only.
!
!
      use p_data
      use t_data
      use basic
      use arrays
      use sparse
      implicit none
      integer numnodes, numvard,iplace, ndofs_sm, numelms,
     &     ibegind, iendplusoned, ibeginsm, iendplusonesm
      real ssterm(3,3)
      real  outarray(numnodes*6*numvard),inarray(numnodes*6),
     &     aterm(18,18), bterm(18,18),fintl(-6:maxi,-6:maxi),
     &     dterm(18,18)
      real, allocatable:: temp(:)
      integer,  save, allocatable :: iboundsm(:), iboundsm2(:)
      real hyp, thimpsm
      integer numvarsm, itri, i, j,izone,izonedim,
     1     l, jone, ier, jtwo, itwo, nbcsm, nbound
!
      numvarsm = 2
      call numfac(numelms)
c$$$      call createvec(temp, numvarsm)
c$$$      call numdofs(numvarsm, ndofs_sm)
      call numdofs(numvarsm, ndofs_sm)
      call createvec(temp, numvarsm)
      temp = 0
      thimpsm = 1.
! 
      if(ifirsts7_lu.ne.0) go to 500
      ifirsts7_lu = 1

      nbound = 0
      do i=1,numnodes
         call zonenod(i,izone,izonedim)
         if(izonedim .ne. 2) nbound = nbound + 1
      enddo
      nbound = 12*numvar*nbound
      allocate(iboundsm(2*nbound),iboundsm2(2*nbound))
!
!
!.....compute LU decomposition only once
!
!.....form matrix
      call zeroarray4solve(s7matrix_sm,numvar2_numbering)
      do itri=1,numelms
!         call getdeex(itri,deex)
         call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
         call calcaterm(itri, aterm, fintl)
         call calcbterm(itri, bterm, fintl)
         call calcdterm(itri, dterm, fintl)
         hyp = hyperc*deex**2*amuc
         do i=1,18
            do j=1,18
              jtwo = isval2(itri,j)
              itwo = isval2(itri,i)
              ssterm(1,1) = dterm(i,j)
              ssterm(1,2) =-aterm(i,j)
              ssterm(2,1) =  +dt*hyp*thimpsm*bterm(i,j)
              ssterm(2,2) = aterm(i,j) - regular*dterm(i,j)
              call insertval(s7matrix_sm,ssterm(1,1),itwo  ,jtwo  ,1)
              call insertval(s7matrix_sm,ssterm(1,2),itwo  ,jtwo+6,1)
              call insertval(s7matrix_sm,ssterm(2,1),itwo+6,jtwo  ,1)
              call insertval(s7matrix_sm,ssterm(2,2),itwo+6,jtwo+6,1)
           enddo
        enddo
      enddo
!
!......define indices for boundary arrays
      call boundarysm(iboundsm,iboundsm2,nbcsm,iplace)
      ! bc stuff...
      do i=1,nbcsm
         call setdiribc(s7matrix_sm, iboundsm(i))
c$$$         if(iboundsm2(i) .ne. 0) then
c$$$            call insertval(s7matrix_sm,1.,iboundsm(i),iboundsm2(i),0)
c$$$         endif
      enddo
      call finalizearray4solve(s7matrix_sm)
!
 500  continue
!
!.....define RHS vector
      do itri=1,numelms
         call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
         call calcdterm(itri, dterm, fintl)
!         call getdeex(itri,deex)
         hyp = hyperc*deex**2*amu
         do i=1,18
            do j=1,18
               jone = isval1(itri,j)
               itwo = isval2(itri,i)

               temp(itwo+6) = temp(itwo+6)+(dterm(i,j)               
     &              - dt*hyp*(1.-thimpsm)*bterm(i,j))*inarray(jone)

            enddo
         enddo
      enddo
      
      do l=1,nbcsm
      if(iboundsm(l).le.0 .or. iboundsm(l).gt.ndofs_sm) then
         write(*,*) "error in smoother"
         call safestop(98)
      endif
      temp(iboundsm(l)) = 0.
      enddo
      call sumshareddofs(temp)
!
!.....perform LU backsubstitution to get outarray solution
      call solve(s7matrix_sm,temp,ier)
      if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "after solve s7"
!
!.....store solution
      do l=1,numnodes
         call entdofs(numvard, l, 0, ibegind, iendplusoned)
         call entdofs(numvarsm, l, 0, ibeginsm, iendplusonesm)
         do i=0,5
            outarray(ibegind+6*(iplace-1)+i) = 
     &           temp(ibeginsm+6+i)
         enddo
      enddo
!
      call deletevec(temp)
      return
!
      end
!============================================================
      subroutine scaleback
!
      use p_data
      use t_data
      use basic
      use arrays
      implicit none
      real alpha, alphas
      integer i
!.....scale solution back if a linear calculation
      alpha = 1.
      if(ekin .gt. .001) alpha = sqrt(.001/ekin)
      vel = alpha*vel
c      do i=301,360
c         write(*,*) 'acbauer is not sure why this is in scaleback'
c         call safestop(334)
c         vel(i) = 1.
c      enddo
      phi = alpha*phi
!
      alphas = alpha*alpha
      ekin = ekin*alphas
      emag = emag*alphas
      write(*,*) 'emag===',emag, alphas
      ekind = ekind*alphas
      emagd = emagd*alphas
!
      ekinp = ekinp*alphas
      emagp = emagp*alphas
      ekinpd = ekinpd*alphas
      emagpd = emagpd*alphas
!
!
      ekint = ekint*alphas
      emagt = emagt*alphas
      ekintd = ekintd*alphas
      emagtd = emagtd*alphas
!
      ekinph = ekinph*alphas
      ekinth = ekinth*alphas
      emagph = emagph*alphas
      emagth = emagth*alphas
!
      return
      end
!===========================
      subroutine fun(x,z,df,alx,alz)
      implicit none
      integer i, j
      real x, z, df, alx, alz, pi, anorm, akx, akz, snx, cox, snz, coz
      dimension df(0:4, 0:4)
!
!.....source term and derivatives for initial perturbation
!     first dimension of df is number of derivatives, second is
!     number of x derivatives
!
      pi = acos(-1.)
      anorm = ((pi/alx)**2+(pi/alz)**2)
!
      akx = pi/alx
      akz = pi/alz
!
      snx = sin(akx*x)
      cox = cos(akx*x)
      snz = sin(akz*z)
      coz = cos(akz*z)
!
      do i=0,4
        do j=0,4
          df(i,j) = 0.
        enddo
      enddo
!.....first index is total derivatives, second is x derivative
      df(0,0) = anorm       *snx*snz
!
      df(1,0) = anorm       *akz   *snx*coz
      df(1,1) = anorm*akx          *cox*snz
!
      df(2,0) =-anorm       *akz**2*snx*snz
      df(2,1) = anorm*akx*   akz   *cox*coz
      df(2,2) =-anorm*akx**2       *snx*snz      
!
      df(3,0) =-anorm*       akz**3*snx*coz   
      df(3,1) =-anorm*akx*   akz**2*cox*snz
      df(3,2) =-anorm*akx**2*akz   *snx*coz
      df(3,3) =-anorm*akx**3       *cox*snz
!
      df(4,0) = anorm*       akz**4*snx*snz
      df(4,1) =-anorm*akx*   akz**3*cox*coz
      df(4,2) = anorm*akx**2*akz**2*snx*snz
      df(4,3) =-anorm*akx**3*akz   *cox*coz
      df(4,4) = anorm*akx**4       *snx*snz
!
      return
      end

!============================================================
      subroutine inverse(inarray,outarray)
!
!     calculates the inverse of the numvar=1 arrary inarray, and puts
!     it into the numvar=1 array outarray
!
      use p_data
      use t_data
      use basic
      use arrays
      implicit none
      integer l, ndofs, ibegin
!
      real inarray(*),outarray(*)
!
      call numdofs(1, ndofs)
      do l=1,ndofs/6
         ibegin = (l-1)*6+1
         outarray(ibegin) = 1./inarray(ibegin)
         outarray(ibegin+1) = -inarray(ibegin+1)/inarray(ibegin)**2
         outarray(ibegin+2) = -inarray(ibegin+2)/inarray(ibegin)**2
         outarray(ibegin+3) = (2.*inarray(ibegin+1)**2-inarray(ibegin)            &
     &                    *inarray(ibegin+3))/inarray(ibegin)**3
         outarray(ibegin+4) = (2.*inarray(ibegin+1)*inarray(ibegin+2)               &
     &         -inarray(ibegin)*inarray(ibegin+4))/inarray(ibegin)**3
         outarray(ibegin+5) = (2.*inarray(ibegin+2)**2-inarray(ibegin)
     &                    *inarray(ibegin+5))/inarray(ibegin)**3
      enddo
!
      return
      end
!
!==================================
      function sech(x)
      sech = 1./cosh(x)
      return
      end
!==================================
      subroutine definerhs(inarray,outarray,numnodes)
!
!      This routine defines a new vector to be equal to the vector
!      on the RHS
!
!
      use p_data
      use t_data
      use basic
      use arrays
      use sparse
      implicit none
      integer numnodes, numnodes1, l, i, numvards, 
     &     nrankds, itri, j, jone, ione, nbcds, 
     &     ier, jer, nbound, numelms, izone, izonedim, ibegin,
     &     iendplusone, ibeginnv, iendplusonenv, ndofs
      real ssterm, ddterm
      dimension ssterm(3,3),ddterm(3,3)
      real, allocatable::phin(:),temp(:)
      integer, save, allocatable::iboundds(:)
      real inarray(:),outarray(:)
      real aterm(18,18), dterm(18,18), fintl(-6:maxi,-6:maxi)
!
      call numfac(numelms)
      call numnod(numnodes1)
      call numdofs(1, ndofs)
      call createvec(phin, 1)
      call createvec(temp, 1)
!.....put psi into a numvar=1 storage array
      do l=1,numnodes1
         call entdofs(1, l, 0, ibegin, iendplusone)
         call entdofs(numvar, l, 0, ibeginnv, iendplusonenv)
         do i=0,iendplusone-ibegin-1
            phin(ibegin+i) = inarray(ibeginnv+i)
         enddo
      enddo
!
      if(ifirsts4_lu.ne.0) go to 500
      ifirsts4_lu = 1
      nbound = 0
      do i=1,numnodes1
         call zonenod(i,izone,izonedim)
         if(izonedim .ne. 2) nbound = nbound + 1
      enddo
      nbound = 12*numvar*nbound
      allocate(iboundds(nbound))
!
      numvards = 1
      nrankds = ndofs*6
!
!.....compute LU decomposition only once
!
!.....form matrices
      call zeroarray4solve(s4matrix_sm,numvar1_numbering)
      call zeroarray4multiply(d4matrix_sm,numvar1_numbering)
      do itri=1,numelms
         call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
         call calcaterm(itri, aterm, fintl)
         call calcdterm(itri, dterm, fintl)
        do j=1,18
        do i=1,18
            jone = isval1(itri,j)
            ione = isval1(itri,i)
            ssterm(1,1) = dterm(i,j)
            ddterm(1,1) = aterm(i,j)
            call insertval(s4matrix_sm, ssterm(1,1),ione,jone,1)
            call insertval(d4matrix_sm, ddterm(1,1),ione,jone,1)
!
        enddo
        enddo
      enddo
!
!......define indices for boundary arrays
      call boundaryds(iboundds,nbcds,1)
!
      ! bc stuff...
      do i=1,nbcds
         call setdiribc(s4matrix_sm, iboundds(i))
      enddo
      call finalizearray4solve(s4matrix_sm)
      call finalizearray4multiply(d4matrix_sm)
!
!.....perform LU decomposition of sparse matrix "s4matrix"
!     store result in opaque object "s4handle"
      jer = 0
 500  continue
      temp = 0
!
!.....multiply operator matrix d4matrix_sm times phin to get RHS temp
      call matrixvectormult(d4matrix_sm,phin,temp)
      do l=1,nbcds
      temp(iboundds(l)) = 0.
      enddo
!
!.....perform LU backsubstitution to get outarray solution

      call solve(s4matrix_sm,temp,ier)
      outarray = temp
!
      call deletevec(phin)
      call deletevec(temp)
      return
!
      end
!============================================================
      subroutine rhsdef
!
      use p_data
      use t_data
      use basic
      use arrays
      implicit none 
      integer itri, i, j, iii, ii, i1, i2, i3, k, l, numelms,
     &     ibegin, iendplusone, ibegin2, iendplusone2, ibeginn,
     &     iendplusonen
      real term, sum, fac1, fac2, fac3, fact, f, alx, alz
!
      call getboundingboxsize(alx, alz)
      b1vecini = 0.
      call numfac(numelms)

      do 505 itri=1,numelms
!     calculate matrix elements of the integration matrix fint
        do i=-4,maxi
        do j=-4,maxi
        fint(i,j) = 0
        if(i.ge.0 .and.j.ge.0)                                          &
     &     fint(i,j) = f(i,j,atri(itri),btri(itri),ctri(itri))
        enddo
        enddo
!.....x,z dependent part follows
!
!     expand x and z in terms of si and eta for this element
      call local(xi,zi,btri(itri),ttri(itri),itri)
!      write(*,400) itri,xi(1),xi(2),xi(3),zi(1),zi(2),zi(3)
! 400  format('=',i3,6(f6.2))
!
!
!     calculate the RHS function and it's first 2 derivatives
      call fun(xi(1),zi(1),df,alx,alz)  
!
!.....start loops to compute the RHS vector for this triangular element,
!
       do iii=1,3
          call entdofs(1, ist(itri,iii)+1, 0, ibegin, iendplusone)
          call entdofs(2, ist(itri,iii)+1, 0, ibegin2, iendplusone2)
          call entdofs(numvar, ist(itri,iii)+1,0,ibeginn, iendplusonen)
         do ii=1,6
           i = (iii-1)*6 + ii
           isval1(itri,i) = ibegin+ii-1 ! 6*ist(itri,iii)+ii
           isval2(itri,i) = ibegin2+ii-1 ! 12*ist(itri,iii)+ii
           isvaln(itri,i) = ibeginn+ii-1  ! numvar*6*ist(itri,iii)+ii
           i1 = isvaln(itri,i)
           i2 = i1 + 6
           i3 = i2 + 6
!
            term = 0
            do k=1,20
!
             sum = 0
             do j=0,4
             do l=0,j
               fac1 = df(j,l)
               do p=0,l
               fac2 = fac1*xi(2)**(l-p)  *xi(3)**p/(fact(p)*fact(l-p))
               do q=0,j-l
               fac3 = fac2*zi(2)**(j-l-q)*zi(3)**q/(fact(q)*fact(j-l-q))
!
               sum = sum + fac3*fint(mi(k)+j-p-q,ni(k)+p+q)
!
               enddo
               enddo
             enddo
             enddo
!
             term =  term + gtri(k,i,itri)*sum
            enddo
!
!............velocity perturbations
             b1vecini(i1) = b1vecini(i1) - term                         &
     &                    *((pi/alx)**2+(pi/alz)**2)
             if(numvar.ge.2) b1vecini(i2) = b1vecini(i2)
!    &              - term*((pi/alx)**2+(pi/alz)**2)
             if(numvar.ge.3) b1vecini(i3) = 0.
!
!............fluxes
             b2vecini(i1) = 0.
             if(numvar.ge.2) b2vecini(i2) = 0.
             if(numvar.ge.3) b2vecini(i3) = 0.
!
           enddo
           enddo
!.....end of do loops on i
!
 505  continue  ! end of loop over numelms
! since a proc is contributing values to parts of the vector
! it does not own, we call sumshareddofs so that these values
! get summed up for all values shared by multiple procs
! and then update these values
      call sumshareddofs(b1vecini)
c      call safestop(1111)
      return
      end

!============================================================
      subroutine init
      use p_data
      use t_data
      use basic
      use arrays
      implicit none
      integer numnodes,izone,izonedim,i
      real eterm, alx, alz

      call getboundingboxsize(alx, alz)
!
      call second(tfirst)  !Get current system time in seconds
      call numnod(numnodes)

      ifirstd1_lu = 0
      ifirsts1_lu = 0
      ifirstd2_lu = 0
      ifirsts2_lu = 0
      ifirstr1_lu = 0
      ifirstr2_lu = 0
      ifirstq2_lu = 0
      ifirsts3_lu = 0
      ifirsts4_lu = 0
      ifirsts5_lu = 0
      ifirsts6_lu = 0
      ifirsts7_lu = 0
      ifirsts8_lu = 0
      ifirstd8_lu = 0
      ifirstq8_lu = 0
      ifirstr8_lu = 0

      gbound = 0.
!
!.....define all input parameters (choose by test problem)
      call input
      ntimemin = min(2,ntimemax-1)
      iboundmax = 0
      do i=1,numnodes
         call zonenod(i,izone,izonedim)
         if(izonedim .ne. 2) iboundmax = iboundmax + 1
      enddo
      iboundmax = 12*numvar*iboundmax
!
      tt = 2.*thimp
      gam = 5./3.
!
!      allocate big arrays
      call space
!
!
!.................................................................
!     PART 1:
!.....initialize needed variables and define geometry and triangles
      if(ntimemax .gt. ntimep) then
            write(*,6665) ntimemax,ntimep
 6665 format(" ntimemax .gt. ntimep", 2i5)
            call safestop(18) 
      endif
!
!
!.....zero out big arrays
      b1vecini = 0.
      b2vecini = 0.
      eterm = 0
      rinv = 0.
!
!.....define properties of triangles
      call tridef
!.....define 1/R for toroidal problem
      if(itor.eq.1) call rinvdef(rinv)
!
!     open all needed output files
      call openf

!
!.....end of section defining geometrical parameters for triangles
      return
      end
!============================================================
