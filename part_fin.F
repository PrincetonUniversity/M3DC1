      subroutine mask(x,z,factor)
      use basic
      implicit none
c$$$      real factor, r1, r2, alpha, rmid, rnorm, r, x, z, alx, alz,
c$$$     &     xmin, zmin
      real, intent(in) :: x, z
      real, intent(out) :: factor

      real :: alphax, alphaz, alx, alz, xmin, zmin
      real :: x_left, x_right, z_bottom, z_top
      factor = 1.
! acbauer -- changed to use xmin and zmin since it apears that 
!            this routine may want to center the domain around (0,0)

      call getboundingboxsize(alx, alz)
      call getmincoord(xmin, zmin)
c$$$!
c$$$!     radius where the factor becomes 0
c$$$      r1 = 2.0
c$$$!
c$$$!     radius where the factor becomes 1
c$$$      r2 = 1.5
c$$$      alpha = 6./(r1 - r2)
c$$$      rmid = (r1+r2)/2.
c$$$!
c$$$      if(iper.eq.0) then
c$$$         rnorm = sqrt(alx*alz)/4.
c$$$         r = sqrt( (x-xmin-alx*.5)**2 + (z-zmin-alz*.5)**2)/rnorm
c$$$      else
c$$$         rnorm = alz/4.
c$$$         r = abs(z-zmin-alz/2.)/rnorm
c$$$      endif
c$$$!
c$$$      factor = (1 + tanh(-alpha*(r-rmid)))*.5
c$$$!
c$$$!     write(66,1001) x,z,r,factor
c$$$ 1001 format(1p4e12.4)
c$$$      return

      x_left = alx / 20.
      x_right = alx - alx / 20.
      z_bottom = alz / 20.
      z_top = alz - alz / 20.
      
      alphax = 40. / alx
      alphaz = 40. / alz

      factor = 1.

      if(iper.eq.0) then
         factor = factor * 0.5*
     &        (tanh(alphax*(x - x_left)) - tanh(alphax*(x - x_right)))
      endif

      if(jper.eq.0) then
         factor = factor * 0.5*
     &        (tanh(alphaz*(z - z_bottom)) - tanh(alphaz*(z-z_top)))
      endif

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
      use sparse
      use nintegrate_mod

      implicit none
      integer numnodes, numelms, i, izone, izonedim 

      if(isplitstep.eq.1) then
         vecsize = numvar
         vecsize1 = 1
      else
         vecsize  = numvar*2 + idens + ipres
         vecsize1 = vecsize
      endif

      ! create numberings
      call createdofnumbering(numvar1_numbering, iper, jper, 
     &     6, 0, 0, 0, maxdofs1)
      call createdofnumbering(numvar2_numbering, iper, jper, 
     &     12, 0, 0, 0, maxdofs2)
      if(vecsize.gt.2) then
         call createdofnumbering(vecsize, iper, jper,
     &        vecsize*6, 0, 0, 0, maxdofsn)
      endif

      call numnod(numnodes)
      call numfac(numelms)
      ntri   = 2*numelms !should this be ntri = numelms???
!
!.....arrays defined at all vertices
c     created with ordering specified by numvar

      if(isplitstep.eq.1) then
         call createvec(vel, vecsize)
         call createvec(vel0, vecsize)
         call createvec(vels, vecsize)
         call createvec(velold, vecsize)
         call createvec(veln, vecsize)
         call createvec(veloldn, vecsize)
         if(iresolve.eq.1) call createvec(phioldn, vecsize)
         call createvec(r4, vecsize)
         call createvec(phip, vecsize)
         call createvec(b2vector, vecsize)
         call createvec(b3vector, vecsize)

         if(ipres.eq.1) then 
            call createvec(qp4, 1)
            call createvec(pres, 1)
            call createvec(press, 1)
            call createvec(pres0, 1)
            call createvec(presold, 1)
            if(iresolve.eq.1) call createvec(presoldn, 1)
         endif
         
         call createvec(den, 1)
         call createvec(den0, 1)
         call createvec(dens, 1)
         call createvec(denold, 1)
         call createvec(qn4, 1)
         if(isources.eq.1) then
            call createvec(sb1, 1)
            call createvec(sb2, 1)
            call createvec(sp1, 1)
         endif
      endif

      call createvec(phi,      vecsize)
      call createvec(phi0,     vecsize)
      call createvec(phis,     vecsize)
      call createvec(phiold,   vecsize)

      call createvec(b1vector, vecsize)
      call createvec(vtemp,    vecsize)
      call createvec(q4,       vecsize)

      call createvec(jphi, 1)
      call createvec(vor, 1)
      call createvec(com, 1)
      call createvec(resistivity, 1)
      call createvec(kappa, 1)
      call createvec(visc, 1)
      call createvec(visc_c, 1)
      call createvec(sigma, 1)
      call createvec(tempvar,1)
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
      real hyp, thimpsm
      integer numvarsm, itri, i, j,izone,izonedim,
     1     l, jone, ier, jtwo, itwo
!
      if(myrank.eq.0 .and. iprint.eq.1) print *, " smoother3"

      numvarsm = 2
      call numfac(numelms)
      call numdofs(numvarsm, ndofs_sm)
      call createvec(temp, numvarsm)
      temp = 0
      thimpsm = 1.
! 
      if(ifirsts7_lu.ne.0) go to 500
      ifirsts7_lu = 1
!
!
!.....compute LU decomposition only once
!
!.....form matrix
      call zerosuperluarray(s7matrix_sm,numvar2_numbering)
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
!.....Insert boundary conditions
      call boundary_com(s7matrix_sm, temp)
      call finalizearray(s7matrix_sm)
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

      call boundary_com(0, temp)
      call sumshareddofs(temp)
!
!.....perform LU backsubstitution to get outarray solution
      call solve(s7matrix_sm,temp,ier)
      if(ier.ne.0) then
         print *, "Error in smoother3 solve"
         call safestop(102)
      endif
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
      use diagnostics

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

c$$$      subroutine inverse(inarray,outarray)
c$$$!
c$$$!     calculates the inverse of the numvar=1 arrary inarray, and puts
c$$$!     it into the numvar=1 array outarray
c$$$!
c$$$      use p_data
c$$$      use t_data
c$$$      use basic
c$$$      use arrays
c$$$      implicit none
c$$$      integer :: l, ndofs, iendplusone, ibegin
c$$$!
c$$$      real inarray(*),outarray(*)
c$$$!
c$$$      call numdofs(1, ndofs)
c$$$      do l=1,ndofs/6
c$$$         ibegin = (l-1)*6+1
c$$$         outarray(ibegin) = 1./inarray(ibegin)
c$$$         outarray(ibegin+1) = -inarray(ibegin+1)/inarray(ibegin)**2
c$$$         outarray(ibegin+2) = -inarray(ibegin+2)/inarray(ibegin)**2
c$$$         outarray(ibegin+3) = (2.*inarray(ibegin+1)**2-inarray(ibegin)         
c$$$     &                    *inarray(ibegin+3))/inarray(ibegin)**3
c$$$         outarray(ibegin+4) = (2.*inarray(ibegin+1)*inarray(ibegin+2)               &
c$$$     &         -inarray(ibegin)*inarray(ibegin+4))/inarray(ibegin)**3
c$$$         outarray(ibegin+5) = (2.*inarray(ibegin+2)**2-inarray(ibegin)
c$$$     &                    *inarray(ibegin+5))/inarray(ibegin)**3
c$$$      enddo
c$$$!
c$$$      return
c$$$      end
c$$$!
!==================================
      function sech(x)
      sech = 1./cosh(x)
      return
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

      real :: f
!
      call numfac(numelms)

      do itri=1,numelms
!        calculate matrix elements of the integration matrix fint
         do i=-6,maxi
            do j=-6,maxi
               fint(i,j) = 0
               if(i.ge.0 .and.j.ge.0)                                          &
     &              fint(i,j) = f(i,j,atri(itri),btri(itri),ctri(itri))
            enddo
         enddo

         do iii=1,3
            call entdofs(1, ist(itri,iii)+1, 0, ibegin,  iendplusone )
            call entdofs(2, ist(itri,iii)+1, 0, ibegin2, iendplusone2)
            call entdofs(vecsize,  ist(itri,iii)+1, 0,                  &
     &           ibeginn, iendplusonen)
            do ii=1,6
               i = (iii-1)*6 + ii
               isval1(itri,i) = ibegin +ii-1 ! 6*ist(itri,iii)+ii
               isval2(itri,i) = ibegin2+ii-1 ! 12*ist(itri,iii)+ii
               isvaln(itri,i) = ibeginn+ii-1 ! numvar*6*ist(itri,iii)+ii
            enddo
         enddo
!.....end of do loops on i
!
      enddo

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
      call numnod(numnodes)

      ifirsts4_lu = 0
      ifirsts5_lu = 0
      ifirsts7_lu = 0

      fbound = 0.
      gbound = 0.
!
!.....define all input parameters (choose by test problem)
      print *, "input..."
      call input
      ntimemin = min(2,ntimemax-1)
!
!      allocate big arrays
      print *, "Calling space..."
      call space
      print *, "Done space."
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
c$$$      b1vecini = 0.
c$$$      b2vecini = 0.
      eterm = 0
      rinv = 0.
!
!.....define properties of triangles
      call tridef
!.....define 1/R for toroidal problem
      if(itor.eq.1) call rinvdef(rinv)
!
!     open all needed output files
      if(maxrank.eq.1 .and. myrank.eq.0) call openf

!
!.....end of section defining geometrical parameters for triangles
      return
      end
