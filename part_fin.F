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
!     ifirstcall = 1 if this is the first call and 0 otherwise
!                  this is usefull for avoiding memory leaks
      subroutine space(ifirstcall)
!
!     allocates space for big arrays
      use p_data
      use t_data
      use basic
      use arrays
      use sparse
      use nintegrate_mod

      implicit none

      integer, intent(in) :: ifirstcall

      integer :: numnodes, numelms, i

      if(isplitstep.eq.1) then
         vecsize = numvar
         vecsize1 = 1
      else
         vecsize  = numvar*2 + idens + ipres
         vecsize1 = vecsize
      endif

!
!.....create numberings
      if(ifirstcall .eq. 1) then
         call createdofnumbering(numvar1_numbering, iper, jper, 
     &        6, 0, 0, 0, maxdofs1)
         call createdofnumbering(numvar2_numbering, iper, jper, 
     &        12, 0, 0, 0, maxdofs2)
         call createdofnumbering(num_fields, iper, jper, 
     &        num_fields*6, 0, 0, 0, maxdofsn)
         if(vecsize.gt.2) then
            call createdofnumbering(vecsize, iper, jper,
     &           vecsize*6, 0, 0, 0, maxdofsn)
         endif
      endif

      call numnod(numnodes)
      call numfac(numelms)
      ntri   = numelms !should this be ntri = numelms??? -- YES, acbauer
!
!.....arrays defined at all vertices
!     createvec will delete the arrays if they have already been allocated

      if(ifirstcall.eq.1) then
      ! Physical Variables
      call createvec(field , num_fields)
      call createvec(field0, num_fields)
      call createvec(fieldi, num_fields)

      ! Auxiliary Variables
      call createvec(jphi, 1)
      call createvec(vor, 1)
      call createvec(com, 1)
      call createvec(resistivity, 1)
      call createvec(kappa, 1)
      call createvec(visc, 1)
      call createvec(visc_c, 1)
      call createvec(sigma, 1)
      call createvec(tempvar,1)
      if(i3d.eq.1) call createvec(bf,1)      

      ! Arrays for implicit time advance
      call createvec(phi,  vecsize)
      call createvec(phiold,   vecsize)

      call createvec(b1vector, vecsize)
      call createvec(vtemp,    vecsize)
      call createvec(q4,       vecsize)

      if(isplitstep.eq.1) then
         call createvec(vel, vecsize)
         call createvec(velold, vecsize)
         call createvec(veln, vecsize)
         call createvec(veloldn, vecsize)
         call createvec(r4, vecsize)
         call createvec(phip, vecsize)
         call createvec(b2vector, vecsize)
         call createvec(b3vector, vecsize)

         if(ipres.eq.1) then
            call createvec(pres, 1)
            call createvec(presold, 1)
            call createvec(qp4, 1)
         endif

         call createvec(den, 1)        
         call createvec(denold, 1)
         call createvec(qn4, 1)

         if(isources.eq.1) then
            call createvec(sb1, 1)
            call createvec(sb2, 1)
            call createvec(sp1, 1)
         endif
      endif

      endif
!
!.....arrays associated with the triangles
      if(ifirstcall.eq.0) then
         deallocate(ist)
         deallocate(isvaln)
         deallocate(isval1)
         deallocate(isval2)
         deallocate(atri)
         deallocate(btri)
         deallocate(ctri)
         deallocate(ttri)
         deallocate(gtri)
      endif
      allocate(ist(ntri,3),isvaln(ntri,18),isval1(ntri,18),             &
     &     isval2(ntri,18))
      allocate(atri(numelms),btri(numelms),ctri(numelms),ttri(numelms),
     &     gtri(20,18,numelms)) 
!
!.....metric terms involved in toroidal geometry
      if(ifirstcall.eq.0) deallocate(rinv)
      allocate(rinv(numnodes*6))
!
!.....assign pointers to proper vectors
      call assign_variables

      return
      end
!============================================================
      subroutine smoother1(inarray,outarray,numnodes,numvard,iplace)
!
!     This routine takes as input the Laplacian of the scalar field in
!     the array inarray, smooths it, and returns the smoothed scalar array
!     in location iplace in the stand-alone array outarray(ordering =numvard).  
!     The LU decomposition takes place the first-time only.
!
!
      use p_data
      use t_data
      use basic
      use arrays
      use sparse

      implicit none

      integer, intent(in) :: numnodes, numvard, iplace
      vectype, intent(out) :: outarray(numnodes*6*numvard)
      vectype, intent(in)  ::  inarray(numnodes*6)

      integer :: itri, ndof_2,
     &     i, j, numvarsm, jtwo, itwo, l, 
     &     ier, jone, nbound, numelms, izone,izonedim,
     &     ibegind, iendplusoned, ibeginsm, iendplusonesm
      real :: hyp, thimpsm
      real, dimension(18,18) :: aterm, bterm, dterm
      real, dimension(-6:maxi,-6:maxi) :: fintl(-6:maxi,-6:maxi)

      logical, save :: first_time = .true.

      real, dimension(2,2) :: ssterm
      vectype, allocatable :: temp(:)
!
      if(myrank.eq.0 .and. iprint.eq.1) print *, " smoother1"

      numvarsm = 2
      call numfac(numelms)
      call numdofs(numvarsm, ndof_2)
      call createvec(temp, numvarsm)
      temp = 0
      thimpsm = 1.
! 
      if(.not.first_time) go to 500
      first_time = .false.
!
!
!.....compute LU decomposition only once
!
!.....form matrix
      call zerosuperlumatrix(s5matrix_sm,icomplex,numvar2_numbering)
      do itri=1,numelms
         call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
         call calcaterm(itri, aterm, fintl)
         call calcbterm(itri, bterm, fintl)
         call calcdterm(itri, dterm, fintl)
!         call getdeex(itri,deex)
         hyp = hyperc*deex**2*amu
         do i=1,18
            do j=1,18
               jtwo = isval2(itri,j)
               itwo = isval2(itri,i)
               ssterm(1,1) = dterm(i,j)
               ssterm(1,2) =-aterm(i,j)
               ssterm(2,1) =  +dt*hyp*thimpsm*bterm(i,j)
               ssterm(2,2) = aterm(i,j)
              call insertval2(s5matrix_sm,ssterm(1,1),0,itwo,  jtwo,  1)
              call insertval2(s5matrix_sm,ssterm(1,2),0,itwo,  jtwo+6,1)
              call insertval2(s5matrix_sm,ssterm(2,1),0,itwo+6,jtwo,  1)
              call insertval2(s5matrix_sm,ssterm(2,2),0,itwo+6,jtwo+6,1)             
            enddo
         enddo
      enddo
!
!.....Insert boundary condtions
      call boundary_vor(s5matrix_sm, temp)
      call finalizematrix(s5matrix_sm)
 500  continue
!
!.....define RHS vector
      do itri=1,numelms
         call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
         call calcbterm(itri, bterm, fintl)
         call calcdterm(itri, dterm, fintl)
        
!         call getdeex(itri,deex)
         hyp = hyperc*deex**2*amu
         do i=1,18
            itwo = isval2(itri,i)
            do j=1,18
               jone = isval1(itri,j)
!
               temp(itwo+6) = temp(itwo+6)+(dterm(i,j)                             &
     &              - dt*hyp*(1.-thimpsm)*bterm(i,j))*inarray(jone)
!
            enddo
         enddo
      enddo

      ! Insert boundary conditions
      call boundary_vor(0, temp)

! since a proc is contributing values to parts of the vector
! it does not own, we call sumsharedppplvecvals so that these values
! get summed up for all values shared by multiple procs
! and then update these values
      call sumsharedppplvecvals(temp)
!
!.....perform LU backsubstitution to get outarray solution

      call solve(s5matrix_sm,temp,ier)
      if(ier.ne.0) then
         print *, "Error in smoother1 solve"
         call safestop(101)
      end if
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

      integer, intent(in) :: numnodes, numvard, iplace
      vectype, intent(out) :: outarray(numnodes*6*numvard)
      vectype, intent(in) :: inarray(numnodes*6)

      integer :: ndofs_sm, numelms,
     &     ibegind, iendplusoned, ibeginsm, iendplusonesm
      integer :: numvarsm, itri, i, j,izone,izonedim,
     &     l, jone, ier, jtwo, itwo
      real :: hyp, thimpsm
      real, dimension(18,18) ::  aterm, bterm, dterm
      real :: fintl(-6:maxi,-6:maxi)

      logical, save :: first_time = .true.

      real, dimension(2,2) :: ssterm
      vectype, allocatable:: temp(:)
!
      if(myrank.eq.0 .and. iprint.eq.1) print *, " smoother3"

      numvarsm = 2
      call numfac(numelms)
      call numdofs(numvarsm, ndofs_sm)
      call createvec(temp, numvarsm)
      temp = 0
      thimpsm = 1.
! 
      if(.not.first_time) go to 500
      first_time = .false.
!
!
!.....compute LU decomposition only once
!
!.....form matrix
      call zerosuperlumatrix(s7matrix_sm,icomplex,numvar2_numbering)
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

              call insertval2(s7matrix_sm,ssterm(1,1),0,itwo  ,jtwo  ,1)
              call insertval2(s7matrix_sm,ssterm(1,2),0,itwo  ,jtwo+6,1)
              call insertval2(s7matrix_sm,ssterm(2,1),0,itwo+6,jtwo  ,1)
              call insertval2(s7matrix_sm,ssterm(2,2),0,itwo+6,jtwo+6,1)
           enddo
        enddo
      enddo
!
!.....Insert boundary conditions
      call boundary_com(s7matrix_sm, temp)
      call finalizematrix(s7matrix_sm)
!
 500  continue
!
!.....define RHS vector
      do itri=1,numelms
         call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
         call calcbterm(itri, bterm, fintl)
         call calcdterm(itri, dterm, fintl)
!         call getdeex(itri,deex)
         hyp = hyperc*deex**2*amuc
         do i=1,18
            itwo = isval2(itri,i)
            do j=1,18
               jone = isval1(itri,j)

               temp(itwo+6) = temp(itwo+6)+(dterm(i,j)               
     &              - dt*hyp*(1.-thimpsm)*bterm(i,j))*inarray(jone)
            enddo
         enddo
      enddo

      call boundary_com(0, temp)
      call sumsharedppplvecvals(temp)
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
      use basic
      use arrays
      use diagnostics

      implicit none

      real, parameter :: scalefac = 1e-10

!.....scale solution back if a linear calculation

      if(ekin .lt. max_ke) return

      field = scalefac*field

      phiold = 1e-10*phiold
      if(isplitstep.eq.1) then
         velold = scalefac*velold
         veloldn = scalefac*veloldn
         denold = scalefac*denold
         if(ipres.eq.1) presold = scalefac*presold
      endif
      
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
      use sparse
      implicit none

      fbound = 0.
      gbound = 0.
!
!.....define all input parameters (choose by test problem)
      print *, "input..."
      call input
!
!      allocate big arrays
      print *, "Calling space..."
      call space(1)
      print *, "Done space."
!
!
!.....zero out big arrays
      rinv = 0.
!
!.....define properties of triangles
      call tridef
!.....define 1/R for toroidal problem
      if(itor.eq.1) call rinvdef(rinv)

      return
      end
