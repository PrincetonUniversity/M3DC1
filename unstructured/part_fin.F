      subroutine mask(x,z,factor)
      use basic
      implicit none

      real, intent(in) :: x, z
      real, intent(out) :: factor

      real :: alphax, alphaz, alx, alz, xmin, zmin
      real :: x_left, x_right, z_bottom, z_top
      factor = 1.

      call getboundingboxsize(alx, alz)
      call getmincoord(xmin, zmin)

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
     &        (tanh(alphaz*(z - z_bottom)) - tanh(alphaz*(z - z_top)))
      endif

      return
      end
!============================================================
!     ifirstcall = 1 if this is the first call and 0 otherwise
!                  this is useful for avoiding memory leaks
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

      integer :: numnodes, numelms

      if(myrank.eq.1 .and. iprint.ge.1) print *, " Entering space..."

      if(isplitstep.eq.1) then
         vecsize_phi = numvar
         vecsize_vel = numvar
         vecsize_n = 1
         vecsize_p = 1
      else
         vecsize_phi  = numvar*2 + idens + ipres
         vecsize_vel  = vecsize_phi
         vecsize_n    = vecsize_phi
         vecsize_p    = vecsize_phi
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
         if(vecsize_phi.gt.2 .and. vecsize_phi.ne.num_fields) then
            call createdofnumbering(vecsize_phi, iper, jper,
     &           vecsize_phi*6, 0, 0, 0, maxdofsn)
         endif
         if(vecsize_vel.gt.2 .and. vecsize_vel.ne.vecsize_phi) then
            call createdofnumbering(vecsize_vel, iper, jper,
     &           vecsize_vel*6, 0, 0, 0, maxdofsn)
         endif
      endif

      call numnod(numnodes)
      call numfac(numelms)
      ntri   = numelms !should this be ntri = numelms??? -- YES, acbauer
!
!.....arrays defined at all vertices
!     createvec will delete the arrays if they have already been allocated

      if(ifirstcall.eq.1) then
         if(myrank.eq.0 .and. iprint.eq.1) print *, 'Allocating...'

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
         if(icomplex.eq.1 .or. ifout.eq.1) call createvec(bf,1)
         if(ibootstrap.gt.0) call createvec(visc_e,1)

         ! Arrays for implicit time advance
         call createvec(phi,      vecsize_phi)
         call createvec(phiold,   vecsize_phi)
         call createvec(phip,     vecsize_phi)
         call createvec(q4,       vecsize_phi)

         call createvec(b1_phi,   vecsize_phi)
         call createvec(b2_phi,   vecsize_phi)

         if(isplitstep.eq.1) then
            call createvec(vel,      vecsize_vel)
            call createvec(velold,   vecsize_vel)
            call createvec(veln,     vecsize_vel)
            call createvec(veloldn,  vecsize_vel)
            call createvec(r4,       vecsize_vel)

            if(ipres.eq.1) then
               call createvec(pres,    vecsize_p)
               call createvec(presold, vecsize_p)
               call createvec(qp4,     vecsize_p)
            endif
            
            call createvec(den,    vecsize_n)        
            call createvec(denold, vecsize_n)
            call createvec(qn4,    vecsize_n)

            call createvec(b1_vel,   vecsize_vel)
            call createvec(b2_vel,   vecsize_vel)

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
         if(myrank.eq.0 .and. iprint.eq.1) print *, ' deallocating...'
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

      if(myrank.eq.0 .and. iprint.eq.1) print *, ' Allocting tri...'
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

      if(myrank.eq.1 .and. iprint.ge.1) print *, " Exiting space."

      return
      end
!======================================================================
      subroutine scaleback
!
      use basic
      use arrays
      use diagnostics

      implicit none

      integer :: i, ndofs, numnodes, ibegin, iendplusone
      integer, allocatable ::  itemp(:)

      real, parameter :: scalefac = 1.e-10

!.....scale solution back if a linear calculation

      if(ekin.lt.max_ke .or. max_ke.eq.0) return
      if(myrank.eq.0) write(*,*) " =>solution scaled back at time", time

      call numnod(numnodes)
      call numdofs(num_fields, ndofs)
      allocate(itemp(ndofs))
      itemp = 1.
      do i=1,numnodes
         call entdofs(num_fields, i, 0, ibegin, iendplusone)
         if(itemp(ibegin) .eq. 0) cycle
         itemp(ibegin) = 0

         field(ibegin:iendplusone-1) = scalefac* 
     &        field(ibegin:iendplusone-1)

         call entdofs(vecsize_phi, i, 0, ibegin, iendplusone)
         phiold(ibegin:iendplusone-1) = scalefac*
     &        phiold(ibegin:iendplusone-1)

         if(isplitstep.eq.1) then
            call entdofs(vecsize_vel, i, 0, ibegin, iendplusone)
            velold(ibegin:iendplusone-1) = scalefac* 
     &           velold(ibegin:iendplusone-1)
            veloldn(ibegin:iendplusone-1) = scalefac*
     &           veloldn(ibegin:iendplusone-1)
            if(idens.eq.1) then 
               call entdofs(vecsize_n, i, 0, ibegin, iendplusone)
               denold(ibegin:iendplusone-1) = scalefac*
     &              denold(ibegin:iendplusone-1)
            endif
            if(ipres.eq.1) then
               call entdofs(vecsize_p, i, 0, ibegin, iendplusone)
               presold(ibegin:iendplusone-1) = scalefac*
     &              presold(ibegin:iendplusone-1)
            endif
         endif
         if(icomplex.eq.1) then
            call entdofs(1, i, 0, ibegin, iendplusone)
            bf(ibegin:iendplusone-1) = scalefac*
     &           bf(ibegin:iendplusone-1)
         endif
      enddo

      deallocate(itemp)

      end
!==================================
      real function sech(x)
      implicit none
      real, intent(in) :: x
      
      sech = 1./cosh(x)
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

      if(myrank.eq.0 .and. iprint.ge.1) print *, " Entering init..."
      
      fbound = 0.
      gbound = 0.
!
!.....define all input parameters (choose by test problem)
      call input
!
!      allocate big arrays
      call space(1)
!
!
!.....zero out big arrays
      rinv = 0.
!
!.....define properties of triangles
      call tridef
!.....define 1/R for toroidal problem
      if(itor.eq.1) call rinvdef(rinv)

      if(myrank.eq.0 .and. iprint.ge.1) print *, " Exiting init."

      return
      end
