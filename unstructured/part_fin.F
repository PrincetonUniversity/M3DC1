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
         if(vecsize_phi.gt.2) then
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
         if(gyro.eq.1) call createvec(gyro_tau,1)
         if(i3d.eq.1 .or. ifout.eq.1) call createvec(bf,1)

         ! Arrays for implicit time advance
         call createvec(phi,      vecsize_phi)
         call createvec(phiold,   vecsize_phi)
         call createvec(q4,       vecsize_phi)

         call createvec(b1_phi,   vecsize_phi)
         call createvec(b2_phi,   vecsize_phi)

         if(isplitstep.eq.1) then
            call createvec(vel,      vecsize_vel)
            call createvec(velold,   vecsize_vel)
            call createvec(veln,     vecsize_vel)
            call createvec(veloldn,  vecsize_vel)
            call createvec(r4,       vecsize_vel)
            call createvec(phip,     vecsize_phi)

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
!============================================================
      subroutine smoother1(inarray,outarray,numvard,iplace)
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
      use nintegrate_mod

      implicit none

#include "finclude/petsc.h"

      integer, intent(in) :: numvard, iplace

      vectype, intent(out) :: outarray(*)
      vectype, intent(in)  ::  inarray(*)

      integer :: itri, 
     &     i, j, numvarsm, jtwo, itwo,
     &     ier, numelms

      real :: hyp, thimpsm
      logical, save :: first_time = .true.

      real, dimension(2,2) :: ssterm
      vectype, allocatable :: temp(:)
      vectype, dimension(20) :: avec

      PetscTruth :: flg
      integer :: ipetsc, isuperlu

!
      if(myrank.eq.0 .and. iprint.eq.1) print *, " smoother1"

      numvarsm = 2
      call numfac(numelms)
      call createvec(temp, numvarsm)
      temp = 0.
      thimpsm = 1.
! 
      if(.not.first_time) go to 500
      first_time = .false.
!
!
!.....compute LU decomposition only once
!
!.....form matrix default linear solver superlu cj-april-09-2008
      ipetsc=0
      isuperlu=1
      call PetscOptionsGetInt
     &(PETSC_NULL_CHARACTER,'-ipetsc',ipetsc,flg,ier)
      if(ipetsc) isuperlu=0
      if(isuperlu) call zerosuperlumatrix(s5matrix_sm,icomplex,numvarsm)
      if(ipetsc) call zeropetscmatrix(s5matrix_sm,icomplex,numvarsm)
         if(isuperlu) 
     &   print *, " smoother1 zerosuperlumatrix", s5matrix_sm
         if(ipetsc) 
     &   print *, "   smoother1 zeropetscmatrix", s5matrix_sm 
      do itri=1,numelms

         call define_fields_79(itri,0,1)

         hyp = hypc*amu

         do i=1,18
            do j=1,18
              jtwo = isval2(itri,j)
              itwo = isval2(itri,i)

              ssterm(1,1) = 
     &             int2(g79(:,OP_1,i),g79(:,OP_1,j),weight_79,79)
              ssterm(1,2) = 
     &              int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j),weight_79,79)
     &             +int2(g79(:,OP_DR,i),g79(:,OP_DR,j),weight_79,79)
              if(itor.eq.1) then
                 ssterm(1,2) = ssterm(1,2) + 2.*
     &                int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),
     &                     weight_79,79)
              endif
              ssterm(2,1) = dt*hyp*thimpsm* 
     &              int2(g79(:,OP_GS,i),g79(:,OP_GS,j),weight_79,79)
              ssterm(2,2) = -ssterm(1,2)
              if(inoslip_pol.eq.0) then
                 ssterm(2,2) = ssterm(2,2) - regular*ssterm(1,1)
              end if

              call insval(s5matrix_sm,ssterm(1,1),0,itwo  ,jtwo  ,1)
              call insval(s5matrix_sm,ssterm(1,2),0,itwo  ,jtwo+6,1)
              call insval(s5matrix_sm,ssterm(2,1),0,itwo+6,jtwo  ,1)
              call insval(s5matrix_sm,ssterm(2,2),0,itwo+6,jtwo+6,1)

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

         call define_fields_79(itri,0,1)

         call calcavector(itri, inarray, 1, 1, avec)
         call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, tm79)

         hyp = hypc*amu

         do i=1,18
            itwo = isval2(itri,i)

            temp(itwo+6) = temp(itwo+6) + 
     &           (int2(g79(:,OP_1,i),tm79(:,OP_1),weight_79,79)
     &           -dt*hyp*(1.-thimpsm)*
     &            int2(g79(:,OP_GS,i),tm79(:,OP_GS),weight_79,79))
         enddo
      enddo

      ! Insert boundary conditions
      call boundary_vor(0, temp)
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
      call copyvec(temp,2,numvarsm,outarray,iplace,numvard)
      
      call deletevec(temp)
      return
!
      end
!============================================================
      subroutine smoother3(inarray,outarray,numvard,iplace)
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

      use nintegrate_mod

      implicit none

#include "finclude/petsc.h"

      integer, intent(in) :: numvard, iplace
      vectype, intent(out) :: outarray(*)
      vectype, intent(in) :: inarray(*)


      integer ::  numelms
      integer :: numvarsm, itri, i, j, ier, jtwo, itwo
      real :: hyp, thimpsm

      logical, save :: first_time = .true.

      real, dimension(2,2) :: ssterm
      vectype, allocatable:: temp(:)
      vectype, dimension(20) :: avec


      PetscTruth :: flg
      integer :: ipetsc, isuperlu
!
      if(myrank.eq.0 .and. iprint.eq.1) print *, " smoother3"

      numvarsm = 2
      call numfac(numelms)
      call createvec(temp, numvarsm)
      temp = 0.
      thimpsm = 1.
! 
      if(.not.first_time) go to 500
      first_time = .false.

      if(myrank.eq.0 .and. iprint.eq.1)
     &     print *, "  defining smoother3 matrix"
!
!
!.....compute LU decomposition only once
!
!.....form matrix default linear solver superlu cj-april-09-2008
      ipetsc=0
      isuperlu=1
      call PetscOptionsGetInt
     &(PETSC_NULL_CHARACTER,'-ipetsc',ipetsc,flg,ier)
      if(ipetsc) isuperlu=0
      if(isuperlu) call zerosuperlumatrix(s7matrix_sm,icomplex,numvarsm)
      if(ipetsc) call zeropetscmatrix(s7matrix_sm,icomplex,numvarsm)
         if(isuperlu) 
     &   print *, " smoother3 zerosuperlumatrix", s7matrix_sm
         if(ipetsc) 
     &   print *, "   smoother3 zeropetscmatrix", s7matrix_sm 
      do itri=1,numelms

         call define_fields_79(itri,0,1)

         hyp = hypc*amuc

         do i=1,18
            do j=1,18
              jtwo = isval2(itri,j)
              itwo = isval2(itri,i)

              ssterm(1,1) = 
     &             int2(g79(:,OP_1,i),g79(:,OP_1,j),weight_79,79)
              ssterm(1,2) = 
     &              int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j),weight_79,79)
     &             +int2(g79(:,OP_DR,i),g79(:,OP_DR,j),weight_79,79)
              ssterm(2,1) = dt*hyp*thimpsm* 
     &              int2(g79(:,OP_LP,i),g79(:,OP_LP,j),weight_79,79)
              ssterm(2,2) = -ssterm(1,2)
              if(inoslip_pol.eq.0) then
                 ssterm(2,2) = ssterm(2,2) - regular*ssterm(1,1)
              end if

              call insval(s7matrix_sm,ssterm(1,1),0,itwo  ,jtwo  ,1)
              call insval(s7matrix_sm,ssterm(1,2),0,itwo  ,jtwo+6,1)
              call insval(s7matrix_sm,ssterm(2,1),0,itwo+6,jtwo  ,1)
              call insval(s7matrix_sm,ssterm(2,2),0,itwo+6,jtwo+6,1)
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

         call define_fields_79(itri,0,1)

         call calcavector(itri, inarray, 1, 1, avec)
         call eval_ops(avec, si_79, eta_79, ttri(itri), ri_79,79, tm79)

         hyp = hypc*amuc

         do i=1,18
            itwo = isval2(itri,i)

            temp(itwo+6) = temp(itwo+6) + 
     &           (int2(g79(:,OP_1,i),tm79(:,OP_1),weight_79,79)
     &           -dt*hyp*(1.-thimpsm)*
     &            int2(g79(:,OP_LP,i),tm79(:,OP_LP),weight_79,79))
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
      call copyvec(temp,2,numvarsm,outarray,iplace,numvard)

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

      if(myrank.eq.1 .and. iprint.ge.1) print *, " Entering init..."

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

      if(myrank.eq.1 .and. iprint.ge.1) print *, " Exiting init."

      return
      end
