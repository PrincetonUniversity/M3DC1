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
      real, intent(out) :: outarray(numnodes*6*numvard)
      real, intent(in)  ::  inarray(numnodes*6)


      integer :: itri, ndof_2,
     &     i, j, numvarsm, jtwo, itwo, l, 
     &     ier, jone, nbound, numelms, izone,izonedim,
     &     ibegind, iendplusoned, ibeginsm, iendplusonesm
      real :: hyp, thimpsm
      real, dimension(18,18) :: aterm, bterm, dterm
      real, dimension(-6:maxi,-6:maxi) :: fintl(-6:maxi,-6:maxi)

      real, dimension(2,2) :: ssterm
      real, allocatable :: temp(:)
!
      if(myrank.eq.0 .and. iprint.eq.1) print *, " smoother1"

      numvarsm = 2
      call numfac(numelms)
      call numdofs(numvarsm, ndof_2)
      call createvec(temp, numvarsm)
      temp = 0
      thimpsm = 1.
! 
      if(ifirsts5_lu.ne.0) go to 500
      ifirsts5_lu = 1
!
!
!.....compute LU decomposition only once
!
!.....form matrix
      call zerosuperluarray(s5matrix_sm,numvar2_numbering)
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
               call insertval(s5matrix_sm,ssterm(1,1),itwo,  jtwo,  1)
               call insertval(s5matrix_sm,ssterm(1,2),itwo,  jtwo+6,1)
               call insertval(s5matrix_sm,ssterm(2,1),itwo+6,jtwo,  1)
               call insertval(s5matrix_sm,ssterm(2,2),itwo+6,jtwo+6,1)             
            enddo
         enddo
      enddo
!
!.....Insert boundary condtions
      call boundary_vor(s5matrix_sm, temp)
      call finalizearray(s5matrix_sm)
 500  continue
!
!.....define RHS vector
      do itri=1,numelms
         call calcfint(fintl, maxi, atri(itri), btri(itri), ctri(itri))
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
! it does not own, we call sumshareddofs so that these values
! get summed up for all values shared by multiple procs
! and then update these values
      call sumshareddofs(temp)
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

