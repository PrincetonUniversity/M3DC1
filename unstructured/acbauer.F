      block data acbblockdata
      integer init1, init2
      double precision p_xmin, p_zmin, p_xmax, p_zmax
      common /acbblock/ p_xmin, p_zmin, p_xmax, p_zmax, init1, init2
      data init1/0/,init2/0/
      end


      subroutine getmincoord(xmin, zmin)
      implicit none
      double precision dp_xmin, dp_zmin
      integer init1,init2
      real p_xmin, p_zmin, p_xmax, p_zmax, xmin, zmin
      common /acbblock/p_xmin, p_zmin, p_xmax, p_zmax, init1, init2
      
      if(init1 .eq. 0) then
         call getmincoord2(dp_xmin, dp_zmin)
         p_xmin = dp_xmin 
         p_zmin = dp_zmin
         init1 = 1
      endif
      xmin = p_xmin
      zmin = p_zmin
      return
      end
!===================
      subroutine getmaxcoord(xmax, zmax)
      implicit none
      real xmax, zmax, p_xmax,p_zmax , p_zmin, p_xmin
      double precision dp_xmax, dp_zmax
      integer init1, init2
      common /acbblock/p_xmin, p_zmin, p_xmax, p_zmax, init1, init2
      
      if(init2.eq. 0) then
         call getmaxcoord2(dp_xmax, dp_zmax)
         p_xmax = dp_xmax 
         p_zmax = dp_zmax
         init2 = 1
      endif

      xmax = p_xmax
      zmax = p_zmax
      return
      end
!===================
      subroutine getboundingboxsize(alx, alz)
      implicit none
      real xmin, xmax, zmin, zmax, alx, alz
      call getmincoord(xmin, zmin)
      call getmaxcoord(xmax, zmax)
!      print *, "getboundingboxsize xmin zmin xmax zmax = ", 
!     1         xmin, zmin, xmax, zmax

      alx = xmax - xmin
      alz = zmax - zmin
      return
      end
!=====
! my routine to print out arrays
      subroutine printarray(aa,n,m, string)
      implicit none
      include "mpif.h"
      double precision aa(*), sum, pos, max
      integer n, i, m, j, myrank
      character*8 string
      
      write(*,*) 'printarray has string ',string
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, j)
      open(10,file='check.txt',form='formatted')
      sum = 0.
      max = 0.

#ifdef REGULAR
      do i=1,n
         write(10,1002) m, i, aa(i)
         pos = aa(i)
         if(pos .lt. 0.) pos = -pos
         if(max .lt. pos) max = pos
         sum = sum + pos
      enddo
#endif
#define NOPER
#ifdef NOPER
      j = 1
      do i=7,24
         write(10,1002) m, j, aa(i)
         pos = aa(i)
         if(pos .lt. 0.) pos = -pos
         if(max .lt. pos) max = pos
         sum = sum + pos
         j = j + 1
      enddo
      do i=31,48
         write(10,1002) m, j, aa(i)
         pos = aa(i)
         if(pos .lt. 0.) pos = -pos
         if(max .lt. pos) max = pos
         sum = sum + pos
         j = j + 1
      enddo
      do i=61,72
         write(10,1002) m, j, aa(i)
         pos = aa(i)
         if(pos .lt. 0.) pos = -pos
         if(max .lt. pos) max = pos
         sum = sum + pos
         j = j + 1
      enddo
      do i=1,6
         write(10,1002) m, j, aa(i)
         pos = aa(i)
         if(pos .lt. 0.) pos = -pos
         if(max .lt. pos) max = pos
         sum = sum + pos
         j = j + 1
      enddo
      do i=25,30
         write(10,1002) m, j, aa(i)
         pos = aa(i)
         if(pos .lt. 0.) pos = -pos
         if(max .lt. pos) max = pos
         sum = sum + pos
         j = j + 1
      enddo
      do i=49,60
         write(10,1002) m, j, aa(i)
         pos = aa(i)
         if(pos .lt. 0.) pos = -pos
         if(max .lt. pos) max = pos
         sum = sum + pos
         j = j + 1
      enddo
#endif

      write(10,1003) sum, max
 1002 format(i3, i6, e12.4)
 1003 format('sum is ',e12.4,' max is ',e12.4)
      close(10)
      if(m.eq.0) call safestop(0)
      return
      end
!==========
!=====
! my routine to print out arrays
      subroutine printlocalarray(aa,numberingid, x,z, string)
      implicit none
      include 'mpif.h'
      double precision aa(*), coords(3), sum, sum2, x, z, tol
      integer i, j, numnodes, ibegin, iendplusone
      integer numberingid, myrank, ndofs
      character*8 string
      
      tol = .001 ! tolerance for how close to node
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, i)
      sum = 0.
      write(*,*) 'printlocalarray has string ',string
      write(*,*) 'checking at',x,z
!      open(10,file='check.txt',form='formatted')
      call numnod(numnodes)
      do i=1,numnodes
         call xyznod(i, coords)
         if(coords(1) .gt. x-tol .and. 
     &        coords(1) .lt. x+tol .and. 
     &        coords(2) .gt. z-tol .and. 
     &        coords(2) .lt. z+tol) then
            call entdofs(numberingid, i, 0, ibegin, iendplusone)
            do j=ibegin,iendplusone-1
!               write(10,*) aa(j)
               write(*,*) myrank, aa(j)
               sum = sum + dabs(aa(j))
            enddo
         endif
      enddo 
!      close(10)
      write(*,*) 'local sum is ',sum, myrank
      call numprocdofs(numberingid, ndofs)
      sum = 0
      do i=1,ndofs
         sum = sum + dabs(aa(i))
      enddo
      call MPI_ALLREDUCE(sum, sum2, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, MPI_COMM_WORLD, i)
      if(myrank .eq. 0) write(*,*) 'global sum is',sum2
      call safestop(8292)
      return
      end
!==========
! my routine to print out arrays
      subroutine iprintarray(aa,n,m)
      implicit none
      integer aa(*), sum, pos
      integer n, i, m
      
      open(10,file='check.txt',form='formatted')
      sum = 0
      do i=1,n
         write(10,*) m, i, aa(i)
         pos = aa(i)
         if(pos .lt. 0.) pos = -pos
         sum = sum + pos
      enddo
      write(10,*) sum
      close(10)
      if(m.eq.0) call safestop(0)
      return
      end
!====================
      subroutine writeglobaldofs(vec, orderingid, myrank, maxrank,
     &     string)
      implicit none
      include 'mpif.h'
      integer i, j, orderingid, myrank, maxrank, icounter,
     &     mpitag, ierr, iprocdof, iarray(2)
      real vec(*)
      double precision dvals(2), pos
      character*8 string

      mpitag = 4998
      icounter = 0
      call numprocdofs(orderingid, iprocdof)
      dvals(1) = 0.
      dvals(2) = 0.
      if(myrank .ne. 0) then
         call mpi_recv(iarray, 2, MPI_INTEGER, myrank-1, mpitag,
     &        MPI_COMM_WORLD, j, ierr)
         call mpi_wait(i,j)
         icounter = iarray(1)
         open(10,file='check1.txt',form='formatted')
      else
         write(*,*) 'printarray has string ',string
         open(10,file='check0.txt',form='formatted')
      endif

      do i=1,iprocdof
         write(10,1002) 0, i+icounter, vec(i)
c         write(*,*) 0, i+icounter, vec(i)
         pos = vec(i)
         if(pos .lt. 0.) pos = -pos
         if(dvals(1) .lt. pos) dvals(1) = pos
         dvals(2) = dvals(2) + pos
      enddo
      write(10,1003) dvals(2), dvals(1)
      write(*,1003) dvals(2), dvals(1)
      close(10)
 1002 format(i3, i6, e12.4)
 1003 format('sum is ',e12.4,' max is ',e12.4)
      icounter = icounter + iprocdof
      if(myrank .ne. maxrank - 1) then
         iarray(1) = icounter
         iarray(2) = icounter
         call mpi_send(iarray, 2, MPI_INTEGER, myrank+1, mpitag,
     &        MPI_COMM_WORLD, ierr)
      endif
      call safestop(4841)
      return 
      end      
!============================================================
      subroutine printmesh()
      implicit none

      integer i, itri, nodeids(4), numfaces, numnodes
      real a, b, c, theta
      double precision coords(3)

      call numnod(numnodes)
      do i=1,numnodes
         call xyznod(i, coords)
         write(*,10)  i, coords(1), coords(2), coords(3)
      enddo
      call numfac(numfaces)
      do itri=1, numfaces
         call nodfac(itri,nodeids)
         write(*,*) itri, nodeids(1), nodeids(2), nodeids(3)
         call getelmparams(itri, a, b, c, theta)
         write(*,12) itri, a, b, c, theta
      enddo

!      call safestop (28)

 10   format(i3, x, 3(x, e8.3))
 12   format('params ',i3, x, 4(x, e10.4))

      return
      end

!============================================================
c     itype is for which numvar field
c     dum is the array of dofs/PPPLVector.val
c     numvare is the numvar for dum
      subroutine exportfield(dum)
      use p_data
      use t_data
      use basic
      implicit none

      real dum(*)

      call createfieldrq(1, dum)
           
      return
      end
!============================================================

      subroutine tecplot(dofs, nv)
      use t_data
      use basic
      implicit none
      integer i, j, nv, mnew, nnew
      real dofs(*), val, temp, x, z, deexnew, deeznew
      
      open(10,file='orig.plt',form='formatted')
      
      mnew = 41
      nnew = 41
      deexnew = 4./40.
      deeznew = 4./40.
      write(10,1)
 1    format('VARIABLES = "X", "Y", "TEMP"')
      write(10,2) nnew*mnew,(nnew-1)*(mnew-1)
 2    format('ZONE N=',i8,' ET=QUADRILATERAL E=',i8,' F=FEPOINT')
      do j=1,mnew
         z = (j-1)*deeznew
         do i=1,nnew
            x = (i-1)*deexnew
            call evaluate(x,z,val, temp, dofs, 1, nv)
            write(10,3) x,z,val
         enddo
      enddo

 3    format(f12.5,1x,f12.5,1xf12.5)
      do j=1,mnew-1
         do i=1,nnew-1
            write(10,4) i+(j-1)*nnew,i+1+(j-1)*nnew,i+1+j*nnew,
     &           i+j*nnew
         enddo
      enddo
 4    format(i8,1x,i8,1x,i8,1x,i8)

      close(10)

      return
      end


!==================
      subroutine printvalatnode(vals, numbering)
      implicit none
      real vals(*)
      double precision coords(3)
      integer numbering, i, j, nnodes, ibegin, iendplusone

      call numnod(nnodes)
      do i=1,nnodes
         call xyznod(i, coords)
         if(coords(1) .lt. .002d0 .and. coords(2) .lt. .002d0) then
            call entdofs(numbering, i, 0, ibegin, iendplusone)
            do j=ibegin,iendplusone-1
               write(*,*) j, vals(j)
            enddo
         endif
      enddo

      return
      end

