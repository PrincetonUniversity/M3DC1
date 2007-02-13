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
	 write(*,*) 'initialized min coordinate'
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
         write(*,*) 'initialized max coordinate'
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

      alx = xmax - xmin
      alz = zmax - zmin
      return
      end
!===================
      subroutine getelmparams(elmid, a, b, c, theta)
      implicit none
      real a, b, c, theta, x1, x2, x3, z1, z2, z3
      real x2p, x3p, z2p, z3p
      integer elmid, nodeids(4)
      double precision coords(3)
      
!      write(*,*) 'elmid ===',elmid
      call nodfac(elmid, nodeids)
      call xyznod(nodeids(1),coords)
      x1 = coords(1)
      z1 = coords(2)
      call xyznod(nodeids(2),coords)
      x2 = coords(1)
      z2 = coords(2)
      call xyznod(nodeids(3),coords)
      x3 = coords(1)
      z3 = coords(2)

      theta = atan2(z2-z1,x2-x1)
      x2 = x2 - x1
      z2 = z2 - z1
      x3 = x3 - x1
      z3 = z3 - z1
      x2p = cos(-theta) * x2 - sin(-theta) * z2
      z2p = sin(-theta) * x2 + cos(-theta) * z2
      x3p = cos(-theta) * x3 - sin(-theta) * z3
      z3p = sin(-theta) * x3 + cos(-theta) * z3
      if(z2p .gt. 0.0001 .or. z2p .lt. -0.0001) then
         write(*,*) "z2p should be 0.d0 but is ", z2p
      endif
      a = x2p-x3p
      b = x3p
      c = z3p
      if(c .le. 0.) then
         write(*,*) 'clockwise node ordering for element',elmid
         call safestop(6894)
      endif

      end subroutine getelmparams
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
      integer n, i, m, j, numnodes, ibegin, iendplusone
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
     &     mpitag, ierr, iprocdof, iarray(2), ndofs
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
!====================
      subroutine global2local(itri, x, z, si, eta)
      implicit none
      integer itri, nodes(4)
      real x, z, si, eta, a, b, c, theta, co, sn, x1, z1
      double precision coord(3)

      call nodfac(itri, nodes)
      call xyznod(nodes(1), coord) 
      x1 = coord(1)
      z1 = coord(2)
      call getelmparams(itri, a, b, c, theta)
      co = cos(theta)
      sn = sin(theta)
      si  = (x-x1)*co + (z-z1)*sn - b
      eta =-(x-x1)*sn + (z-z1)*co
      
      return
      end
!====================
      subroutine local2global(itri, si, eta, x, z)
      implicit none
      integer itri, nodes(4)
      real x, z, si, eta, a, b, c, theta, co, sn, x1, z1
      double precision coord(3)

      call nodfac(itri, nodes)
      call xyznod(itri, coord)
      x1 = coord(1)
      z1 = coord(2)
      call getelmparams(itri, a, b, c, theta)
      co = cos(theta)
      sn = sin(theta)
      x = x1 + co*(b+si) - sn*eta
      z = z1 + sn*(b+si) + co*eta
      
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
      subroutine getdeex(itri, deex)
      use t_data
      implicit none
      
      integer, intent(in) :: itri
      real, intent(out) :: deex
      
c$$$      real area, length
c$$$      double precision coords(9)
c$$$      integer nodeids(4)

      real :: alx, alz
      integer :: numnodes, numfaces, numedges, numregions

c$$$      area = .5*(atri(itri)+btri(itri))*ctri(itri)
c$$$      call nodfac(itri, nodeids)
c$$$      call xyznod(nodeids(1),coords)
c$$$      call xyznod(nodeids(2),coords(4))
c$$$      call xyznod(nodeids(3),coords(7))
c$$$      
c$$$      deex = (coords(4)-coords(1))*(coords(4)-coords(1)) +
c$$$     &       (coords(5)-coords(2))*(coords(5)-coords(2))
c$$$      length =  (coords(4)-coords(7))*(coords(4)-coords(7)) +
c$$$     &          (coords(5)-coords(8))*(coords(5)-coords(8))
c$$$      if(deex .gt. length) deex = length
c$$$      length = (coords(1)-coords(7))*(coords(1)-coords(7)) +
c$$$     &         (coords(2)-coords(8))*(coords(2)-coords(8))
c$$$
c$$$      if(deex .gt. length) deex = length
c$$$
c$$$      deex = 2.*area/sqrt(deex)

      call getboundingboxsize(alx,alz)
      call numglobalents(numnodes, numedges, numfaces, numregions)

      deex = sqrt(alx*alz/2.*numfaces)

      deex = alx/30.
      print *, "shouldn't be here!", deex

      return
      end

!============================================================
      subroutine exportfield2(itype,numvare, dum, iwhich)
      use p_data
      use t_data
      use basic
      implicit none
      integer itype, numvare, numelms, itri, idof, iii, ii, ndof,
     &     i, index, k, isforder, itensororder, iwhich
    
      real si, eta, x, z, wlocal(18), dum(*), rlocal(18), sum, 
     &     avector(20), temp
      double precision coords(3)
      double precision, allocatable:: vals(:), dx(:), dz(:)
      
!     
!.....itype=1 for first velocity (flux) variable
!
      itensororder = 1          ! the order of the tensor
      isforder = 5              ! the order of lagrange shape functions used
      call createfield("phi", itensororder, isforder, iwhich) 
      call getnumlagrangetriinterppoint(isforder, ndof)
      allocate(dx(ndof), dz(ndof), vals(ndof)) ! vals should take into account the tensor order for non-scalars!
      dx = 0.d0
      dz = 0.d0
      call numfac(numelms)
      do itri=1,numelms
         call getlagrangetriinterppoint(itri, isforder, dx, dz) ! interpolation points in global CS
         do idof=1,ndof
            x = dx(idof)
            z = dz(idof)
            call evaluate(x,z,vals(idof), temp, dum,itype,numvare, itri)
         enddo ! idof loop
         call setfieldontri(itri, itensororder, isforder, vals, iwhich)
      enddo ! itri loop

!.....update the field
      call outputfield(iwhich, iwhich)
      deallocate(dx,dz, vals)
           
      return
      end
!============================================================

      subroutine tecplot(dofs, nv)
      use t_data
      use basic
      implicit none
      integer i, j, nv, mnew, nnew
      real dofs(*), val, temp, dum, x, z, deexnew, deeznew
      
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

