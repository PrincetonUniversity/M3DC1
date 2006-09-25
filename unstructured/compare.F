      module arrayscomp
!.....arrays defined at all vertices
      real ,allocatable::                                           
     &     vel(:), vels(:), vel0(:), phi(:), phis(:), phi0(:),
     &     jphi(:),vor(:),com(:),den(:),den0(:)
      integer ,allocatable:: ist(:,:)

      integer  ntimer, ioddm, maxplots, ntimep, numvar, m, n,
     &     mi(20), ni(20)
      parameter(ioddm=2)
      real timer, gtri(20,18,ioddm), deex, deez, alx, alz

      data mi /0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0/
      data ni /0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5/
      end module arrayscomp
c==================================
      subroutine getparams(nout, mout, deexout, deezout)
      use arrayscomp
      implicit none
      integer nout, mout
      real deexout, deezout

      nout = n
      mout = m
      deexout = deex
      deezout = deez
      
      return
      end

c==================================
c#define DRIVER
#ifdef DRIVER
      program driver

      implicit none
      integer numvar

      numvar = 1
      call compload()

      call comptecplot(numvar)
      call compfree()

      stop
      end
#endif

c===================

      subroutine compload()
      use arrayscomp
      implicit none
      integer iodd, k, jjj, jj, j, q, ll, itri, inc, i
      real ti(20,20),rot(18,18), sum, atri, btri, ctri, ttri

      call comprdrestart()
      do iodd=1,ioddm
         if(iodd .eq. 1) then
            ttri = 0.
            atri = 0.
            btri = deex
            ctri = deez
         else if(iodd .eq. 2) then
            ttri = atan2(deez,deex)
            atri = deex*cos(ttri)
            btri = deez*sin(ttri)
            ctri = deez*cos(ttri)
         endif
c         write(*,*) 'abc', atri,btri, ctri, deex, deez
         call comptmatrix(ti,20,atri,btri,ctri)
                               
!     calculate the rotation matrix rot
         call comprotation(rot,18,ttri)
!
!     form the matrix g by multiplying ti and rot
         do k=1,20
            do jjj=1,3
               do jj=1,6
                  j = (jjj-1)*6 + jj
                  sum = 0
                  do q = 1,18
                     sum = sum + ti(k,q)*rot(q,j)
                  enddo
                  gtri(k,j,iodd) = sum
               enddo
            enddo
         enddo
      enddo

      itri = 0
      do j=1,m-1
         do i=1,n-1
            itri = itri + 1
            ist(itri,1) = i+n*(j-1)-1
            ist(itri,2) = i+1+n*(j-1)-1
            ist(itri,3) = i+1+j*n-1
            ist(itri+(n-1)*(m-1),1) = i+n*(j-1)-1
            ist(itri+(n-1)*(m-1),2) = i+1+j*n-1
            ist(itri+(n-1)*(m-1),3) = i+n*j-1
         enddo
      enddo

      return
      end

!============================================================
      subroutine compevaluate(x, z, val, itype)
      use arrayscomp
      implicit none
      real rval, x, z, x1, z1, val
      integer itri, iodd, iii, ii, i, numvare, index, itype, k
      real wlocal(18),avector(20), atri, btri, ctri, ttri,
     &     si, eta, sum
!     
!.....itype=1 for first velocity (flux) variable
!
!.....evaluate the solution to get the value [ans] at one point (x,z)
!
!
      numvare = numvar          !will likely need to change this later
!     first find out what triangle x,z is in
      call compwhattri(x,z,itri,iodd,x1,z1)
!
!.....construct the 18 vector corresponding to this triangle
!.....calculate the index and local coordinates for this triangle
      do iii=1,3     
        do ii=1,6
           i = (iii-1)*6 + ii
           index = numvare*6*ist(itri,iii)+ii+(itype-1)*6
           wlocal(i) = phi(index) ! later change this to be for any dofs
        enddo
      enddo
!
!.....calculate local coordinates
      if(iodd .eq. 1) then
         ttri = 0.
         atri = 0.
         btri = deex
         ctri = deez
      else if(iodd .eq. 2) then
         ttri = atan2(deez,deex)
         atri = deex*cos(ttri)
         btri = deez*sin(ttri)
         ctri = deez*cos(ttri)
      else
         write(*,*) 'bad iodd'
         stop
      endif
      si  = (x-x1)*cos(ttri) + (z-z1)*sin(ttri) - btri
      eta =-(x-x1)*sin(ttri) + (z-z1)*cos(ttri)
      if(iodd .eq. 2) then
c         write(*,*) 'si', si,eta,iodd, x, z, x1, z1
c         write(*,*) atri, btri,ctri, ttri
      endif
!
!
!     calculate the function value corresponding to this point
      do i=1,20
         sum = 0.
         do k=1,18
            sum = sum + gtri(i,k,iodd)*wlocal(k)
         enddo
         avector(i) = sum
      enddo
 
!
!.....evaluate the polynomial and second derivative
      rval = 0.
      do i=1,20
         rval = rval + avector(i)*si**mi(i)*eta**ni(i)
      enddo
      
      val = rval

      return
      end


!============================================================

      subroutine comprdrestart()
      use arrayscomp
      implicit none
      integer mmnn18, mmnn6, itime, j1, k

      open(56,file='truC1restartout20ts',form='formatted',status='old')
      read(56,101) n,m,numvar,alx, alz
 101  format(I5, 1X, I5, 1X, I5, 1X, E12.4, 1X, E12.4)
      deex = alx/(n-1.)
      deez = alz/(m-1.)
      call compspace()

      mmnn18 = (m+1)*n*numvar*6
      do j1=1,mmnn18
         read(56,1204) k, vel(j1)
      enddo
      do j1=1,mmnn18
         read(56,1204) k, vel0(j1)
      enddo
      do j1=1,mmnn18
         read(56,1204) k, vels(j1)
      enddo
      do j1=1,mmnn18
         read(56,1204) k, phi(j1)
      enddo
      do j1=1,mmnn18
         read(56,1204) k, phi0(j1)
      enddo
      do j1=1,mmnn18
         read(56,1204) k, phis(j1)
      enddo
!
      mmnn6 = 6*m*n
      do j1=1,mmnn6
         read(56,1204) k, jphi(j1)
      enddo
      do j1=1,mmnn6
         read(56,1204) k, vor(j1)
      enddo
      do j1=1,mmnn6
         read(56,1204) k, com(j1)
      enddo
      do j1=1,mmnn6
         read(56,1204) k, den(j1)
      enddo
      do j1=1,mmnn6
         read(56,1204) k, den0(j1)
      enddo

      read(56,1204) ntimer,timer
      close(56)
 1204 format(I5, 1X, E12.4)
      write(*,1201) ntimer,timer
 1201 format("* * * restart file read at cycle,time = ",i5,1pe12.4,     &
     &     "  * * *")
      return
      end
!============================================================
      subroutine compspace()
!
!     allocates space for big arrays
      use arrayscomp
      implicit none
      integer mmnn, mmnn18, mmnn6
      mmnn = m*n
      mmnn18 = (m+1)*n*numvar*6
      mmnn6  = m*n*6
!
!.....arrays defined at all vertices
      allocate(vel(mmnn18), vels(mmnn18), vel0(mmnn18), 
     &     phi(mmnn18), phis(mmnn18), phi0(mmnn18),
     &     jphi(mmnn6), vor(mmnn6), com(mmnn6),den(mmnn6),
     &     den0(mmnn6)) 
      allocate(ist(2*m*(n+1),3))
      return
      end
!
!============================================================
      subroutine compfree()
      use arrayscomp
      implicit none

      deallocate(vel, vels, vel0, phi, phis, phi0,
     &     jphi, vor, com,den, den0)
      deallocate(ist)
      return
      end
!============================================================
      subroutine compwhattri(x,z,itri,iodd,xref,zref)
      use arrayscomp
      implicit none
      integer nmax, mmax, irect, jrect, ll, itri, iodd
      real x, z, xref, zref, x1, z1
      nmax = n-1
      mmax = m-1
!.....first determine what rectangle x lies within
      irect = ifix(x/deex) + 1
      jrect = ifix(z/deez) + 1
      if(irect.gt.nmax) irect = nmax
      if(jrect.gt.mmax) jrect = mmax
!
!     next determine what region number
      ll = irect + (jrect-1)*(nmax)
!
!     determine if odd or even
      iodd = 1
      x1 = (irect-1)*deex
      z1 = (jrect-1)*deez
      if((x- x1      )**2 + (z-(z1+deez))**2 .lt.                       &
     &   (x-(x1+deex))**2 + (z- z1      )**2) iodd = 2
      itri = (iodd-1)*nmax*mmax + ll
      xref = x1
      zref = z1
      if(iodd.eq.2) zref = z1! + deez
!
      return
      end
!============================================================
      subroutine comptecplot(nv)
      use arrayscomp
      implicit none
      integer i, j, nv, mnew, nnew
      real val, temp, dum, x, z, deexnew, deeznew
      integer idivide
      open(10,file='compare.plt',form='formatted')
      
      idivide = 4
      mnew = idivide*(m-1)+1
      nnew = idivide*(n-1)+1
      deexnew = deex/idivide
      deeznew = deez/idivide
      write(10,1)
 1    format('VARIABLES = "X", "Y", "TEMP"')
      write(10,2) nnew*mnew,(nnew-1)*(mnew-1)
 2    format('ZONE N=',i8,' ET=QUADRILATERAL E=',i8,' F=FEPOINT')
      do j=1,mnew
         z = (j-1)*deeznew
         do i=1,nnew
            x = (i-1)*deexnew
            call compevaluate(x,z,val, 1)
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
!============================================================
      subroutine comprotation(rot,ndim,theta)
      dimension rot(ndim,*),r1(6,6)

!
!     calculate the rotation matrix R from Table 2 of Ref 2
!
      co = cos(theta)
      sn = sin(theta)
      do i=1,6
        do j=1,6
          r1(i,j) = 0
        enddo
      enddo
!
      r1(1,1) = 1.
!
      r1(2,2) = co
      r1(2,3) = sn
!
      r1(3,2) = -sn
      r1(3,3) = co
!
      r1(4,4) = co**2
      r1(4,5) = 2.*sn*co
      r1(4,6) = sn**2
!
      r1(5,4) = -sn*co
      r1(5,5) = co**2-sn**2
      r1(5,6) = sn*co
!
      r1(6,4) = sn**2
      r1(6,5) = -2*sn*co
      r1(6,6) = co**2
! 
      do i=1,18
        do j=1,18
          rot(i,j) = 0
        enddo
      enddo
!
      do i=1,6
        do j=1,6
          rot(i,j)       = r1(i,j)
          rot(i+6,j+6)   = r1(i,j)
          rot(i+12,j+12) = r1(i,j)
        enddo
      enddo
!
      return
      end
!============================================================
      subroutine comptmatrix(ti,ndim,a,b,c)
      dimension t(20,20),ti(ndim,*),wkspce(9400),ipiv(20)
      ierrorchk = 0

!   define the 20 x 20 Transformation Matrix that enforces the condition that
!   the nomal slope between triangles has only cubic variation..
!...first initialize to zero
!
!
      do i=1,20
        do j=1,20
          t(i,j) = 0.
        enddo
      enddo
!
!...Table 1 of Ref. [2]
!
      t(1,1)   = 1.
      t(1,2)   = -b
      t(1,4)   = b**2
      t(1,7)   = -b**3
      t(1,11)  = b**4
      t(1,16)  = -b**5
!
      t(2,2)   = 1
      t(2,4)   = -2*b
      t(2,7)   = 3*b**2
      t(2,11)  = -4*b**3
      t(2,16)  = 5*b**4
!
      t(3,3)   = 1
      t(3,5)   = -b
      t(3,8)   = b**2
      t(3,12)  = -b**3
!
      t(4,4)   = 2.
      t(4,7)   = -6.*b
      t(4,11)  = 12*b**2
      t(4,16)  = -20*b**3
!
      t(5,5)   = 1.
      t(5,8)   = -2.*b
      t(5,12)  = 3*b**2
!
      t(6,6)   = 2.
      t(6,9)   = -2*b
      t(6,13)  = 2*b**2
      t(6,17)  = -2*b**3
!
      t(7,1)   = 1.
      t(7,2)   = a
      t(7,4)   = a**2
      t(7,7)   = a**3
      t(7,11)  = a**4
      t(7,16)  = a**5
!
      t(8,2)   = 1.
      t(8,4)   = 2*a
      t(8,7)   = 3*a**2
      t(8,11)  = 4*a**3
      t(8,16)  = 5*a**4
!
      t(9,3)   = 1.
      t(9,5)   = a
      t(9,8)   = a**2
      t(9,12)  = a**3
!
      t(10,4)  = 2
      t(10,7)  = 6*a
      t(10,11) = 12*a**2
      t(10,16) = 20*a**3
!
      t(11,5)  = 1.
      t(11,8)  = 2.*a
      t(11,12) = 3*a**2
!
      t(12,6)  = 2.
      t(12,9)  = 2*a
      t(12,13) = 2*a**2
      t(12,17) = 2*a**3
!
      t(13,1)  = 1
      t(13,3)  = c
      t(13,6)  = c**2
      t(13,10) = c**3
      t(13,15) = c**4
      t(13,20) = c**5
!
      t(14,2)  = 1.
      t(14,5)  = c
      t(14,9)  = c**2
      t(14,14) = c**3
      t(14,19) = c**4
!
      t(15,3)  = 1.
      t(15,6)  = 2*c
      t(15,10) = 3*c**2
      t(15,15) = 4*c**3
      t(15,20) = 5*c**4
!
      t(16,4)  = 2.
      t(16,8)  = 2*c
      t(16,13) = 2*c**2
      t(16,18) = 2*c**3
!
      t(17,5)  = 1.
      t(17,9)  = 2*c
      t(17,14) = 3*c**2
      t(17,19) = 4*c**3
!
      t(18,6)  = 2.
      t(18,10) = 6*c
      t(18,15) = 12*c**2
      t(18,20) = 20*c**3
!
      t(19,16) = 5*a**4*c
      t(19,17) = 3*a**2*c**3 - 2*a**4*c
      t(19,18) = -2*a*c**4+3*a**3*c**2
      t(19,19) = c**5-4*a**2*c**3
      t(19,20) = 5*a*c**4
!
      t(20,16) = 5*b**4*c
      t(20,17) = 3*b**2*c**3 - 2*b**4*c
      t(20,18) = 2*b*c**4 - 3*b**3*c**2
      t(20,19) = c**5 - 4*b**2*c**3
      t(20,20) = -5*b*c**4
!
      if(ierrorchk.eq.0) go to 100
!
!.....analytic formula for determinant
      danaly = -64*(a+b)**17*c**20*(a**2+c**2)*(b**2+c**2)
!
!.....calculate determinant using nag
      ifail = 0
      do i=1,20
        do j=1,20
          ti(i,j) = t(i,j)
        enddo
      enddo
      det = 0.
!     call f03aaf(ti,20,20,det,wkspce,ifail)
!
      diff = det - danaly
      percent = 100* diff / danaly
!
      if(abs(percent) .gt. 1.e-12)                                      &
     &write(*,1001) percent
 1001 format("percent error in determinant",1pe12.4)
!
 100  continue
!
!
!.....calculate the inverse of t using NAG library routines
      info1 = 0
      info2 = 0
      do i=1,20
        do j=1,20
          ti(i,j) = t(i,j)
        enddo
      enddo
      call f07adf(20,20,ti,20,ipiv,info1)
      call f07ajf(20,ti,20,ipiv,wkspce,400,info2)
      if(info1.ne.0.or.info2.ne.0)write(*,1002) info1,info2
 1002 format('infos from f07a... ',3i5)
!
! error checking

      return
      end
