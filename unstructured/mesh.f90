!==========================================================
! calcavector
! ~~~~~~~~~~~
!
! calculates the 20 polynomial coefficients avector
! of field inarr in element itri.
!==========================================================
subroutine calcavector(itri, inarr, itype, numvare, avector)

  use t_data

  implicit none

  integer, intent(in) :: itri, itype, numvare
  vectype, dimension(*), intent(in) :: inarr
  vectype, dimension(20), intent(out) :: avector

  integer :: ibegin, iendplusone
    
  integer :: i, ii, iii, k
  vectype, dimension(18) :: wlocal

  ! construct the 18 vector corresponding to this triangle
  ! calculate the index and local coordinates for this triangle
  do iii=1,3  
     call entdofs(numvare, ist(itri,iii)+1, 0, ibegin, iendplusone)
     do ii=1,6
        i = (iii-1)*6 + ii
        wlocal(i) = inarr(ibegin+ii-1+(itype-1)*6)
     enddo
  enddo

  ! calculate the function value corresponding to this point
  do i=1,20
     avector(i) = 0.
     do k=1,18
        avector(i) = avector(i) + gtri(i,k,itri)*wlocal(k)
     enddo
  enddo

end subroutine calcavector


!============================================================
! f 
! ~
! calculates the double integral over a triangle
! with base a + b and height c  [see ref 2, figure 1]
!
! f = Int{si**m*eta**n}d(si)d(eta)
!============================================================
real function f(m,n,a,b,c)
  use math

  implicit none

  integer, intent(in) :: m,n
  real, intent(in) :: a, b, c
  real :: anum, denom
  
  anum = c**(n+1)*(a**(m+1)-(-b)**(m+1))*factorial(m)*factorial(n)
  denom = factorial(m+n+2)
  f = anum/denom
end function f


!============================================================
! rotation
! ~~~~~~~~
! calculates the rotation matrix rot given angle theta
!============================================================
subroutine rotation(rot,ndim,theta)
  implicit none

  integer, intent(in) :: ndim
  real, intent(in) :: theta
  real, intent(out) :: rot(ndim,*)

  integer :: i, j
  real :: r1(6,6), co, sn

  co = cos(theta)
  sn = sin(theta)
  do i=1,6
     do j=1,6
        r1(i,j) = 0.
     enddo
  enddo

  r1(1,1) = 1.

  r1(2,2) = co
  r1(2,3) = sn

  r1(3,2) = -sn
  r1(3,3) = co

  r1(4,4) = co**2
  r1(4,5) = 2.*sn*co
  r1(4,6) = sn**2

  r1(5,4) = -sn*co
  r1(5,5) = co**2-sn**2
  r1(5,6) = sn*co

  r1(6,4) = sn**2
  r1(6,5) = -2.*sn*co
  r1(6,6) = co**2

  do i=1,18
     do j=1,18
        rot(i,j) = 0.
     enddo
  enddo

  do i=1,6
     do j=1,6
        rot(i,j)       = r1(i,j)
        rot(i+6,j+6)   = r1(i,j)
        rot(i+12,j+12) = r1(i,j)
     enddo
  enddo

  return
end subroutine rotation

!============================================================
! tmatrix
! ~~~~~~~
! define the 20 x 20 Transformation Matrix that enforces the condition that
! the nomal slope between triangles has only cubic variation..
!============================================================
subroutine tmatrix(ti,ndim,a,b,c)
  implicit none
  integer, intent(in) :: ndim
  real, intent(out) :: ti(ndim,*)
  real, intent(in) :: a, b, c

  integer :: i, j, ierrorchk, ifail, info1, info2
  real :: danaly, det, percent, diff
  real :: t(20,20), wkspce(9400)
  integer :: ipiv(20)
 
  ierrorchk = 0

  ! first initialize to zero

  do i=1,20
     do j=1,20
        t(i,j) = 0.
     enddo
  enddo

  ! Table 1 of Ref. [2]
  t(1,1)   = 1.
  t(1,2)   = -b
  t(1,4)   = b**2
  t(1,7)   = -b**3
  t(1,11)  = b**4
  t(1,16)  = -b**5

  t(2,2)   = 1
  t(2,4)   = -2*b
  t(2,7)   = 3*b**2
  t(2,11)  = -4*b**3
  t(2,16)  = 5*b**4

  t(3,3)   = 1
  t(3,5)   = -b
  t(3,8)   = b**2
  t(3,12)  = -b**3

  t(4,4)   = 2.
  t(4,7)   = -6.*b
  t(4,11)  = 12*b**2
  t(4,16)  = -20*b**3

  t(5,5)   = 1.
  t(5,8)   = -2.*b
  t(5,12)  = 3*b**2

  t(6,6)   = 2.
  t(6,9)   = -2*b
  t(6,13)  = 2*b**2
  t(6,17)  = -2*b**3

  t(7,1)   = 1.
  t(7,2)   = a
  t(7,4)   = a**2
  t(7,7)   = a**3
  t(7,11)  = a**4
  t(7,16)  = a**5

  t(8,2)   = 1.
  t(8,4)   = 2*a
  t(8,7)   = 3*a**2
  t(8,11)  = 4*a**3
  t(8,16)  = 5*a**4

  t(9,3)   = 1.
  t(9,5)   = a
  t(9,8)   = a**2
  t(9,12)  = a**3
  
  t(10,4)  = 2
  t(10,7)  = 6*a
  t(10,11) = 12*a**2
  t(10,16) = 20*a**3

  t(11,5)  = 1.
  t(11,8)  = 2.*a
  t(11,12) = 3*a**2

  t(12,6)  = 2.
  t(12,9)  = 2*a
  t(12,13) = 2*a**2
  t(12,17) = 2*a**3

  t(13,1)  = 1
  t(13,3)  = c
  t(13,6)  = c**2
  t(13,10) = c**3
  t(13,15) = c**4
  t(13,20) = c**5

  t(14,2)  = 1.
  t(14,5)  = c
  t(14,9)  = c**2
  t(14,14) = c**3
  t(14,19) = c**4

  t(15,3)  = 1.
  t(15,6)  = 2*c
  t(15,10) = 3*c**2
  t(15,15) = 4*c**3
  t(15,20) = 5*c**4

  t(16,4)  = 2.
  t(16,8)  = 2*c
  t(16,13) = 2*c**2
  t(16,18) = 2*c**3

  t(17,5)  = 1.
  t(17,9)  = 2*c
  t(17,14) = 3*c**2
  t(17,19) = 4*c**3

  t(18,6)  = 2.
  t(18,10) = 6*c
  t(18,15) = 12*c**2
  t(18,20) = 20*c**3

  t(19,16) = 5*a**4*c
  t(19,17) = 3*a**2*c**3 - 2*a**4*c
  t(19,18) = -2*a*c**4+3*a**3*c**2
  t(19,19) = c**5-4*a**2*c**3
  t(19,20) = 5*a*c**4
  
  t(20,16) = 5*b**4*c
  t(20,17) = 3*b**2*c**3 - 2*b**4*c
  t(20,18) = 2*b*c**4 - 3*b**3*c**2
  t(20,19) = c**5 - 4*b**2*c**3
  t(20,20) = -5*b*c**4

  if(ierrorchk.eq.0) go to 100

  ! analytic formula for determinant
  danaly = -64*(a+b)**17*c**20*(a**2+c**2)*(b**2+c**2)

  ! calculate determinant using nag
  ifail = 0
  do i=1,20
     do j=1,20
        ti(i,j) = t(i,j)
     enddo
  enddo
  det = 0.
  !     call f03aaf(ti,20,20,det,wkspce,ifail)

  diff = det - danaly
  percent = 100* diff / danaly
  
  if(abs(percent) .gt. 1.e-12)                                      &
       &write(*,1001) percent
1001 format("percent error in determinant",1pe12.4)
  
100 continue

  ! calculate the inverse of t using NAG library routines
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
 1002 format(3i5)
  
end subroutine tmatrix


!============================================================
! getelmparams
! ~~~~~~~~~~~~
! calculates element parameters a, b, c, and theta from
! node coordinates
!============================================================
subroutine getelmparams(elmid, a, b, c, theta)
  implicit none
  integer, intent(in) :: elmid
  real, intent(out) :: a, b, c, theta
  real :: x1, x2, x3, z1, z2, z3, x2p, x3p, z2p, z3p
  integer :: nodeids(4)
  double precision :: coords(3)
      
  call nodfac(elmid, nodeids)
  call xyznod(nodeids(1),coords)
  x1 = coords(1)
  z1 = coords(2)
  call xyznod(nodeids(2),coords)
  x2 = coords(1) - x1
  z2 = coords(2) - z1
  call xyznod(nodeids(3),coords)
  x3 = coords(1) - x1
  z3 = coords(2) - z1

  theta = atan2(z2,x2)
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
     write(*,*) 'ERROR: clockwise node ordering for element',elmid
     call safestop(6894)
  endif
  
end subroutine getelmparams



!============================================================
! tridef
! ~~~~~~
! populates the *tri, fint, and isval arrays
!============================================================
subroutine tridef
  use p_data
  use t_data
  use basic
  use arrays

  implicit none
  
  integer :: itri, i, j, k, ii, iii, numelms
  integer :: ibegin, iendplusone
  integer :: ibegin2, iendplusone2
  integer :: ibeginn, iendplusonen
  integer, dimension(4) :: nodeids
  real :: ti(20,20),rot(18,18), sum
  real :: f
  
  if(myrank.eq.0 .and. iprint.ge.1) print *, " Entering tridef..."

  call numfac(numelms)
     
  ! start the loop over triangles within a rectangular region
  do itri=1,numelms

     call nodfac(itri, nodeids)
     ist(itri,1) = nodeids(1)-1
     ist(itri,2) = nodeids(2)-1
     ist(itri,3) = nodeids(3)-1    

     ! define a,b,c and theta
     call getelmparams(itri,atri(itri),btri(itri),ctri(itri),ttri(itri))  

     ! define the Inverse Transformation Matrix that enforces the 
     ! condition that the normal slope between triangles has only 
     ! cubic variation
     call tmatrix(ti,20,atri(itri),btri(itri),ctri(itri))

     ! calculate the rotation matrix rot
     call rotation(rot,18,ttri(itri))

     ! form the matrix g by multiplying ti and rot
     do k=1,20
        do j=1,18
           sum = 0.
           do ii = 1,18
              sum = sum + ti(k,ii)*rot(ii,j)
           enddo
           gtri(k,j,itri) = sum
        enddo
     enddo

     ! calculate fint (used for analytic integrations)
     do i=-6,maxi
        do j=-6,maxi
           fint(i,j) = 0.
           if(i.ge.0 .and.j.ge.0) &
                fint(i,j) = f(i,j,atri(itri),btri(itri),ctri(itri))
        enddo
     enddo

     !  populate isval* arrays
     do iii=1,3
        call entdofs(1, ist(itri,iii)+1, 0, ibegin,  iendplusone )
        call entdofs(2, ist(itri,iii)+1, 0, ibegin2, iendplusone2)
        call entdofs(vecsize_vel,  ist(itri,iii)+1, 0, ibeginn, iendplusonen)
        do ii=1,6
           i = (iii-1)*6 + ii
           isval1(itri,i) = ibegin +ii-1 ! 6*ist(itri,iii)+ii
           isval2(itri,i) = ibegin2+ii-1 ! 12*ist(itri,iii)+ii
           isvaln(itri,i) = ibeginn+ii-1 ! numvar*6*ist(itri,iii)+ii
        enddo
     enddo
  end do

  if(myrank.eq.0 .and. iprint.ge.1) print *, " Exiting tridef."
end subroutine tridef


!============================================================
! nodcoord
! ~~~~~~~~
! get the global coordinates (x,z) of node inode
!============================================================
subroutine nodcoord(inode,x,z)
  use basic
  implicit none
  
  integer, intent(in) :: inode
  real, intent(out) :: x, z
  
  double precision :: coords(3)
  real :: xmin, zmin

  call xyznod(inode,coords)
      
  select case(nonrect)
  case(0)
     call getmincoord(xmin,zmin)
     x = coords(1) - xmin + xzero
     z = coords(2) - zmin + zzero
  case(1)
     x = coords(1)
     z = coords(2)
  end select
end subroutine nodcoord

!============================================================
! whattri
! ~~~~~~~
! Gets the element itri in which global coordinates (x,z)
! fall.  If no element on this process contains (x,z) then
! itri = -1.  Otherwise, xref and zref are global coordinates
! of the first node of element itri.
!============================================================
subroutine whattri(x,z,itri,xref,zref)
  use basic
  implicit none

  real, intent(in) :: x, z
  integer, intent(out) :: itri
  real, intent(out) :: xref, zref

  double precision :: x2, z2
  integer :: nodeids(4)
  real :: xmin, zmin

  call getmincoord(xmin,zmin)
  select case(nonrect)
  case(0)
     x2 = x - xzero + xmin
     z2 = z - zzero - zmin
  case(1)
     x2 = x
     z2 = z
  end select
  
  call usesearchstructure(x2,z2,itri)
      
  if(itri.lt.0) return

  call nodfac(itri,nodeids)
  call nodcoord(nodeids(1), xref, zref)
end subroutine whattri

!============================================================
! evaluate
! ~~~~~~~~
! calculates the value ans of field dum at global coordinates
! (x,z).  itri is the element containing (x,z).  (If this
! element does not reside on this process, itri=-1).
!============================================================
subroutine evaluate(x,z,ans,ans2,dum,itype,numvare,itri)
  
  use p_data
  use t_data
  use basic

  use nintegrate_mod

  implicit none

  include 'mpif.h'

  integer, intent(in) :: itype, numvare
  integer, intent(inout) :: itri
  real, intent(in) :: x, z
  vectype, intent(in) :: dum(*)
  real, intent(out) :: ans, ans2

  integer :: p, nodeids(4), ier
  real :: x1, z1
  vectype, dimension(20) :: avector
  real :: ri, si, eta, co, sn
  real :: term1, term2
  real, dimension(2) :: temp1, temp2
  integer :: hasval, tothasval

  ! evaluate the solution to get the value [ans] at one point (x,z)

  ! first find out what triangle x,z is in.  whattri
  ! returns itri, x1, and z1 with x1 and z1 being
  ! the coordinates of the first node/vertex

  if(itri.eq.0) then
     call whattri(x,z,itri,x1,z1)
  else
     call nodfac(itri,nodeids)
     call nodcoord(nodeids(1), x1, z1)
  endif

  ans = 0.
  ans2 = 0.


  ! if this process contains the point, evaluate the field at that point.
  if(itri.gt.0) then

     ! calculate local coordinates
     co = cos(ttri(itri))
     sn = sin(ttri(itri))
  
     si  = (x-x1)*co + (z-z1)*sn - btri(itri)
     eta =-(x-x1)*sn + (z-z1)*co

     ! calculate the inverse radius
     if(itor.eq.1) then
        ri = 1./x
     else
        ri = 1.
     endif

     ! calculate the value of the function
     call calcavector(itri, dum, itype, numvare, avector)
     
     do p=1,20
     
        term1 = si**mi(p)*eta**ni(p)
        term2 = 0.
        
        if(mi(p).ge.1) then
           if(itor.eq.1) then
              term2 = term2 - 2.*co*(mi(p)*si**(mi(p)-1) * eta**ni(p))*ri
           endif
           
           if(mi(p).ge.2) then
              term2 = term2 + si**(mi(p)-2)*(mi(p)-1)*mi(p) * eta**ni(p)
           endif
        endif
     
        if(ni(p).ge.1) then
           if(itor.eq.1) then
              term2 = term2 + 2.*sn*(si**mi(p) * eta**(ni(p)-1)*ni(p))*ri
           endif
           
           if(ni(p).ge.2) then
              term2 = term2 + si**mi(p) * eta**(ni(p)-2)*(ni(p)-1)*ni(p)
           endif
        endif
     
        ans = ans + avector(p)*term1
        ans2 = ans2 + avector(p)*term2
        hasval = 1
     enddo
  else
     hasval = 0
  endif
     

  ! Distribute the result if necessary
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(maxrank.gt.1) then
     ! Determine number of processes whose domains contain this point
     call mpi_allreduce(hasval, tothasval, 1, MPI_INTEGER, MPI_SUM, &
          MPI_COMM_WORLD, ier)

     ! Find the average value at this point over all processes containing
     ! the point.  (Each value should be identical.)
     temp1(1) = ans
     temp1(2) = ans2
     call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
          MPI_COMM_WORLD, ier)
     ans = temp2(1)/tothasval
     ans2 = temp2(2)/tothasval
  endif

end subroutine evaluate
