module mesh_mod

! subroutine incbdef
! subroutine tridef
! subroutine offset
! subroutine tmatrix
! subroutine local
! subroutine whattri

implicit none

contains

!==========================================================
subroutine calcavector(itri, iodd, inarr, itype, numvare, avector)

  use t_data

  implicit none

  integer, intent(in) :: itri, iodd, itype, numvare
  real, dimension(*), intent(in) :: inarr
  real, dimension(20), intent(out) :: avector
    
  integer :: i, ii, iii, index, k
  real, dimension(18) :: wlocal

  ! construct the 18 vector corresponding to this triangle
  ! calculate the index and local coordinates for this triangle
  do iii=1,3
     do ii=1,6
        i = (iii-1)*6 + ii
        index = numvare*6*ist(itri,iii)+ii+(itype-1)*6
        wlocal(i) = inarr(index)
     enddo
  enddo

  ! calculate the function value corresponding to this point
  do i=1,20
     avector(i) = 0.
     do k=1,18
        avector(i) = avector(i) + gtri(i,k,iodd)*wlocal(k)
     enddo
  enddo

end subroutine calcavector
!===================================================
subroutine incbdef
  use basic
  use t_data

  implicit none

  integer :: ibign, ibigm, ibig

  ! jbig = ibig + incb(ibig,l),l=1,itypemax(ibig)
  
  ! region #1...interior
  do ibign=2,n-1
     do ibigm=2,m-1
        ibig = (ibigm-1)*n + ibign
        itypemax(ibig) = 9
        incb(ibig,1) = -n-1
        incb(ibig,2) = -n
        incb(ibig,3) = -n+1
        incb(ibig,4) =   -1
        incb(ibig,5) =    0
        incb(ibig,6) =   +1
        incb(ibig,7) =  n-1
        incb(ibig,8) =  n
        incb(ibig,9) =  n+1
     enddo
  enddo

  ! region #2...bottom
  ibigm=1
  do ibign=2,n-1
     ibig = (ibigm-1)*n + ibign
     itypemax(ibig) = 6
     incb(ibig,1) =   -1
     incb(ibig,2) =    0
     incb(ibig,3) =   +1
     incb(ibig,4) =  n-1
     incb(ibig,5) =  n
     incb(ibig,6) =  n+1
     if(jper.eq.1) then
        itypemax(ibig) = 9
        incb(ibig,7) = (m-1)*n-1
        incb(ibig,8) = (m-1)*n
        incb(ibig,9) = (m-1)*n+1
     endif
  enddo

  ! region #3...top
  ibigm = m
  do ibign=2,n-1
     ibig = (ibigm-1)*n + ibign
     itypemax(ibig) = 6
     incb(ibig,1) = -n-1
     incb(ibig,2) = -n
     incb(ibig,3) = -n+1
     incb(ibig,4) =   -1
     incb(ibig,5) =    0
     incb(ibig,6) =   +1
     if(jper.eq.1) then
        itypemax(ibig) = 9
        incb(ibig,7) = -(m-1)*n-1
        incb(ibig,8) = -(m-1)*n
        incb(ibig,9) = -(m-1)*n+1
     endif
  enddo
  
  ! region #4...right
  ibign = n
  do ibigm=2,m-1
     ibig = (ibigm-1)*n + ibign
     itypemax(ibig) = 6
     incb(ibig,1) = -n-1
     incb(ibig,2) = -n
     incb(ibig,3) =   -1
     incb(ibig,4) =    0
     incb(ibig,5) =  n-1
     incb(ibig,6) =  n
     if(iper.eq.1) then
        itypemax(ibig) = 9
        incb(ibig,7) = -2*n+1
        incb(ibig,8) =   -n+1
        incb(ibig,9) =     +1
     endif
  enddo

  ! region #5...left
  ibign= 1
  do ibigm=2,m-1
     ibig = (ibigm-1)*n + ibign
     itypemax(ibig) = 6
     incb(ibig,1) = -n
     incb(ibig,2) = -n+1
     incb(ibig,3) =    0
     incb(ibig,4) =   +1
     incb(ibig,5) =  n
     incb(ibig,6) =  n+1
     if(iper.eq.1) then
        itypemax(ibig) = 9
        incb(ibig,7) =    -1
        incb(ibig,8) =   n-1
        incb(ibig,9) = 2*n-1
     endif
  enddo

  ! region #6...bottom right
  ibign=n
  ibigm=1
  ibig = (ibigm-1)*n + ibign
  itypemax(ibig) = 4
  incb(ibig,1) =   -1
  incb(ibig,2) =    0
  incb(ibig,3) =  n-1
  incb(ibig,4) =  n
  if(iper.eq.1) then
     itypemax(ibig) = 6
     incb(ibig,5) =  -n+1
     incb(ibig,6) =    +1
     if(jper.eq.1) then
        itypemax(ibig) = 9
        incb(ibig,7) = n*(m-1)-1
        incb(ibig,8) = n*(m-1)
        incb(ibig,9) = n*(m-2)+1
     endif
  else
     if(jper.eq.1) then
        itypemax(ibig) = 6
        incb(ibig,5) = n*(m-1)-1
        incb(ibig,6) = n*(m-1)
     endif
  endif

  ! region #7...bottom left
  ibign=1
  ibigm=1
  ibig = (ibigm-1)*n + ibign
  itypemax(ibig) = 4
  incb(ibig,1) =    0
  incb(ibig,2) =   +1
  incb(ibig,3) =  n
  incb(ibig,4) =  n+1
  if(iper.eq.1) then
     itypemax(ibig) = 6
     incb(ibig,5) =   n-1
     incb(ibig,6) = 2*n-1
     if(jper.eq.1) then
        itypemax(ibig) = 9
        incb(ibig,7) = (m-1)*n
        incb(ibig,8) = (m-1)*n+1
        incb(ibig,9) = m*n-1
     endif
  else
     if(jper.eq.1) then
        itypemax(ibig) = 6
        incb(ibig,5) = (m-1)*n
        incb(ibig,6) = (m-1)*n+1
     endif
  endif

  ! region #8...top left
  ibign=1
  ibigm=m
  ibig = (ibigm-1)*n + ibign
  itypemax(ibig) = 4
  incb(ibig,1) = -n
  incb(ibig,2) = -n+1
  incb(ibig,3) =    0
  incb(ibig,4) =   +1
  if(iper.eq.1) then
     itypemax(ibig) = 6
     incb(ibig,5) =  n-1
     incb(ibig,6) =   -1
     if(jper.eq.1) then
        itypemax(ibig) = 9
        incb(ibig,7) = -(m-1)*n
        incb(ibig,8) = -(m-1)*n+1
        incb(ibig,9) = -(m-2)*n-1
     endif
  else
     if(jper.eq.1) then
        itypemax(ibig) = 6
        incb(ibig,5) = -(m-1)*n
        incb(ibig,6) = -(m-1)*n+1
     endif
  endif

  ! region #9...top right
  ibign=n
  ibigm=m
  ibig = (ibigm-1)*n + ibign
  itypemax(ibig) = 4
  incb(ibig,1) = -n-1
  incb(ibig,2) = -n
  incb(ibig,3) =   -1
  incb(ibig,4) =    0
  if(iper.eq.1) then
     itypemax(ibig) = 6
     incb(ibig,5) =   -n+1
     incb(ibig,6) = -2*n+1
     if(jper.eq.1) then
        itypemax(ibig) = 9
        incb(ibig,7) = -(m-1)*n
        incb(ibig,8) = -(m-1)*n-1
        incb(ibig,9) = -m*n+1
     endif
  else
     if(jper.eq.1) then
        itypemax(ibig) = 6
        incb(ibig,5) = -(m-1)*n
        incb(ibig,6) = -(m-1)*n-1
     endif
  endif

  return
end subroutine incbdef


subroutine tridef
  use p_data
  use t_data
  use basic

  implicit none

  real, dimension(20,20) :: ti
  real, dimension(18,18) :: rot
  integer :: itri, iodd, irect, jrect, k, j, jj, jjj, q, ll
  real :: ve, sum

  itri = 0
      
  !     start the loop over triangles within a rectangular region
  do iodd=1,ioddmx
     
     ! define a,b,c and theta
     go to (101,102,104,103,102,101,103,104),iodd
101  ttri(iodd) = 0.
     atri(iodd) = 0.
     btri(iodd) = deex
     ctri(iodd) = deez
     go to 105
102  ttri(iodd) = atan2(deez,deex)
     atri(iodd) = deex*cos(ttri(iodd))
     btri(iodd) = deez*sin(ttri(iodd))
     ctri(iodd) = deez*cos(ttri(iodd))
     ! testing new orientation for iodd=2
     ttri(iodd) = 3*pi/2.
     atri(iodd) = deez
     btri(iodd) = 0.
     ctri(iodd) = deex
     go to 105
103  ttri(iodd) = 0.
     atri(iodd) = deex
     btri(iodd) = 0
     ctri(iodd) = deez
     go to 105
104  ttri(iodd) = -atan2(deez,deex)
     ve = -ttri(iodd)
     atri(iodd) = deez*sin(ve)
     btri(iodd) = deex*cos(ve)
     ctri(iodd) = deex*sin(ve)
105  continue
      

     ! define the Inverse Transformation Matrix that enforces the 
     ! condition that the normal slope between triangles has only
     ! cubic variation

     call tmatrix(ti,20,atri(iodd),btri(iodd),ctri(iodd))
     

     ! calculate the rotation matrix rot
     call rotation(rot,18,ttri(iodd))
     
     ! form the matrix g by multiplying ti and rot
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

     ! start the loop over the rectangular regions to define vertex points
     ! that correspond to triangles
     do ll = 1,nreg  
        itri = itri + 1  
        call offset(itri,ll,iodd)
     enddo

  enddo ! on iodd
  
  do jrect = 1,m
     do irect = 1,n
        ll = irect + (jrect-1)*n
        xcord(ll) = deex*(irect-1)
        zcord(ll) = deez*(jrect-1)
     enddo
  enddo
     
  return
end subroutine tridef

!============================================================
subroutine tmatrix(ti,ndim,a,b,c)

  implicit none

  real, intent(out), dimension(ndim, *) :: ti
  integer, intent(in) :: ndim
  real, intent(in) :: a, b, c

  real, dimension(20,20) :: t
  real, dimension(9400) :: wkspce
  integer, dimension(20) :: ipiv(20)
  integer :: ierrorchk, i, j, info1, info2
  real :: danaly, det, diff, percent, ifail

  ierrorchk = 0

  ! define the 20 x 20 Transformation Matrix that enforces the condition that
  ! the nomal slope between triangles has only cubic variation..

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
!  call f03aaf(ti,20,20,det,wkspce,ifail)

  diff = det - danaly
  percent = 100* diff / danaly
  
  if(abs(percent) .gt. 1.e-12)                                      &
       write(*,1001) percent
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
  
  ! error checking
  return
end subroutine tmatrix

!============================================================
subroutine offset(itri,ll,iodd)
  use basic
  use t_data
  
  implicit none
  
  integer, intent(in) :: itri, ll, iodd
  
  integer :: inc, inc2, index, i1, i2, i3
  
  if(ioddmx.eq.2) go to 100

  ! ioddmx = 8
  inc = ifix(2.*(ll-1)/(n-1))
  index = 2*ll+n+(n+1)*inc - 1
  go to (201,202,203,204,205,206,207,208),iodd
  
  ! note:  this needs to be patched up for iper=1

201 ist(itri,2) = (index + 1)
  ist(itri,3) = (index + 1 + n)
  ist(itri,1) = index
  return
202 ist(itri,2) = (index + 1 + n)
  ist(itri,3) = (index + n)
  ist(itri,1) = index
  return
203 ist(itri,3) = (index + n)
  ist(itri,1) = (index + n -1)
  ist(itri,2) = index
  return
204 ist(itri,3) = (index + n -1)
  ist(itri,1) = (index - 1)
  ist(itri,2) = index
  return
205 ist(itri,3) = (index - 1)
  ist(itri,1) = (index - n - 1)
  ist(itri,2) = index
  return
206 ist(itri,1) = (index - n - 1)
  ist(itri,2) = (index - n)
  ist(itri,3) = index
  return
207 ist(itri,1) = (index - n)
  ist(itri,2) = (index - n + 1)
  ist(itri,3) = index
  return
208 ist(itri,2) = (index - n + 1)
  ist(itri,3) = (index + 1)
  ist(itri,1) = index
  return
  ! ioddmx = 2
100 continue

  inc = 0
  if(iper.eq.0) inc = ifix((ll-1.+1.e-8)/(n-1.))
  inc2 = 0
  if(iper.eq.1) inc2 = (mod(ll-1,n)+1)/n
  go to (301,302),iodd
301 continue
  i1 = ll-1+inc
  ist0(itri,1) = i1
  if(jper.eq.1 .and. i1.ge.m*n) i1 = i1-m*n
  i2 = ll  +inc-inc2*n
  ist0(itri,2) = i2
  if(jper.eq.1 .and. i2.ge.m*n) i2 = i2-m*n
  i3 = ll+n+inc-inc2*n
  ist0(itri,3) = i3
  if(jper.eq.1 .and. i3.ge.m*n) i3 = i3-m*n
  ist(itri,1) = i1
  ist(itri,2) = i2
  ist(itri,3) = i3
  return
302 continue
  i1 = ll-1+inc
  ist0(itri,2) = i1
  if(jper.eq.1 .and. i1.ge.m*n) i1 = i1-m*n
  i2 = ll+n+inc-inc2*n
  ist0(itri,3) = i2
  if(jper.eq.1 .and. i2.ge.m*n) i2 = i2-m*n
  i3 = ll+n-1+inc                
  ist0(itri,1) = i3
  if(jper.eq.1 .and. i3.ge.m*n) i3 = i3-m*n
  ist(itri,2) = i1
  ist(itri,3) = i2
  ist(itri,1) = i3
  return
  
end subroutine offset

!============================================================
subroutine local(xi,zi,ll,b,theta,iodd)
  use basic
  
  implicit none
  
  real, intent(out), dimension(*) :: xi,zi
  real, intent(in) :: b, theta
  integer, intent(in) :: ll, iodd
  
  real :: small, co, sn, x1, z1
  integer :: ix, iz
  
  small = 1.e-8

  ! calculates the coefficients xi and zi for region ll in the local
  ! coordinate expansion:  x = Sum{ xi(i)*si**mi(i)*eta**ni(i)}; {i=1,3}
  !                        z = Sum{ zi(i)*si**mi(i)*eta**ni(i)}; {i=1,3}
  if(iper.eq.0) then
     iz = ifix((ll-1.+small)/(n-1.))+1
     ix = ll - (n-1)*(iz-1)
  else
     iz = ifix((ll-1.+small)/n)+1
     ix = ll -      n*(iz-1)
  endif
      
  x1 = (ix-1)*deex
  z1 = (iz-1)*deez
  co = cos(theta)
  sn = sin(theta)
  xi(1) = x1 + co*b
  xi(2) =  co
  xi(3) = -sn
  zi(1) = z1 + sn*b
  if(iodd.eq.2) zi(1) = zi(1) + deez
  zi(2) =  sn
  zi(3) =  co
  return
end subroutine local
  
subroutine whattri(x,z,itri,iodd,xref,zref)
  use basic

  implicit none

  real, intent(in) :: x, z
  integer, intent(out) :: itri, iodd 
  real, intent(out) :: xref, zref

  real :: angle, x1, z1
  integer :: nmax, mmax, irect, jrect, ll
      
  if(ioddmx.eq.2) go to 100
  ! this will not work for iper=1
  
  ! ioddmx = 8
  irect = ifix(x/(2*deex)) + 1
  jrect = ifix(z/(2*deez)) + 1
  if(irect.gt.(n-1)/2) irect = (n-1)/2
  if(jrect.gt.(m-1)/2) jrect = (m-1)/2 
  x1 = ((irect-1)*2+1)*deex
  z1 = ((jrect-1)*2+1)*deez


  ! next determine what region number
  ll = irect + (jrect-1)*(n-1)/2

  angle = atan2(z-z1,x-x1) 
  if(angle.lt.0) angle = angle + 2.*pi
  iodd = angle/(pi/4.) + 1
  if(iodd.gt.8) iodd = 8
  if(iodd.lt.1) iodd = 1
  itri = (iodd-1)*(n-1)*(m-1)/4 + ll

  ! coordinates of point p1
  go to (1,1,2,3,4,4,5,1),iodd
1 xref = x1
  zref = z1
  go to 6
2 xref = x1 - deex
  zref = z1 + deez
  go to 6
3 xref = x1 - deex
  zref = z1
  go to 6
4 xref = x1 - deex
  zref = z1 - deez
  go to 6
5 xref = x1
  zref = z1 - deez
  go to 6
6 continue
  return
  
  ! ioddmx = 2
100 continue

  
  nmax = n-1
  if(iper.eq.1) nmax = n
  mmax = m-1
  if(jper.eq.1) mmax = m
  !.....first determine what rectangle x lies within
  irect = ifix(x/deex) + 1
  jrect = ifix(z/deez) + 1
  if(irect.gt.nmax) irect = nmax
  if(jrect.gt.mmax) jrect = mmax
  
  ! next determine what region number
  ll = irect + (jrect-1)*(nmax)

  ! determine if odd or even
  iodd = 1
  x1 = (irect-1)*deex
  z1 = (jrect-1)*deez
  if((x- x1      )**2 + (z-(z1+deez))**2 .lt.                       &
       (x-(x1+deex))**2 + (z- z1      )**2) iodd = 2
  itri = (iodd-1)*nmax*mmax + ll
  xref = x1
  zref = z1
  if(iodd.eq.2) zref = z1 + deez

  return
end subroutine whattri

end module mesh_mod
