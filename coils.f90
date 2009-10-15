module coils
  implicit none

  integer, parameter :: maxbundles = 200
  integer, parameter :: subcoils = 1
  integer, parameter :: maxcoils = maxbundles*subcoils**2

!!$  character*8, parameter :: coil_filename = "coil.dat"
!!$  character*11, parameter :: current_filename = "current.dat"

contains

  !======================================================
  ! load_coils
  ! ~~~~~~~~~~
  !
  ! reads coil and current data from file
  !======================================================
 subroutine load_coils(xc, zc, ic, numcoils, coil_filename, current_filename)
   use basic

   implicit none

   include 'mpif.h'

   real, intent(out), dimension(maxcoils) :: xc, zc, ic
   character*(*) :: coil_filename, current_filename
   integer, intent(out) :: numcoils

   real :: x, z, w, h, a1, a2, c
   real, dimension(maxcoils) :: vbuff

   integer :: fcoil, fcurr, i, j, k, s, ier, ibuff
   real :: dx, dz

   xc = 0.
   zc = 0.
   ic = 0.
   numcoils = 0

   if(myrank.eq.0) then
      ! Read coil data
      if(iprint.ge.1) print *, "Reading coil data..", coil_filename

      fcoil = 10
      open(unit=fcoil,file=coil_filename,status="old",err=200,action='read')
      fcurr = 20
      open(unit=fcurr,file=current_filename,status="old",err=201,action='read')


      s = 0
      do i=1, maxbundles
         read(fcoil,'(6f12.4)',end=100) x, z, w, h, a1, a2
         read(fcurr,'(1f12.4)',end=100) c
         
         a1 = a1*pi/180.
         a2 = a2*pi/180.
         if(a2.ne.0.) a2 = a2 + pi/2.
         a1 = tan(a1)
         a2 = tan(a2)
         
         c = pi*4.e-7 * 1000. * c / (subcoils**2) / (2.*pi)
         
         do j=1, subcoils
            do k=1, subcoils
               s = s + 1
               dx = w*(2.*j/subcoils - 1.)/2.
               dz = h*(2.*k/subcoils - 1.)/2.
               xc(s) = x + dx - dz*a2
               zc(s) = z + dz + dx*a1
               
               ic(s) = c
            end do
         end do
      end do

      print *, "Error: too many coils!"

100   if(iprint.ge.1) &
           print *, "Read ", i-1, " coils."

      ! divide coils into sub-coils
      numcoils = (i-1)*subcoils**2
   
      goto 300

200   if(myrank.eq.0 .and. iprint.ge.1) &
           print *, "Error reading ", coil_filename
      goto 300
201   if(myrank.eq.0 .and. iprint.ge.1) &
           print *, "Error reading ", coil_filename
      
300   close(fcoil)
      close(fcurr)
   end if

   ! share coil data
   call mpi_allreduce(xc, vbuff, maxcoils, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, ier)
   xc = vbuff
   call mpi_allreduce(zc, vbuff, maxcoils, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, ier)
   zc = vbuff
   call mpi_allreduce(ic, vbuff, maxcoils, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, ier)
   ic = vbuff
   call mpi_allreduce(numcoils, ibuff, 1, MPI_INTEGER, &
        MPI_SUM, MPI_COMM_WORLD, ier)
   numcoils = ibuff
 end subroutine load_coils


 subroutine field_from_coils(xc, zc, ic, nc, field, isize, iplace, ipole)
   implicit none

   real, intent(in), dimension(nc) :: xc, zc, ic
   integer, intent(in) :: nc
   vectype, intent(inout), dimension(*) :: field
   integer, intent(in) :: isize, iplace, ipole

   integer :: i, numnodes, ineg, ibegin, iendplusone, iend, k
   real :: x, z
   real, dimension(nc) :: xp, zp
   real, dimension(6,maxcoils) :: g

   call numnod(numnodes)
   do i=1,numnodes
     
      call entdofs(isize, i, 0, ibegin, iendplusone)
      ibegin = ibegin + (iplace-1)*6
      iend = ibegin+5

      call nodcoord(i,x,z)
      xp = x
      zp = z
     
      ! Field due to coil currents
      call gvect(xp,zp,xc,zc,nc,g,ipole,ineg)
      do k=1,nc
         field(ibegin:iend) = field(ibegin:iend) - g(:,k)*ic(k)
      end do
   end do
 end subroutine field_from_coils


 subroutine load_field_from_coils(coil_filename, current_filename, &
      field, isize, iplace)
   
   implicit none

   character*(*) :: coil_filename, current_filename
   vectype, dimension(*), intent(inout) :: field
   integer, intent(in) :: isize, iplace

   real, dimension(maxcoils) :: xc, zc, ic
   integer :: nc

   call load_coils(xc, zc, ic, nc, coil_filename, current_filename)
!   call field_from_coils(xc, zc, ic, nc, field, isize, iplace, 0)
   call field_from_coils_2(xc, zc, ic, nc, field, isize, iplace)

 end subroutine load_field_from_coils

!============================================================
subroutine gvect(r,z,xi,zi,n,g,nmult,ineg)
  ! calculates derivatives wrt first argument

  implicit none

  integer, intent(in) :: n, nmult
  integer, intent(out) :: ineg
  real, dimension(n), intent(in) :: r, z, xi, zi
  real, dimension(6,n), intent(out) :: g
  
  real :: a0,a1,a2,a3,a4
  real :: b0,b1,b2,b3,b4
  real :: c1,c2,c3,c4
  real :: d1,d2,d3,d4
  real :: pi,tpi

  real :: rpxi, rxi, zmzi, rk, ce, ck, term1, term2, rz, co
  real :: rksq, sqrxi, x

  integer :: i, imult, imultp

  data a0,a1,a2,a3,a4/1.38629436112,9.666344259e-2,                 &
       3.590092383e-2,3.742563713e-2,1.451196212e-2/
  data b0,b1,b2,b3,b4/.5,.12498593597,6.880248576e-2,               &
       3.328355346e-2,4.41787012e-3/
  data c1,c2,c3,c4/.44325141463,6.260601220e-2,                     &
       4.757383546e-2,1.736506451e-2/
  data d1,d2,d3,d4/.24998368310,9.200180037e-2,                     &
       4.069697526e-2,5.26449639e-3/
  data pi,tpi/3.1415926535,6.283185308/
  
  if(nmult.gt.0) go to 101
  do i=1,n
     rpxi=r(i)+xi(i)
     rxi=r(i)*xi(i)
     zmzi=z(i)-zi(i)
     rksq=4.*rxi/(rpxi**2+zmzi**2)
     rk=sqrt(rksq)
     sqrxi=sqrt(rxi)
     x=1.-rksq
     ce=1.+x*(c1+x*(c2+x*(c3+x*c4)))+                                  &
          x*(d1+x*(d2+x*(d3+x*d4)))*(-alog(x))
     ck=a0+x*(a1+x*(a2+x*(a3+x*a4)))+                                  &
          (b0+x*(b1+x*(b2+x*(b3+x*b4))))*(-alog(x))
     
     term1=2.*ck-2.*ce-ce*rksq/x
     term2=2.*xi(i)-rksq*rpxi
     
     g(1,i) =- sqrxi*(2.*ck-2.*ce-ck*rksq)/rk
     g(2,i)=-rk*0.25/sqrxi*(rpxi*term1                                  &
          +2.*xi(i)*(ce/x-ck))
     g(3,i)=-rk*0.25*zmzi/sqrxi*term1
     g(5,i)=0.0625*zmzi*(rk/sqrxi)**3*(rpxi*term1                      &
          +(ce-ck+2.*ce*rksq/x)*                                            &
          (term2)/x)
     g(6,i)=-rk*0.25/sqrxi*(term1*                                     &
          (1.-rksq*zmzi**2/(4.*rxi))+zmzi**2*rksq**2/(4.*rxi*x)             &
          *(ce-ck+2.*ce*rksq/x))
     g(4,i)=-rk*0.25/sqrxi*(-rksq*rpxi/(4.*rxi)*                       &
          (rpxi*term1+2.*xi(i)*(ce/x-ck))+term1-                            &
          rksq*rpxi/(4.*rxi*x)*(ce-ck+2.*ce*rksq/x)*                        &
          (term2)+rksq/(2.*r(i)*x)*(2.*ce/x-ck)*term2)
100 end do
 
  return

101 continue
     
     ! check for multipolar coils
  do i=1,n
     if(xi(i) .lt. 100.) go to 200
     rz = zi(i)
     imult = ifix(xi(i) - 100.)
     if(imult .lt. 0 .or. imult.gt.10) go to 250
     imultp = imult + 1
     go to(10,11,12,13,14,15,16,17,18,19,20),imultp
10   continue
        
     ! even nullapole
     g(1,i) = tpi*rz**2
     g(2,i) = 0.
     g(3,i) = 0.
     g(4,i) = 0.
     g(5,i) = 0.
     g(6,i) = 0.
     go to 200
11   continue

     ! odd nullapole
     g(1,i) = 0.
     g(2,i) = 0.
     g(3,i) = 0.
     g(4,i) = 0.
     g(5,i) = 0.
     g(6,i) = 0.
     go to 200
12   continue

     ! even dipole
     g(1,i) = tpi*(r(i)**2 - rz**2)/2.
     g(2,i) = tpi*r(i)
     g(3,i) = 0.
     g(5,i) = 0.
     g(6,i) = 0.
     g(4,i) = tpi
     go to 200
13   continue
     
     ! odd dipole
     co=tpi/rz
     g(1,i) = co*(r(i)**2*z(i))
     g(2,i) = co*(2.*r(i)*z(i))
     g(3,i) = co*(r(i)**2)
     g(5,i) = co*2.*r(i)
     g(6,i) = 0.
     g(4,i) = co*2*z(i)
     go to 200
14   continue
        
     ! even quadrapole
     co=pi/(4.*rz**2)
     g(1,i) = co*(r(i)**4-4.*r(i)**2*z(i)**2 - 2.*r(i)**2*rz**2+rz**4)
     g(2,i) = co*(4.*r(i)**3-8.*r(i)*z(i)**2-4.*r(i)*rz**2)
     g(3,i) = co*(-8.*r(i)**2*z(i))
     g(5,i) = co*(-16.*r(i)*z(i))
     g(6,i) = co*(-8.*r(i)**2)
     g(4,i) = co*(12.*r(i)**2 - 8.*z(i)**2 - 4.*rz**2)
     go to 200
15   continue

     ! odd quadrapole
     co=pi/(3.*rz**3)
     g(1,i) = co*r(i)**2*z(i)*(3.*r(i)**2-4.*z(i)**2-3.*rz**2)
     g(2,i) = co*(12.*r(i)**3*z(i)-8.*r(i)*z(i)**3-6.*r(i)*z(i)*rz**2)
     g(3,i) = co*(3.*r(i)**4 - 12.*r(i)**2*z(i)**2 - 3.*r(i)**2*rz**2)
     g(5,i) = co*(12.*r(i)**3-24.*r(i)*z(i)**2 - 6.*r(i)*rz**2)
     g(6,i) = co*(-24.*r(i)**2*z(i))
     g(4,i) = co*(36.*r(i)**2*z(i)-8.*z(i)**3-6.*z(i)*rz**2)
     go to 200
16   continue

     ! even hexapole
     co=pi/(12.*rz**4)
     g(1,i) = co*(r(i)**6 - 12.*r(i)**4*z(i)**2 - 3.*r(i)**4*rz**2       &
          + 8.*r(i)**2*z(i)**4 + 12.*r(i)**2*z(i)**2*rz**2           &
          + 3.*r(i)**2*rz**4 - rz**6 )
     g(2,i)= co*(6.*r(i)**5 - 48.*r(i)**3*z(i)**2                       &
          - 12.*r(i)**3*rz**2 + 16.*r(i)*z(i)**4                     &
          + 24.*r(i)*z(i)**2*rz**2 + 6.*r(i)*rz**4 )
     g(3,i)= co*(-24.*r(i)**4*z(i) + 32.*r(i)**2*z(i)**3                &
          + 24.*r(i)**2*z(i)*rz**2)
     g(5,i)=co*(-96.*r(i)**3*z(i)+64.*r(i)*z(i)**3                     &
          + 48.*r(i)*z(i)*rz**2 )
     g(4,i)=co*(30.*r(i)**4-144.*r(i)**2*z(i)**2-36.*r(i)**2*rz**2     &
          + 16.*z(i)**4 + 24.*z(i)**2*rz**2 + 6.*rz**4)
     g(6,i)=co*(-24.*r(i)**4 + 96.*r(i)**2*z(i)**2 + 24.*r(i)**2*rz**2)
     go to 200
17   continue
        
     ! odd hexapole
     co=pi/(30.*rz**5)
     g(1,i) = co*(15.*r(i)**6*z(i) - 60.*r(i)**4*z(i)**3                 &
          - 30.*r(i)**4*z(i)*rz**2 + 24.*r(i)**2*z(i)**5             &
          + 40.*r(i)**2*z(i)**3*rz**2 + 15.*r(i)**2*z(i)*rz**4)
     g(2,i)= co*(90.*r(i)**5*z(i) - 240.*r(i)**3*z(i)**3                &
          - 120.*r(i)**3*z(i)*rz**2 + 48.*r(i)*z(i)**5               &
          + 80.*r(i)*z(i)**3*rz**2 + 30.*r(i)*z(i)*rz**4)
     g(3,i)= co*(15.*r(i)**6 - 180.*r(i)**4*z(i)**2                     &
          - 30.*r(i)**4*rz**2 + 120.*r(i)**2*z(i)**4                 &
          +120.*r(i)**2*z(i)**2*rz**2 + 15.*r(i)**2*rz**4)
     g(5,i)=co*(90.*r(i)**5 - 720.*r(i)**3*z(i)**2                     &
          - 120.*r(i)**3*rz**2 + 240.*r(i)*z(i)**4                   &
          + 240.*r(i)*z(i)**2*rz**2 + 30.*r(i)*rz**4)
     g(6,i)=co*(-360.*r(i)**4*z(i) + 480.*r(i)**2*z(i)**3              &
          + 240.*r(i)**2*z(i)*rz**2)
     g(4,i)=co*(450.*r(i)**4*z(i) - 720.*r(i)**2*z(i)**3               &
          - 360.*r(i)**2*z(i)*rz**2 + 48.*z(i)**5                    &
          + 80.*z(i)**3*rz**2 + 30.*z(i)*rz**4)
     go to 200
18   continue
        
     ! even octapole
     co=pi/(160.*rz**6)
     g(1,i) = co*(5.*r(i)**8 - 120.*r(i)**6*z(i)**2 - 20.*r(i)**6*rz**2  &
          + 240.*r(i)**4*z(i)**4 + 240.*r(i)**4*z(i)**2*rz**2        &
          + 30.*r(i)**4*rz**4 - 64.*r(i)**2*z(i)**6                  &
          - 160.*r(i)**2*z(i)**4*rz**2 - 120.*r(i)**2*z(i)**2*rz**4  &
          - 20.*r(i)**2*rz**6 + 5.*rz**8)
     g(2,i)= co*(40.*r(i)**7 - 720.*r(i)**5*z(i)**2 - 120.*r(i)**5*rz**2 &
          + 960.*r(i)**3*z(i)**4 + 960.*r(i)**3*z(i)**2*rz**2        &
          + 120.*r(i)**3*rz**4 - 128.*r(i)*z(i)**6                   &
          - 320.*r(i)*z(i)**4*rz**2 - 240.*r(i)*z(i)**2*rz**4        &
          - 40.*r(i)*rz**6)
     g(3,i)= co*(-240.*r(i)**6*z(i) + 960.*r(i)**4*z(i)**3              &
          + 480.*r(i)**4*z(i)*rz**2 - 384.*r(i)**2*z(i)**5           &
          - 640.*r(i)**2*z(i)**3*rz**2 - 240.*r(i)**2*z(i)*rz**4)
     g(5,i)=co*(-1440.*r(i)**5*z(i) + 3840.*r(i)**3*z(i)**3            &
          + 1920.*r(i)**3*z(i)*rz**2 - 768.*r(i)*z(i)**5             &
          - 1280.*r(i)*z(i)**3*rz**2 - 480.*r(i)*z(i)*rz**4)
     g(6,i)=co*(-240.*r(i)**6 + 2880.*r(i)**4*z(i)**2                  &
          + 480.*r(i)**4*rz**2 - 1920.*r(i)**2*z(i)**4               &
          - 1920.*r(i)**2*z(i)**2*rz**2 - 240.*r(i)**2*rz**4)
     g(4,i)=co*(280.*r(i)**6 - 3600.*r(i)**4*z(i)**2                   &
          - 600.*r(i)**4*rz**2 + 2880.*r(i)**2*z(i)**4               &
          + 2880.*r(i)**2*z(i)**2*rz**2                              &
          + 360.*r(i)**2*rz**4 - 128.*z(i)**6                        &
          - 320.*z(i)**4*rz**2 - 240.*z(i)**2*rz**4 - 40.*rz**6)
     go to 200
19   continue

     ! odd octapole
     co=pi/(140.*rz**7)
     g(1,i) = co*r(i)**2*z(i)*(35.*r(i)**6 - 280.*r(i)**4*z(i)**2        &
          - 105.*r(i)**4*rz**2 + 336.*r(i)**2*z(i)**4                &
          + 420.*r(i)**2*z(i)**2*rz**2 + 105.*r(i)**2*rz**4          &
          - 64.*z(i)**6 - 168.*z(i)**4*rz**2                         &
          - 140.*z(i)**2*rz**4 - 35.*rz**6)
     g(2,i)= co*(280.*r(i)**7*z(i) - 1680.*r(i)**5*z(i)**3              &
          - 630.*r(i)**5*z(i)*rz**2 + 1344.*r(i)**3*z(i)**5          &
          + 1680.*r(i)**3*z(i)**3*rz**2 + 420.*r(i)**3*z(i)*rz**4    &
          - 128.*r(i)*z(i)**7 - 336.*r(i)*z(i)**5*rz**2              &
          - 280.*r(i)*z(i)**3*rz**4 - 70.*r(i)*z(i)*rz**6)
     g(3,i)= co*(35.*r(i)**8-840.*r(i)**6*z(i)**2-105.*r(i)**6*rz**2    &
          + 1680.*r(i)**4*z(i)**4 + 1260.*r(i)**4*z(i)**2*rz**2      &
          + 105.*r(i)**4*rz**4 - 448.*r(i)**2*z(i)**6                &
          - 840.*r(i)**2*z(i)**4*rz**2 - 420.*r(i)**2*z(i)**2*rz**4  &
          - 35.*r(i)**2*rz**6)
     g(5,i)=co*(280.*r(i)**7 - 5040.*r(i)**5*z(i)**2                   &
          - 630.*r(i)**5*rz**2 + 6720.*r(i)**3*z(i)**4               &
          + 5040.*r(i)**3*z(i)**2*rz**2 + 420.*r(i)**3*rz**4         &
          - 896.*r(i)*z(i)**6 - 1680.*r(i)*z(i)**4*rz**2             &
          - 840.*r(i)*z(i)**2*rz**4 - 70.*r(i)*rz**6)
     g(6,i)=co*(-1680.*r(i)**6*z(i) + 6720.*r(i)**4*z(i)**3            &
          + 2520.*r(i)**4*z(i)*rz**2 - 2688.*r(i)**2*z(i)**5         &
          - 3360.*r(i)**2*z(i)**3*rz**2 - 840.*r(i)**2*z(i)*rz**4)
     g(4,i)=co*(1960.*r(i)**6*z(i) - 8400.*r(i)**4*z(i)**3             &
          - 3150.*r(i)**4*z(i)*rz**2 + 4032.*r(i)**2*z(i)**5         &
          + 5040.*r(i)**2*z(i)**3*rz**2 + 1260.*r(i)**2*z(i)*rz**4   &
          - 128.*z(i)**7 - 336.*z(i)**5*rz**2                        &
          - 280.*z(i)**3*rz**4 - 70.*z(i)*rz**6)
     go to 200
20   continue
        
     ! even decapole
     co=pi/(560.*rz**8)
     g(1,i) = co*(7.*r(i)**10 - 280.*r(i)**8*z(i)**2 - 35.*r(i)**8*rz**2 &
          + 1120.*r(i)**6*z(i)**4 + 840.*r(i)**6*z(i)**2*rz**2       &
          + 70.*r(i)**6*rz**4 - 896.*r(i)**4*z(i)**6                 &
          - 1680.*r(i)**4*z(i)**4*rz**2 - 840.*r(i)**4*z(i)**2*rz**4 &
          - 70.*r(i)**4*rz**6 + 128.*r(i)**2*z(i)**8                 &
          + 448.*r(i)**2*z(i)**6*rz**2 + 560.*r(i)**2*z(i)**4*rz**4  &
          + 280.*r(i)**2*z(i)**2*rz**6 + 35.*r(i)**2*rz**8           &
          - 7.*rz**10)
     g(2,i)= co*(70.*r(i)**9 - 2240.*r(i)**7*z(i)**2                    &
          - 280.*r(i)**7*rz**2 + 6720.*r(i)**5*z(i)**4               &
          + 5040.*r(i)**5*z(i)**2*rz**2 + 420.*r(i)**5*rz**4         &
          - 3584.*r(i)**3*z(i)**6 - 6720.*r(i)**3*z(i)**4*rz**2      &
          - 3360.*r(i)**3*z(i)**2*rz**4 - 280.*r(i)**3*rz**6         &
          + 256.*r(i)*z(i)**8 + 896.*r(i)*z(i)**6*rz**2              &
          + 1120.*r(i)*z(i)**4*rz**4 + 560.*r(i)*z(i)**2*rz**6       &
          + 70.*r(i)*rz**8)
     g(3,i)= co*(-560.*r(i)**8*z(i) + 4480.*r(i)**6*z(i)**3             &
          + 1680.*r(i)**6*z(i)*rz**2 - 5376.*r(i)**4*z(i)**5         &
          - 6720.*r(i)**4*z(i)**3*rz**2 - 1680.*r(i)**4*z(i)*rz**4   &
          + 1024.*r(i)**2*z(i)**7 + 2688.*r(i)**2*z(i)**5*rz**2      &
          + 2240.*r(i)**2*z(i)**3*rz**4 + 560.*r(i)**2*z(i)*rz**6)
     g(5,i)=co*(-4480.*r(i)**7*z(i) + 26880.*r(i)**5*z(i)**3           &
          + 10080.*r(i)**5*z(i)*rz**2 - 21504.*r(i)**3*z(i)**5       &
          - 26880.*r(i)**3*z(i)**3*rz**2 - 6720.*r(i)**3*z(i)*rz**4  &
          + 2048.*r(i)*z(i)**7 + 5376.*r(i)*z(i)**5*rz**2            &
          + 4480.*r(i)*z(i)**3*rz**4 + 1120.*r(i)*z(i)*rz**6)
     g(6,i)=co*(-560.*r(i)**8 + 13440.*r(i)**6*z(i)**2                 &
          + 1680.*r(i)**6*rz**2 - 26880.*r(i)**4*z(i)**4             &
          - 20160.*r(i)**4*z(i)**2*rz**2 - 1680.*r(i)**4*rz**4       &
          + 7168.*r(i)**2*z(i)**6 + 13440.*r(i)**2*z(i)**4*rz**2     &
          + 6720.*r(i)**2*z(i)**2*rz**4 + 560.*r(i)**2*rz**6)
     g(4,i)=co*(630.*r(i)**8 - 15680.*r(i)**6*z(i)**2                  &
          - 1960.*r(i)**6*rz**2 + 33600*r(i)**4*z(i)**4              &
          + 25200.*r(i)**4*z(i)**2*rz**2 + 2100.*r(i)**4*rz**4       &
          - 10752.*r(i)**2*z(i)**6 - 20160.*r(i)**2*z(i)**4*rz**2    &
          - 10080.*r(i)**2*z(i)**2*rz**4 - 840.*r(i)**2*rz**6        &
          + 256.*z(i)**8 + 896.*z(i)**6*rz**2                        &
          + 1120.*z(i)**4*rz**4 + 560.*z(i)**2*rz**6                 &
          + 70.*rz**8)
     go to 200
200 end do

  return
250 continue
   
  ! error
  ineg=39
  
  return
end subroutine gvect

! Int(cos(2 ntor x)/(1 + k2*sin(x)^2)^(3/2) dx)/pi     0 < x < pi/2
subroutine integral(npts,k2,ntor,f)
  implicit none

  real, intent(in) :: k2(npts)
  real, intent(out) :: f(npts)
  integer, intent(in) :: ntor, npts

  real, parameter :: pi = 3.141592653589793
  real :: dx, dx2, x
  real, dimension(npts) :: f0, f1, f2

  ! If k2 is large, contribution will be small
  ! and simple asymptotic solution suffices
!!$  if(real(k2(1)).gt.10.) then
!!$     f = 4./(pi * sqrt(k2))
!!$     return
!!$  endif

  dx = (pi/2.) / 100.
  dx2 = dx/2.

  ! Integrate using Simpson's rule
  f = 0.
  f0 = 1.
  do x=0., pi/2., dx
     f1 = cos(2.*ntor*(x+dx2))/sqrt((1.+k2*sin(x+dx2)**2)**3)
     f2 = cos(2.*ntor*(x+dx))/sqrt((1.+k2*sin(x+dx)**2)**3)
     f = f + (f0 + 4.*f1 + f2)
     f0 = f2
  end do
  f = f * (dx/6.) / pi
end subroutine integral


subroutine coil(curr, r1, z1, npts, r0, z0, ntor, fr, fphi, fz)
  implicit none

  real, intent(in) :: curr
  real, intent(in) :: r1, z1
  integer, intent(in) :: ntor, npts
  real, intent(in), dimension(npts) :: r0, z0
  vectype, intent(inout), dimension(npts) :: fr, fphi, fz

  real, dimension(npts) :: dr2, k2, temp, co, fac

  dr2 = (r0 - r1)**2 + (z0 - z1)**2
  k2 = 4.*r1*r0/dr2

  fac = curr*r1/sqrt(dr2**3)

  ! calculate contribution from coils
  call integral(npts,k2,ntor,temp)
  fz = fz + fac*r1*temp
  
  call integral(npts,k2,ntor+1,co)
  call integral(npts,k2,ntor-1,temp)

  fr = fr + fac*(z0 - z1)*(co + temp)/2.
  fz = fz - fac*r0*(co + temp)/2.
  fphi = fphi - fac*(0,1)*(z0 - z1)*(co - temp)/2.
  
end subroutine coil

subroutine surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr, fphi, fz)
  implicit none
  
  real, intent(in) :: r, z, dldZ, dldR
  integer, intent(in) :: npts, ntor
  real, intent(in), dimension(npts) :: r0, z0
  vectype, intent(out), dimension(npts) :: fr, fphi, fz
  
  real, dimension(npts) :: dr2, fac, k2, temp1, temp2

  dr2 = (r - r0)**2 + (z - z0)**2
  k2 = 4.*r*r0/dr2
  fac = -ntor/sqrt(dr2**3)

  call integral(npts, k2, ntor+1, temp1)
  call integral(npts, k2, ntor-1, temp2)

  fr =  fac*(dldZ*r+dldR*(z0-z))*(temp1 - temp2)/2.
  fz = -fac*dldR*r0*(temp1 - temp2)/2.
  fphi =  -fac*(0,1)*(dldZ*r + dldR*(z0 - z))*(temp1 + temp2)/2.

  call integral(npts, k2, ntor, temp1)
  fphi = fphi + fac*(0,1)*dldZ*r0*temp1
end subroutine surface

!======================================================================
! pane
! ~~~~
! Adds 'window pane' contribution to magnetic field (fr, fphi, fz)
! at observation point r=(r0, z0) from two axisymmetric coils
! at (r1,z1) and (r2,z2) carrying currents curr and -curr
!======================================================================
subroutine pane(curr, r1, r2, z1, z2, npts, r0, z0, ntor, fr, fphi, fz)
  implicit none

  real, intent(in) :: curr
  real, intent(in) :: r1, r2, z1, z2
  integer, intent(in) :: ntor, npts
  real, intent(in), dimension(npts) :: r0, z0
  vectype, intent(inout), dimension(npts) :: fr, fphi, fz

  vectype, dimension(npts) :: fr0, fphi0, fz0
  vectype, dimension(npts) :: fr1, fphi1, fz1
  vectype, dimension(npts) :: fr2, fphi2, fz2
  real :: l, dl, lmax, dldZ, dldR, r, z, dl2

  ! calculate conitrbution from coils
  call coil( curr,r1,z1,npts,r0,z0,ntor,fr,fphi,fz)
  call coil(-curr,r2,z2,npts,r0,z0,ntor,fr,fphi,fz)

  lmax = sqrt((r2-r1)**2 + (z2-z1)**2)
  dl = lmax/100.
  dl2 = dl/2.
  
  dldR = (r2 - r1)/lmax
  dldZ = (z2 - z1)/lmax

  ! calculate contribution from surface currents 
  r = r1
  z = z1
  call surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr0, fphi0, fz0)
  do l=0., lmax, dl
     r = (r2*(l+dl2) + r1*(lmax-(l+dl2)))/lmax
     z = (z2*(l+dl2) + z1*(lmax-(l+dl2)))/lmax
     call surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr1, fphi1, fz1)
       
     r = (r2*(l+dl) + r1*(lmax-(l+dl)))/lmax
     z = (z2*(l+dl) + z1*(lmax-(l+dl)))/lmax
     call surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr2, fphi2, fz2)

     fr = fr + curr*dl*(fr0 + 4.*fr1 + fr2)/6.
     fphi = fphi + curr*dl*(fphi0 + 4.*fphi1 + fphi2)/6.
     fz = fz + curr*dl*(fz0 + 4.*fz1 + fz2)/6.

     fr0 = fr2
     fphi0 = fphi2
     fz0 = fz2
  end do
  
end subroutine pane


subroutine field_from_coils_2(xc, zc, ic, nc, fin, isize, iplace)
  use t_data
  use basic
  use sparse
  use arrays
  use nintegrate_mod
  use newvar_mod

  implicit none

  include 'mpif.h'

  real, intent(in), dimension(nc) :: xc, zc, ic
  integer, intent(in) :: nc
  vectype, intent(inout), dimension(*) :: fin
  integer, intent(in) :: isize, iplace

  vectype, allocatable :: psi(:)
  integer :: i, ii, iii, j, jj, jjj, k, l, itri, nelms, ier, i1, j1
  integer :: ibegin, iendplusone, jbegin, jendplusone
  vectype :: temp(3,3)

  integer, parameter :: imethod = 3
  integer :: inumb

  integer :: numnodes
  real :: x, z
  vectype :: divb, norm2, buff1(5), buff2(5), divbr, divbz, divbphi
  real :: norm

  if(myrank.eq.0) print *, 'IN FIELD_FROM_COILS_2'
  
  select case(imethod)
  case(0)               ! set psi=BR, p=BZ, i=bphi
     inumb = 3
  case(1)               ! find least-squares solution to 
     inumb = 1          ! dpsi/dR = -R*BR;  dpsi/dR = R*BZ
  case(2)
     inumb = 2
  case(3)               ! find least-squares solution to 
     inumb = 2          ! BR = -dpsi/DZ/R - d^2f/dRdphi
                        ! BZ =  dpsi/DR/R - d^2f/dZdphi
                        ! Bphi = R del_perp^2(f)
  end select


  call createvec(psi,inumb)
  psi = 0.

  call zerosuperlumatrix(brmatrix_sm, icomplex, inumb)

  call numfac(nelms)

  do itri=1,nelms
        
     call define_triangle_quadrature(itri,25)
     call define_fields(itri,0,1,0)

     temp79a = 0.    ! B_R
     temp79b = 0.    ! B_phi
     temp79c = 0.    ! B_Z

     do i=1, nc, 2
        call pane(ic(i),xc(i),xc(i+1),zc(i),zc(i+1),npoints,x_79,z_79,ntor,&
             temp79a,temp79b,temp79c)
     end do

     temp79a = -temp79a*2.*pi
     temp79b = -temp79b*2.*pi
     temp79c = -temp79c*2.*pi
     
     do iii=1,3
     call entdofs(inumb,  ist(itri,iii)+1, 0, ibegin, iendplusone)
     do ii=1,6
        i = (iii-1)*6 + ii
        i1 = ibegin + ii - 1

        do jjj=1,3
        call entdofs(inumb,  ist(itri,jjj)+1, 0, jbegin, jendplusone)
        do jj=1,6
           j = (jjj-1)*6 + jj
           j1 = jbegin + jj - 1

           temp = 0.
           select case(imethod)
           case(0)
              temp(1,1) = int2(g79(:,OP_1,i),g79(:,OP_1,j))
              temp(2,2) = temp(1,1)
              temp(3,3) = temp(1,1)

           case(1)
              temp(1,1) = int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j)) &
                        + int2(g79(:,OP_DR,i),g79(:,OP_DR,j))

           case(2)
#ifdef USECOMPLEX
              temp(1,1) = -int3(ri_79,g79(:,OP_1,i),g79(:,OP_DZ,j))
              temp(1,2) = int2(g79(:,OP_1,i),g79(:,OP_DRP,j))
              temp(2,1) = int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j))
              temp(2,2) = int2(g79(:,OP_1,i),g79(:,OP_DZP,j))
#endif            
           case(3)
#ifdef USECOMPLEX
              temp(1,1) = int3(ri2_79,g79(:,OP_DR,i),g79(:,OP_DR,j)) &
                   +      int3(ri2_79,g79(:,OP_DZ,i),g79(:,OP_DZ,j)) 
              temp(1,2) = int3(ri_79, g79(:,OP_DZ,i),g79(:,OP_DRP,j)) &
                   -      int3(ri_79, g79(:,OP_DR,i),g79(:,OP_DZP,j))
              temp(2,1) = int3(ri_79, g79(:,OP_DRP,i),g79(:,OP_DZ,j)) &
                   -      int3(ri_79, g79(:,OP_DZP,i),g79(:,OP_DR,j))
              temp(2,2) = int2(g79(:,OP_DRP,i),g79(:,OP_DRP,j)) &
                   +      int2(g79(:,OP_DZP,i),g79(:,OP_DZP,j)) &
                   +      int3(r2_79,g79(:,OP_LP,i),g79(:,OP_LP,j))
#endif
           end select
                    
           do k=1,inumb
              do l=1,inumb
                 call insertval(brmatrix_sm, temp(k,l), icomplex, &
                      i1+6*(k-1),j1+6*(l-1),1)
              end do
           end do
        end do
        end do

        select case(imethod)
        case(0)
           psi(i1   ) = psi(i1   ) + int2(g79(:,OP_1,i),temp79a)
           psi(i1+6 ) = psi(i1+6 ) + int3(r_79,g79(:,OP_1,i),temp79b)
           psi(i1+12) = psi(i1+12) + int2(g79(:,OP_1,i),temp79c)

        case(1)
           psi(i1) = psi(i1) &
                + int3(r_79,g79(:,OP_DR,i),temp79c) &
                - int3(r_79,g79(:,OP_DZ,i),temp79a)

        case(2)
           psi(i1  ) = psi(i1  ) + int2(g79(:,OP_1,i),temp79a)
           psi(i1+6) = psi(i1+6) + int2(g79(:,OP_1,i),temp79c)

        case(3)
#ifdef USECOMPLEX           
           psi(i1   ) = psi(i1   ) &
                + int3(ri_79,g79(:,OP_DR,i),temp79c) &
                - int3(ri_79,g79(:,OP_DZ,i),temp79a)
           psi(i1+6 ) = psi(i1+6 ) &
                - int2(g79(:,OP_DRP,i),temp79a) &
                - int2(g79(:,OP_DZP,i),temp79c) &
                + int3(r_79,g79(:,OP_LP,i),temp79b)
#endif
        end select
     end do
     end do
  end do

  call sumsharedppplvecvals(psi)
  call finalizematrix(brmatrix_sm)

  call solve(brmatrix_sm,psi,ier)
  
  select case(imethod)
  case(0)
     call copyvec(psi,1,inumb,field,psi_g,num_fields)
     call copyvec(psi,2,inumb,field,bz_g,num_fields)
     call copyvec(psi,3,inumb,field,chi_g,num_fields)

#ifdef USECOMPLEX
#define CONJUGATE(x) conjg(x)
#else
#define CONJUGATE(x) x
#endif

     divb = 0.
     divbr = 0.
     divbz = 0.
     divbz = 0.
     norm2 = 0.
     call numnod(numnodes)
     do i=1, numnodes
        call assign_local_pointers(i)
        call nodcoord(i,x,z)
        divbr = divbr + psi1_l(2) + psi1_l(1)/x
        divbz = divbz + chi1_l(3)
        divbphi = divbphi + (0,1)*ntor*bz1_l(1)/x**2
        divb = divb + psi1_l(2) + psi1_l(1)/x + chi1_l(3) &
             + (0,1)*ntor * bz1_l(1)/x**2
        norm2 = norm2 &
             + psi1_l(2)*CONJUGATE(psi1_l(2)) &
             + chi1_l(3)*CONJUGATE(chi1_l(3)) &
             + (ntor/x)**2*bz1_l(1)*CONJUGATE(bz1_l(1))
     end do
     if(maxrank.gt.1) then
        buff1(1) = divb
        buff1(2) = norm2
        buff1(3) = divbr
        buff1(4) = divbz
        buff1(5) = divbphi
        call mpi_allreduce(buff1, buff2, 2*5, MPI_DOUBLE_PRECISION, &
             MPI_SUM, MPI_COMM_WORLD, ier)
        divb = buff2(1)
        norm2 = buff2(2)
        divbr = buff2(3)
        divbz = buff2(4)
        divbphi = buff2(5)
     end if
     norm = real(sqrt(norm2))
     divb = divb / norm
     if(myrank.eq.0) then
        print *, 'Div(B), norm2 = ', divb, norm2
        print *, 'R, Z, phi = ', divbr, divbz, divbphi
     end if

  case(1)
     call copyvec(psi,1,inumb,field,psi_g,num_fields)
  case(2)
     call copyvec(psi,1,inumb,field,psi_g,num_fields)
     call copyvec(psi,2,inumb,field,bz_g,num_fields)
  case(3)
     call copyvec(psi,1,inumb,field,psi_g,num_fields)
     call copyvec(psi,2,inumb,bf,1,1)

     call newvarn(mass_matrix_lhs,field,bz_g,num_fields, &
          psi,2,inumb,bf_matrix_rhs,NV_NOBOUND,field)
  end select

  call deletevec(psi)

end subroutine field_from_coils_2


end module coils
