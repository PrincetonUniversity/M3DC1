module coils
  implicit none

  integer, parameter :: maxbundles = 20
  integer, parameter :: subcoils = 5
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
      iend = ibegin+5

      call nodcoord(i,x,z)
      xp = x
      zp = z
     
      ! Field due to coil currents
      call gvect(xp,zp,xc,zc,nc,g,ipole,ineg)
      do k=1,nc

         field(ibegin+(iplace-1)*6:iend+(iplace-1)*6) = &
              field(ibegin+(iplace-1)*6:iend+(iplace-1)*6) + g(:,k)*ic(k)
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
   call field_from_coils(xc, zc, ic, nc, field, isize, iplace, 0)

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



end module coils
