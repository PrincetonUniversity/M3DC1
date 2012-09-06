module coils
  implicit none

  integer, parameter :: maxbundles = 200
  integer, parameter :: subcoils = 1
  integer, parameter :: maxcoils = maxbundles*subcoils**2

contains

  !======================================================
  ! load_coils
  ! ~~~~~~~~~~
  !
  ! reads coil and current data from file
  !======================================================
 subroutine load_coils(xc, zc, ic, numcoils, coil_filename, current_filename, &
      ntor)
   use math

   implicit none

   include 'mpif.h'

   real, intent(out), dimension(maxcoils) :: xc, zc  ! coordinates of each coil
   complex, intent(out), dimension(maxcoils) :: ic   ! current in each coil
   integer, intent(out) :: numcoils                  ! number of coils read
   character*(*) :: coil_filename, current_filename  ! input files
   integer, intent(in) :: ntor                       ! toroidal mode number

   real :: x, z, w, h, a1, a2, c, phase
   real, dimension(maxcoils) :: vbuff
   complex, dimension(maxcoils) :: cbuff

   integer :: fcoil, fcurr, i, j, k, s, ier, ibuff, rank
   real :: dx, dz

   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ier)

   xc = 0.
   zc = 0.
   ic = 0.
   numcoils = 0

   if(rank.eq.0) then
      ! Read coil data
      print *, "Reading coil data..", coil_filename

      fcoil = 10
      open(unit=fcoil,file=coil_filename,status="old",err=200,action='read')
      fcurr = 20
      open(unit=fcurr,file=current_filename,status="old",err=201,action='read')


      s = 0
      do i=1, maxbundles
         read(fcoil,'(6F12.4)',end=100) x, z, w, h, a1, a2
         write(*,'(A,2F12.4)') "x, z ", x, z
         read(fcurr,'(2F12.4)',end=100) c, phase
         write(*,'(A,2F12.4)') "current (kA), phase (deg): ", c, phase
         
         a1 = a1*pi/180.
         a2 = a2*pi/180.
         if(a2.ne.0.) a2 = a2 + pi/2.
         a1 = tan(a1)
         a2 = tan(a2)
         
         c = amu0 * 1000. * c / (subcoils**2) / twopi
         
         ! divide coils into sub-coils
         do j=1, subcoils
            do k=1, subcoils
               s = s + 1
               dx = w*(2.*j/subcoils - 1.)/2.
               dz = h*(2.*k/subcoils - 1.)/2.
               xc(s) = x + dx - dz*a2
               zc(s) = z + dz + dx*a1
               
               ic(s) = c*(cos(pi*phase/180.) &
                    - (0.,1.)*sin(pi*phase/180.))
            end do
         end do
      end do

      print *, "Error: too many coils!"

100   print *, "Read ", i-1, " coils."

      ! divide coils into sub-coils
      numcoils = (i-1)*subcoils**2
   
      goto 300

200   if(rank.eq.0) print *, "Error reading ", coil_filename
      goto 300
201   if(rank.eq.0) print *, "Error reading ", current_filename
      
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
   call mpi_allreduce(ic, cbuff, 2*maxcoils, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, ier)
   ic = cbuff
   call mpi_allreduce(numcoils, ibuff, 1, MPI_INTEGER, &
        MPI_SUM, MPI_COMM_WORLD, ier)
   numcoils = ibuff
 end subroutine load_coils


 subroutine field_from_coils(xc, zc, ic, nc, f, ipole, ierr)
   use field
   use mesh_mod

   implicit none

   include 'mpif.h'

   real, intent(in), dimension(nc) :: xc, zc   ! array of coil positions
   complex, intent(in), dimension(nc) :: ic    ! array of coil currents
   integer, intent(in) :: nc                   ! number of coils
   type(field_type), intent(inout) :: f        ! poloidal flux field
   integer, intent(in) :: ipole                ! type of field to add
   integer, intent(out) :: ierr

   integer :: i, numnodes, k, ier, itmp
   real :: x, phi, z
   real, dimension(nc) :: xp, zp
   real, dimension(dofs_per_node,maxcoils) :: g
   vectype, dimension(dofs_per_node) :: data

   numnodes = owned_nodes()
   g = 0.
   do i=1,numnodes
     
      call get_node_pos(i,x,phi,z)
      xp = x
      zp = z
     
      ! Field due to coil currents
      call gvect(xp,zp,xc,zc,nc,g(1:6,:),ipole,ierr)
      if(ierr.ne.0) exit

      call get_node_data(f, i, data)
      do k=1,nc
         data = data - g(:,k)*ic(k)
      end do
      call set_node_data(f, i, data)      
   end do

   call mpi_allreduce(ierr, itmp, 1, &
        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ier)
   ierr = itmp
   if(ierr.ne.0) return

   call finalize(f%vec)

 end subroutine field_from_coils

!============================================================
subroutine gvect(r,z,xi,zi,n,g,nmult,ierr)
  ! calculates derivatives wrt first argument
  use math

  implicit none

  integer, intent(in) :: n, nmult
  integer, intent(out) :: ierr
  real, dimension(n), intent(in) :: r, z, xi, zi
  real, dimension(6,n), intent(out) :: g
  
  real :: a0,a1,a2,a3,a4
  real :: b0,b1,b2,b3,b4
  real :: c1,c2,c3,c4
  real :: d1,d2,d3,d4

  real :: rpxi, rxi, zmzi, rk, ce, ck, term1, term2, rz, co
  real :: rksq, sqrxi, x

  integer :: i, imult

  data a0,a1,a2,a3,a4/1.38629436112,9.666344259e-2,                 &
       3.590092383e-2,3.742563713e-2,1.451196212e-2/
  data b0,b1,b2,b3,b4/.5,.12498593597,6.880248576e-2,               &
       3.328355346e-2,4.41787012e-3/
  data c1,c2,c3,c4/.44325141463,6.260601220e-2,                     &
       4.757383546e-2,1.736506451e-2/
  data d1,d2,d3,d4/.24998368310,9.200180037e-2,                     &
       4.069697526e-2,5.26449639e-3/

  ierr = 0
  if(nmult.le.0) then
     do i=1,n
        rpxi=r(i)+xi(i)
        rxi=r(i)*xi(i)
        zmzi=z(i)-zi(i)
        rksq=4.*rxi/(rpxi**2+zmzi**2)
        rk=sqrt(rksq)
        sqrxi=sqrt(rxi)
        x=1.-rksq
        if(x.eq.0) then
           print *, 'Error: current is precisely on node'
           ierr = 1
           return
        end if

        ce=1.+x*(c1+x*(c2+x*(c3+x*c4)))+                                  &
             x*(d1+x*(d2+x*(d3+x*d4)))*(-alog(x))
        ck=a0+x*(a1+x*(a2+x*(a3+x*a4)))+                                  &
             (b0+x*(b1+x*(b2+x*(b3+x*b4))))*(-alog(x))

        term1=2.*ck-2.*ce-ce*rksq/x
        term2=2.*xi(i)-rksq*rpxi
        
        g(1,i) =- sqrxi*(2.*ck-2.*ce-ck*rksq)/rk
        g(2,i)=-rk*0.25/sqrxi*(rpxi*term1                                 &
             +2.*xi(i)*(ce/x-ck))
        g(3,i)=-rk*0.25*zmzi/sqrxi*term1
        g(5,i)=0.0625*zmzi*(rk/sqrxi)**3*(rpxi*term1                      &
             +(ce-ck+2.*ce*rksq/x)*                                       &
             (term2)/x)
        g(6,i)=-rk*0.25/sqrxi*(term1*                                     &
             (1.-rksq*zmzi**2/(4.*rxi))+zmzi**2*rksq**2/(4.*rxi*x)        &
             *(ce-ck+2.*ce*rksq/x))
        g(4,i)=-rk*0.25/sqrxi*(-rksq*rpxi/(4.*rxi)*                       &
             (rpxi*term1+2.*xi(i)*(ce/x-ck))+term1-                       &
             rksq*rpxi/(4.*rxi*x)*(ce-ck+2.*ce*rksq/x)*                   &
             (term2)+rksq/(2.*r(i)*x)*(2.*ce/x-ck)*term2)
     end do
 
     return
  endif
     
  ! check for multipolar coils
  do i=1,n
     if(xi(i) .lt. 100.) cycle
     rz = zi(i)
     imult = int(xi(i) - 100.)
     if(imult .lt. 0 .or. imult.gt.10) then
        ! error
        ierr=39
        return
     endif
        
     select case(imult)
     case(0)
        ! even nullapole
        g(1,i) = twopi*rz**2
        g(2,i) = 0.
        g(3,i) = 0.
        g(4,i) = 0.
        g(5,i) = 0.
        g(6,i) = 0.

     case(1)
        ! odd nullapole
        g(1,i) = 0.
        g(2,i) = 0.
        g(3,i) = 0.
        g(4,i) = 0.
        g(5,i) = 0.
        g(6,i) = 0.

     case(2)
        ! even dipole
        g(1,i) = twopi*(r(i)**2 - rz**2)/2.
        g(2,i) = twopi*r(i)
        g(3,i) = 0.
        g(5,i) = 0.
        g(6,i) = 0.
        g(4,i) = twopi

     case(3)
        ! odd dipole
        co=twopi/rz
        g(1,i) = co*(r(i)**2*z(i))
        g(2,i) = co*(2.*r(i)*z(i))
        g(3,i) = co*(r(i)**2)
        g(5,i) = co*2.*r(i)
        g(6,i) = 0.
        g(4,i) = co*2*z(i)

     case(4)
        ! even quadrapole
        co=pi/(4.*rz**2)
        g(1,i) = co*(r(i)**4-4.*r(i)**2*z(i)**2 - 2.*r(i)**2*rz**2+rz**4)
        g(2,i) = co*(4.*r(i)**3-8.*r(i)*z(i)**2-4.*r(i)*rz**2)
        g(3,i) = co*(-8.*r(i)**2*z(i))
        g(5,i) = co*(-16.*r(i)*z(i))
        g(6,i) = co*(-8.*r(i)**2)
        g(4,i) = co*(12.*r(i)**2 - 8.*z(i)**2 - 4.*rz**2)
        
     case(5)
        ! odd quadrapole
        co=pi/(3.*rz**3)
        g(1,i) = co*r(i)**2*z(i)*(3.*r(i)**2-4.*z(i)**2-3.*rz**2)
        g(2,i) = co*(12.*r(i)**3*z(i)-8.*r(i)*z(i)**3-6.*r(i)*z(i)*rz**2)
        g(3,i) = co*(3.*r(i)**4 - 12.*r(i)**2*z(i)**2 - 3.*r(i)**2*rz**2)
        g(5,i) = co*(12.*r(i)**3-24.*r(i)*z(i)**2 - 6.*r(i)*rz**2)
        g(6,i) = co*(-24.*r(i)**2*z(i))
        g(4,i) = co*(36.*r(i)**2*z(i)-8.*z(i)**3-6.*z(i)*rz**2)

     case(6)
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

     case(7)
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

     case(8)
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

     case(9)
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

     case(10)
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

     case default
        print *, 'Error: unknown multipole ', imult
     end select
  end do
end subroutine gvect

! Int(cos(2 ntor x)/(1 + k2*sin(x)^2)^(3/2) dx)/pi     0 < x < pi/2
subroutine integral(npts,k2,ntor,f)
  use math

  implicit none

  real, intent(in) :: k2(npts)
  real, intent(out) :: f(npts)
  integer, intent(in) :: ntor, npts

  real :: dx, dx2, x
  real, dimension(npts) :: f0, f1, f2
  integer :: ix
  integer, parameter :: ixmax = 100

  dx = (pi/2.) / ixmax
  dx2 = dx/2.

  ! Integrate using Simpson's rule
  f = 0.
  f0 = 1.
  x = 0.
  do ix=1, ixmax
     f1 = cos(2.*ntor*(x+dx2))/sqrt((1.+k2*sin(x+dx2)**2)**3)
     f2 = cos(2.*ntor*(x+dx))/sqrt((1.+k2*sin(x+dx)**2)**3)
     f = f + (f0 + 4.*f1 + f2)
     f0 = f2
     x = x + dx
  end do
  f = f * (dx/6.) / pi
end subroutine integral


subroutine coil(curr, r1, z1, npts, r0, z0, ntor, fr, fphi, fz)
  implicit none

  complex, intent(in) :: curr
  real, intent(in) :: r1, z1
  integer, intent(in) :: ntor, npts
  real, intent(in), dimension(npts) :: r0, z0
  complex, intent(inout), dimension(npts) :: fr, fphi, fz

  real, dimension(npts) :: dr2, k2, temp, co
  complex, dimension(npts) :: fac

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
  fphi = fphi - fac*(0.,1.)*(z0 - z1)*(co - temp)/2.
  
end subroutine coil

subroutine surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr, fphi, fz)
  implicit none
  
  real, intent(in) :: r, z, dldZ, dldR
  integer, intent(in) :: npts, ntor
  real, intent(in), dimension(npts) :: r0, z0
  complex, intent(out), dimension(npts) :: fr, fphi, fz
  
  real, dimension(npts) :: dr2, fac, k2, temp1, temp2

  dr2 = (r - r0)**2 + (z - z0)**2
  k2 = 4.*r*r0/dr2
  fac = -ntor/sqrt(dr2**3)

  call integral(npts, k2, ntor+1, temp1)
  call integral(npts, k2, ntor-1, temp2)

  fr =  fac*(dldZ*r+dldR*(z0-z))*(temp1 - temp2)/2.
  fz = -fac*dldR*r0*(temp1 - temp2)/2.
  fphi =  -fac*(0.,1.)*(dldZ*r + dldR*(z0 - z))*(temp1 + temp2)/2.

  call integral(npts, k2, ntor, temp1)
  fphi = fphi + fac*(0.,1.)*dldZ*r0*temp1
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

  complex, intent(in) :: curr
  real, intent(in) :: r1, r2, z1, z2
  integer, intent(in) :: ntor, npts
  real, intent(in), dimension(npts) :: r0, z0
  complex, intent(inout), dimension(npts) :: fr, fphi, fz

  complex, dimension(npts) :: fr0, fphi0, fz0
  complex, dimension(npts) :: fr1, fphi1, fz1
  complex, dimension(npts) :: fr2, fphi2, fz2
  real :: l, dl, lmax, dldZ, dldR, r, z, dl2
  integer :: il
  integer, parameter :: ilmax = 100

  ! calculate conitrbution from coils
  call coil( curr,r1,z1,npts,r0,z0,ntor,fr,fphi,fz)
  call coil(-curr,r2,z2,npts,r0,z0,ntor,fr,fphi,fz)

  lmax = sqrt((r2-r1)**2 + (z2-z1)**2)
  dl = lmax/ilmax
  dl2 = dl/2.
  
  dldR = (r2 - r1)/lmax
  dldZ = (z2 - z1)/lmax

  ! calculate contribution from surface currents 
  r = r1
  z = z1
  call surface(r, z, dldR, dldZ, npts, r0, z0, ntor, fr0, fphi0, fz0)
  l = 0.
  do il=1, ilmax
     
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
     
     l = l + dl
  end do
  
end subroutine pane

end module coils
