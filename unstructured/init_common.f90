module init_common

  implicit none

contains

subroutine init_planewave(outarr,x,phi,z,kx,kphi,kz,amp,phasex,phasephi,phasez)
  use element

  implicit none

  vectype, intent(out), dimension(MAX_PTS) :: outarr
  real, intent(in), dimension(MAX_PTS)  :: x, phi, z
  real, intent(in) :: kx, kphi, kz, amp, phasex, phasephi, phasez

  real, dimension(MAX_PTS) :: argx,argp,argz
  argx = kx*x     + phasex
  argp = kphi*phi + phasephi
  argz = kz*z     + phasez

  outarr = amp*cos(argp)*sin(argx)*sin(argz)

#ifdef USECOMPLEX
  outarr = outarr + amp*(0,1)*sin(argp)*sin(argx)*sin(argz)
#endif
end subroutine init_planewave


subroutine init_random(x,phi,z,outarr)
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  real, intent(in), dimension(MAX_PTS) :: x, phi, z
  vectype, intent(out), dimension(MAX_PTS) :: outarr
  integer, allocatable :: seed(:)
  integer :: i, j, n
  real, dimension(MAX_PTS) :: xx, zz, rsq, r, ri, ri3, rexp, co, sn
  vectype, dimension(MAX_PTS) :: temp
  real :: alx, alz, kx, kp, kz, random, roundoff

  call get_bounding_box_size(alx, alz)
!
! changed to be consistent with fortran95 7/12/2011 ... scj
  call random_seed(SIZE = n)
  allocate(seed(n))
  seed = 23
  call random_seed(PUT = seed)
  deallocate(seed)

  roundoff = 1.e-12

  xx = x - xzero
  zz = z - zzero
  rsq = xx**2 + zz**2 + roundoff
  r = sqrt(rsq)
  ri = 1./sqrt(rsq + roundoff)
  ri3 = ri/rsq
  rexp = exp(-rsq/ln)
  co = cos(phi)
  sn = sin(phi)
  outarr = 0.

  select case (icsym)

  case (0)   !  original option...no symmetry imposed
     do i=1,maxn
        kx = pi*i/alx
        do j=1, maxn
           kz = j*pi/alz
           kp = j
           call random_number(random)
           call init_planewave(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
                0.,0.,0.)
           outarr = outarr + temp
        end do
     end do

  case (1)  !   odd symmetry about midplane
     do i=1,maxn
        kx = pi*i/alx
        do j=1, maxn/2
           kz = 2.*j*pi/alz
           kp = j
           call random_number(random)
           call init_planewave(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
                0.,0.,0.)
           outarr = outarr + temp
        end do
     end do

  case (2)  !   even symmetry about midplane
     do i=1,maxn
        kx = pi*i/alx
        do j=1, maxn/2
           kz = (2.*j-1)*pi/alz
           kp = j
           call random_number(random)
           call init_planewave(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
                0.,0.,0.)
           outarr = outarr + temp
        end do
     end do

  case (3)  !   NOT RANDOM....start in (1,1) eigenfunction
     outarr = eps* r * rexp*(zz*co - xx*sn)
     
  end select
end subroutine init_random


end module init_common
