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
  real, dimension(MAX_PTS) :: xx, zz, rsq, r, ri, ri3, rexp, theta
  vectype, dimension(MAX_PTS) :: temp, phase
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
  theta = atan2(zz,xx)
#ifdef USECOMPLEX
  phase = exp((0,1)*(ntor*phi - mpol*theta))
#else
  phase = cos(ntor*phi - mpol*theta)
#endif
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
     outarr = eps*r*rexp*phase
     
  end select
end subroutine init_random

subroutine init_perturbations
  use basic
  use arrays
  use field
  use m3dc1_nint
  use newvar_mod

  implicit none

  type(field_type) :: psi_vec, phi_vec
  integer :: itri, numelms, i
  vectype, dimension(dofs_per_element) :: dofs

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining initial perturbations'

  call create_field(psi_vec)
  call create_field(phi_vec)

  psi_vec = 0.
  phi_vec = 0.

  numelms = local_elements()
  
  do itri=1,numelms
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,0,1,0)

     call eval_ops(itri, p_field(0), p079)
     
     temp79a = (p079(:,OP_1) - pedge)/p0

     ps179 = 0.
     ph179 = 0.

     ! calculate perturbed fields
     call init_random(x_79-xmag, phi_79, z_79, ph179(:,OP_1))

     ph179(:,OP_1) = ph179(:,OP_1) + r_79*verzero
     
     ph179(:,OP_1) = ph179(:,OP_1)*temp79a

     ! populate vectors for solves
        
     ! psi
     do i=1, dofs_per_element
        dofs(i) = int2(mu79(:,OP_1,i),ps179(:,OP_1))
     end do
     call vector_insert_block(psi_vec%vec,itri,1,dofs,VEC_ADD)
     
     ! phi
     do i=1, dofs_per_element
        dofs(i) = int2(mu79(:,OP_1,i),ph179(:,OP_1))
     end do
     call vector_insert_block(phi_vec%vec,itri,1,dofs,VEC_ADD)
  end do
     
  ! do solves
  call newvar_solve(psi_vec%vec,mass_mat_lhs)
  psi_field(1) = psi_vec
  
  call newvar_solve(phi_vec%vec,mass_mat_lhs)
  u_field(1) = phi_vec

  call destroy_field(psi_vec)
  call destroy_field(phi_vec)
end subroutine init_perturbations

end module init_common
