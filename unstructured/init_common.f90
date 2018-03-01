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
  real, dimension(MAX_PTS) :: xx, zz, rsq, r, theta
  vectype, dimension(MAX_PTS) :: temp, phase
  real :: alx, alz, kx, kp, kz, px, pp, pz, random, roundoff

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
           call random_number(px)
           call random_number(pp)
           call random_number(pz)
           px = 2.*pi*px
           pp = 2.*pi*pp
           pz = 2.*pi*pz
           call random_number(random)
           call init_planewave(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
                px,pp,pz)
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
     rsq = xx**2 + zz**2 + roundoff
     r = sqrt(rsq)

     outarr = eps*r*phase
     if(ln.gt.0.) outarr = outarr*exp(-rsq/ln)
             
  end select
end subroutine init_random

subroutine init_perturbations
  use basic
  use arrays
  use field
  use m3dc1_nint
  use newvar_mod
  use diagnostics

  implicit none

  type(field_type) :: psi_vec, phi_vec
  integer :: itri, numelms, i, izone, imr
  vectype, dimension(dofs_per_element) :: dofs

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining initial perturbations'

  call create_field(psi_vec)
  call create_field(phi_vec)

  psi_vec = 0.
  phi_vec = 0.

  numelms = local_elements()
  do itri=1,numelms
     call get_zone(itri, izone)

     if(izone.ne.1) cycle

     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,0,1,0)

     call eval_ops(itri, p_field(0), p079)
     call eval_ops(itri, psi_field(0), ps079)

     ps179 = 0.
     ph179 = 0.

     ! calculate perturbed fields
     call init_random(x_79-xmag, phi_79, z_79, ph179(:,OP_1))

     ph179(:,OP_1) = ph179(:,OP_1) + r_79*verzero

     ! apply mask
     if(p0 .gt. 0.) then 
        do i=1, npoints
           imr = magnetic_region(ps079(i,OP_1),ps079(i,OP_DR),ps079(i,OP_DZ), &
                x_79(i), z_79(i))
           if(imr.eq.0) then
              if(real(p079(i,OP_1)).gt.pedge) then
                 temp79a(i) = (p079(i,OP_1) - pedge)/p0
              else
                 temp79a(i) = 0.
              end if
           else
              temp79a(i) = 0.
           end if
        end do
        ph179(:,OP_1) = ph179(:,OP_1)*temp79a
     end if

     ! populate vectors for solves
        
     ! psi
     dofs = intx2(mu79(:,:,OP_1),ps179(:,OP_1))
     call vector_insert_block(psi_vec%vec,itri,1,dofs,VEC_ADD)
     
     ! phi
     dofs = intx2(mu79(:,:,OP_1),ph179(:,OP_1))
     call vector_insert_block(phi_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  ! do solves
  call newvar_solve(psi_vec%vec,mass_mat_lhs)
  psi_field(1) = psi_vec

  ! use dirichlet boundary conditions to avoid perturbations on boundary
!  call newvar_solve(phi_vec%vec,mass_mat_lhs_dc)
  call newvar_solve(phi_vec%vec,mass_mat_lhs)
  u_field(1) = phi_vec

  call destroy_field(psi_vec)
  call destroy_field(phi_vec)
end subroutine init_perturbations

subroutine den_eq
  use basic
  use arrays
  use diagnostics
  use math
  use mesh_mod
  use m3dc1_nint
  use newvar_mod
  use pellet

  implicit none

  type(field_type) :: den_vec
  integer :: itri, numelms, def_fields
  vectype, dimension(dofs_per_element) :: dofs
  real, dimension(MAX_PTS) :: n, p
  
  if((idenfunc.eq.0 .or. idenfunc.eq.4) .and. .not.(ipellet.gt.0 .and. linear.eq.1)) return

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Defining density equilibrium'
  call create_field(den_vec)
  
  def_fields = FIELD_PSI + FIELD_N

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,def_fields,1,0)

     select case(idenfunc)
     case(1)
        n079(:,OP_1) = den0*.5* &
             (1. + &
             tanh((real(ps079(:,OP_1))-(psibound+denoff*(psibound-psimin)))&
             /(dendelt*(psibound-psimin))))
        
     case(2)        
        n079(:,OP_1) = den0
        if(den0.ne.den_edge) then
           temp79a = ((real(ps079(:,OP_1))-psimin)/(psibound-psimin) - denoff)&
                /dendelt
           n079(:,OP_1) = n079(:,OP_1) &
                + .5*(den_edge-den0)*(1. + tanh(real(temp79a)))
        end if
     end select

     if(ipellet.gt.0 .and. linear.eq.1) then
        n = 0.
        p = 0.
        n079(:,OP_1) = n079(:,OP_1) + &
             pellet_deposition(x_79, phi_79, z_79, p, n, pellet_rate)
     end if

     dofs = intx2(mu79(:,:,OP_1),n079(:,OP_1))
     call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(den_vec%vec,mass_mat_lhs)
  den_field(0) = den_vec

  call destroy_field(den_vec)

end subroutine den_eq

subroutine den_per
  use basic
  use arrays
  use pellet
  use diagnostics
  use field
  use m3dc1_nint
  use newvar_mod
  implicit none

  type(field_type) :: den_vec
  integer :: numelms, itri
  vectype, dimension(dofs_per_element) :: dofs
  real, dimension(MAX_PTS) :: n, p

  if(ipellet.ge.0) return

  call create_field(den_vec)
  den_vec = 0.

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,0,1,0)

     n179(:,OP_1) = 0.
     if(ipellet.lt.0) then
        n = 0.
        p = 0.
        n179(:,OP_1) = n179(:,OP_1) + &
             pellet_deposition(x_79, phi_79, z_79, p, n, pellet_rate)
     end if

     dofs = intx2(mu79(:,:,OP_1),n179(:,OP_1))

     call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(den_vec%vec,mass_mat_lhs)
  den_field(1) = den_vec

  call destroy_field(den_vec)

end subroutine den_per


end module init_common
