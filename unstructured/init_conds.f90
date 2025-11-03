!==============================
subroutine constant_field(outarr, val)
  use element

  implicit none

  vectype, dimension(dofs_per_node), intent(out) :: outarr
  real, intent(in) :: val

  outarr(1) = val
  outarr(2:6) = 0.
#ifdef USE3D
  outarr(7:12) = 0.
#endif

end subroutine constant_field
!==============================
subroutine plane_wave(outarr, x, z, kx, kz, amp, phase)
  implicit none

  real, intent(in) :: x, z
  vectype, dimension(6), intent(out) :: outarr
  real, intent(in)  :: kx, kz, amp, phase

  real :: arg
  arg = kx*x + kz*z + phase

  outarr(1) =  amp*cos(arg)
  outarr(2) = -amp*sin(arg)*kx
  outarr(3) = -amp*sin(arg)*kz
  outarr(4) = -amp*cos(arg)*kx*kx
  outarr(5) = -amp*cos(arg)*kx*kz
  outarr(6) = -amp*cos(arg)*kz*kz
end subroutine plane_wave
!==============================
subroutine plane_wave2(outarr,x,phi,z,kx,kphi,kz,amp,phasex,phasephi,phasez)
  use element

  implicit none

  vectype, dimension(dofs_per_node), intent(out) :: outarr
  real, intent(in)  :: x, phi, z, kx, kphi, kz, amp, phasex, phasephi, phasez

  real :: argx,argp,argz,cox,cop,coz,six,sip,siz
  argx = kx*x     + phasex
  argp = kphi*phi + phasephi
  argz = kz*z     + phasez

  six = sin(argx)
  cox = cos(argx)
  sip = sin(argp)
  cop = cos(argp)
  siz = sin(argz)
  coz = cos(argz)

  outarr(1) =  amp*cop*six*siz
  outarr(2) =  amp*cop*cox*siz*kx
  outarr(3) =  amp*cop*six*coz*kz
  outarr(4) = -amp*cop*six*siz*kx*kx
  outarr(5) =  amp*cop*cox*coz*kx*kz
  outarr(6) = -amp*cop*six*siz*kz*kz
#ifdef USE3D
  outarr(7 ) = -amp*sip*six*siz*kphi
  outarr(8 ) = -amp*sip*cox*siz*kx*kphi
  outarr(9 ) = -amp*sip*six*coz*kz*kphi
  outarr(10) =  amp*sip*six*siz*kx*kx*kphi
  outarr(11) = -amp*sip*cox*coz*kx*kz*kphi
  outarr(12) =  amp*sip*six*siz*kz*kz*kphi
#endif
end subroutine plane_wave2
!==============================
subroutine add_angular_velocity(outarr, x,omega)
  use arrays

  implicit none

  vectype, dimension(6), intent(inout) :: outarr
  vectype, dimension(6), intent(in) :: omega
  real, intent(in) :: x

  vectype, dimension(6) :: temp

  temp(1) = x**2 * omega(1)
  temp(2) = x**2 * omega(2) + 2.*x*omega(1)
  temp(3) = x**2 * omega(3)
  temp(4) = x**2 * omega(4) + 4.*x*omega(2) + 2.*omega(1)
  temp(5) = x**2 * omega(5) + 2.*x*omega(3)
  temp(6) = x**2 * omega(6)

  outarr = outarr + temp

end subroutine add_angular_velocity
!=============================
subroutine add_product(x,a,b)
  vectype, dimension(6), intent(in) :: a,b
  vectype, dimension(6), intent(inout) :: x

  x(1) = x(1) + a(1)*b(1)
  x(2) = x(2) + a(1)*b(2) + a(2)*b(1)
  x(3) = x(3) + a(1)*b(3) + a(3)*b(1)
  x(4) = x(4) + a(1)*b(4) + a(4)*b(1) + 2.*a(2)*b(2)
  x(5) = x(5) + a(1)*b(5) + a(5)*b(1) + a(2)*b(3) + a(3)*b(2)
  x(6) = x(6) + a(1)*b(6) + a(6)*b(1) + 2.*a(3)*b(3)
end subroutine add_product
!=============================
subroutine random_per(x,phi,z,fac)
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  real, intent(in) :: x, phi, z
  vectype, intent(in), dimension(dofs_per_node) :: fac
  integer, allocatable :: seed(:)
  integer :: i, j, n
  real :: alx, alz, kx, kp, kz, xx, zz, random, rsq, r,ri,roundoff,ri3,rexp,co,sn
  vectype, dimension(dofs_per_node) :: temp

  call get_bounding_box_size(alx, alz)
!
! changed to be consistent with fortran95 7/12/2011 ... scj
  call random_seed(SIZE = n)
  allocate(seed(n))
  seed = 23
  call random_seed(PUT = seed)
  deallocate(seed)

  temp = 0.
  roundoff = 1.e-12

  xx = x - xzero
  zz = z - zzero
  rsq = xx**2 + zz**2 + roundoff
  r = sqrt(rsq)
  ri = 1./sqrt(rsq + roundoff)
  ri3 = ri/rsq
  rexp = exp(-rsq/ln)
  if(itor.eq.1) then
     co = cos(phi)
     sn = sin(phi)
  else  
     co = cos(phi/rzero)
     sn = sin(phi/rzero)
  end if

  do i=1,maxn
     kx = pi*i/alx
     select case (icsym)

     case (0)   !  original option...no symmetry imposed
     do j=1, maxn
        kz = j*pi/alz
        kp = j
        call random_number(random)
        call plane_wave2(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
             0.,0.,0.)
        call add_product(psi1_l,fac,temp)
        call random_number(random)
        call plane_wave2(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
             0.,0.,0.)
        call add_product(u1_l,fac,temp)
     end do

     case (1)  !   make U odd symmetry about midplane:  perturb only U
     do j=1, maxn/2
        kz = 2.*j*pi/alz
        kp = j
        call random_number(random)
        call plane_wave2(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
             0.,0.,0.)
        call add_product(u1_l,fac,temp)
     end do

     case (2)  !   make U even  symmetry about midplane:  perturb only U
     do j=1, maxn/2
        kz = (2.*j-1)*pi/alz
        kp = j
        call random_number(random)
        call plane_wave2(temp,xx,phi,zz,kx,kp,kz,2.*eps*(random-.5), &
             0.,0.,0.)
        call add_product(u1_l,fac,temp)
     end do

     case (3)  !   NOT RANDOM....start in (1,1) eigenfunction
     temp(1) = eps* r * rexp*(zz*co - xx*sn)
     temp(2) = eps* ri * rexp*(zz*xx*co - xx*xx*sn)   &
             - eps*(2./ln)* r * rexp*(zz*xx*co - xx*xx*sn)    &
             - eps* r * rexp*sn 
     temp(3) = eps* ri * rexp*(zz*zz*co - xx*zz*sn)   &
             - eps*(2./ln)* r * rexp*(zz*zz*co - xx*zz*sn)    &
             + eps* r * rexp*co
     temp(4) = -eps* ri3 * rexp*(zz*xx*xx*co - xx*xx*xx*sn)   &
             + eps*(2./ln)* r * rexp*((2.*zz*xx*xx/ln - zz)*co - (2*xx*xx*xx/ln - 3.*xx)*sn)    &
             + eps* ri * rexp*((zz - 4.*zz*xx*xx/ln)*co - (3*xx-4*xx**3/ln)*sn)               
          
     temp(5) = -eps* ri3 * rexp*(zz*zz*xx*co - zz*xx*xx*sn)   &
             + eps*(2./ln)* r * rexp*((2.*zz*zz*xx/ln - xx)*co - (2*zz*xx*xx/ln - zz)*sn)    &
             + eps* ri * rexp*((xx - 4.*zz*zz*xx/ln)*co - (zz-4*zz*xx**2/ln)*sn)               
         
     temp(6) = -eps* ri3 * rexp*(zz*zz*zz*co - xx*zz*zz*sn)   &
             + eps*(2./ln)* r * rexp*((2.*zz*zz*zz/ln - 3*zz)*co - (2*xx*zz*zz/ln - xx)*sn)    &
             + eps* ri * rexp*((3*zz - 4.*zz*zz*zz/ln)*co - (xx-4*xx*zz*zz/ln)*sn)   
#ifdef USE3D
     temp(7) = eps* r * rexp*(-zz*sn - xx*co)
     temp(8) = eps* ri * rexp*(-zz*xx*sn - xx*xx*co)   &
             - eps*(2./ln)* r * rexp*(-zz*xx*sn - xx*xx*co)    &
             - eps* r * rexp*co 
     temp(9) = eps* ri * rexp*(-zz*zz*sn - xx*zz*co)   &
             - eps*(2./ln)* r * rexp*(-zz*zz*sn - xx*zz*co)    &
             - eps* r * rexp*sn
     temp(10) = -eps* ri3 * rexp*(-zz*xx*xx*sn - xx*xx*xx*co)   &
             + eps*(2./ln)* r * rexp*(-(2.*zz*xx*xx/ln - zz)*sn - (2*xx*xx*xx/ln - 3.*xx)*co)    &
             + eps* ri * rexp*(-(zz - 4.*zz*xx*xx/ln)*sn - (3*xx-4*xx**3/ln)*co)               
          
     temp(11) = -eps* ri3 * rexp*(-zz*zz*xx*sn - zz*xx*xx*co)   &
             + eps*(2./ln)* r * rexp*(-(2.*zz*zz*xx/ln - xx)*sn - (2*zz*xx*xx/ln - zz)*co)    &
             + eps* ri * rexp*(-(xx - 4.*zz*zz*xx/ln)*sn - (zz-4*zz*xx**2/ln)*co)               
         
     temp(12) = -eps* ri3 * rexp*(-zz*zz*zz*sn - xx*zz*zz*co)   &
             + eps*(2./ln)* r * rexp*(-(2.*zz*zz*zz/ln - 3*zz)*sn - (2*xx*zz*zz/ln - xx)*co)    &
             + eps* ri * rexp*(-(3*zz - 4.*zz*zz*zz/ln)*sn - (xx-4*xx*zz*zz/ln)*co)   
#endif               
 
     call add_product(u1_l,fac,temp)

     end select
  end do

end subroutine random_per
!===========================
subroutine cartesian_to_cylindrical(x,vec)
  implicit none

  vectype, dimension(6), intent(inout) :: vec
  real, intent(in) :: x

  vec(6) = vec(6) * x
  vec(5) = vec(5) * x +    vec(3)
  vec(4) = vec(4) * x + 2.*vec(2)
  vec(3) = vec(3) * x
  vec(2) = vec(2) * x +    vec(1)
  vec(1) = vec(1) * x 
end subroutine cartesian_to_cylindrical
!===========================
subroutine cartesian_to_cylindrical_all()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: inode, numnodes, icounter_tt
  real :: x, phi, z

  numnodes = owned_nodes()

   do icounter_tt=1,numnodes
     inode = nodes_owned(icounter_tt)

     call get_node_pos(inode, x, phi, z)

     call get_local_vals(inode)

     call cartesian_to_cylindrical(x,psi0_l)
     call cartesian_to_cylindrical(x,psi1_l)
     call cartesian_to_cylindrical(x,  u0_l)
     call cartesian_to_cylindrical(x,  u1_l)
     
     if(numvar.ge.2) then
        call cartesian_to_cylindrical(x,bz0_l)
        call cartesian_to_cylindrical(x,bz1_l)
        call cartesian_to_cylindrical(x,vz0_l)
        call cartesian_to_cylindrical(x,vz1_l)
     endif
     
     call set_local_vals(inode)
  end do
end subroutine cartesian_to_cylindrical_all


!=========================================
! set_neo_vel
!=========================================
subroutine set_neo_vel
  use basic
  use arrays
  use m3dc1_nint
  use newvar_mod
  use neo
  use math
  use sparse
  use model
  use diagnostics

  implicit none

  integer :: i, j, itri, nelms, ier

  type(field_type) :: vz_vec, u_f, chi_f, diamag
  type(vector_type) :: vp_vec
  type(matrix_type) :: vpol_mat
  vectype, dimension(dofs_per_element, dofs_per_element, 2, 2) :: temp
  vectype, dimension(dofs_per_element) :: temp2, temp4
  vectype, dimension(dofs_per_element, 2) :: temp3

  real, dimension(MAX_PTS) :: theta, vtor, vpol, psival
  vectype, dimension(MAX_PTS) :: vz, vp, dia

  integer, dimension(dofs_per_element) :: imask_vor, imask_chi
  integer, dimension(MAX_PTS) :: iout
  integer :: imag

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Setting velocity from NEO data"

     if(db.ne.0. .and. ineo_subtract_diamag.eq.1) call create_field(diamag)
  call create_field(vz_vec)
  call create_vector(vp_vec, 2)
  call associate_field(u_f, vp_vec, 1)
  call associate_field(chi_f, vp_vec, 2)

  call set_matrix_index(vpol_mat, vpol_mat_index)
  call create_mat(vpol_mat, 2, 2, icomplex, 1)

  nelms = local_elements()
  do itri=1,nelms
        
     call define_element_quadrature(itri,int_pts_main,5)
     call define_fields(itri,0,1,0)
     call eval_ops(itri, psi_field(0), ps079)
     if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then 
        call eval_ops(itri, den_field(0), n079)
        call eval_ops(itri, p_field(0), p079)
        call eval_ops(itri, pe_field(0), pe079)
        pi079 = p079 - pe079

        ! zeff does not appear here because db includes zeff
        dia = db*(pi079(:,OP_DR)*ps079(:,OP_DR)+pi079(:,OP_DZ)*ps079(:,OP_DZ))&
             / n079(:,OP_1)
     end if

     theta = atan2(z_79-zmag,x_79-xmag)
     psival = -(ps079(:,OP_1) - psimin)
     call neo_eval_vel(int_pts_main, psival, theta, vpol, vtor, iout)
     vz = vtor / (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)
     vp = vpol / (b0_norm/sqrt(4.*pi*1.6726e-24*ion_mass*n0_norm)/l0_norm)

     ! NEO coordinates are (r, theta, phi_neo) = (r, phi, theta)
     ! --> theta = grad(r)xgrad(phi)
     ! if grad(psi).grad(r) < 0 then vpol is opposite direction to Bpol
     if(psimin.gt.psibound) vp = -vp
     ! phi_neo = -phi
     vz = -vz

     temp79e = sqrt((ps079(:,OP_DR)**2 + ps079(:,OP_DZ)**2)*ri2_79)

     do i=1, int_pts_main
        call magnetic_region(ps079(i,OP_1),ps079(i,OP_DR),ps079(i,OP_DZ), &
             x_79(i), z_79(i), imag)
        if(imag.ne.REGION_PLASMA) then
           vz(i) = 0.
           vp(i) = 0.
           iout(i) = 1
        endif
     end do

     where(iout.eq.1 .or. abs(temp79e).lt.1e-2)
        temp79f = 0.
        dia = 0.
     elsewhere
        temp79f = 1./temp79e
        dia = dia / &
             (ps079(:,OP_DR)*ps079(:,OP_DR) + ps079(:,OP_DZ)*ps079(:,OP_DZ))
     end where

     call get_vor_mask(itri, imask_vor)
     call get_chi_mask(itri, imask_chi)

     do i=1,dofs_per_element          
        ! assemble matrix
        do j=1, dofs_per_element

           ! vorticity equation
           if(imask_vor(i).eq.0) then
              temp(i,j,1,:) = 0.
           else
              temp(i,j,1,1) = &
                   int3(r2_79,mu79(i,:,OP_DR),nu79(j,:,OP_DR)) + &
                   int3(r2_79,mu79(i,:,OP_DZ),nu79(j,:,OP_DZ))
              temp(i,j,1,2) = &
                   int3(ri_79,mu79(i,:,OP_DR),nu79(j,:,OP_DZ)) - &
                   int3(ri_79,mu79(i,:,OP_DZ),nu79(j,:,OP_DR))
           end if
           
           ! compression equation
           if(imask_chi(i).eq.0) then
              temp(i,j,2,:) = 0.
           else
              temp(i,j,2,1) = &
                   int3(ri_79,mu79(i,:,OP_DZ),nu79(j,:,OP_DR)) - &
                   int3(ri_79,mu79(i,:,OP_DR),nu79(j,:,OP_DZ))
              temp(i,j,2,2) = &
                   int3(ri4_79,mu79(i,:,OP_DR),nu79(j,:,OP_DR)) + &
                   int3(ri4_79,mu79(i,:,OP_DZ),nu79(j,:,OP_DZ))
           end if
        end do

        ! assemble RHS
        ! toroidal rotation
        temp2(i) = int3(ri_79,mu79(i,:,OP_1),vz)

        ! vorticity
        if(imask_vor(i).eq.0) then
           temp3(i,1) = 0.
        else
           temp3(i,1) = &
                int4(temp79f,vp,mu79(i,:,OP_DR),ps079(:,OP_DR)) + &
                int4(temp79f,vp,mu79(i,:,OP_DZ),ps079(:,OP_DZ))
        endif

        ! compression
        if(imask_vor(i).eq.0) then
           temp3(i,2) = 0.
        else
           temp3(i,2) = &
                int5(ri3_79,temp79f,vp,mu79(i,:,OP_DZ),ps079(:,OP_DR)) - &
                int5(ri3_79,temp79f,vp,mu79(i,:,OP_DR),ps079(:,OP_DZ))
        endif

        ! diamagnetic term
        if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then
           temp4(i) = int2(mu79(i,:,OP_1), dia)
        end if
     end do

     call insert_block(vpol_mat, itri, 1, 1, temp(:,:,1,1), MAT_ADD)
     call insert_block(vpol_mat, itri, 1, 2, temp(:,:,1,2), MAT_ADD)
     call insert_block(vpol_mat, itri, 2, 1, temp(:,:,2,1), MAT_ADD)
     call insert_block(vpol_mat, itri, 2, 2, temp(:,:,2,2), MAT_ADD)
     call vector_insert_block(vp_vec, itri, 1, temp3(:,1), VEC_ADD)
     call vector_insert_block(vp_vec, itri, 2, temp3(:,2), VEC_ADD)
     call vector_insert_block(vz_vec%vec, itri, 1, temp2(:), VEC_ADD)

     if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then
        call vector_insert_block(diamag%vec, itri, 1, temp4(:), VEC_ADD)
     end if
  end do
  call sum_shared(vz_vec%vec)
  call sum_shared(vp_vec)
  if(db.ne.0. .and. ineo_subtract_diamag.eq.1) call sum_shared(diamag%vec)

  ! solve vz
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving vz..."
  call newsolve(mass_mat_lhs%mat,vz_vec%vec,ier)

  ! add neoclassical rotation to base rotation
  call add_field_to_field(vz_field(0), vz_vec)

  ! subtract diamagnetic term from rotation (because neo is adding this in)
  if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then
     ! solve diamagnetic term
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving diamag..."
     call newsolve(mass_mat_lhs%mat,diamag%vec,ier)
     call add_field_to_field(vz_field(0), diamag, -1.)
  end if

  ! solve vpol
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving vpol..."
  call boundary_vpol(vp_vec, u_f, chi_f, vpol_mat)
  call finalize(vpol_mat)
  call newsolve(vpol_mat, vp_vec, ier)
  u_field(0) = u_f
  chi_field(0) = chi_f

  call destroy_field(vz_vec)
  call destroy_vector(vp_vec)
  call destroy_mat(vpol_mat)
  if(db.ne.0. .and. ineo_subtract_diamag.eq.1) call destroy_field(diamag)

  if(myrank.eq.0 .and. iprint.ge.1) print *, " Done calculating velocity"
end subroutine set_neo_vel


!=====================================
subroutine initial_conditions()
  use basic
  use vector_mod
  use arrays

  use tilting_cylinder
  use taylor_reconnection
  use force_free_state
  use gem_reconnection
  use wave_propagation
  use gradshafranov
  use mri
  use grav
  use strauss
  use circular_field
  use rotate
  use eqdsk_eq
  use dskbal_eq
  use jsolver_eq
  use biharmonic
  use threed_wave_test
  use threed_diffusion_test
  use frs
  use ftz
  use eigen
  use neo
  use int_kink
  use rwm
  use solovev
  use auxiliary_fields
  use basicq
  use basicj
  use rmp
  use init_common
  use init_vmec
  use kprad_m3dc1
  use pellet
  use diagnostics
  use cylinder

  implicit none

  integer :: ierr

  if(iread_neo.eq.1) then
     call read_neo(ierr)
     if(ierr.ne.0) return
  end if

  if(iread_eqdsk.ge.1) then
     call eqdsk_init()
  else if(iread_dskbal.ge.1) then
     call dskbal_init()
  else if(iread_jsolver.ge.1) then
     call jsolver_init()
  else
     if(itor.eq.0) then
        ! slab equilibria
        select case(itaylor)
        case(0)
           call tilting_cylinder_init()
        case(1)
           call taylor_reconnection_init()
        case(2)
           call force_free_init()
        case(3)
           call gem_reconnection_init()
        case(4)
           call wave_init()
        case(5)
           call grav_init()
        case(6)
           call strauss_init()
        case(7)
           call circular_field_init()
        case(8)
           call biharmonic_init(1)
        case(9)
           call biharmonic_init(0)
        case(10,11)
           print *, "Circ_shell_only equilibrium has been deprecated"
           call safestop(1)
        case(12,13)
           print *, "Resistive wall test equilibrium has been deprecated"
           call safestop(1)
        case(14)
           call threed_wave_test_init()
        case(15)
           call threed_diffusion_test_init()
        case(16)
           call frs_init()
        case(17)
           call ftz_init()
        case(18)
           call eigen_init()
        case(19)
           call int_kink_init()
        case(20)
           call kstar_profiles()
        case(21,22,25,26,27,28,30,32,34)
           call fixed_q_profiles()
        case(23)
           call frs1_init()
        case(24)
           call rwm_init()
        case(29,31)
           call basicj_init()
        case(33)
           call cyl_init()
        end select
     else
        ! toroidal equilibria
        select case(itaylor)
        case(-1)
           bz_field(0) = bzero*rzero
           vz_field(0) = vzero
           den_field(0) = den0
           p_field(0) = p0
           pe_field(0) = p0*pefac
        case(0)
           call tilting_cylinder_init()
           call cartesian_to_cylindrical_all()
        case(1)
           call gradshafranov_init()
        case(2)
           call mri_init()
        case(3)
           call rotate_init()
        case(7)
           call circular_field_init()
           call cartesian_to_cylindrical_all()
        case(10,11)
           print *, "Circ_shell_only equilibrium has been deprecated"
           call safestop(1)
        case(12,13)
           print *, "Resistive_wall_test equilibrium has been deprecated"
           call safestop(1)
        case(14)
           call threed_wave_test_init()
        case(15)
           call threed_diffusion_test_init()
        case(16)
           call frs_init()
        case(17)
           call ftz_init()
        case(18)
           call eigen_init()
        case(19)
           call solovev_init()
        case(24)
           call rwm_init()
        case(29,31)
           call basicj_init()
#ifdef USEST
        case(40) ! Fixed boundary stellarator
           if (igeometry.eq.1 .and. iread_vmec.eq.1 .and. bloat_factor.eq.0) then
              call vmec_init()
           else
              if(myrank.eq.0) print *, &
                'VMEC equilibrium needs igeometry=1, iread_vmec=1, and bloat_factor=0!'
              call safestop(1)
           end if
#endif
        case(41) ! Free boundary stellarator or 3D fields
           if(iread_ext_field.eq.0) then  
              if(myrank.eq.0) print *, &
                "Invalid input: Free boundary stellarator needs external field."
              call safestop(1)
           end if
           if(type_ext_field.ge.1) &
                call load_stellarator_field
        end select
     end if
  end if

  if(iread_neo.eq.1) then
     call set_neo_vel
     call unload_neo
  end if

  if(ipellet.ne.0) then
     ! need to calculate norm for pellet_distribution
     call calculate_Lor_vol
  end if
  call den_eq
  call den_per
  call kprad_init_conds

#if defined(USEST) && defined(USEPARTICLES)
  if(kinetic .eq. 1) then
     call kinetic_eq
  endif
#endif
  if(irunaway .gt. 0) then
     call nre_eq
     call nre_per
   endif

  ! For RMP, 3D vacuum, and error fields
  if(irmp.ge.1 .or. iread_ext_field.ge.1 .or. &
       tf_tilt.ne.0. .or. tf_shift.ne.0. .or. &
       any(pf_tilt.ne.0.) .or. any(pf_shift.ne.0.)) then
     ! External fields already loaded for itaylor = 41
     if(itaylor.eq.41 .and. extsubtract.eq.0) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, &
           "Skipping: RMP specification not currently implemented for ST."
     else
        call rmp_per(1)
     end if
  end if
#ifdef USEST
  if(igeometry.eq.1 .and. iread_vmec.ge.1) then
     call destroy_vmec
  end if   
#endif

  ! calculate equilibrium and perturbed ne and temperature profiles
  call calculate_ne(0, den_field(0), ne_field(0), 1)
  call calculate_ne(1, den_field(1), ne_field(1), 1)
  call calculate_temperatures(0, te_field(0), ti_field(0), &
       pe_field(0), p_field(0), ne_field(0), den_field(0), &
       1)
  ! the first '0' below is intentional and not a typo! --YZ
  call calculate_temperatures(0, te_field(1), ti_field(1), &
       pe_field(1), p_field(1), ne_field(1), den_field(1), &
       1)

  if(iflip_b.eq.1) call mult(bz_field(0), -1.)
  if(iflip_j.gt.0) then 
     call mult(psi_field(0), -1.)
     if(icsubtract.eq.1) call mult(psi_coil_field, -1.)
     psilim = -psilim
     psibound = -psibound
  end if
  if(iflip_v.eq.1) call mult(vz_field(0), -1.)
  if(iflip_v.eq.-1) call mult(vz_field(0), 0.)
end subroutine initial_conditions
                                                                     
subroutine kstar_profiles()

  use basic
  use math
  use mesh_mod
  use sparse
  use arrays
  use m3dc1_nint
  use newvar_mod
  use boundary_conditions
  use model
  use gradshafranov
  use int_kink

  vectype, dimension (dofs_per_element,dofs_per_element) :: temp
  vectype, dimension (dofs_per_element) :: temp2
  vectype, dimension (MAX_PTS) :: co, sn, r, theta, rdpsidr
  real :: x, phi, z
  real, parameter :: e=2.7183
  real :: a, r1, r2, u, fa, ra
  real :: b0,r0
  integer :: m,n
  integer :: numnodes, nelms, l, itri, i, j, ier, icounter_tt
  integer :: imask(dofs_per_element)
  type (field_type) :: psi_vec
  type(matrix_type) :: psi_mat
  
  !call create_field(dpsi_dr)
  call create_field(psi_vec)
  
  call set_matrix_index(psi_mat, psi_mat_index)
  call create_mat(psi_mat,1,1,icomplex, 1)

  psi_vec = 0.

  !input variables: B0,fA,r1,r2,q0,R0,m,n,u,rA
  b0 = bzero
  r0 = rzero
  m = mpol
  n = igs
  a = alpha0
  r1 = alpha1
  r2 = alpha2
  u = alpha3
  fa = p1
  ra = p2
  
  if(myrank.eq.0 .and. iprint.ge.1) write(*,*)   &
                                    'b0,r0,m,n,e,a,r1,r2,u,fa,ra'   &
                                    ,b0,r0,m,n,e,a,r1,r2,u,fa,ra

  numnodes = owned_nodes()

  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l,x,phi,z)
     
     call constant_field(den0_l,1.)
     call constant_field(bz0_l,bzero)
     call constant_field(p0_l,p0)
     call constant_field(pe0_l,p0/2.)
     call int_kink_per(x,phi,z)
     
     call set_node_data(den_field(0),l,den0_l)
     call set_node_data(bz_field(0),l,bz0_l)
     call set_node_data(p_field(0),l,p0_l)
     call set_node_data(pe_field(0),l,pe0_l)
     call set_node_data(u_field(1),l,u1_l)
  enddo

  nelms = local_elements()
  do itri=1,nelms
     
     call define_element_quadrature(itri,int_pts_diag, int_pts_tor)
     call define_fields(itri,0,1,0) ! defines x_79,z_79,mu,nu
     !   call eval_ops(itri,dpsi_dr,dpt79)
     
     r = sqrt((x_79-xmag)**2 + (z_79-zmag)**2)
     theta = atan2(z_79-zmag,x_79-xmag)
     co = cos(theta)
     sn = sin(theta)

     rdpsidr = (B0*r**2)/ &
          ((1+fA/E**((-r1+r/a)**2/r2**2))*q0*R0* &
          (1+((r*Abs(-1+(m/(n*q0))**u)**(1/(2.*u)))/(a*rA))**(2*u))**(1/u))
      
     call get_boundary_mask(itri, BOUNDARY_DIRICHLET, imask, BOUND_DOMAIN)

     !  assemble matrix    
     do i=1,dofs_per_element
        if(imask(i).eq.0) then
           temp(i,:) = 0.
        else
           do j=1,dofs_per_element
              temp(i,j) = int4(mu79(i,:,OP_1),nu79(j,:,OP_DR),co,r) &
                   +      int4(mu79(i,:,OP_1),nu79(j,:,OP_DZ),sn,r)
           enddo
        end if
        !  assemble rhs
        temp2(i) = int2(mu79(i,:,OP_1),rdpsidr)
     enddo
     
     call insert_block(psi_mat, itri, 1,1, temp(:,:), MAT_ADD)
     
     call vector_insert_block(psi_vec%vec, itri, 1, temp2(:), MAT_ADD)
  enddo
  
  call sum_shared(psi_vec%vec)
  call flush(psi_mat)
  call boundary_gs(psi_vec%vec, 0., psi_mat)
  call finalize(psi_mat)

! solve for psi
  if(myrank.eq.0 .and. iprint.ge.1) print *, "solving psi"
 
  call newsolve(psi_mat,psi_vec%vec,ier)
  psi_field(0) = psi_vec
  
  call destroy_mat(psi_mat)
  call destroy_field(psi_vec)
  !call destroy_field(dpsi_dr)
  
end subroutine kstar_profiles

