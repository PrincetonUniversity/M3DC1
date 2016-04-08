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
  co = cos(phi)
  sn = sin(phi)

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
!==============================================================================
subroutine den_eq
  use basic
  use arrays
  use diagnostics
  use math
  use mesh_mod
  use m3dc1_nint
  use newvar_mod

  type(field_type) :: den_vec
  integer :: itri, numelms, i, def_fields
  vectype, dimension(dofs_per_element) :: dofs

  real :: k, kx
  
  if(idenfunc.eq.0) return

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
        k = 1./dendelt
        kx = k*(real((psi0_l(1)-psimin)/(psibound-psimin)) - denoff)
        
        n079(:,OP_1) = 1. + tanh(kx)
     end select

     do i=1, dofs_per_element
        dofs(i) = int2(mu79(:,OP_1,i),n079(:,OP_1))
     end do
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
  use field
  use m3dc1_nint
  use newvar_mod
  implicit none

  type(field_type) :: den_vec
  integer :: numelms, itri, i
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
             pellet_deposition(x_79, phi_79, z_79, p, n)
     end if

     do i=1, dofs_per_element
        dofs(i) = int2(mu79(:,OP_1,i),n179(:,OP_1))
     end do
     call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(den_vec%vec,mass_mat_lhs)
  den_field(1) = den_vec

  call destroy_field(den_vec)

end subroutine den_per

!=========================================================================
subroutine calculate_external_fields(sf)
  use basic
  use math
  use mesh_mod
  use sparse
  use arrays
  use coils
  use m3dc1_nint
  use newvar_mod
  use boundary_conditions
  use read_schaffer_field

  implicit none

  include 'mpif.h'

  type(schaffer_field), dimension(*) :: sf
  type(matrix_type) :: br_mat, bf_mat
  type(vector_type) :: psi_vec, bz_vec, p_vec
  integer :: i, j, itri, nelms, ier

  complex, dimension(int_pts_main) :: fr, fphi, fz
#ifdef USE3D
  real, dimension(MAX_PTS) :: gr, gphi, gz, p
#endif

  real, dimension(maxfilaments) :: xc_na, zc_na
  complex, dimension(maxfilaments) :: ic_na
  integer :: nc_na
  logical :: read_p

  vectype, dimension(dofs_per_element,dofs_per_element) :: temp, temp_bf
  vectype, dimension(dofs_per_element) :: temp2, temp3, temp4

  type(field_type) :: psi_f, bz_f, p_f

#ifdef USECOMPLEX
  complex :: sfac
#else
  real, dimension(int_pts_main) :: co, sn
#endif

  if(myrank.eq.0 .and. iprint.ge.2) print *, "Calculating error fields"

  if(irmp .ne. 0) then
     call load_coils(xc_na, zc_na, ic_na, nc_na, &
          'rmp_coil.dat', 'rmp_current.dat', ntor)
  end if

  call create_vector(p_vec,1)
  call associate_field(p_f,p_vec,1)

  call create_vector(psi_vec,1)
  call associate_field(psi_f,psi_vec,1)

  call create_vector(bz_vec,1)
  call associate_field(bz_f,bz_vec,1)

  call set_matrix_index(br_mat, br_mat_index)
  call create_mat(br_mat, 1, 1, icomplex, 1)
#ifdef CJ_MATRIX_DUMP
  print *, "create_mat coils br_mat", br_mat%imatrix 
#endif

  call set_matrix_index(bf_mat, bf_mat_index)
  call create_mat(bf_mat, 1, 1, icomplex, 0)
#ifdef CJ_MATRIX_DUMP
  print *, "create_mat coils br_mat", br_mat%imatrix 
#endif

  read_p = .false.
  if(iread_ext_field.ne.0) then
     do i=1, iread_ext_field 
        if(sf(i)%vmec) read_p = .true.
     end do
  end if

  if(myrank.eq.0 .and. iprint.ge.2) print *, 'calculating field values...'
  nelms = local_elements()
  do itri=1,nelms
        
     call define_element_quadrature(itri,int_pts_main,5)
     call define_fields(itri,0,1,0)

     fr   = 0.    ! B_R
     fphi = 0.    ! B_phi
     fz   = 0.    ! B_Z

     if(irmp.eq.0) then
        temp79a = 0.
        temp79b = 0.
        temp79c = 0.
        temp79d = 0.
     else
        do i=1, nc_na, 2
           call pane(ic_na(i),xc_na(i),xc_na(i+1),zc_na(i),zc_na(i+1), &
                npoints_pol,x_79,z_79,ntor,fr,fphi,fz)
        end do

#ifdef USECOMPLEX
        temp79a = fr
        temp79b = fphi
        temp79c = fz
        temp79d = 0.
#else
        do i=1, npoints_tor
           co(1:npoints_pol) = &
                cos(ntor*phi_79((i-1)*npoints_pol+1:i*npoints_pol))
           sn(1:npoints_pol) = &
                sin(ntor*phi_79((i-1)*npoints_pol+1:i*npoints_pol))
           temp79a((i-1)*npoints_pol+1:i*npoints_pol) = &
                real(fr(1:npoints_pol))*co(1:npoints_pol) + &
                aimag(fr(1:npoints_pol))*sn(1:npoints_pol)
           temp79b((i-1)*npoints_pol+1:i*npoints_pol) = &
                real(fphi(1:npoints_pol))*co(1:npoints_pol) + &
                aimag(fphi(1:npoints_pol))*sn(1:npoints_pol)
           temp79c((i-1)*npoints_pol+1:i*npoints_pol) = &
                real(fz(1:npoints_pol))*co(1:npoints_pol) + &
                aimag(fz(1:npoints_pol))*sn(1:npoints_pol)
        end do
#endif

        temp79a = -twopi*temp79a
        temp79b = -twopi*temp79b
        temp79c = -twopi*temp79c
     end if

     if(iread_ext_field.ne.0) then
        do i=1, iread_ext_field 
#if defined(USECOMPLEX)
           sfac = exp(-cmplx(0.,1.)*ntor*shift_ext_field(i)*pi/180.)
           call get_external_field_ft(sf(i),x_79,z_79,fr,fphi,fz,npoints)
           temp79a = temp79a + (1e4/b0_norm)*fr  *scale_ext_field*sfac
           temp79b = temp79b + (1e4/b0_norm)*fphi*scale_ext_field*sfac
           temp79c = temp79c + (1e4/b0_norm)*fz  *scale_ext_field*sfac
#elif defined(USE3D)
           call get_external_field(sf(i),&
                x_79,phi_79-shift_ext_field(i)*pi/180.,z_79,&
                gr,gphi,gz,p,npoints)
           temp79a = temp79a + (1e4/b0_norm)*gr  *scale_ext_field
           temp79b = temp79b + (1e4/b0_norm)*gphi*scale_ext_field
           temp79c = temp79c + (1e4/b0_norm)*gz  *scale_ext_field
           if(sf(i)%vmec) then
              temp79d = temp79d + p/p0_norm
           end if
#endif
        end do
     end if

     ! assemble matrix
     do i=1,dofs_per_element
        do j=1,dofs_per_element
           temp(i,j) = int3(ri2_79,mu79(:,OP_DR,i),nu79(:,OP_DR,j)) &
                +      int3(ri2_79,mu79(:,OP_DZ,i),nu79(:,OP_DZ,j)) &
                + regular*int3(ri4_79,mu79(:,OP_1,i),nu79(:,OP_1,j))
#if defined(USECOMPLEX) || defined(USE3D)
           temp_bf(i,j) = &
             + int3(ri_79,mu79(:,OP_DR,i),nu79(:,OP_DZP,j)) &
             - int3(ri_79,mu79(:,OP_DZ,i),nu79(:,OP_DRP,j))
#endif
        end do

        ! assemble RHS
        temp2(i) = &
             + int3(ri_79,mu79(:,OP_DR,i),temp79c) &
             - int3(ri_79,mu79(:,OP_DZ,i),temp79a)

        temp3(i) = int3(r_79,mu79(:,OP_1,i),temp79b)

        if(read_p) temp4(i) = int2(mu79(:,OP_1,i),temp79d)
     end do

     call insert_block(br_mat, itri, 1, 1, temp(:,:), MAT_ADD)
     call insert_block(bf_mat, itri, 1, 1, temp_bf(:,:), MAT_ADD)

     call vector_insert_block(psi_vec, itri, 1, temp2(:), MAT_ADD)
     call vector_insert_block(bz_vec, itri, 1, temp3(:), MAT_ADD)
     if(read_p) call vector_insert_block(p_vec, itri, 1, temp4(:), MAT_ADD)
  end do

  if(myrank.eq.0 .and. iprint.ge.2) print *, 'Finalizing...'
  call finalize(br_mat)
  call finalize(bf_mat)
  call sum_shared(psi_vec)
  call sum_shared(bz_vec)
  if(read_p) call sum_shared(p_vec)
 

  ! create external fields
  if(extsubtract.eq.1) then
     call create_field(psi_ext)
     call create_field(bz_ext)
     call create_field(bf_ext)
     use_external_fields = .true.
  end if

  !solve p
  if(read_p) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving p..."

     call newsolve(mass_mat_lhs%mat,p_vec,ier)
     p_field(1) = p_f
  end if

  ! solve bz
  if(numvar.ge.2) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving bz..."
     call newsolve(mass_mat_lhs%mat,bz_vec,ier)
     if(extsubtract.eq.1) then
        bz_ext = bz_f     
     else
        bz_field(1) = bz_f
     end if

#if defined(USECOMPLEX) || defined(USE3D)
     ! calculate f and add Grad_perp(f') to RHS
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving f..."
     if(extsubtract.eq.1) then
        bf_ext = 0.
        call solve_newvar1(bf_mat_lhs,bf_ext,mass_mat_rhs_bf, &
             bz_ext, bf_ext)
        call matvecmult(bf_mat, bf_ext, bz_vec)
     else
        bf_field(1) = 0.
        call solve_newvar1(bf_mat_lhs,bf_field(1),mass_mat_rhs_bf, &
             bz_field(1), bf_field(1))
        call matvecmult(bf_mat, bf_field(1), bz_vec)
     end if

     call add(psi_vec,bz_vec)
#endif
  end if

  ! solve -Del*(psi)/R^2 = Curl(B).Grad(phi)
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving psi..."
  call newsolve(br_mat,psi_vec, ier)
  if(extsubtract.eq.1) then
     psi_ext = psi_f
  else
     psi_field(1) = psi_f
  end if

  call destroy_vector(psi_vec)
  call destroy_vector(bz_vec)
  call destroy_vector(p_vec)
  call destroy_mat(br_mat)
  call destroy_mat(bf_mat)

  if(myrank.eq.0 .and. iprint.ge.2) print *, "Done calculating error fields"
end subroutine calculate_external_fields


!==============================================================================
subroutine rmp_per
  use basic
  use arrays
  use coils
  use boundary_conditions
  use read_schaffer_field

  implicit none

  logical :: is_boundary
  integer :: izone, izonedim, numnodes, l, ierr, icounter_tt
  real :: normal(2), curv, x, z, r2, dx, dz
  character(len=13) :: ext_field_name
  type(schaffer_field), allocatable :: sf(:)
#ifdef USECOMPLEX
  vectype :: ii
#endif

  if(irmp.eq.3) then
     numnodes = owned_nodes()
     do icounter_tt=1,numnodes
        l = nodes_owned(icounter_tt)

        call boundary_node(l,is_boundary,izone,izonedim,normal,curv,x,z)
        if(.not.is_boundary) cycle

        call get_local_vals(l)

        dx = x - xmag
        dz = z - zmag
        r2 = dx**2 + dz**2

        ! psi = 0.5*bx0*exp( i*2*theta)*exp(i*ntor*phi)
        !     = 0.5*bx0*(cos(2*theta) + i*sin(2*theta))*exp(i*ntor*phi)

        ! cos(2*theta)
        psi1_l(1) = (dx**2 - dz**2)/r2
        psi1_l(2) = 4.*dx*dz**2/r2**2
        psi1_l(3) =-4.*dz*dx**2/r2**2
        psi1_l(4) = 4.*dz**2*(1.-4.*dx**2/r2)/r2**2
        psi1_l(5) = 8.*dx*dz*(dx**2 - dz**2)/r2**3
        psi1_l(6) =-4.*dx**2*(1.-4.*dz**2/r2)/r2**2

#ifdef USECOMPLEX
        if(ntor.lt.0) then
           ii = -(0,1.)
        else
           ii =  (0,1.)
        endif
        ! sin(2*theta)
        psi1_l(1) = psi1_l(1) + ii*2.*dx*dz/r2
        psi1_l(2) = psi1_l(2) + ii*2.*dz*(1. - 2.*dx**2/r2)/r2
        psi1_l(3) = psi1_l(3) + ii*2.*dx*(1. - 2.*dz**2/r2)/r2
        psi1_l(4) = psi1_l(4) + ii*4.*dx*dz*(1. - 4.*dz**2/r2)/r2**2
        psi1_l(5) = psi1_l(5) - ii*2.*(1. - 8.*dx**2*dz**2/r2**2)/r2
        psi1_l(6) = psi1_l(6) + ii*4.*dx*dz*(1. - 4.*dx**2/r2)/r2**2
#endif
        psi1_l = psi1_l*0.5*bx0

        call set_local_vals(l)
     end do
     return
  endif

  ! load external field data from schaffer file
  if(iread_ext_field.ge.1) then
     allocate(sf(iread_ext_field))
     do l=1, iread_ext_field
        if(iread_ext_field.eq.1) then
           ext_field_name = 'error_field'
        else
           write(ext_field_name, '("error_field",I2.2)') l
        end if
        call load_schaffer_field(sf(l),ext_field_name,isample_ext_field, &
             isample_ext_field_pol,ierr)
        if(ierr.ne.0) call safestop(6)

#ifdef USECOMPLEX
        if(myrank.eq.0 .and. iprint.gt.2) then
           print *, 'Calculating field FT'
        end if
        call calculate_external_field_ft(sf(l), ntor)
#endif
     end do
  end if

  ! calculate external fields from non-axisymmetric coils and external field
  call calculate_external_fields(sf)

  ! unload data
  if(iread_ext_field.ge.1) then
     do l=1, iread_ext_field
        call unload_schaffer_field(sf(l))
     end do
     deallocate(sf)
  end if

  ! leave perturbation only on the boundary
  if(irmp.eq.2) then
     numnodes = owned_nodes()
     do icounter_tt=1,numnodes
        l = nodes_owned(icounter_tt)
        call boundary_node(l,is_boundary,izone,izonedim,normal,curv,x,z)
        if(.not.is_boundary) then
           call get_local_vals(l)
           psi1_l = 0.
           call set_local_vals(l)
        endif
     end do
  endif
end subroutine rmp_per

!==============================================================================
! Tilting Cylinder (itaylor = 0)
!==============================================================================
module tilting_cylinder

  implicit none

  real, parameter, private :: k = 3.8317059702

contains

subroutine tilting_cylinder_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, alx, alz

  call get_bounding_box_size(alx, alz)

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

     call get_local_vals(l)

     call cylinder_equ(x, z)
     call cylinder_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine tilting_cylinder_init

subroutine cylinder_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: rr, ri, arg, befo, ff, fp, fpp, j0, j1, kb
  real :: dbesj0, dbesj1
  
  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  !.....Use the equilibrium given in R. Richard, et al,
  !     Phys. Fluids B, 2 (3) 1990 p 489
  !.....(also Strauss and Longcope)
  !
  rr = sqrt(x**2 + z**2)
  ri = 0.
  if(rr.ne.0) ri = 1./rr
  arg = k*rr
  if(rr.le.1) then
     befo = 2./dbesj0(dble(k))
     j0 = dbesj0(dble(arg))
     j1 = dbesj1(dble(arg))
     ff = .5*befo
     fp = 0.
     fpp = -.125*k**2*befo
     if(arg.ne.0)  then
        ff   = befo*j1/arg
        fp  = befo*(j0 - 2.*j1/arg)/rr
        fpp = befo*(-3*j0 + 6.*j1/arg - arg*j1)/rr**2
     endif
  else
     ff   = (1.-1./rr**2)
     fp  = 2./rr**3
     fpp = -6/rr**4
  endif
   
  psi0_l(1) = ff* z
  psi0_l(2) = fp*z*x *ri
  psi0_l(3) = fp*z**2*ri + ff
  psi0_l(4) = fpp*x**2*z*ri**2 + fp*(  z*ri - x**2*z*ri**3)
  psi0_l(5) = fpp*z**2*x*ri**2 + fp*(  x*ri - z**2*x*ri**3)
  psi0_l(6) = fpp*z**3  *ri**2 + fp*(3*z*ri - z**3  *ri**3)

  if(rr.gt.1) then
     call constant_field(bz0_l, bzero   )
     call constant_field(pe0_l, p0-pi0)
     call constant_field( p0_l, p0)
  else 
     if(numvar.ge.2) then
        kb = k**2*(1.-beta)
        bz0_l(1) = sqrt(kb*psi0_l(1)**2+bzero**2)
        bz0_l(2) = kb/bz0_l(1)*psi0_l(1)*psi0_l(2)
        bz0_l(3) = kb/bz0_l(1)*psi0_l(1)*psi0_l(3)
        bz0_l(4) = kb/bz0_l(1)                              &
             *(-kb/bz0_l(1)**2*(psi0_l(1)*psi0_l(2))**2     &
             + psi0_l(2)**2+psi0_l(1)*psi0_l(4))
        bz0_l(5) = kb/bz0_l(1)*(-kb/bz0_l(1)**2*            &
             psi0_l(1)**2*psi0_l(2)*psi0_l(3)               &
             + psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
        bz0_l(6) = kb/bz0_l(1)                              &
             *(-kb/bz0_l(1)**2*(psi0_l(1)*psi0_l(3))**2     &
             + psi0_l(3)**2+psi0_l(1)*psi0_l(6))
     else 
       call constant_field(bz0_l, bzero)
     end if

     if(numvar.ge.3) then
        kb = k**2*beta*pefac
        pe0_l(1) = 0.5*kb*psi0_l(1)**2 + p0 - pi0
        pe0_l(2) = kb*psi0_l(1)*psi0_l(2)
        pe0_l(3) = kb*psi0_l(1)*psi0_l(3)
        pe0_l(4) = kb*(psi0_l(2)**2+psi0_l(1)*psi0_l(4))
        pe0_l(5) = kb*(psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
        pe0_l(6) = kb*(psi0_l(3)**2+psi0_l(1)*psi0_l(6))
        
        kb = k**2*beta
        p0_l(1) = 0.5*kb*psi0_l(1)**2 + p0
        p0_l(2) = kb*psi0_l(1)*psi0_l(2)
        p0_l(3) = kb*psi0_l(1)*psi0_l(3)
        p0_l(4) = kb*(psi0_l(2)**2+psi0_l(1)*psi0_l(4))
        p0_l(5) = kb*(psi0_l(2)*psi0_l(3)+psi0_l(1)*psi0_l(5))
        p0_l(6) = kb*(psi0_l(3)**2+psi0_l(1)*psi0_l(6))
     else
        call constant_field(pe0_l, p0-pi0)
        call constant_field( p0_l, p0)
     end if
  endif

  call constant_field(den0_l, 1.)

end subroutine cylinder_equ

subroutine cylinder_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: befo, uzero, czero
  
  befo = eps*exp(-(x**2+z**2))
  uzero = befo*cos(z)
  czero = befo*sin(z)
  u1_l(1) = uzero
  u1_l(2) = - 2.*x*uzero
  u1_l(3) = - 2.*z*uzero - czero
  u1_l(4) = (-2.+4.*x**2)*uzero
  u1_l(5) = 4.*x*z*uzero + 2.*x*czero
  u1_l(6) = (-3.+4.*z**2)*uzero + 4.*z*czero
  
  vz1_l =  0.

  chi1_l(1) = czero
  chi1_l(2) = -2*x*czero
  chi1_l(3) = -2*z*czero + uzero
  chi1_l(4) = (-2.+4*x**2)*czero
  chi1_l(5) = 4.*x*z*czero - 2.*x*uzero
  chi1_l(6) = (-3.+4.*z**2)*czero - 4*z*uzero

  psi1_l =  0.
  bz1_l = 0.
  pe1_l = 0.
  p1_l = 0.
  den1_l = 0.

end subroutine cylinder_per

end module tilting_cylinder



!==============================================================================
! Taylor Reconnection (itaylor = 1)
!==============================================================================

module taylor_reconnection

contains


subroutine taylor_reconnection_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, alx, alz

  call get_bounding_box_size(alx, alz)

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

     call get_local_vals(l)

     call taylor_reconnection_equ(x, z)
     call taylor_reconnection_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine taylor_reconnection_init



subroutine taylor_reconnection_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  psi0_l(1) = -0.5*z**2
  psi0_l(2) =  0.
  psi0_l(3) = -z
  psi0_l(4) =  0.
  psi0_l(5) =  0.
  psi0_l(6) = -1.

  call constant_field( bz0_l, bzero)
  call constant_field( pe0_l, p0-pi0)
  call constant_field(den0_l, 1.)
  
end subroutine taylor_reconnection_equ

subroutine taylor_reconnection_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.
  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.
  p1_l = 0.
  den1_l = 0.

end subroutine taylor_reconnection_per


end module taylor_reconnection


!==============================================================================
! Force-Free Taylor State (itaylor = 2)
!==============================================================================
module force_free_state

contains

subroutine force_free_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, alx, alz

  call get_bounding_box_size(alx, alz)

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)

     call get_node_pos(l, x, phi, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

     call get_local_vals(l)

     call force_free_equ(x, z)
     call force_free_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine force_free_init

subroutine force_free_equ(x, z)
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  real, intent(in) :: x, z
  
  real :: kx, kz, alx, alz, alam
  
  call get_bounding_box_size(alx, alz)

  kx = pi/alx
  kz = pi/alz

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  alam = sqrt(kx**2+kz**2)

  psi0_l(1) = sin(kx*x)*sin(kz*z)
  psi0_l(2) = kx*cos(kx*x)*sin(kz*z)
  psi0_l(3) = kz*sin(kx*x)*cos(kz*z)
  psi0_l(4) = -kx**2*sin(kx*x)*sin(kz*z)
  psi0_l(5) =  kx*kz*cos(kx*x)*cos(kz*z)
  psi0_l(6) = -kz**2*sin(kz*x)*sin(kz*z)

  bz0_l = alam*psi0_l

  call constant_field( pe0_l, p0-pi0)
  call constant_field(  p0_l, p0)
  call constant_field(den0_l, 1.)
  
end subroutine force_free_equ

subroutine force_free_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.

  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.

  p1_l = 0.
  den1_l = 0.

end subroutine force_free_per

end module force_free_state

!===========================================================================
! GEM Reconnection (itaylor = 3)
!===========================================================================
module gem_reconnection

  real, private :: akx, akz

contains

subroutine gem_reconnection_init()
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, alx, alz

  call get_bounding_box_size(alx, alz)

  akx = twopi/alx
  akz = pi/alz

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     call get_local_vals(l)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

     call gem_reconnection_equ(x, z)
     call gem_reconnection_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine gem_reconnection_init

subroutine gem_reconnection_equ(x, z)
  use basic
  use arrays
  use math

  implicit none

  real, intent(in) :: x, z

  call constant_field(  u0_l, 0.)
  if(numvar.ge.2) call constant_field( vz0_l, 0.)
  if(numvar.ge.3) call constant_field(chi0_l, 0.)

  psi0_l(1) = 0.5*alog(cosh(2.*z))
  psi0_l(2) = 0.0
  psi0_l(3) = tanh(2.*z)  
  psi0_l(4) = 0.0
  psi0_l(5) = 0.0
  psi0_l(6) = 2.*sech(2.*z)**2

  ! if numvar = 2, then use Bz to satisfy force balance
  if(numvar.eq.2) then
     bz0_l(1) = sqrt(bzero**2 + sech(2.*z)**2)
     bz0_l(2) = 0.
     bz0_l(3) = -2.*tanh(2.*z)*sech(2.*z)**2/bz0_l(1)
     bz0_l(4) = 0.
     bz0_l(5) = 0.
     bz0_l(6) = (2.*sech(2.*z)**2/bz0_l(1))                      &
          *(bz0_l(3)*tanh(2.*z)/bz0_l(1)                          &
          - 2.*sech(2.*z)**2 + 4.*tanh(2.*z)**2)
  endif

  ! if numvar >= 3, then use pressure to satisfy force balance
  if(numvar.ge.3) then

     if(bzero.eq.0) then
        bz0_l(1) = sqrt(1.-2.*p0)*sech(2.*z)
        bz0_l(2) = 0.
        bz0_l(3) = -2.*bz0_l(1)*tanh(2.*z)
        bz0_l(4) = 0.
        bz0_l(5) = 0.
        bz0_l(6) = -2.*bz0_l(3)*tanh(2.*z) &
             -4.*bz0_l(1)*(1.-tanh(2.*z)**2)
     else
        bz0_l(1) = sqrt(bzero**2 + (1.-2.*p0)*sech(2.*z)**2)
        bz0_l(2) = 0.
        bz0_l(3) = -(1.-2.*p0)*2.*tanh(2.*z)*sech(2.*z)**2/bz0_l(1)
        bz0_l(4) = 0.
        bz0_l(5) = 0.
        bz0_l(6) = (1.-2.*p0)*(2.*sech(2.*z)**2/bz0_l(1))            &
             *(bz0_l(3)*tanh(2.*z)/bz0_l(1)                          &
             - 2.*sech(2.*z)**2 + 4.*tanh(2.*z)**2)
     endif
  endif

  p0_l(1) = p0*(sech(2.*z)**2 + 0.2)
  p0_l(2) = 0.
  p0_l(3) = p0*(-4.*sech(2.*z)**2*tanh(2.*z))
  p0_l(4) = 0.
  p0_l(5) = 0.
  p0_l(6) = p0*(-8.*sech(2.*z)**2*(sech(2.*z)**2-2.*tanh(2.*z)**2)) 

  pe0_l = p0_l*pefac

  if(idens.eq.1) then
     den0_l(1) = sech(2*z)**2 + 0.2
     den0_l(2) = 0.
     den0_l(3) = -4.*sech(2*z)**2*tanh(2*z)
     den0_l(4) = 0.
     den0_l(5) = 0.
     den0_l(6) = -8.*sech(2*z)**2*(sech(2*z)**2-2.*tanh(2*z)**2)
  else
     call constant_field(den0_l, 1.)
  endif

end subroutine gem_reconnection_equ

subroutine gem_reconnection_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.
  
  psi1_l(1) =  eps*cos(akx*x)*cos(akz*z)
  psi1_l(2) = -eps*sin(akx*x)*cos(akz*z)*akx
  psi1_l(3) = -eps*cos(akx*x)*sin(akz*z)*akz
  psi1_l(4) = -eps*cos(akx*x)*cos(akz*z)*akx**2
  psi1_l(5) =  eps*sin(akx*x)*sin(akz*z)*akx*akz
  psi1_l(6) = -eps*cos(akx*x)*cos(akz*z)*akz**2

  bz1_l = 0.
  pe1_l = 0. 
  p1_l = 0.
  den1_l = 0.

end subroutine gem_reconnection_per

end module gem_reconnection


!=============================================================================
! Wave Propagation (itaylor = 4)
!=============================================================================
module wave_propagation

  implicit none
  real, private :: alx, alz, akx, akx2, omega
  real, private :: psiper, phiper, bzper, vzper, peper, chiper, nper, pper

contains

subroutine wave_init()
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z
  real :: b2,a2
  real :: kp,km,t1,t2,t3
  real :: coef(4)
  real :: root(3)
  real :: error(3)
  real :: bi

  call get_bounding_box_size(alx, alz)

  ! for itaylorw=3, set up a phi perturbation only
  if(iwave.eq.3) then
     phiper = eps
     vzper = 0.
     chiper = 0.
     psiper = 0.
     bzper = 0.
     nper = 0.
     pper = 0.
     peper = 0.
     goto 1
  endif

  akx = twopi/alx      
  akx2 = akx**2
  b2 = bzero*bzero + bx0*bx0
  a2 = bx0**2/b2
  bi = 2.*pi0*gyro
  
  ! numvar=2 =================
  if(numvar.eq.2) then
     kp = akx * (b2 + bi*(3.*a2-1.)/4.)
     km = akx * (b2 - bi*(3.*a2-1.)/4.)

     if(iwave.eq.0) then ! fast wave
        omega = (akx/2.)*sqrt(a2/b2)*(sqrt(4.*b2**2+km**2)+kp)
        psiper = eps
        phiper = -eps*2.*b2 / (sqrt(4.*b2**2+km**2)+km)
        bzper = akx*psiper
        vzper = akx*phiper
        
     else ! slow wave
        omega = (akx/2.)*sqrt(a2/b2)*(sqrt(4.*b2**2+km**2)-kp)
        psiper = eps
        phiper = -eps*2.*b2 / (sqrt(4.*b2**2+km**2)-km)
        bzper = -akx*psiper
        vzper = -akx*phiper
     endif
     chiper = 0.
     nper = 0.
     peper = 0.
     pper = 0.
     
  ! numvar=3 =================
  elseif(numvar.ge.3) then
     
     coef(4) = 96.*b2**4
     coef(3) = -6.*akx2*b2**3*(16.*b2**2*(1.+a2+akx2*a2)            &
          + akx2*bi**2*(1.+6.*a2-3.*a2**2)                          &
          + 16.*gam*p0*b2)
     coef(2) = 6.*a2*akx2**2*b2**3                                  &
          *(b2*((4.*b2-bi*akx2*(1.+a2))**2                          &
          +4.*bi**2*akx2*(1.-a2)*(1.+akx2*a2))                      &
          + gam*p0*(32.*b2**2+16.*akx2*b2**2                        &
          +bi**2*akx2*(1.-3.*a2)**2))
     coef(1) = -6.*gam*p0*a2**2*akx2**3*b2**4                       &
          *(4.*b2+akx2*bi*(1.-3.*a2))**2

     if(myrank.eq.0 .and. iprint.ge.1) then
        write(*,*) "Coefs: ", coef(1), coef(2), coef(3), coef(4)
     endif

     call cubic_roots(coef, root, error)

     if(myrank.eq.0 .and. iprint.ge.1) then
        write(*,*) "Coefs: ", coef(1), coef(2), coef(3), coef(4)
        write(*,*) "Roots: ", root(1), root(2), root(3)
        write(*,*) "Error: ", error(1), error(2), error(3)
     endif

     select case(iwave)
     case(1)
        omega=sqrt(root(1))
     case(2)
        omega=sqrt(root(2))
     case default 
        omega=sqrt(root(3))
     end select

     t1 = akx2*bi/(4.*b2)
     t2 = (akx2/omega**2)*                                          &
          (bzero**2/(1.-gam*p0*akx2/(omega**2)))
     t3 = bx0 * akx / omega

     psiper = eps
     phiper = psiper*(akx2*t3*(1.+(t3**2/(1.-t2-t3**2)))-(1.-t2)/t3)&
          / (1. - t2 + t1*t2*(1.+3.*a2)                             &
          + t1*t3**2*(2.*t2-(1.-3.*a2))/(1-t2-t3**2))
     vzper = phiper*t1*t3*(2.*t2-(1.-3.*a2))/(1-t2-t3**2)           &
          - psiper*akx2*t3**2/(1.-t2-t3**2)
     bzper = (-phiper*t1*t2*(1.+3.*a2) - vzper*t3 + psiper*akx2*t3) &
          / (1.-t2)
     chiper = bzero / (1.-gam*p0*akx2/(omega**2)) / omega           &
          * ((1.+3.*a2)*t1*phiper - bzper)
     peper = -chiper*gam*(p0-pi0)*akx2 / omega
     nper = -chiper*akx2/omega
     pper = -chiper*gam*p0*akx2 / omega
  endif

  if(myrank.eq.0) then
     print *, "Wave angular frequency: ", omega
     print *, "Wave phase velocity: ", omega/akx
     print *, "Wave transit time: ", 2.*pi/omega
  end if

1 continue

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     x = x - alx*.5 - xzero
     z = z - alz*.5 - zzero

     call get_local_vals(l)

     call wave_equ(x, z)
     call wave_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine wave_init

subroutine wave_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.
  
  psi0_l(1) = z*bx0
  psi0_l(2) = 0.0
  psi0_l(3) = bx0
  psi0_l(4) = 0.0
  psi0_l(5) = 0.0
  psi0_l(6) = 0.0

  call constant_field( bz0_l, bzero)
  call constant_field( pe0_l, p0-pi0)
  call constant_field(  p0_l, p0)
  call constant_field(den0_l, 1.)

end subroutine wave_equ


subroutine wave_per(x, z)
  use math
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  call plane_wave(u1_l, x, z,  akx, 0., phiper, pi/2.)
  call plane_wave(vz1_l, x, z, akx, 0., vzper, pi/2.)
  call plane_wave(chi1_l, x, z, akx, 0., chiper, pi)
  
  call plane_wave(psi1_l, x, z, akx, 0., psiper, pi/2.)
  call plane_wave( bz1_l, x, z, akx, 0., bzper, pi/2.)

  call plane_wave(pe1_l, x, z, akx, 0., peper, pi/2.)
  call plane_wave( p1_l, x, z, akx, 0.,  pper, pi/2.)

  if(idens.eq.1) then
     call plane_wave(den1_l, x, z, akx, 0., nper, pi/2.)    
  endif

end subroutine wave_per

end module wave_propagation



!============================================================================
! Gravitational Instability Equilibrium (itor = 0, itaylor = 5)
!============================================================================
module grav

contains

subroutine grav_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_local_vals(l)

     call get_node_pos(l, x, phi, z)

     call grav_equ(x, z)
     call grav_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine grav_init

subroutine grav_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: fac1, n0, pn0, ppn0

  n0 = exp(z/ln)
  pn0 = n0/ln
  ppn0 = pn0/ln

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  psi0_l = 0.

  fac1 = gravz*ln+gam*p0
  bz0_l(1) = sqrt(bzero**2 - 2.*fac1*(n0-1.))
  bz0_l(2) = 0.
  bz0_l(3) = -pn0*fac1/bz0_l(1)
  bz0_l(4) = 0.
  bz0_l(5) = 0.
  bz0_l(6) = (fac1/bz0_l(1))*(pn0*bz0_l(3)/bz0_l(1) - ppn0)

  
  fac1 = p0 - pi0
  pe0_l(1) = fac1+gam*fac1*(n0-1.)
  pe0_l(2) = 0.
  pe0_l(3) = gam*fac1*pn0
  pe0_l(4) = 0.
  pe0_l(5) = 0.
  pe0_l(6) = gam*fac1*ppn0

  den0_l(1) = n0
  den0_l(2) = 0.
  den0_l(3) = pn0
  den0_l(4) = 0.
  den0_l(5) = 0.
  den0_l(6) = ppn0

  fac1 = p0
  p0_l(1) = fac1+gam*fac1*(n0-1.)
  p0_l(2) = 0.
  p0_l(3) = gam*fac1*pn0
  p0_l(4) = 0.
  p0_l(5) = 0.
  p0_l(6) = gam*fac1*ppn0

end subroutine grav_equ


subroutine grav_per(x, z)
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  real, intent(in) :: x, z

  real :: kx, kz, alx, alz

  call get_bounding_box_size(alx, alz)
  kx = twopi/alx
  kz = pi/alz

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.

  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.
  p1_l = 0.

  den1_l(1) =  eps*sin(kx*(x-xzero))*sin(kz*(z-zzero))
  den1_l(2) =  eps*cos(kx*(x-xzero))*sin(kz*(z-zzero))*kx
  den1_l(3) =  eps*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kz
  den1_l(4) = -eps*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kx**2
  den1_l(5) =  eps*cos(kx*(x-xzero))*cos(kz*(z-zzero))*kx*kz
  den1_l(6) = -eps*sin(kx*(x-xzero))*sin(kz*(z-zzero))*kz**2

end subroutine grav_per

end module grav


!==============================================================================
! Strauss Equilibrium (itor = 0, itaylor = 6)
!
! H.R. Strauss, Phys. Fluids 19(1):134 (1976)
!==============================================================================
module strauss

  implicit none

  real, private :: alx,alz

contains

subroutine strauss_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z

  call get_bounding_box_size(alx, alz)

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     x = x - alx/2.
     z = z - alz/2.

     call get_local_vals(l)

     call strauss_equ(x, z)
     call strauss_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine strauss_init

subroutine strauss_equ(x, z)
  use math
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: kx, kz, a0

    u0_l = 0.
  if(numvar.ge.2) vz0_l = 0.
  if(numvar.ge.3) chi0_l = 0.

  kx = pi/alx
  kz = pi/alz

  a0 = 2.*alz*alx*bzero/(pi*ln*q0)

  psi0_l(1) = -a0*cos(kx*x)*cos(kz*z)
  psi0_l(2) =  a0*sin(kx*x)*cos(kz*z)*kx
  psi0_l(3) =  a0*cos(kx*x)*sin(kz*z)*kz
  psi0_l(4) =  a0*cos(kx*x)*cos(kz*z)*kx**2
  psi0_l(5) = -a0*sin(kx*x)*sin(kz*z)*kx*kz
  psi0_l(6) =  a0*cos(kx*x)*cos(kz*z)*kz**2

  call constant_field(bz0_l, bzero)
  call constant_field(pe0_l, p0 - pi0)
  call constant_field(den0_l, 1.)
  call constant_field(p0_l, p0)

end subroutine strauss_equ


subroutine strauss_per(x, z)
  use math
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: akx, akz

  akx = pi/alx
  akz = twopi/alz

  psi1_l(1) = -eps*cos(akx*x)*sin(akz*z)
  psi1_l(2) =  eps*sin(akx*x)*sin(akz*z)*akx
  psi1_l(3) = -eps*cos(akx*x)*cos(akz*z)*akz
  psi1_l(4) =  eps*cos(akx*x)*sin(akz*z)*akx**2
  psi1_l(5) =  eps*sin(akx*x)*cos(akz*z)*akx*akz
  psi1_l(6) =  eps*cos(akx*x)*sin(akz*z)*akz**2
  u1_l = 0.
  
  bz1_l = 0.
  vz1_l = 0.
  pe1_l = 0.
  chi1_l = 0.
  den1_l = 0.
  p1_l = 0.

end subroutine strauss_per

end module strauss


!==============================================================================
! Circular magnetic field equilibrium (for parallel conduction tests)
!==============================================================================
module circular_field

  real, private :: alx, alz, x0, z0

contains

subroutine circular_field_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z

  call get_bounding_box_size(alx, alz)

  x0 = alx/4.
  z0 = 0.

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)
     x = x - xmag
     z = z - zmag

     call get_local_vals(l)

     call circular_field_equ(x, z)
     call circular_field_per(x, phi, z)

     call set_local_vals(l)
  enddo

  call finalize(field0_vec)
  call finalize(field_vec)

end subroutine circular_field_init

subroutine circular_field_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: j0, r2, fac, l2, jf
  

  j0 = 2.*bzero/q0
  l2 = ln**2
  r2 = (x**2 + z**2)/l2
  fac = exp(-r2)
  jf = j0*fac

  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  psi0_l(1) = jf*l2/4.
  psi0_l(2) = -(2./l2)*psi0_l(1)*x
  psi0_l(3) = -(2./l2)*psi0_l(1)*z
  psi0_l(4) = -(2./l2)*(psi0_l(1) + psi0_l(2)*x)
  psi0_l(5) =  (2./l2)**2 * psi0_l(1)*x*z
  psi0_l(6) = -(2./l2)*(psi0_l(1) + psi0_l(3)*z)

  bz0_l(1) = sqrt(bzero**2 - (j0**2*l2/8.)*(1. - (1.-2.*r2)*fac**2))
  bz0_l(2) = -jf**2*x*(1.-r2)/(2.*bz0_l(1))
  bz0_l(3) = -jf**2*z*(1.-r2)/(2.*bz0_l(1))
  bz0_l(4) = jf**2/(2.*bz0_l(1)) * &
       ((x*bz0_l(2)/bz0_l(1) - 1.) * (1.-r2) + 2.*(x**2/l2)*(3.-2.*r2))
  bz0_l(5) = -jf**2*x*z/bz0_l(1) * &
       (((1.-r2)*jf/(2.*bz0_l(1)))**2 - (3.-2.*r2)/l2)
  bz0_l(6) = jf**2/(2.*bz0_l(1)) * &
       ((z*bz0_l(3)/bz0_l(1) - 1.) * (1.-r2) + 2.*(z**2/l2)*(3.-2.*r2))

  call constant_field(p0_l, p0)
  call constant_field(pe0_l, p0-pi0)
  call constant_field(den0_l, 1.)

end subroutine circular_field_equ


subroutine circular_field_per(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z
  real :: r2, l2, fac

  l2 = ln**2
  r2 = (x**2 + z**2)/l2
  fac = eps*exp(-r2)

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.

  psi1_l = 0.
  bz1_l = 0.
  p1_l = 0.

  psi1_l(1) = fac* sin(z/ln)
  psi1_l(2) = fac* sin(z/ln)*(-2.*x/l2)
  psi1_l(3) = fac*(sin(z/ln)*(-2.*z/l2) + cos(z/ln)/ln)
  psi1_l(4) = fac* sin(z/ln)*(-2.  /l2 + (2.*x/l2)**2)
  psi1_l(5) = fac*(sin(z/ln)*(-2.*z/l2) + cos(z/ln)/ln)*(-2.*x/l2)
  psi1_l(6) = fac*(sin(z/ln)*(-2.*z/l2) + cos(z/ln)/ln)*(-2.*z/l2) &
       - fac*(2.*cos(z/ln)*z/ln + 3.*sin(z/ln))/l2

#ifdef USE3D
  psi1_l(7:12) = -ntor*psi1_l(1:6)*sin(ntor*phi)
  psi1_l(1:6) = psi1_l(1:6)*cos(ntor*phi)
#endif

!!$  p1_l(1) = eps*exp(-((x-x0)**2+z**2)/(2.*ln**2))
!!$  p1_l(2) = -(x-x0)*p1_l(1)/ln**2
!!$  p1_l(3) = -(z-z0)*p1_l(1)/ln**2
!!$  p1_l(4) = (((x-x0)/ln)**2 - 1.)*p1_l(1)/ln**2
!!$  p1_l(5) =  (x-x0)*(z-z0)*p1_l(1)/ln**4
!!$  p1_l(6) = (((z-z0)/ln)**2 - 1.)*p1_l(1)/ln**2

  ! for viscosity test..
!  vz1_l = p1_l
!  p1_l = 0.

  ! for parallel viscosity test...

!!$  u1_l(1) = eps*exp(-(x**2+z**2)/(2.*ss**2))
!!$  u1_l(2) = -x*u1_l(1)/ss**2
!!$  u1_l(3) = -z*u1_l(1)/ss**2
!!$  u1_l(4) = ((x/ss)**2 - 1.)*u1_l(1)/ss**2
!!$  u1_l(5) =  x*z*u1_l(1)/ss**4
!!$  u1_l(6) = ((z/ss)**2 - 1.)*u1_l(1)/ss**2
!!$  vz1_l(1) = vzero*exp(-(x**2+z**2)/(2.*ss**2))
!!$  vz1_l(2) = -x*vz1_l(1)/ss**2
!!$  vz1_l(3) = -z*vz1_l(1)/ss**2
!!$  vz1_l(4) = ((x/ss)**2 - 1.)*vz1_l(1)/ss**2
!!$  vz1_l(5) =  x*z*psi0_l(1)/ss**4
!!$  vz1_l(6) = ((z/ss)**2 - 1.)*vz1_l(1)/ss**2
!!$  p1_l = 0.

  pe1_l = pefac*p1_l

  den1_l = 0.

end subroutine circular_field_per

end module circular_field


!==============================================================================
! Magnetorotational Equilibrium (itor = 1, itaylor = 2)
!==============================================================================
module mri

  real, private :: kx, kz

contains

subroutine mri_init()
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, alx, alz

  call get_bounding_box_size(alx, alz)

  kx = pi/alx
  kz = twopi/alz

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     z = z - alz*.5

     call get_local_vals(l)

     call mri_equ(x, z)
     call mri_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine mri_init

subroutine mri_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: nu, fac1, fac2

  nu = db*0.5*gyro*pi0/bzero
  fac1 = sqrt(gravr)
  fac2 = (9./128.)*nu**2/fac1

  call constant_field(  u0_l,0.)

  psi0_l(1) = bzero * x**2 / 2.
  psi0_l(2) = bzero * x
  psi0_l(3) = 0.
  psi0_l(4) = bzero
  psi0_l(5) = 0.
  psi0_l(6) = 0.

  if(numvar.ge.2) then
     bz0_l = 0.

     vz0_l(1)  = fac1*sqrt(x) + (3./8.)*nu - fac2/sqrt(x)
     vz0_l(2)  = 0.5*fac1/sqrt(x) + 0.5*fac2/(sqrt(x)**3)
     vz0_l(3)  = 0.
     vz0_l(4)  = -0.25*fac1/(sqrt(x)**3) - 0.75*fac2/(sqrt(x)**5)
     vz0_l(5) = 0.
     vz0_l(6) = 0.
  end if

  if(numvar.ge.3) then
     call constant_field( pe0_l,p0 - pi0)
     call constant_field(chi0_l,0.)
  end if

  call constant_field(den0_l, 1.)
  call constant_field(p0_l,p0)

end subroutine mri_equ


subroutine mri_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  real :: fac1

  u1_l = 0.
  if(numvar.ge.2)  then
     fac1 = eps*sqrt(gravr*x)
     vz1_l(1) =  fac1*sin(kx*(x-xzero))*sin(kz*(z-zzero))
     vz1_l(2) =  fac1*cos(kx*(x-xzero))*sin(kz*(z-zzero))*kx
     vz1_l(3) =  fac1*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kz
     vz1_l(4) = -fac1*sin(kx*(x-xzero))*cos(kz*(z-zzero))*kx**2
     vz1_l(5) =  fac1*cos(kx*(x-xzero))*cos(kz*(z-zzero))*kx*kz
     vz1_l(6) = -fac1*sin(kx*(x-xzero))*sin(kz*(z-zzero))*kz**2
  endif
  chi1_l = 0.
  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.

  den1_l = 0.
  p1_l = 0.

end subroutine mri_per

end module mri


!==============================================================================
! Rotating cylinder
! ~~~~~~~~~~~~~~~~~
!
! This is a rotating equilibrium with a radially increasing density to offset
! the centrifugal force.  In the isothermal case, with do density diffusion or
! thermal diffusion, this is a static equilibrium 
! (with or without gyroviscosity).
!==============================================================================
module rotate

  real, private :: kx, kz

contains

subroutine rotate_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, alx, alz

  call get_bounding_box_size(alx, alz)

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     z = z - alz*.5

     call get_local_vals(l)

     call rotate_equ(x, z)
     call rotate_per(x, z)

     call set_local_vals(l)
  enddo

end subroutine rotate_init

subroutine rotate_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z


  u0_l = 0.
  vz0_l(1) =    vzero*x**2
  vz0_l(2) = 2.*vzero*x
  vz0_l(3) = 0.
  vz0_l(4) = 2.*vzero
  vz0_l(5) = 0.
  vz0_l(6) = 0.
  chi0_l = 0.

  psi0_l = 0.
  bz0_l(1) =    bzero*x**2
  bz0_l(2) = 2.*bzero*x
  bz0_l(3) = 0.
  bz0_l(4) = 2.*bzero
  bz0_l(5) = 0.
  bz0_l(6) = 0.
  pe0_l(1) = p0 - pi0
  pe0_l(2:6) = 0.
  p0_l(1) = p0
  p0_l(2:6) = 0.
  den0_l(1) =  2.*(bzero/vzero)**2
  den0_l(2:6) = 0.
  

end subroutine rotate_equ


subroutine rotate_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.
  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.
  den1_l = 0.
  p1_l = 0.

end subroutine rotate_per

end module rotate


!==============================================================================
! Eqdsk_eq
! ~~~~~~~~
!
! Loads eqdsk equilibrium
!==============================================================================
module eqdsk_eq

contains

subroutine eqdsk_init()
  use math
  use basic
  use arrays
  use eqdsk
  use gradshafranov
  use newvar_mod
  use sparse
  use diagnostics
  use mesh_mod
  use m3dc1_nint

  implicit none

  integer :: l, ll, ierr, itri, k, numelms, i
  real :: x, phi, z , dpsi, ffp2, pp2
  vectype, parameter ::  negone = -1
  vectype, dimension(dofs_per_element) :: dofs
  type(field_type) :: psi_vec, bz_vec, den_vec, p_vec

  real, allocatable :: flux(:), nflux(:)

!!$  numnodes = owned_nodes()

  if(myrank.eq.0 .and. iprint.gt.0) print *, "before load_eqdsk", iread_eqdsk
  call load_eqdsk(ierr)
  if(ierr.ne.0) call safestop(1)

  press = press*amu0
  pprime = pprime*amu0
  current = current*amu0

  if(myrank.eq.0 .and. iprint.ge.1) then
     print *, 'normalized current ', current
  end if

  ! create spline containing q profile as function of psi_N
  if( .not. isnan(qpsi(nw))) then
    allocate(nflux(nw))
    do l=1, nw
       nflux(l) = l / (nw-1.)
    end do
    call create_spline(q_spline,nw,nflux,qpsi)
    deallocate(nflux)
  endif

  if(current.lt.0 .and. iflip_j.eq.1) then
     if(myrank.eq.0) then 
        print *, 'WARNING: iflip_j = 1 for case with negative current.'
        print *, 'Default behavior has been changed to take negative current'
        print *, 'into account automatically.  To override this, set'
        print *, 'iflip_j = 2.'
     end if
     call safestop(1)
  end if

  if(iflip_z.eq.1) zmaxis = -zmaxis

  tcuro = current
  xmag = rmaxis
  zmag = zmaxis
  if(xmag0.eq.0.) then
     xmag0 = rmaxis
     zmag0 = zmaxis
  end if
  rzero = rmaxis

  if(ifixedb.eq.0) then 
     if(iread_eqdsk.eq.3 .and. ifixedb.eq.0) then
        igs_calculate_ip_fields = .true.
     end if
     if(iread_eqdsk.eq.3 .or. idevice.eq.-1 .or. imulti_region.eq.1) then
        igs_calculate_pf_fields = .true.
     end if
     call vacuum_field
  end if

  if(iread_eqdsk.eq.3) then 
     
 ! define initial field associated with delta-function or gaussian source
  !     corresponding to current tcuro at location (xmag,zmag)
     if(myrank.eq.0 .and. iprint.gt.0) write(*,1001) xmag,zmag,tcuro,sigma0
1001 format(' in gradshafranov_init',/,'   xmag,zmag,tcuro,sigma0',1p4e12.4)
     if(sigma0 .eq.0) then
        call deltafun(xmag,zmag,tcuro,jphi_field)
     else
        call gaussianfun(xmag,zmag,tcuro,sigma0,jphi_field)
     endif
   
  else
     
    if(myrank.eq.0 .and. iprint.ge.1) &
         print *, "Interpolating geqdsk Equilibrium"

    call create_field(psi_vec)
    call create_field(bz_vec)
    call create_field(p_vec)
    call create_field(den_vec)

    numelms = local_elements()

    do k=0,1
       psi_vec = 0.
       bz_vec = 0.
       p_vec = 0.
       den_vec = 0.
       
       do itri=1,numelms
          call define_element_quadrature(itri,int_pts_main,int_pts_tor)
          call define_fields(itri,0,1,0)

          if(k.eq.0) then 
             ! calculate equilibrium fields
             call eqdsk_equ
          else
             ! calculate perturbed fields
             call eqdsk_per
          end if

          ! populate vectors for solves

          ! psi
          do i=1, dofs_per_element
             dofs(i) = int2(mu79(:,OP_1,i),ps079(:,OP_1))
          end do
          call vector_insert_block(psi_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! bz
          do i=1, dofs_per_element
             dofs(i) = int2(mu79(:,OP_1,i),bz079(:,OP_1))
          end do
          call vector_insert_block(bz_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! p
          do i=1, dofs_per_element
             dofs(i) = int2(mu79(:,OP_1,i),p079(:,OP_1))
          end do
          call vector_insert_block(p_vec%vec,itri,1,dofs,VEC_ADD)
          
          ! den
          do i=1, dofs_per_element
             dofs(i) = int2(mu79(:,OP_1,i),n079(:,OP_1))
          end do
          call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
       end do

       ! do solves
       call newvar_solve(psi_vec%vec,mass_mat_lhs)
       psi_field(k) = psi_vec
       
       call newvar_solve(bz_vec%vec,mass_mat_lhs)
       bz_field(k) = bz_vec
       
       call newvar_solve(p_vec%vec,mass_mat_lhs)
       p_field(k) = p_vec
       pe_field(k) = p_vec
       call mult(pe_field(k), pefac)
       
       call newvar_solve(den_vec%vec,mass_mat_lhs)
       den_field(k) = den_vec
    end do

    call destroy_field(psi_vec)
    call destroy_field(bz_vec)
    call destroy_field(den_vec)
    call destroy_field(p_vec)

  end if
!
! Bateman scaling parameter reintroduced
  if(igs_pp_ffp_rescale.ne.1) fpol(nw) = fpol(nw)*batemanscale
!
  bzero = fpol(nw)/rzero
  if(iprint.ge.1 .and. myrank.eq.0) then 
     write(*,'(A,2F12.4)') 'Setting bzero, rzero = ', bzero, rzero
  end if

  if(igs.gt.0) then
     if(iread_eqdsk.eq.2) then
        call default_profiles
     else
        allocate(flux(nw))
        dpsi = (sibry-simag)/(nw-1.)
        
        do l=1,nw
           flux(l) = (l-1)*dpsi
           ll = nw - l
           if(batemanscale.eq.1.0 .or. igs_pp_ffp_rescale.eq.1) cycle
! ...Apply Bateman scaling --- redefine fpol keeping ffprim fixed
           if(ll.gt.0) fpol(ll) = sign(1.0,fpol(nw)) &
                *sqrt(fpol(ll+1)**2 - dpsi*(ffprim(ll)+ffprim(ll+1)))
        end do
        call create_profile(nw,press,pprime,fpol,ffprim,flux)
        
        if(.not. isnan(qpsi(nw))) then
           call create_rho_from_q(nw,flux,qpsi)
        else
           if(myrank.eq.0) print *, "toroidal flux not defined"
           flux=0
        endif
        if(myrank.eq.0 .and. iprint.ge.1) then
           open(unit=77,file="debug-out",status="unknown")
           write(77,2010) sibry,simag,tcuro,xmag,zmag
2010       format("sibry,simag,tcuro,xmag,zmag =",1p5e12.4,/,  &
                "  l   press       pprime      fpol        ffprim      ffp2        pp2         flux        qpsi")
           do l=1,nw
              if(l.gt.1 .and. l.lt.nw)  then
                ffp2 = fpol(l)*(fpol(l+1)-fpol(l-1))/(2*dpsi)
                pp2 = (press(l+1)-press(l-1))/(2*dpsi)
              else
                ffp2 = 0
                pp2 = 0
              endif
              write(77,2011) l,press(l),pprime(l),fpol(l),ffprim(l),ffp2,pp2,flux(l),qpsi(l)
           enddo
2011       format(i3,1p8e12.4)
           close(77)
        endif
        deallocate(flux)
     end if

     psibound = sibry
     psimin = simag

     call gradshafranov_solve
     call gradshafranov_per
  else
     psibound = sibry
     psimin = simag
  endif

  call unload_eqdsk

  ! flip psi sign convention
  ! (iread_eqdsk==3 does not use the eqdsk psi)
  if(iread_eqdsk.eq.1 .or. iread_eqdsk.eq.2) then
     call mult(psi_field(0), negone)
     call mult(psi_field(1), negone)
     if(icsubtract.eq.1) call mult(psi_coil_field, negone)
     psibound = -psibound
     psimin = -psimin
     psilim = -psilim
  end if

  if(iprint.ge.1 .and. myrank.eq.0) &
       write(*,2012) sibry,simag,psimin,psilim,psibound
2012 format(" sibry, simag, psimin, psilim,psibound =",1p5e12.4)

end subroutine eqdsk_init

subroutine eqdsk_equ()
  use basic
  use arrays
  use m3dc1_nint

  use eqdsk

  implicit none

  real :: p, q
  integer :: i, j, n, m, k

  real :: dx, dz, temp, dpsi
  real, dimension(4,4) :: a
  real, dimension(4) :: b, c
  
  dx = rdim/(nw - 1.)
  dz = zdim/(nh - 1.)

  ps079(:,OP_1) = 0.
  p079(:,OP_1) = 0.
  n079(:,OP_1) = 0.
  bz079(:,OP_1) = 0.

  do k=1, npoints
     p = (x_79(k)-rleft)/dx + 1.
     q = (z_79(k)-zmid)/dz + nh/2. + .5
     i = p
     j = q

     call bicubic_interpolation_coeffs(psirz,nw,nh,i,j,a)

     do n=1, 4
        do m=1, 4
           temp = a(n,m)
           if(n.gt.1) temp = temp*(p-i)**(n-1)
           if(m.gt.1) temp = temp*(q-j)**(m-1)
           ps079(k,OP_1) = ps079(k,OP_1) + temp
        end do
     end do

     ! calculation of p
     ! ~~~~~~~~~~~~~~~~
     dpsi = (sibry - simag)/(nw - 1.)
     p = (ps079(k,OP_1) - simag)/dpsi + 1.
     i = p

     if(i.gt.nw) then
        p079(k,OP_1) = press(nw)
        bz079(k,OP_1) = fpol(nw)
     else
        ! use press and fpol to calculate values of p and I
        call cubic_interpolation_coeffs(press,nw,i,b)
        call cubic_interpolation_coeffs(fpol,nw,i,c)

        do n=1,4
           temp = b(n)
           if(n.gt.1) temp = temp*(p-i)**(n-1)
           p079(k,OP_1) = p079(k,OP_1) + temp
           
           temp = c(n)
           if(n.gt.1) temp = temp*(p-i)**(n-1)
           bz079(k,OP_1) = bz079(k,OP_1) + temp
        end do
     endif
  end do

  if(pedge.ge.0.) p079(:,OP_1) = p079(:,OP_1) + pedge

  where(real(p079(:,OP_1)).lt.0.) p079(:,OP_1) = 0.

  ! Set density
  if(expn.eq.0.) then
     n079(:,OP_1) = 1.
  else
     n079(:,OP_1) = (p079(:,OP_1)/p0)**expn
  end if
end subroutine eqdsk_equ

subroutine eqdsk_per
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  p079(:,OP_1) = 0.
  ps079(:,OP_1) = 0.
  bz079(:,OP_1) = 0.
  n079(:,OP_1) = 0.

end subroutine eqdsk_per
  
end module eqdsk_eq


!==============================================================================
! Dskbal_eq
! ~~~~~~~~~
!
! Loads dskbal equilibrium
!==============================================================================
module dskbal_eq

  implicit none

contains

subroutine dskbal_init()
  use math
  use basic
  use arrays
  use dskbal
  use gradshafranov
  use newvar_mod
  use sparse
  use diagnostics

  implicit none

  integer :: i, inode, numnodes, icounter_tt
  real :: dp, a(4), minden, minte, minti

  print *, "dskbal_init called"

  call load_dskbal
  p_bal = p_bal*amu0
  pprime_bal = pprime_bal*amu0

  ! find "vacuum" values
  minden = 0.
  minte = 0.
  minti = 0.
  do i=1,npsi_bal
     if(ne_bal(i) .le. 0.) exit
     minden = ne_bal(i)
     minte = te_bal(i)
     minti = ti_bal(i)
  end do

  ! replace zeroed-out values with "vacuum" values
  do i=1,npsi_bal
     if(ne_bal(i) .le. 0.) ne_bal(i) = minden
     if(ti_bal(i) .le. 0.) ti_bal(i) = minti
     if(te_bal(i) .le. 0.) te_bal(i) = minte
     if(p_bal(i) .le. 0.) p_bal(i) = minden*(minte + minti)*1.6022e-12*amu0/10.
  end do

  xmag = x_bal(1,1)
  zmag = z_bal(1,1)
  rzero = xmag
  bzero = f_bal(npsi_bal)/rzero

  ifixedb = 1

  if(iread_dskbal.eq.2) then
     call default_profiles
  else
     call create_profile(npsi_bal,p_bal,pprime_bal,f_bal,ffprime_bal,psi_bal)
  end if
  
  ! initial plasma current filament
  call deltafun(xmag,zmag,tcuro,jphi_field)

  call gradshafranov_solve
  call gradshafranov_per  

  ! set density profile
  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     inode = nodes_owned(icounter_tt)
     call get_local_vals(inode)
     do i=1, npsi_bal-1
        if((psi_bal(i+1)-psi_bal(1))/(psi_bal(npsi_bal)-psi_bal(1)) &
             .gt. (real(psi0_l(1))-psimin)/(psilim-psimin)) exit
     end do

     call cubic_interpolation_coeffs(ne_bal,npsi_bal,i,a)

     dp = ((real(psi0_l(1))-psimin) &
          *(psi_bal(npsi_bal) - psi_bal(1))/(psilim - psimin) &
          -(psi_bal(i)-psi_bal(1))) / (psi_bal(i+1) - psi_bal(i))

     den0_l(1) = a(1) + a(2)*dp + a(3)*dp**2 + a(4)*dp**3
     den0_l(2) = (a(2) + 2.*a(3)*dp + 3.*a(4)*dp**2)*psi0_l(2)
     den0_l(3) = (a(2) + 2.*a(3)*dp + 3.*a(4)*dp**2)*psi0_l(3)
     den0_l(4) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(4)
     den0_l(5) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(5)
     den0_l(6) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(6)
     den0_l = den0_l / n0_norm
     call set_local_vals(inode)
  enddo

  call unload_dskbal

end subroutine dskbal_init
  
end module dskbal_eq


!==============================================================================
! Jsolver_eq
! ~~~~~~~~~~
!
! Loads jsolver equilibrium
!==============================================================================
module jsolver_eq

contains

subroutine jsolver_init()
  use basic
  use arrays
  use gradshafranov
  use newvar_mod
  use sparse
  use coils
  use diagnostics
  use jsolver
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, pzero_jsv, gzero_jsv
  real, allocatable :: ffprime(:),ppxx_jsv2(:),gpx_jsv2(:)

  if(myrank.eq.0 .and. iprint .ge. 1) print *, "jsolver_init called"

  call load_jsolver

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     call get_local_vals(l)

     call jsolver_equ(x, z)
     call jsolver_per(x, z)

     call set_local_vals(l)
  enddo

  xmag = x_jsv(1,1)
  zmag = z_jsv(1,1)

if(myrank.eq.0 .and. iprint .ge. 1) print *, "xmag,zmag=",xmag,zmag

  ! remove factor of 1/xzero_jsv in definition of g
  gxx_jsv = gxx_jsv*xzero_jsv
  gpx_jsv = gpx_jsv*xzero_jsv
  pzero_jsv = 1.5*p_jsv(1)   - .5*p_jsv(2)
  gzero_jsv = 1.5*gxx_jsv(1) - .5*gxx_jsv(2)
!
!.calculate ppxx_jsv and gpx_jsv another way for comparison
  allocate(ppxx_jsv2(npsi_jsv),gpx_jsv2(npsi_jsv))
  do l=2,npsi_jsv-1
     ppxx_jsv2(l) =  (p_jsv(l) - p_jsv(l-1)) /(psival_jsv(l) - psival_jsv(l-1))
     gpx_jsv2(l) = (gxx_jsv(l) - gxx_jsv(l-1))/(psival_jsv(l) - psival_jsv(l-1))
  enddo
! shift gxx and p to be defined at psi_j locations
! instead of at psi_{j+1/2} locations
  do l=npsi_jsv,2,-1
     p_jsv(l) = (p_jsv(l-1) + p_jsv(l))/2.
     gxx_jsv(l) = (gxx_jsv(l-1) + gxx_jsv(l))/2.
  end do
  p_jsv(1)   = pzero_jsv
  gxx_jsv(1) = gzero_jsv

!...Adopt JSOLVER convention that R*B_T = xzero_jsv
  rzero = xzero_jsv
  bzero = 1.0
  if(myrank.eq.0 .and. iprint .ge. 1) print *, "rzero,bzero=",rzero,bzero

  if(igs.gt.0) then
     if(iread_jsolver.eq.2) then
        call default_profiles
     else
        allocate(ffprime(npsi_jsv))
        ffprime = gpx_jsv*gxx_jsv
        call create_profile(npsi_jsv,p_jsv(1:npsi_jsv),ppxx_jsv(1:npsi_jsv), &
             gxx_jsv(1:npsi_jsv),ffprime,psival_jsv(1:npsi_jsv))
!
      if(myrank.eq.0 .and. iprint.ge.1) then
        open(unit=77,file="debug-out",status="unknown")
        write(77,2010) xmag,zmag,tcuro
  2010 format("xmag,zmag,tcuro =",1p3e12.4,/,  &
     "i   press       pprime      fpol        ffprim      flux")
        do l=1,npsi_jsv
        write(77,2011) l,p_jsv(l),ppxx_jsv(l),gxx_jsv(l),ffprime(l),psival_jsv(l)
        enddo
  2011 format(i3,1p5e12.4)
       close(77)
      endif
!
        deallocate(ffprime)
     end if
     if(sigma0.eq.0) then
        call deltafun(xmag,zmag,tcuro,jphi_field)
     else
        call gaussianfun(xmag,zmag,tcuro,sigma0,jphi_field)
     endif
     psibound = psival_jsv(npsi_jsv)
     psimin = psival_jsv(1)

     call gradshafranov_solve
     call gradshafranov_per
  else
     psibound = psival_jsv(npsi_jsv)
     psimin = psival_jsv(1)
  endif

  call unload_jsolver

end subroutine jsolver_init

subroutine jsolver_equ(x, z)
  use basic
  use arrays

  use jsolver

  implicit none

  real, intent(in) :: x, z

  integer :: i, j
  real :: i0, j0, d, dmin, dj
  real, dimension(4) :: a

  ! determine jsolver index (i0,j0) of this node
  i0 = 3
  j0 = 1
  dmin = (x - x_jsv(3,1))**2 + (z - z_jsv(3,1))**2
  do i=3, nthe+2
     do j=1, npsi_jsv
        d = (x - x_jsv(i,j))**2 + (z - z_jsv(i,j))**2
        if(d.lt.dmin) then
           dmin = d
           i0 = i
           j0 = j
        end if
     end do
  end do
 
  j = j0
  dj = j - j0
  call cubic_interpolation_coeffs(psival_jsv,npsi_jsv,j,a)
  psi0_l(1) = a(1) + a(2)*dj + a(3)*dj**2 + a(4)*dj**3

  call cubic_interpolation_coeffs(gxx_jsv,npsi_jsv,j,a)
  bz0_l(1) = a(1) + a(2)*dj + a(3)*dj**2 + a(4)*dj**3

  call cubic_interpolation_coeffs(p_jsv,npsi_jsv,j,a)
  p0_l(1) = a(1) + a(2)*dj + a(3)*dj**2 + a(4)*dj**3
  
  u0_l = 0.
  vz0_l = 0.
  chi0_l = 0.

  if(pedge.ge.0.) p0_l = p0_l + pedge

  ! Set electron pressure and density
  pe0_l = pefac*p0_l


end subroutine jsolver_equ


subroutine jsolver_per(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z

  u1_l = 0.
  vz1_l = 0.
  chi1_l = 0.
  psi1_l = 0.
  bz1_l = 0.
  pe1_l = 0.
  den1_l = 0.
  p1_l = 0.

end subroutine jsolver_per
  
end module jsolver_eq


!==============================================================================
! 3D wave Test
! ~~~~~~~~~~~~
! This is a test case which initializes a 2-dimensional equilibrium
! With a 3-dimensional initial perturbation
!==============================================================================
module threed_wave_test

  real, private :: kx, kz, kphi, omega, k2, kp2, kB, b2

contains

subroutine threed_wave_test_init()
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, x1, x2, z1, z2

  call get_bounding_box(x1, z1, x2, z2)
  kx = pi/(x2-x1)
  kz = pi/(z2-z1)
  kphi = ntor/rzero
  kp2 = kx**2 + kz**2
  k2 = kp2 + kphi**2
  kB = kx*Bx0 + kphi*bzero
  b2 = bzero**2 + Bx0**2

  select case(numvar)
  case(1)
     select case(iwave)
     case(0)
        omega = kB

     case default
        if(myrank.eq.0) then
           print *, 'Error: iwave must be one of the following:'
           print *, ' iwave = 0: Shear Alfven wave'
        end if
        call safestop(4)
     end select
     
  case(2)
     select case(iwave)
     case(0)
        if(myrank.eq.0) print *, 'Shear Alfven wave'
        
        omega = bx0*kx*kB + 0.5*b2*kphi**2 &
             + 0.5*sqrt((2.*bx0*kx*kB + B2*kphi**2)**2 &
             -4.*(bx0*kx*kB)**2*k2/kp2)
        omega = sqrt(omega)

     case(1)
        if(myrank.eq.0) print *, 'Fast magnetosonic wave'

        omega = bx0*kx*kB + 0.5*b2*kphi**2 &
             - 0.5*sqrt((2.*bx0*kx*kB + B2*kphi**2)**2 &
             -4.*(bx0*kx*kB)**2*k2/kp2)
        omega = sqrt(omega)

     case default
        if(myrank.eq.0) then
           print *, 'Error: iwave must be one of the following:'
           print *, ' iwave = 0: Shear Alfven wave'
           print *, ' iwave = 1: Fast wave'
        end if
        call safestop(4)
     end select

  case(3)
     if(bx0 .ne. 0.) then
        if(myrank.eq.0) print *, 'Only Bx0=0 supported for numvar=3'
        call safestop(4)
     endif

     select case(iwave)
     case(0)
        if(myrank.eq.0) print *, 'Shear Alfven wave'
        omega = bzero*kphi

     case(1)
        if(myrank.eq.0) print *, 'Fast magnetosonic wave'
        omega = 0.5*k2*(bzero**2 + gam*p0) &
             + 0.5*sqrt((k2*(bzero**2 + gam*p0))**2 &
             -4.*(bzero*kphi)**2*k2*gam*p0)
        omega = sqrt(omega)

     case(2)
        if(myrank.eq.0) print *, 'Slow magnetosonic wave'
        omega = 0.5*k2*(bzero**2 + gam*p0) &
             - 0.5*sqrt((k2*(bzero**2 + gam*p0))**2 &
             -4.*(bzero*kphi)**2*k2*gam*p0)
        omega = sqrt(omega)

     case default
        if(myrank.eq.0) then
           print *, 'Error: iwave must be one of the following:'
           print *, ' iwave = 0: Shear Alfven wave'
           print *, ' iwave = 1: Fast wave'
           print *, ' iwave = 2: Slow wave'
        end if
        call safestop(4)
     end select
  end select

  if(myrank.eq.0) then
     print *, ' kphi = ', myrank, kphi
     print *, ' k2 = ', myrank, k2
     print *, ' kp2 = ', myrank, kp2
     print *, ' omega = ', myrank, omega
     print *, ' rzero = ', myrank, rzero
     print *, ' bzero = ', myrank, bzero
     print *, ' p0 = ', myrank, p0
  end if

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     x = x - x1
     z = z - z1

     call get_local_vals(l)

     call threed_wave_test_equ(x, phi, z)
     call threed_wave_test_per(x, phi, z)

     call set_local_vals(l)
  enddo
end subroutine threed_wave_test_init

subroutine threed_wave_test_equ(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z

  psi0_l = 0.

  call constant_field(bz0_l, bzero)
  call constant_field(den0_l, 1.)
  call constant_field(p0_l, p0)
  call constant_field(pe0_l, p0-pi0)

  psi0_l(:) = 0.
  psi0_l(1) = -Bx0*z
  psi0_l(3) = -Bx0

end subroutine threed_wave_test_equ


subroutine threed_wave_test_per(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z
  vectype, dimension(dofs_per_node) :: val

  val(1) = eps*sin(kx*x)*sin(kz*z)*cos(kphi*phi)
  val(2) = eps*cos(kx*x)*sin(kz*z)*cos(kphi*phi)*kx
  val(3) = eps*sin(kx*x)*cos(kz*z)*cos(kphi*phi)*kz
  val(4) =-eps*sin(kx*x)*sin(kz*z)*cos(kphi*phi)*kx**2
  val(5) = eps*cos(kx*x)*cos(kz*z)*cos(kphi*phi)*kx*kz
  val(6) =-eps*sin(kx*x)*sin(kz*z)*cos(kphi*phi)*kz**2
#ifdef USE3D
  val(7)  =-eps*sin(kx*x)*sin(kz*z)*sin(kphi*phi)*kphi
  val(8)  =-eps*cos(kx*x)*sin(kz*z)*sin(kphi*phi)*kphi*kx
  val(9)  =-eps*sin(kx*x)*cos(kz*z)*sin(kphi*phi)*kphi*kz
  val(10) = eps*sin(kx*x)*sin(kz*z)*sin(kphi*phi)*kphi*kx**2
  val(11) =-eps*cos(kx*x)*cos(kz*z)*sin(kphi*phi)*kphi*kx*kz
  val(12) = eps*sin(kx*x)*sin(kz*z)*sin(kphi*phi)*kphi*kz**2
#endif

  select case(numvar)
  case(1)
     psi1_l = val
     u1_l = -val*kB/omega

  case(2)
     psi1_l = val
     u1_l = -val*kB/omega
     vz1_l = -val*(0,1)*Bx0*kz*kphi*kp2 &
          / (k2*(Bx0*kx)**2 - kp2*omega**2)
     bz1_l = (-kp2)* (val*(0,1)*Bx0**2*kx*kz*kphi) &
          / (k2*(Bx0*kx)**2 - kp2*omega**2)
     
  case(3)
     if(iwave.eq.0) then
        psi1_l = val
        u1_l = -val*bzero*kphi/omega
     else
        bz1_l =-val*kp2
        vz1_l = val*bzero*kphi*kp2*k2*p0*gam / &
             (omega*(k2*p0*gam - omega**2))
#if defined(USE3D)
        ! unlike other fields, chi has -sin(phi) dependence
        chi1_l(1:6) = (val(7:12)/kphi)* &
             (0,1)*bzero*k2*(kphi**2*p0*gam - omega**2) / &
             (omega*(k2*p0*gam - omega**2))
        ! toroidal derivative has -cos(phi) dependence
        chi1_l(7:12) = (-val(1:6)*kphi)* &
             (0,1)*bzero*k2*(kphi**2*p0*gam - omega**2) / &
             (omega*(k2*p0*gam - omega**2))
#elif defined(USECOMPLEX)
        chi1_l = val*(0,1)*bzero*k2*(kphi**2*p0*gam - omega**2) / &
             (omega*(k2*p0*gam - omega**2))
#endif
        p1_l = val*bzero*kp2*k2*p0*gam / &
             (k2*p0*gam - omega**2)
     endif
  end select

end subroutine threed_wave_test_per

end module threed_wave_test


!==============================================================================
! 3D diffusion Test
! ~~~~~~~~~~~~
! This is a test case which initializes a 2-dimensional equilibrium
! With a 3-dimensional initial perturbation
!==============================================================================
module threed_diffusion_test

  real, private :: kx, kz

contains

subroutine threed_diffusion_test_init()
  use math
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes, icounter_tt
  real :: x, phi, z, x1, x2, z1, z2
  real :: phi0, x0, z0

  call get_bounding_box(x1, z1, x2, z2)
  kx = pi/(x2-x1)
  kz = pi/(z2-z1)

#ifdef USE3D
  phi0 = pi
#else
  phi0 = 0
#endif
  z0 = (z1 + z2)/2.
  x0 = (x1 + 2.*x2)/3.

  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

     call get_local_vals(l)

     call threed_diffusion_test_equ(x-x1, phi, z-z1)
     call threed_diffusion_test_per(x-x0, phi-phi0, z-z0)

     call set_local_vals(l)
  enddo
end subroutine threed_diffusion_test_init

subroutine threed_diffusion_test_equ(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z

  psi0_l(1) = bx0*sin(kx*x)*sin(kz*z)
  psi0_l(2) = bx0*cos(kx*x)*sin(kz*z)*kx
  psi0_l(3) = bx0*sin(kx*x)*cos(kz*z)*kz
  psi0_l(4) =-bx0*sin(kx*x)*sin(kz*z)*kx**2
  psi0_l(5) = bx0*cos(kx*x)*cos(kz*z)*kx*kz
  psi0_l(6) =-bx0*sin(kx*x)*sin(kz*z)*kz**2

  call constant_field(bz0_l, bzero)
  call constant_field(den0_l, 1.)
  call constant_field(p0_l, p0)
  call constant_field(pe0_l, p0-pi0)

end subroutine threed_diffusion_test_equ


subroutine threed_diffusion_test_per(x, phi, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, phi, z
  real :: ex

  ex = eps*exp(-(rzero*phi/ln)**2-(x/ln)**2-(z/ln)**2)

  p1_l(1) = ex
  p1_l(2) = -2.*ex*x/ln**2
  p1_l(3) = -2.*ex*z/ln**2
  p1_l(4) = 4.*ex*x**2/ln**4 - 2.*ex*ln**2
  p1_l(5) = 4.*ex*x*z/ln**4
  p1_l(6) = 4.*ex*z**2/ln**4 - 2.*ex*ln**2
#ifdef USE3D
  p1_l(7:12) = -2.*phi*p1_l(1:6)*(rzero/ln)**2
#endif

  if(idens.eq.1) den1_l = p1_l
  pe1_l = p1_l*pefac

end subroutine threed_diffusion_test_per

end module threed_diffusion_test



!==============================================================================
! frs
! ~~~~~~~~~~~~
!==============================================================================
module frs

!    real, private :: kx, kz

contains

!========================================================
! init
!========================================================
subroutine frs_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes,m, icounter_tt
  real :: x, phi, z, alx, alz, Bp0,r0,rs,Bz_edge


  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

!     z = z - alz*.5

     call get_local_vals(l)

     call frs_equ(x-rzero, z)
     call frs_per(x-rzero, phi, z)

     call set_local_vals(l)
  enddo
  call finalize(field_vec)

end subroutine frs_init


!========================================================
! equ
!========================================================
subroutine frs_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: r,Bp0,dpsidr,d2psidr,integral,rs, Bz_edge,Bz,dBzdx,dBzdz,beta0,r0
  integer :: m,n


  r0=xlim ! Current channel with
  rs=zlim! Position of singular surface

  m=mpol 
  n=ntor
  Bz_edge=1. 


  Bp0=sqrt((Bz_edge**2 -2.*p0*(1.-rs**2)) &
       /((real(m)/real(n)*rzero*(rs/r0)/(1+(rs/r0)**2)/rs )**2 &
       +2.*(-1./2.)*((1.+(rs/r0)**2)**(-2)-(1.+(1./r0)**2)**(-2))))

  call constant_field(den0_l, 1.)

  r=sqrt(x**2 + z**2)

  
  p0_l(1)=p0*(1.- r**2)
  p0_l(2)=p0*(- 2.*x)
  p0_l(3)=p0*(- 2.*z)
  p0_l(4)=p0*(- 2.)
  p0_l(5)=0.
  p0_l(6)=p0*(- 2.)



  pe0_l(1)=0.
  pe0_l(2)=0.
  pe0_l(3)=0.
  pe0_l(4)=0.
  pe0_l(5)=0.
  pe0_l(6)=0.


  psi0_l(1) = (-rzero*Bp0*r0/2.*log(1.+(r/r0)**2) &
       +rzero*Bp0*r0/2.*log(1.+(1./r0)**2))
  psi0_l(2) = -rzero*Bp0*(x/r0)/(1.+(r/r0)**2)
  psi0_l(3) = -rzero*Bp0*(z/r0)/(1.+(r/r0)**2)
  psi0_l(4) = -rzero*Bp0*r0/2.*(-(2.*x/r0**2)**2/(1+(r/r0)**2)**2 &
       +2./r0**2/(1+(r/r0)**2))
  psi0_l(5) = -rzero*Bp0*r0/2.*(-4.*z*x)/r0**4/(1+(r/r0)**2)**2
  psi0_l(6) = -rzero*Bp0*r0/2.*(-(2.*z/r0**2)**2/(1+(r/r0)**2)**2 &
       +2./r0**2/(1+(r/r0)**2))
  


  integral=(-1./2.)*((1.+(r/r0)**2)**(-2)-(1.+(1./r0)**2)**(-2))

  Bz=sqrt(Bz_edge**2-2.*p0_l(1)-2.*Bp0**2*integral)

  dBzdx=(1./Bz)*(-Bp0**2*2.*x/r0**2/(1.+(r/r0)**2)**3+2.*p0*x)
  dBzdz=(1./Bz)*(-Bp0**2*2.*z/r0**2/(1.+(r/r0)**2)**3+2.*p0*z)
  bz0_l(1)=rzero*Bz
  bz0_l(2)=rzero*dBzdx
  bz0_l(3)=rzero*dBzdz
  bz0_l(4)=rzero* &
       (-1./Bz**2*dBzdx*(-Bp0**2*2.*x/r0**2/(1.+(r/r0)**2)**3+2.*p0*x) &
       +1./Bz*(-Bp0**2*(2./r0**2/(1.+(r/r0)**2)**3 &
                        - 3.*(2.*x/r0**2)**2/(1.+(r/r0)**2)**4  )+2.*p0))
  bz0_l(5)=rzero* &
       (-1./Bz**2*dBzdz*(-Bp0**2*2.*x/r0**2/(1.+(r/r0)**2)**3+2.*p0*x) &
       +1./Bz*(-Bp0**2*(-3.)*4.*x*z/r0**4/(1.+(r/r0)**2)**4))
  bz0_l(6)=rzero* &
       (-1./Bz**2*dBzdz*(-Bp0**2*2.*z/r0**2/(1.+(r/r0)**2)**3+2.*p0*z) &
       +1./Bz*(-Bp0**2*(2./r0**2/(1.+(r/r0)**2)**3 &
                        - 3.*(2.*z/r0**2)**2/(1.+(r/r0)**2)**4  )+2.*p0))

end subroutine frs_equ

subroutine frs1_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes,m
  real :: x, phi, z, alx, alz, Bp0,r0,rs,Bz_edge


  numnodes = owned_nodes()
  do l=1, numnodes
     call get_node_pos(l, x, phi, z)

!     z = z - alz*.5

     call get_local_vals(l)

     call frs1_equ(x-rzero, z)
     call frs_per(x-rzero, phi, z)

     call set_local_vals(l)
  enddo

end subroutine frs1_init

!========================================================
! equ
!========================================================
subroutine frs1_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: r,Bp0,dpsidr,d2psidr,integral,rs, Bz_edge,Bz,dBzdx,dBzdz,beta0,r0
  integer :: m,n

  Bz_edge=rzero

  call constant_field(den0_l, 1.)

  r=sqrt(x**2 + z**2)

  
  call constant_field(p0_l,p0)

  psi0_l(1) = -.5/q0*(1.-.5*r**2)**2
  psi0_l(2) =   x/q0*(1.-.5*r**2)
  psi0_l(3) =   z/q0*(1.-.5*r**2)
  psi0_l(4) =  1./q0*(1.-.5*r**2-x**2)
  psi0_l(5) = -1./q0*x*z
  psi0_l(6) =  1./q0*(1.-.5*r**2-z**2)
  


  Bz=sqrt(Bz_edge**2+4./q0**2*(5./24. - r**2/2. + 3.*r**4/8. - r**6/12))

  bz0_l(1)= Bz
  bz0_l(2)= 2*x/(bz*q0**2)*(-1+3*r**2/2 -r**4/2.)
  bz0_l(3)= 2*z/(bz*q0**2)*(-1+3*r**2/2 -r**4/2.)
  bz0_l(4)= -4*x**2/(bz**3*q0**4)*(-1+3*r**2/2 -r**4/2.)**2   &
          + 2/(bz*q0**2)*(-1+3*r**2/2 -r**4/2.)   &
          + 2*x**2/(bz*q0**2)*(3-2*r)
  bz0_l(5)=  -4*x*z/(bz**3*q0**4)*(-1+3*r**2/2 -r**4/2.)**2   &
          + 2*x*z/(bz*q0**2)*(3-2*r)
  bz0_l(6)= -4*z**2/(bz**3*q0**4)*(-1+3*r**2/2 -r**4/2.)**2   &
          + 2/(bz*q0**2)*(-1+3*r**2/2 -r**4/2.)   &
          + 2*z**2/(bz*q0**2)*(3-2*r)

end subroutine frs1_equ

!========================================================
! per
!========================================================
subroutine frs_per(x, phi, z)

  use basic
  use arrays
  use diagnostics
  use mesh_mod

  implicit none

  integer :: i, numnodes
  real :: x, phi, z
  vectype, dimension(dofs_per_node) :: vmask

     vmask = 1.
     vmask(1:6) = p0_l(1:6)/0.05

     
     ! initial parallel rotation
     u1_l = phizero*vmask

     ! allow for initial toroidal rotation
     vz1_l = 0.
!     if(vzero.ne.0) call add_angular_velocity(vz1_l, x+xzero, vzero*vmask)

     ! add random perturbations
     if(nonrect.eq.0) then
        vmask(1) = 1.
        vmask(2:6) = 0.
     endif
     call random_per(x,phi,z,vmask)
end subroutine frs_per

end module frs

!==============================================================================
! ftz
! ~~~~~~~~~~~~
!==============================================================================
module ftz

!    real, private :: kx, kz

contains

!========================================================
! init
!========================================================
subroutine ftz_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes,m, icounter_tt
  real :: x, phi, z, alx, alz, Bp0,r0,rs,Bz_edge



  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

!     z = z - alz*.5

     call get_local_vals(l)

     call ftz_equ(x-rzero, z)
     call ftz_per(x-rzero, phi, z)

     call set_local_vals(l)
  enddo
  call finalize(field_vec)
end subroutine ftz_init


!========================================================
! equ
!========================================================
subroutine ftz_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: r,Bp0,dpsidr,d2psidr,integral,rs, Bz_edge,Bz,dBzdx,dBzdz,beta0,r0,j0
  integer :: m,n


  call constant_field(den0_l, 1.)
  
  
  j0=4./rzero/2.8

  r=sqrt(x**2 + z**2)

  
  p0_l(1)=p0*(1.- r**2)
  p0_l(2)=p0*(- 2.*x)
  p0_l(3)=p0*(- 2.*z)
  p0_l(4)=p0*(- 2.)
  p0_l(5)=0.
  p0_l(6)=p0*(- 2.)

  pe0_l(1)=0.5*p0_l(1)
  pe0_l(2)=0.5*p0_l(2)
  pe0_l(3)=0.5*p0_l(3)
  pe0_l(4)=0.5*p0_l(4)
  pe0_l(5)=0.5*p0_l(5)
  pe0_l(6)=0.5*p0_l(6)

  psi0_l(1) = -rzero*j0*(r**2/4.-r**4/16.-3./16.)
  psi0_l(2) = -rzero*j0*(1./2.-r**2/4.)*x
  psi0_l(3) = -rzero*j0*(1./2.-r**2/4.)*z
  psi0_l(4) = -rzero*j0*(-1./2.*x**2+1./2.-r**2/4.)
  psi0_l(5) = 1./2.*rzero*j0*x*z
  psi0_l(6) =  -rzero*j0*(-1./2.*z**2+1./2.-r**2/4.)
  




  bz0_l(1)=rzero
  bz0_l(2)=0.
  bz0_l(3)=0.
  bz0_l(4)=0.
  bz0_l(5)=0.
  bz0_l(6)=0.

  

end subroutine ftz_equ


!========================================================
! per
!========================================================
subroutine ftz_per(x, phi, z)

  use basic
  use arrays
  use diagnostics
  use mesh_mod

  implicit none

  integer :: i, numnodes
  real :: x, phi, z
  vectype, dimension(dofs_per_node) :: vmask
!!$
!!$
!!$
!!$
     vmask = 1.
!     vmask(1:6) = /0.05
!     vmask(1) = vmask(1) !- pedge/p0
     
     ! initial parallel rotation
     u1_l = phizero*vmask

     ! allow for initial toroidal rotation
     vz1_l = 0.
!     if(vzero.ne.0) call add_angular_velocity(vz1_l, x+xzero, vzero*vmask)

     ! add random perturbations
     if(nonrect.eq.0) then
        vmask(1) = 1.
        vmask(2:6) = 0.
     endif
     call random_per(x,phi,z,vmask)
end subroutine ftz_per

end module ftz

!==============================================================================
! eigen
! ~~~~~~~~~~~~
!==============================================================================
module eigen

!    real, private :: kx, kz

contains

!========================================================
! init
!========================================================
subroutine eigen_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes,m, icounter_tt
  real :: x, phi, z, alx, alz, Bp0,r0,rs,Bz_edge


  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)
     call get_node_pos(l, x, phi, z)

!     z = z - alz*.5

     call get_local_vals(l)

     call eigen_equ(x-rzero, z)
     call eigen_per(x-rzero, phi, z)

     call set_local_vals(l)
  enddo
  call finalize(field_vec)
end subroutine eigen_init


!========================================================
! equ
!========================================================
subroutine eigen_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
 


  call constant_field(den0_l, 1.)
  call constant_field(p0_l, 0.1)
  call constant_field(pe0_l, 0.05)
  call constant_field(psi0_l, 0.)
  call constant_field(bz0_l, rzero)
  


end subroutine eigen_equ


!========================================================
! per
!========================================================
subroutine eigen_per(x, phi, z)

  use basic
  use arrays
  use diagnostics
  use mesh_mod

  implicit none

  integer :: i, numnodes
  real :: x, phi, z
  vectype, dimension(dofs_per_node) :: vmask
!!$
!!$
!!$
!!$
     vmask = 1.
!     vmask(1:6) = /0.05
!     vmask(1) = vmask(1) !- pedge/p0
     
     ! initial parallel rotation
     u1_l = phizero*vmask

     ! allow for initial toroidal rotation
     vz1_l = 0.
!     if(vzero.ne.0) call add_angular_velocity(vz1_l, x+xzero, vzero*vmask)

     ! add random perturbations
     if(nonrect.eq.0) then
        vmask(1) = 1.
        vmask(2:6) = 0.
     endif
     call random_per(x,phi,z,vmask)
end subroutine eigen_per

end module eigen

!==============================================================================
! int_kink
! ~~~~~~~~~~~~
!==============================================================================
module int_kink


contains

!========================================================
! init
!========================================================
subroutine int_kink_init()
  use basic
  use arrays
  use mesh_mod

  implicit none

  integer :: l, numnodes,m, icounter_tt
  real :: x, phi, z, alx, alz, Bp0,r0,rs,Bz_edge


  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     l = nodes_owned(icounter_tt)

     call get_node_pos(l, x, phi, z)

     call get_local_vals(l)

     call int_kink_equ(x, z)
     call int_kink_per(x, phi, z)

     call set_local_vals(l)
  enddo
  call finalize(field_vec)

end subroutine int_kink_init


!========================================================
! equ
!========================================================
subroutine int_kink_equ(x, z)
  use basic
  use arrays

  implicit none

  real, intent(in) :: x, z
  real :: r, preszero
  real, parameter :: alpha=0.6579, rmin=1., pi=3.1415926536
  
  preszero=p0
  r=sqrt(x**2 + z**2)

  if (NUMVAR.eq.1) then
    call constant_field(den0_l, 1.)
    call constant_field(p0_l,p0)
    call constant_field(pe0_l,p0/2.)
    call constant_field(bz0_l,bzero)
  else ! only for rmin=1.

    bz0_l(1) = 0.225079* &
         Sqrt(1.80383 - 0.957984*preszero - 0.108208*(x**2 + z**2) &
         + 3.22684*preszero*(x**2 + z**2) + 0.162312*(x**2 + z**2)**2 &
         - 4.40951*preszero*(x**2 + z**2)**2 - 0.120231*(x**2 + z**2)**3 &
         + 3.28718*preszero*(x**2 + z**2)**3 + 0.0450867*(x**2 + z**2)**4 &
         - 1.5192*preszero*(x**2 + z**2)**4 - 0.00721387*(x**2 + z**2)**5 &
         + 0.45427*preszero*(x**2 + z**2)**5 &
         - 0.0952705*preszero*(x**2 + z**2)**6 &
         + 0.0159547*preszero*(x**2 + z**2)**7 &
         - 0.00245458*preszero*(x**2 + z**2)**8 &
         + 0.00018182*preszero*(x**2 + z**2)**9)

    bz0_l(2) = (0.11253953951963827* &
         (0. - 0.21641620500000003*x + 6.45367005*preszero*x &
         + 0.6492486150000001*x*(x**2 + z**2) &
         - 17.638059524400003*preszero*x*(x**2 + z**2) &
         - 0.7213873500000001*x*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*x*(x**2 + z**2)**2 &
         + 0.36069367500000005*x*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*x*(x**2 + z**2)**3 &
         - 0.07213873500000001*x*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*x*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*x*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*x*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*x*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*x*(x**2 + z**2)**8)) &
         / Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    bz0_l(3) = (0.11253953951963827* &
         (0. - 0.21641620500000003*z + 6.45367005*preszero*z &
         + 0.6492486150000001*z*(x**2 + z**2) &
         - 17.638059524400003*preszero*z*(x**2 + z**2) &
         - 0.7213873500000001*z*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*z*(x**2 + z**2)**2 &
         + 0.36069367500000005*z*(x**2 + z**2)**3 & 
         - 12.153599550456903*preszero*z*(x**2 + z**2)**3 &
         - 0.07213873500000001*z*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*z*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*z*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*z*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*z*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*z*(x**2 + z**2)**8))/ &
         Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    bz0_l(4) = (-0.056269769759819135* &
         (0. - 0.21641620500000003*x + 6.45367005*preszero*x &
         + 0.6492486150000001*x*(x**2 + z**2) &
         - 17.638059524400003*preszero*x*(x**2 + z**2) &
         - 0.7213873500000001*x*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*x*(x**2 + z**2)**2 &
         + 0.36069367500000005*x*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*x*(x**2 + z**2)**3 &
         - 0.07213873500000001*x*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*x*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*x*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*x*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*x*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*x*(x**2 + z**2)**8)**2)/ &
         (1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 & 
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)**1.5 &
         + (0.11253953951963827*&
         (-0.21641620500000003 + 6.45367005*preszero &
         + 1.2984972300000002*x**2 - 35.276119048800005*preszero*x**2 &
         + 0.6492486150000001*(x**2 + z**2) &
         - 17.638059524400003*preszero*(x**2 + z**2) &
         - 2.8855494000000004*x**2*(x**2 + z**2) &
         + 78.89236034011381*preszero*x**2*(x**2 + z**2) &
         - 0.7213873500000001*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*(x**2 + z**2)**2 &
         + 2.1641620500000003*x**2*(x**2 + z**2)**2 &
         - 72.92159730274142*preszero*x**2*(x**2 + z**2)**2 &
         + 0.36069367500000005*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*(x**2 + z**2)**3 &
         - 0.5771098800000001*x**2*(x**2 + z**2)**3 &
         + 36.341613059299846*preszero*x**2*(x**2 + z**2)**3 &
         - 0.07213873500000001*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*(x**2 + z**2)**4 &
         - 11.432456270094134*preszero*x**2*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*(x**2 + z**2)**5 &
         + 2.680396653036631*preszero*x**2*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*(x**2 + z**2)**6 &
         - 0.5498249544690526*preszero*x**2*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*(x**2 + z**2)**7 &
         + 0.052364281378004994*preszero*x**2*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*(x**2 + z**2)**8))/ &
         Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    bz0_l(5) = (-0.056269769759819135* &
         (0. - 0.21641620500000003*x + 6.45367005*preszero*x &
         + 0.6492486150000001*x*(x**2 + z**2) &
         - 17.638059524400003*preszero*x*(x**2 + z**2) &
         - 0.7213873500000001*x*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*x*(x**2 + z**2)**2 &
         + 0.36069367500000005*x*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*x*(x**2 + z**2)**3 &
         - 0.07213873500000001*x*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*x*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*x*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*x*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*x*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*x*(x**2 + z**2)**8)* &
         (0. - 0.21641620500000003*z + 6.45367005*preszero*z &
         + 0.6492486150000001*z*(x**2 + z**2) &
         - 17.638059524400003*preszero*z*(x**2 + z**2) &
         - 0.7213873500000001*z*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*z*(x**2 + z**2)**2 &
         + 0.36069367500000005*z*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*z*(x**2 + z**2)**3 &
         - 0.07213873500000001*z*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*z*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*z*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*z*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*z*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*z*(x**2 + z**2)**8))/ &
         (1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)**1.5 &
         +(0.11253953951963827* &
         (1.2984972300000002*x*z - 35.276119048800005*preszero*x*z &
         - 2.8855494000000004*x*z*(x**2 + z**2) &
         + 78.89236034011381*preszero*x*z*(x**2 + z**2) &
         + 2.1641620500000003*x*z*(x**2 + z**2)**2 &
         - 72.92159730274142*preszero*x*z*(x**2 + z**2)**2 &
         - 0.5771098800000001*x*z*(x**2 + z**2)**3 &
         + 36.341613059299846*preszero*x*z*(x**2 + z**2)**3 &
         - 11.432456270094134*preszero*x*z*(x**2 + z**2)**4 &
         + 2.680396653036631*preszero*x*z*(x**2 + z**2)**5 &
         - 0.5498249544690526*preszero*x*z*(x**2 + z**2)**6 &
         + 0.052364281378004994*preszero*x*z*(x**2 + z**2)**7))/ &
         Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    bz0_l(6) = (-0.056269769759819135* &
         (0. - 0.21641620500000003*z + 6.45367005*preszero*z &
         + 0.6492486150000001*z*(x**2 + z**2) &
         - 17.638059524400003*preszero*z*(x**2 + z**2) &
         - 0.7213873500000001*z*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*z*(x**2 + z**2)**2 &
         + 0.36069367500000005*z*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*z*(x**2 + z**2)**3 &
         - 0.07213873500000001*z*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*z*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*z*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*z*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*z*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*z*(x**2 + z**2)**8)**2)/ &
         (1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 &
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)**1.5 &
         + (0.11253953951963827* &
         (-0.21641620500000003 + 6.45367005*preszero &
         + 1.2984972300000002*z**2 - 35.276119048800005*preszero*z**2 &
         + 0.6492486150000001*(x**2 + z**2) &
         - 17.638059524400003*preszero*(x**2 + z**2) &
         - 2.8855494000000004*z**2*(x**2 + z**2) &
         + 78.89236034011381*preszero*z**2*(x**2 + z**2) &
         - 0.7213873500000001*(x**2 + z**2)**2 &
         + 19.723090085028453*preszero*(x**2 + z**2)**2 &
         + 2.1641620500000003*z**2*(x**2 + z**2)**2 &
         - 72.92159730274142*preszero*z**2*(x**2 + z**2)**2 &
         + 0.36069367500000005*(x**2 + z**2)**3 &
         - 12.153599550456903*preszero*(x**2 + z**2)**3 &
         - 0.5771098800000001*z**2*(x**2 + z**2)**3 &
         + 36.341613059299846*preszero*z**2*(x**2 + z**2)**3 &
         - 0.07213873500000001*(x**2 + z**2)**4 &
         + 4.542701632412481*preszero*(x**2 + z**2)**4 &
         - 11.432456270094134*preszero*z**2*(x**2 + z**2)**4 &
         - 1.1432456270094133*preszero*(x**2 + z**2)**5 &
         + 2.680396653036631*preszero*z**2*(x**2 + z**2)**5 &
         + 0.2233663877530526*preszero*(x**2 + z**2)**6 &
         - 0.5498249544690526*preszero*z**2*(x**2 + z**2)**6 &
         - 0.03927321103350375*preszero*(x**2 + z**2)**7 &
         + 0.052364281378004994*preszero*z**2*(x**2 + z**2)**7 &
         + 0.003272767586125312*preszero*(x**2 + z**2)**8))/ &
         Sqrt(1.8038251465091153 - 0.9579835619689772*preszero &
         - 0.10820810250000001*(x**2 + z**2) &
         + 3.226835025*preszero*(x**2 + z**2) &
         + 0.16231215375000002*(x**2 + z**2)**2 &
         - 4.409514881100001*preszero*(x**2 + z**2)**2 &
         - 0.12023122500000001*(x**2 + z**2)**3 &
         + 3.287181680838075*preszero*(x**2 + z**2)**3 &
         + 0.045086709375000006*(x**2 + z**2)**4 &
         - 1.5191999438071129*preszero*(x**2 + z**2)**4 &
         - 0.0072138735000000015*(x**2 + z**2)**5 &
         + 0.4542701632412481*preszero*(x**2 + z**2)**5 &
         - 0.09527046891745111*preszero*(x**2 + z**2)**6 &
         + 0.0159547419823609*preszero*(x**2 + z**2)**7 & 
         - 0.0024545756895939844*preszero*(x**2 + z**2)**8 &
         + 0.00018182042145140624*preszero*(x**2 + z**2)**9)

    p0_l(1) =preszero*(1. - 3.106*(x**2 + z**2) + 3.54*(x**2 + z**2)**2 - 1.408*(x**2 + z**2)**3)                         
    p0_l(2) =preszero*(-6.212*x + 14.16*x*(x**2 + z**2) - 8.448*x*(x**2 + z**2)**2)                                       
    p0_l(3) =preszero*(-6.212*z + 14.16*z*(x**2 + z**2) - 8.448*z*(x**2 + z**2)**2)                                       
    p0_l(4) =preszero*(-6.212 + 28.32*(x**2) + 14.16*(x**2 + z**2) - 33.792*(x**2)*(x**2 + z**2) - 8.448*(x**2 + z**2)**2) 
    p0_l(5) =preszero*(28.32*x*z - 33.792*x*z*(x**2 + z**2))
    p0_l(6) =preszero*(-6.212 + 28.32*(z**2) + 14.16*(x**2 + z**2) - 33.792*(z**2)*(x**2 + z**2) - 8.448*(x**2 + z**2)**2) 

    pe0_l(1) = p0_l(1)/2.
    pe0_l(2) = p0_l(2)/2.
    pe0_l(3) = p0_l(3)/2.
    pe0_l(4) = p0_l(4)/2.
    pe0_l(5) = p0_l(5)/2.
    pe0_l(6) = p0_l(6)/2.

    den0_l(1) =1. - 1.6*(x**2 + z**2) + 0.8*(x**2 + z**2)**2 - 3.677e-8*(x**2 + z**2)**3
    den0_l(2) =-3.2*x + 3.2*x*(x**2 + z**2) - 2.2062e-7*x*(x**2 + z**2)**2
    den0_l(3) =-3.2*z + 3.2*z*(x**2 + z**2) - 2.2062e-7*z*(x**2 + z**2)**2
    den0_l(4) =-3.2 + 6.4*(x**2) + 3.2*(x**2 + z**2) - 8.8248e-7*(x**2)*(x**2 + z**2) - 2.2062e-7*(x**2 + z**2)**2
    den0_l(5) =6.4*x*z - 8.8248e-7*x*z*(x**2 + z**2)
    den0_l(6) =-3.2 + 6.4*(z**2) + 3.2*(x**2 + z**2) - 8.8248e-7*(z**2)*(x**2 + z**2) - 2.2062e-7*(x**2 + z**2)**2

  end if

  call constant_field(u0_l, 0.)
  call constant_field(vz0_l, 0.)
  call constant_field(chi0_l, 0.)

  psi0_l(1) = alpha * ( r**2/4. - r**4/(8. * rmin**2)               + r**6/(36. * rmin**4)    )
  psi0_l(2) = alpha * ( x/2.    - r**2 * x/(2. * rmin**2)           + r**4 * x/(6. * rmin**4) )
  psi0_l(3) = alpha * ( z/2.    - r**2 * z/(2. * rmin**2)           + r**4 * z/(6. * rmin**4) )
  psi0_l(4) = alpha * ( 1./2.   - (3. * x**2 + z**2)/(2. * rmin**2) + (5. * x**4 + 6. * x**2 * z**2 + z**4)/(6. * rmin**4) )
  psi0_l(5) = alpha * (         - x * z / rmin**2                   + (x**3 * z + x * z**3) * 2./(3. * rmin**4) )
  psi0_l(6) = alpha * ( 1./2.   - (x**2 + 3. * z**2)/(2. * rmin**2) + (x**4 + 6. * x**2 * z**2 + 5. * z**4)/(6. * rmin**4) )
 

end subroutine int_kink_equ


!========================================================
! per
!========================================================
subroutine int_kink_per(x, phi, z)

  use basic
  use arrays
  use diagnostics
  use mesh_mod

  implicit none

  integer :: i, numnodes
  real :: x, phi, z
  vectype, dimension(dofs_per_node) :: vmask

     ! add random perturbations
        vmask(1) = 1.
        vmask(2:dofs_per_node) = 0.
     call random_per(x,phi,z,vmask)





end subroutine int_kink_per

end module int_kink


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
  integer :: imag, magnetic_region

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
        imag = magnetic_region(ps079(i,:), x_79(i), z_79(i))
        if(imag.ne.0) then
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
                   int3(r2_79,mu79(:,OP_DR,i),nu79(:,OP_DR,j)) + &
                   int3(r2_79,mu79(:,OP_DZ,i),nu79(:,OP_DZ,j))
              temp(i,j,1,2) = &
                   int3(ri_79,mu79(:,OP_DR,i),nu79(:,OP_DZ,j)) - &
                   int3(ri_79,mu79(:,OP_DZ,i),nu79(:,OP_DR,j))
           end if
           
           ! compression equation
           if(imask_chi(i).eq.0) then
              temp(i,j,2,:) = 0.
           else
              temp(i,j,2,1) = &
                   int3(ri_79,mu79(:,OP_DZ,i),nu79(:,OP_DR,j)) - &
                   int3(ri_79,mu79(:,OP_DR,i),nu79(:,OP_DZ,j))
              temp(i,j,2,2) = &
                   int3(ri4_79,mu79(:,OP_DR,i),nu79(:,OP_DR,j)) + &
                   int3(ri4_79,mu79(:,OP_DZ,i),nu79(:,OP_DZ,j))
           end if
        end do

        ! assemble RHS
        ! toroidal rotation
        select case(ivform)
        case(0)
           temp2(i) = int3(r_79,mu79(:,OP_1,i),vz)
        case(1)
           temp2(i) = int3(ri_79,mu79(:,OP_1,i),vz)
        end select

        ! vorticity
        if(imask_vor(i).eq.0) then
           temp3(i,1) = 0.
        else
           temp3(i,1) = &
                int4(temp79f,vp,mu79(:,OP_DR,i),ps079(:,OP_DR)) + &
                int4(temp79f,vp,mu79(:,OP_DZ,i),ps079(:,OP_DZ))
        endif

        ! compression
        if(imask_vor(i).eq.0) then
           temp3(i,2) = 0.
        else
           temp3(i,2) = &
                int5(ri3_79,temp79f,vp,mu79(:,OP_DZ,i),ps079(:,OP_DR)) - &
                int5(ri3_79,temp79f,vp,mu79(:,OP_DR,i),ps079(:,OP_DZ))
        endif

        ! diamagnetic term
        if(db.ne.0. .and. ineo_subtract_diamag.eq.1) then
           temp4(i) = int2(mu79(:,OP_1,i), dia)
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
  use circ_shell_only
  use resistive_wall_test
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
           call circ_shell_only_init()
        case(12,13)
           call resistive_wall_test_init()
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
        case(21,22,25,26,27,28)
           call fixed_q_profiles()
        case(23)
           call frs1_init()
        case(24)
           call rwm_init()
        end select
     else
        ! toroidal equilibria
        select case(itaylor)
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
           call circ_shell_only_init()
           call cartesian_to_cylindrical_all()
        case(12,13)
           call resistive_wall_test_init()
           call cartesian_to_cylindrical_all()
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
        end select
     endif
  end if

  if(iread_neo.eq.1) then
     call set_neo_vel
     call unload_neo
  end if
     
  call den_eq()
  call den_per()

  if(irmp.ge.1 .or. iread_ext_field.ge.1) call rmp_per()

  ! calculate equilibrium and perturbed temperature profiles
  call calculate_temperatures(0, te_field(0),ti_field(0), 1)
  call calculate_temperatures(1, te_field(1),ti_field(1), 1)

  if(iflip_b.eq.1) call mult(bz_field(0), -1.)
  if(iflip_j.gt.0) call mult(psi_field(0), -1.)
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
      
     call get_boundary_mask(itri, BOUNDARY_DIRICHLET, imask, domain_boundary)

     !  assemble matrix    
     do i=1,dofs_per_element
        if(imask(i).eq.0) then
           temp(i,:) = 0.
        else
           do j=1,dofs_per_element
              temp(i,j) = int4(mu79(:,OP_1,i),nu79(:,OP_DR,j),co(:),r(:)) &
                   +      int4(mu79(:,OP_1,i),nu79(:,OP_DZ,j),sn(:),r(:))
           enddo
        end if
        !  assemble rhs
        temp2(i) = int2(mu79(:,OP_1,i),rdpsidr)
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
  if(eqsubtract.eq.1) then
     psi_field(0) = psi_vec
  else
     psi_field(1) = psi_vec
  endif
  
  call destroy_mat(psi_mat)
  call destroy_field(psi_vec)
  !call destroy_field(dpsi_dr)
  
end subroutine kstar_profiles

module basicq
implicit none
real :: q0_qp, rzero_qp, p0_qp, bz_qp, r0_qp, r1_qp, q2_qp, q4_qp, pedge_qp
real :: q6_qp, q8_qp, q10_qp, q12_qp, q14_qp
real :: kappa_qp, kappae_qp, coolrate_qp
integer :: myrank_qp, iprint_qp, itaylor_qp
end module basicq
module LZeqbm
  IMPLICIT NONE

  !Parameters
  REAL, PARAMETER :: aspect_ratio = 18.0  ! R_major / r_wall_b
  REAL, PARAMETER :: B0 = 1.0             ! Toroidal field at r=0
  REAL, PARAMETER :: dfrac = 0.1          ! Bndry layer width / r_plas_a
  REAL, PARAMETER :: Trat = 100.0         ! Ratio of plas to vac temperature
  REAL, PARAMETER :: nrat = 1.0           ! Ratio of plas to vac density
  REAL, PARAMETER :: q_interior = 1.0     ! Safety factor across bulk plasma
  REAL, PARAMETER :: q_a = 0.75           ! Safety factor on plasma surface
  REAL, PARAMETER :: r_plas_a = 0.6       ! Plasma minor radius
  REAL, PARAMETER :: r_wall_b = 1.0       ! Minor radius of ideal wall
  REAL, PARAMETER :: r_tile_c = 0.7       ! Minor radius of tile surface

  !Derived quantities
  REAL R_major                            ! Major radius of magnetic axis
  REAL delta                              ! Width of boundary layer
  REAL r_interior                         ! Radius to interior of bndry layer
  REAL prat                               ! Ratio pf plas to vac pressure
  REAL p_interior                         ! Pressure in bulk plasma

  REAL p_vacuum                           ! Pressure in vacuum region
  REAL q0R0                               ! Central safety factor * major radius
  REAL Jdenom                             ! Intermediate derived quantity
  REAL Bbz                                ! Toroidal field in bndry layer
  REAL gamma                              ! Pressure gradient in bndry layer
end module LZeqbm

 subroutine init_qp
use basic
use basicq
implicit none
 !input variables:
  bz_qp = bzero
  r0_qp = alpha0
  r1_qp = th_gs
  q0_qp = q0
  q2_qp = alpha1
  q4_qp = alpha2
  rzero_qp = rzero
  p0_qp = p0
  pedge_qp = pedge
  kappa_qp = kappa0
  kappae_qp = alpha3
  iprint_qp = iprint
  myrank_qp = myrank
  itaylor_qp = itaylor
  coolrate_qp = coolrate
  q6_qp = libetap
  q8_qp = p1
  q10_qp = p2
  q12_qp = djdpsi
  q14_qp = divcur
end subroutine init_qp
subroutine fixed_q_profiles()

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
use basicq
implicit none


vectype, dimension (dofs_per_node) :: vec_l
vectype, dimension (dofs_per_element) :: dofsps, dofsbz, dofspr
real , dimension(npoints) :: rtemp79a, rtemp79b, rtemp79c
real :: x, phi, z, r
integer :: numnodes, nelms, l, itri, i, j, ier, icounter_tt
type (field_type) :: psi_vec, bz_vec, p_vec

call create_field(psi_vec)
call create_field(bz_vec)
call create_field(p_vec)


 if(myrank.eq.0 .and. iprint.ge.1) write (*,2000) bz_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2001) r0_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2002) q0_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2003) q2_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2004) q4_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2005) rzero_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2006) p0_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2007) pedge_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2008) kappa_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2009) kappae_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2010) iprint_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2011) myrank_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2012) itaylor_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2013) coolrate_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2014) q6_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2015) q8_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2016) q10_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2017) q12_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2018) q14_qp 
 if(myrank.eq.0 .and. iprint.ge.1) write (*,2019) r1_qp

 2000 format( 'bz_qp =', 1pe12.4)
 2001 format( 'r0_qp =', 1pe12.4)
 2002 format( 'q0_qp =', 1pe12.4)
 2003 format( 'q2_qp =', 1pe12.4)
 2004 format( 'q4_qp =', 1pe12.4)
 2005 format( 'rzero_qp =', 1pe12.4)
 2006 format( 'p0_qp =', 1pe12.4)
 2007 format( 'pedge_qp =', 1pe12.4)
 2008 format( 'kappa_qp =', 1pe12.4)
 2009 format( 'kappae_qp =', 1pe12.4)
 2010 format( 'iprint_qp =', i5)
 2011 format( 'myrank_qp =', i5)
 2012 format( 'itaylor_qp =', i5)
 2013 format( 'coolrate_qp =', 1pe12.4)
 2014 format( 'q6_qp =', 1pe12.4)
 2015 format( 'q8_qp =', 1pe12.4)
 2016 format( 'q10_qp =', 1pe12.4)
 2017 format( 'q12_qp =', 1pe12.4)
 2018 format( 'q14_qp =', 1pe12.4)
 2019 format( 'r1_qp  =', 1pe12.4)

if(itaylor.eq.22) call setupLZeqbm

numnodes = owned_nodes()

if(myrank.eq.0 .and. iprint.eq.1) write(*,*) "numnodes = ", numnodes

do icounter_tt=1,numnodes
   l = nodes_owned(icounter_tt)

   call get_node_pos(l,x,phi,z)

         call constant_field(den0_l,1.)
!   call constant_field(bz0_l,bz_qp)
!   call constant_field(p0_l,p0)

         call set_node_data(den_field(0),l,den0_l)
!   call set_node_data(bz_field(0),l,bz0_l)
!   call set_node_data(p_field(0),l,p0_l)
   
         call constant_field(u1_l,0.)
           call int_kink_per(x, phi, z)
         call set_node_data(u_field(1),l,u1_l)
enddo

if(myrank.eq.0 .and. iprint.ge.1) write(*,*) "before loop over elements"

nelms = local_elements()
do itri=1,nelms

   call define_element_quadrature(itri,int_pts_diag, int_pts_tor)
   call define_fields(itri,0,1,0) ! defines x_79,z_79,mu,nu

!  assemble matrix

   do i=1,dofs_per_element
      do j=1,npoints
         r = (sqrt((x_79(j)-xmag)**2 + (z_79(j)-zmag)**2))/r0_qp  ! normalized radius
         call getvals_qsolver(r,rtemp79a(j),rtemp79b(j),rtemp79c(j))
      enddo
#ifdef USECOMPLEX
      temp79a = cmplx(rtemp79a)
      temp79b = cmplx(rtemp79b)
      temp79c = cmplx(rtemp79c)
#else
      temp79a = rtemp79a
      temp79b = rtemp79b
      temp79c = rtemp79c
#endif
      dofsps(i) = int2(mu79(:,OP_1,i),temp79a)
      dofsbz(i) = int2(mu79(:,OP_1,i),temp79b)
      dofspr(i) = int2(mu79(:,OP_1,i),temp79c)
   enddo
   call vector_insert_block(psi_vec%vec,itri,1,dofsps,VEC_ADD)
   call vector_insert_block(bz_vec%vec ,itri,1,dofsbz,VEC_ADD)
   call vector_insert_block(p_vec%vec  ,itri,1,dofspr,VEC_ADD)
enddo

! solve for psi
 if(myrank.eq.0 .and. iprint.ge.1) print *, "solving psi"
 
  call newvar_solve(psi_vec%vec,mass_mat_lhs)
  call newvar_solve(bz_vec%vec ,mass_mat_lhs)
  call newvar_solve(p_vec%vec  ,mass_mat_lhs)
 if(eqsubtract.eq.1) then
   psi_field(0) = psi_vec
   bz_field(0)  = bz_vec
   p_field(0)   = p_vec
 else
   psi_field(1) = psi_vec
   bz_field(1)  = bz_vec
   p_field(1)   = p_vec
 endif
   pe_field(0) = p_field(0)
   pe_field(1) = p_field(1)
   call mult(pe_field(0),pefac)
   call mult(pe_field(1),pefac)

 call destroy_field(psi_vec)
 call destroy_field(bz_vec)
 call destroy_field(p_vec)

  call finalize(field_vec)
end subroutine fixed_q_profiles

subroutine getvals_qsolver(rval,bpsival,ival,pval)
use basicq
implicit none

integer, parameter :: N=1000  !  (number of intervals)
integer :: ifirstq = 1
real :: dpsi,rval,psival,bpsival,ival,pval,cubicinterp,qfunc, qpfunc, pfunc, ppfunc
real :: psimid, qmid, qpmid, ppmid, denom,fterm,gterm,aquad,bquad,cquad,disc,A_qp
integer :: j
real, dimension (0:N) :: bpsi, btor, bpolor, psi, jphi, jthor, gradpor, equor, pary

if(ifirstq.eq.1) then
   ifirstq = 0
   A_qp = rzero_qp/r0_qp

!           small psi is normalized r**2
   dpsi = 1./N
   do j=0,N
      psi(j) = j*dpsi
!  DEBUG
   if(iprint_qp .ge.1 .and. myrank_qp .eq.0) write(*,4000) j, psi(j), qfunc(psi(j))
 4000 format('j   psi   qfunc(psi)', i5, 1p2e12.4)
   enddo

!  boundary condition at edge
   btor(N) = bz_qp
   bpsi(N) = 0.
   bpolor(N) = btor(N)/(2.*A_qp*qfunc(psi(N)))
   if(myrank_qp.eq.0 .and. iprint_qp.ge.1) write(*,3000) btor(N),bpsi(N),bpolor(N)
 3000 format( 'btor(N), bpsi(N), bpolor(N) =',1p3e12.4)
   if(myrank_qp.eq.0 .and. iprint_qp.ge.1) write(*,*) "diagnostics to follow"
!  integrate first order ode from boundary in
   do j=N,1,-1
      psimid = (j-.5)*dpsi
      qmid = A_qp*qfunc(psimid)
      qpmid= A_qp*qpfunc(psimid)
      ppmid= ppfunc(psimid)
      denom = psimid + qmid**2
      fterm = -(1 -psimid*qpmid/qmid)/denom
      gterm = -ppmid*qmid**2/denom

      aquad = 1. + .5*dpsi*fterm
      bquad = dpsi*fterm*btor(j)
      cquad = -btor(j)**2*(1.-.5*dpsi*fterm) + 2.*dpsi*gterm

      disc = bquad**2 - 4.*aquad*cquad
      btor(j-1) = (-bquad + sqrt(disc))/(2.*aquad)
      bpolor(j-1) =btor(j-1)*r0_qp**2/(2.*rzero_qp*qfunc(psi(j-1)))
      bpsi(j-1)   = bpsi(j) - .5*dpsi*(bpolor(j)+bpolor(j-1)) 
      if(myrank_qp.eq.0 .and. iprint_qp.ge.1) &
      write(*,1002) qmid,denom,fterm,gterm,aquad,bquad,cquad,disc,btor(j-1)
   enddo
 1002 format(1p9e9.1)
  
!  calculate poloidal and toroidal fields in cell centers
   do j=1,N-1
     jphi(j) = 4.*((j+.5)*(bpsi(j+1)-bpsi(j))-(j-.5)*(bpsi(j)-bpsi(j-1)))/(dpsi*r0_qp**2) 
     jthor(j) =  2*((bpsi(j+1)-bpsi(j))*rzero_qp*qfunc((j+0.5)*dpsi) &
                  - (bpsi(j)-bpsi(j-1))*rzero_qp*qfunc((j-0.5)*dpsi))/(dpsi**2*r0_qp**2)
     gradpor(j) =  ppfunc(j*dpsi)
     pary(j) = pfunc(psi(j))
!    error in equilibrium equation
     equor(j) = (jphi(j)*bpolor(j)+jthor(j)*btor(j)+gradpor(j))*sqrt(j*dpsi)
   enddo
     pary(0) =  pfunc(psi(0))
     pary(N) =  pfunc(psi(N))
!
if(myrank_qp .eq. 0 .and. iprint_qp .ge. 1) then
   write(6,1001)
 1001 format(" j       r**2       bpsi        btor         p          equil")
   do j=0,N
     write(6,1000) j,psi(j),bpsi(j),btor(j),pary(j),equor(j)
   enddo
1000 format(i3,1p7e12.4)
endif

endif !   end of initialization

psival = rval*rval
bpsival = cubicinterp(psival,psi,bpsi,N)
   ival = cubicinterp(psival,psi,btor,N)
   pval = cubicinterp(psival,psi,pary,N)

end subroutine getvals_qsolver

function cubicinterp(x,xary,yary,N)
implicit none
real :: x,xt,del,m1,m2,a,b,c,d,cubicinterp
real, dimension(0:N) :: xary, yary
integer :: i,N
  xt = 0
  a = 0
  b = 0
  c = 0
  d = 0
if      (x .le. xary(0))   then
  a = yary(0)
else if (x .ge. xary(N))   then
  a = yary(N)
else if (x .le. xary(1))   then
  xt = x - xary(0)
  del = xary(1) - xary(0)
  m2 =  (yary(2)-yary(0))/(xary(2)-xary(0))
  a = yary(0)
  b = 2.*(yary(1) - yary(0))/del    - m2
  c =   -(yary(1) - yary(0))/del**2 + m2/del
else if (x .ge. xary(N-1)) then
  xt = x - xary(N-1)
  del = xary(N) - xary(N-1)
  m1 =  (yary(N)-yary(N-2))/(xary(N)-xary(N-2))
  a = yary(N-1)
  b = m1
  c = (yary(N) - yary(N-1) - m1*del)/del**2
else

  do i=1,N-2
     if(x.ge.xary(i) .and. x.le.xary(i+1)) then
       xt = x - xary(i)
       del = xary(i+1) - xary(i)
       m1 = (yary(i+1)-yary(i-1))/(xary(i+1)-xary(i-1))
       m2 = (yary(i+2)-yary(i  ))/(xary(i+2)-xary(i  ))
       a = yary(i)
       b = m1
       c =  3.*(yary(i+1) - yary(i) - m1*del)/del**2 - (m2 - m1)/del
       d = -2.*(yary(i+1) - yary(i) - m1*del)/del**3 + (m2 - m1)/del**2
       exit
    endif
  enddo

endif
cubicinterp = a + b*xt + c*xt**2 + d*xt**3
return
end function cubicinterp

function qfunc(psi)    !   q  (safety factor)
use basicq
real :: psi,qfunc,q_LZ  !  note:  psi = r**2
real :: c0,c1,c2,c3,c4 
real :: asq, bigA, bigB,psis  
real ra0

select case(itaylor_qp)

case(21)
   c0 = 4.179343
   c1 = -0.080417
   c2=-8.659146
   c3 = 10.668674
   c4 = -4.108323
   qfunc = (q0_qp) + psi**2*(c0+c1*psi+c2*psi**2+c3*psi**3+c4*psi**4)

case(22)
    qfunc = q_LZ(psi)

case(25)
   qfunc = (q0_qp) + psi*(q2_qp + q4_qp*psi)

case(26)
   qfunc = q0_qp*(1. + (psi/q2_qp)**q4_qp )**(1./q4_qp)
 
case(27)
!new coding
   asq = q4_qp**2
   bigA = (-2. + 3.*q0_qp/q2_qp)/asq 
   bigB = (1. -  2.*q0_qp/q2_qp)/asq**2
   if(psi .le. asq) then
     qfunc = q0_qp/(1. + bigA*psi + bigB*psi**2)
   else
     qfunc = q2_qp*psi/asq
   endif

case(28)
   ra0 = q4_qp*abs(((q8_qp/q10_qp/q0_qp)**(q12_qp+q14_qp*q4_qp**2)-(1.,0.))**(-0.5/(q12_qp+q14_qp*q4_qp**2)))
   qfunc = (1+(psi/ra0**2)**(q12_qp+psi*q14_qp))**(1/(q12_qp+psi*q14_qp))*q0_qp*(1+q6_qp/exp((sqrt(psi)-r1_qp)**2/q2_qp**2))

end select
return
end function qfunc

function qpfunc(psi)   !   derivative of q wrt psi
use basicq
real :: psi,qpfunc,qprime_LZ   !  note:  psi=r^2
real :: c0,c1,c2,c3,c4   
real :: asq, bigA, bigB  
real ra0

select case (itaylor_qp)

case(21)
   c0 = 4.179343
   c1 = -0.080417
   c2=-8.659146
   c3 = 10.668674
   c4 = -4.108323
   qpfunc =  psi*(2.*c0+3.*c1*psi+4.*c2*psi**2+5.*c3*psi**3+6.*c4*psi**4)

case(22)
   qpfunc = qprime_LZ(psi)

case(25)
   qpfunc = (q2_qp + 2.*q4_qp*psi)

case(26)
   qpfunc = q0_qp*(1. + (psi/q2_qp)**q4_qp )**((1.-q4_qp)/q4_qp)       &
                *(1./q2_qp)**q4_qp*q4_qp*psi**(q4_qp-1)
case(27)
!new coding
   asq = q4_qp**2
   bigA = (-2. + 3.*q0_qp/q2_qp)/asq 
   bigB = (1. -  2.*q0_qp/q2_qp)/asq**2
   if(psi .le. asq) then
     qpfunc = -q0_qp*(bigA + 2.*bigB*psi)/(1. + bigA*psi + bigB*psi**2)**2
   else
     qpfunc = q2_qp/asq
   endif

case(28)
   psis = max(1.e-5,psi)
   ra0 = q4_qp*abs(((q8_qp/q10_qp/q0_qp)**(q12_qp+q14_qp*q4_qp**2)-(1.,0.))**(-0.5/(q12_qp+q14_qp*q4_qp**2)))
   qpfunc = (1+(psis/ra0**2)**(q12_qp+psis*q14_qp))**(1/(q12_qp+psis*q14_qp))*q0_qp*(-((log(1+(psis/ra0**2)**(q12_qp+psis*q14_qp))*q14_qp)/(q12_qp+psis*q14_qp)**2)+((psis/ra0**2)**(q12_qp+psis*q14_qp)*(log(psis/ra0**2)*q14_qp+(q12_qp+psis*q14_qp)/psis))/((1+(psis/ra0**2)**(q12_qp+psis*q14_qp))*(q12_qp+psis*q14_qp)))*(1+q6_qp/exp((sqrt(psis)-r1_qp)**2/q2_qp**2))-((1+(psis/ra0**2)**(q12_qp+psis*q14_qp))**(1/(q12_qp+psis*q14_qp))*q0_qp*q6_qp*(sqrt(psis)-r1_qp))/(exp((sqrt(psis)-r1_qp)**2/q2_qp**2)*sqrt(psis)*q2_qp**2)

end select
return

end function qpfunc

function pfunc(psi)    !   p  (pressure)
use basicq
real :: psi,pfunc,p_LZ   !  note:  psi=r^2
real :: asq, bigA, bigB 
select case(itaylor_qp)

case(21,25,26)
 pfunc = p0_qp * (1. - 3.2*psi + 4.16*psi**2 - 2.56*psi**3 + 0.64*psi**4)

case(22)
  pfunc = p_LZ(psi)

case(27)
!new coding
   asq = q4_qp**2
   bigA = (-4. + 6.*q0_qp/q2_qp)/asq 
   bigB = (3. -  6.*q0_qp/q2_qp)/asq**2
   if(psi .le. asq) then
     pfunc = p0_qp*(1+ bigA*psi + bigB*psi**2)**(2./3.) + pedge_qp
   else
     pfunc = pedge_qp
   endif

end select
return
end function pfunc

function ppfunc(psi)    !  derivative of p wrt psi
use basicq
real :: psi,ppfunc,pprime_LZ   !  note:  psi=r^2
real :: asq, bigA, bigB 
select case(itaylor_qp)

case(21,25,26)
  ppfunc = p0_qp * (-3.2 + 8.32*psi - 7.68*psi**2 + 2.56*psi**3)
case(22)
  ppfunc = pprime_LZ(psi)
case(27)
   asq = q4_qp**2
   bigA = (-4. + 6.*q0_qp/q2_qp)/asq 
   bigB = ( 3. - 6.*q0_qp/q2_qp)/asq**2
   if(psi .lt. asq) then
     ppfunc = p0_qp*(2./3.)*(1+ bigA*psi + bigB*psi**2)**(-1./3.)    &
                           *(bigA + 2.*bigB*psi)
   else
     ppfunc = 0.
   endif

end select
return
end function ppfunc




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Profiles for a pseudo-VMEC equilibrium structure resembling L. Zakharov's    !
! circular cross-section, straight-cylinder, flat q, ideal equilibrium with a  !
! surface current.                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setupLZeqbm
    use lzeqbm
    IMPLICIT NONE

    delta = dfrac * r_plas_a
    r_interior = r_plas_a - delta
    R_major = aspect_ratio * r_wall_b
    q0R0 = q_interior * R_major
    Jdenom = q0R0**2 + r_interior**2
    Bbz = q0R0**2 * B0 / Jdenom
    gamma = (1.5 * (q0R0**2 * B0**2 * r_interior**4 / Jdenom**2 - &
         r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / &
         (r_interior**3 - r_plas_a**3)) ! /aspect_ratio
    prat = nrat * Trat
    p_vacuum = gamma*delta/(prat - 1.0)
    p_interior = p_vacuum + gamma*delta
  end subroutine setupLZeqbm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function p_LZ(psi)
    use lzeqbm
    implicit none
    REAL p_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       p_LZ = p_interior
    ELSE IF (rmin < r_plas_a) THEN
       p_LZ = p_interior - gamma*(rmin - r_interior)
    ELSE
       p_LZ = p_vacuum
    endIF
  end function p_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pprime_LZ(psi)
    use lzeqbm
    implicit none
    REAL pprime_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, dpdr

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       pprime_LZ = 0.0
       RETURN
    endIF

    IF (rmin < r_plas_a) THEN
       dpdr = -gamma
    ELSE
       dpdr = 0.0
    endIF

    pprime_LZ = 0.5*dpdr/rmin
  end function pprime_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function q_LZ(psi)
    use lzeqbm
    implicit none
    intrinsic sqrt
    REAL q_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, Btheta


    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       q_LZ = q_interior
    ELSE IF (rmin < r_plas_a) THEN
       Btheta = SQRT(2.0*gamma*(rmin**3 - r_plas_a**3)/3.0 + &
            r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / rmin
       q_LZ = rmin * Bbz / (R_major * Btheta)
    ELSE
       q_LZ = q_a * (rmin/r_plas_a)**2
    endIF
  end function q_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function qprime_LZ(psi)
    use lzeqbm
    implicit none
    REAL qprime_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, Btheta, Bthetaprime, dqdr

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       qprime_LZ = 0.0
       RETURN
    endIF

    IF (rmin < r_plas_a) THEN
       Btheta = SQRT(2.0*gamma*(rmin**3 - r_plas_a**3)/3.0 + &
            r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / rmin
       dqdr = Bbz * (2.0 - gamma*rmin/Btheta**2) / (R_major * Btheta)
    ELSE
       dqdr = 2.0*rmin*q_a / r_plas_a**2
    endIF

    qprime_LZ = 0.5*dqdr/rmin
  end function qprime_LZ
function get_kappa(psi)  ! thermal conductivity for itaylor=27, ikappafunc=12
use basicq
real :: psi,get_kappa   !  note:  psi=r^2
real :: asq, bigA, bigB, num1, num2, denom, jedge

   asq = q4_qp**2
   bigA = (1. - 1.5*q0_qp/q2_qp)/asq 
   bigB = (1. -  2.*q0_qp/q2_qp)/asq**2
   jedge = (pedge_qp/p0_qp)**1.5
   if(psi .gt. asq) psi = asq    !     temporary fix
   if(psi .le. asq) then
     num1 = 1. - 2*psi*bigA + psi**2*bigB + jedge
     num2 = (1. - 4.*psi*bigA + 3.*psi**2*bigB + jedge )**(1./3.)
     denom =  asq*(bigA - 1.5*psi*bigB)
     get_kappa = kappa_qp*num1*num2/denom
   else
     get_kappa = kappae_qp
   endif

return
end function get_kappa

function hsink_qp(psi)
use basicq
real :: psi, hsink_qp, pfunc
    hsink_qp = coolrate_qp*pfunc(psi)

return

end function hsink_qp
