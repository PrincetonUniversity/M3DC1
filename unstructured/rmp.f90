module rmp

contains

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

  type(schaffer_field), allocatable :: sf(:)
  type(matrix_type) :: br_mat, bf_mat
  type(vector_type) :: psi_vec, bz_vec, p_vec
  integer :: i, j, itri, nelms, ier

  complex, dimension(int_pts_main) :: fr, fphi, fz
  complex, dimension(MAX_PTS) :: br, btheta, bz, phase
  real, dimension(MAX_PTS) :: r, theta, arg
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

  if(irmp.eq.1) then
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
  print *, "create_mat coils bf_mat", bf_mat%imatrix 
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

     select case(irmp)
     case(1)
        fr   = 0.    ! B_R
        fphi = 0.    ! B_phi
        fz   = 0.    ! B_Z

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
        temp79d = 0.


     ! m/n Vacuum fields
     case(2)

        if(itor.eq.1) then 
           temp79a = 0.
           temp79b = 0.
           temp79c = 0.
        else
           r = sqrt((x_79 - xzero)**2 + (z_79 - zzero)**2 + regular**2)
           theta = -atan2(z_79 - zzero, x_79 - xzero)
           arg = ntor*r/rzero
           phase = exp((0,1)*(ntor*phi_79/rzero - mpol*theta))
           do i=1, npoints
              br(i)     = 0.5*(Bessel_I(mpol-1, arg(i)) + Bessel_I(mpol+1, arg(i)))
              btheta(i) = Bessel_I(mpol, arg(i)) / r(i)
              bz(i)     = Bessel_I(mpol, arg(i))
           end do
           br     =        (ntor/rzero)*phase*br     *eps
           btheta = -(0,1)*mpol        *phase*btheta *eps
           bz     =  (0,1)*(ntor/rzero)*phase*bz     *eps
#ifdef USECOMPLEX
           temp79a =  br*cos(theta) - btheta*sin(theta)
           temp79b =  bz
           temp79c = -br*sin(theta) - btheta*cos(theta)
#else
           temp79a =  real(br)*cos(theta) - real(btheta)*sin(theta)
           temp79b =  real(bz)
           temp79c = -real(br)*sin(theta) - real(btheta)*cos(theta)
#endif
        end if
        temp79d = 0.


     case default
        print *, 'ZERO'
        temp79a = 0.
        temp79b = 0.
        temp79c = 0.
        temp79d = 0.
     end select


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
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving psi: ier = ", ier
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

!!$  if(irmp.eq.3) then
!!$     numnodes = owned_nodes()
!!$     do icounter_tt=1,numnodes
!!$        l = nodes_owned(icounter_tt)
!!$
!!$        call boundary_node(l,is_boundary,izone,izonedim,normal,curv,x,z)
!!$        if(.not.is_boundary) cycle
!!$
!!$        call get_local_vals(l)
!!$
!!$        dx = x - xmag
!!$        dz = z - zmag
!!$        r2 = dx**2 + dz**2
!!$
!!$        ! psi = 0.5*bx0*exp( i*2*theta)*exp(i*ntor*phi)
!!$        !     = 0.5*bx0*(cos(2*theta) + i*sin(2*theta))*exp(i*ntor*phi)
!!$
!!$        ! cos(2*theta)
!!$        psi1_l(1) = (dx**2 - dz**2)/r2
!!$        psi1_l(2) = 4.*dx*dz**2/r2**2
!!$        psi1_l(3) =-4.*dz*dx**2/r2**2
!!$        psi1_l(4) = 4.*dz**2*(1.-4.*dx**2/r2)/r2**2
!!$        psi1_l(5) = 8.*dx*dz*(dx**2 - dz**2)/r2**3
!!$        psi1_l(6) =-4.*dx**2*(1.-4.*dz**2/r2)/r2**2
!!$
!!$#ifdef USECOMPLEX
!!$        if(ntor.lt.0) then
!!$           ii = -(0,1.)
!!$        else
!!$           ii =  (0,1.)
!!$        endif
!!$        ! sin(2*theta)
!!$        psi1_l(1) = psi1_l(1) + ii*2.*dx*dz/r2
!!$        psi1_l(2) = psi1_l(2) + ii*2.*dz*(1. - 2.*dx**2/r2)/r2
!!$        psi1_l(3) = psi1_l(3) + ii*2.*dx*(1. - 2.*dz**2/r2)/r2
!!$        psi1_l(4) = psi1_l(4) + ii*4.*dx*dz*(1. - 4.*dz**2/r2)/r2**2
!!$        psi1_l(5) = psi1_l(5) - ii*2.*(1. - 8.*dx**2*dz**2/r2**2)/r2
!!$        psi1_l(6) = psi1_l(6) + ii*4.*dx*dz*(1. - 4.*dx**2/r2)/r2**2
!!$#endif
!!$        psi1_l = psi1_l*0.5*bx0
!!$
!!$        call set_local_vals(l)
!!$     end do
!!$     return
!!$  endif

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

!!$  ! leave perturbation only on the boundary
!!$  if(irmp.eq.2) then
!!$     numnodes = owned_nodes()
!!$     do icounter_tt=1,numnodes
!!$        l = nodes_owned(icounter_tt)
!!$        call boundary_node(l,is_boundary,izone,izonedim,normal,curv,x,z)
!!$        if(.not.is_boundary) then
!!$           call get_local_vals(l)
!!$           psi1_l = 0.
!!$           call set_local_vals(l)
!!$        endif
!!$     end do
!!$  endif
end subroutine rmp_per

end module rmp
