module rmp

  use coils

  real :: tf_tilt, tf_tilt_angle
  real :: tf_shift, tf_shift_angle
  real, dimension(maxcoils) :: pf_tilt, pf_tilt_angle
  real, dimension(maxcoils) :: pf_shift, pf_shift_angle

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
  use gradshafranov

  implicit none

  include 'mpif.h'

  type(schaffer_field), allocatable :: sf(:)
  type(matrix_type) :: br_mat
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

  vectype, dimension(dofs_per_element,dofs_per_element) :: temp
  vectype, dimension(dofs_per_element) :: temp2, temp3, temp4
  vectype, dimension(MAX_PTS) :: bbr, bbz, drbr, dzbr, drbz, dzbz

  type(field_type) :: psi_f, bz_f, p_f
  real, dimension(MAX_PTS, 1, 6) :: g

#ifdef USECOMPLEX
  complex :: sfac
  complex :: tilt_co, tilt_sn
  complex :: shift_co, shift_sn
#else
  real, dimension(int_pts_main) :: co, sn
  real, dimension(int_pts_main) :: tilt_co, tilt_sn
  real, dimension(int_pts_main) :: shift_co, shift_sn
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

     temp79a = 0.
     temp79b = 0.
     temp79c = 0.
     temp79d = 0.

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

     end select

     if(tf_tilt.ne.0. .or. tf_shift.ne.0.) then
#ifdef USECOMPLEX
        tilt_co  = tf_tilt*deg2rad*exp(-(0,1)*tf_tilt_angle*deg2rad )
        tilt_sn  = tf_tilt*deg2rad*exp(-(0,1)*tf_tilt_angle*deg2rad )*(0,-1)
        shift_co = tf_shift*exp(-(0,1)*tf_shift_angle*deg2rad)
        shift_sn = tf_shift*exp(-(0,1)*tf_shift_angle*deg2rad)*(0,-1)
#else
        tilt_co  = tf_tilt*deg2rad*cos(phi_79 - tf_tilt_angle*deg2rad)
        tilt_sn  = tf_tilt*deg2rad*sin(phi_79 - tf_tilt_angle*deg2rad)
        shift_co = tf_shift*cos(phi_79 - tf_shift_angle*deg2rad)
        shift_sn = tf_shift*sin(phi_79 - tf_shift_angle*deg2rad)
#endif
        temp79a = temp79a + bzero*rzero * ( &
             - (z_79/x_79**2)*tilt_co    &
             - (  1./x_79**2)*shift_sn)
        temp79b = temp79b + bzero*rzero * ( &
             - (z_79/x_79**2)*tilt_sn    &
             + (  1./x_79**2)*shift_co)
        temp79c = temp79c + bzero*rzero * ( &
             + (  1./x_79   )*tilt_co)
     end if

     do i=1, numcoils_vac
        j = coil_mask(i)
        if(pf_tilt(j).ne.0. .or. pf_shift(j).ne.0.) then
#ifdef USECOMPLEX
           tilt_co  = pf_tilt(j)*deg2rad* &
                exp(-(0,1)*pf_tilt_angle(j)*deg2rad )
           tilt_sn  = pf_tilt(j)*deg2rad* &
                exp(-(0,1)*pf_tilt_angle(j)*deg2rad )*(0,-1)
           shift_co = pf_shift(j)* &
                exp(-(0,1)*pf_shift_angle(j)*deg2rad)
           shift_sn = pf_shift(j)* &
                exp(-(0,1)*pf_shift_angle(j)*deg2rad)*(0,-1)
#else
           tilt_co  = pf_tilt(j)*deg2rad* &
                cos(phi_79 - pf_tilt_angle(j)*deg2rad)
           tilt_sn  = pf_tilt(j)*deg2rad* &
                sin(phi_79 - pf_tilt_angle(j)*deg2rad)
           shift_co = pf_shift(j)* &
                cos(phi_79 - pf_shift_angle(j)*deg2rad)
           shift_sn = pf_shift(j)* &
                sin(phi_79 - pf_shift_angle(j)*deg2rad)
#endif
           call gvect(x_79, z_79, npoints, xc_vac(i), zc_vac(i), 1, g, 0, ier)
           bbr =   ic_vac(i)* g(:,1,3)/x_79
           bbz =  -ic_vac(i)* g(:,1,2)/x_79
           drbr =  ic_vac(i)*(g(:,1,5)/x_79 - g(:,1,3)/x_79**2)
           dzbr =  ic_vac(i)* g(:,1,6)/x_79
           drbz = -ic_vac(i)*(g(:,1,4)/x_79 - g(:,1,2)/x_79**2)
           dzbz = -ic_vac(i)* g(:,1,5)/x_79

           temp79a = temp79a + &
                (-shift_co*drbr &
                + tilt_sn*(-bbz - x_79*dzbr + z_79*drbr))
           temp79b = temp79b + &
                (shift_sn*bbr/x_79 &
                +tilt_co*(z_79/x_79*bbr - bbz))
           temp79c = temp79c + &
                (-shift_co*drbz &
                + tilt_sn*(bbr - x_79*dzbz + z_79*drbz))
        end if
     end do


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
              temp79d = temp79d + p/p0_norm + pedge
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
        end do

        ! assemble RHS
        temp2(i) = &
             + int3(ri_79,mu79(:,OP_DR,i),temp79c) &
             - int3(ri_79,mu79(:,OP_DZ,i),temp79a)

        temp3(i) = int3(r_79,mu79(:,OP_1,i),temp79b)

        if(read_p) temp4(i) = int2(mu79(:,OP_1,i),temp79d)
     end do

     call insert_block(br_mat, itri, 1, 1, temp(:,:), MAT_ADD)

     call vector_insert_block(psi_vec, itri, 1, temp2(:), MAT_ADD)
     call vector_insert_block(bz_vec, itri, 1, temp3(:), MAT_ADD)
     if(read_p) call vector_insert_block(p_vec, itri, 1, temp4(:), MAT_ADD)
  end do

  if(myrank.eq.0 .and. iprint.ge.2) print *, 'Finalizing...'
  call finalize(br_mat)
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

  if(numvar.ge.2) then
     ! Solve F = R B_phi
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving bz..."
     call newsolve(mass_mat_lhs%mat,bz_vec,ier)
     if(extsubtract.eq.1) then
        bz_ext = bz_f     
     else
        bz_field(1) = bz_f
     end if

#if defined(USECOMPLEX) || defined(USE3D)
     ! Solve R^2 del_perp^2 f = F
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving f..."
     if(extsubtract.eq.1) then
        bf_ext = 0.
        call solve_newvar1(bf_mat_lhs,bf_ext,mass_mat_rhs_bf, &
             bz_ext, bf_ext)
     else
        bf_field(1) = 0.
        call solve_newvar1(bf_mat_lhs,bf_field(1),mass_mat_rhs_bf, &
             bz_field(1), bf_field(1))
     end if
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

  integer :: l, ierr
  character(len=13) :: ext_field_name
  type(schaffer_field), allocatable :: sf(:)

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
end subroutine rmp_per

end module rmp
