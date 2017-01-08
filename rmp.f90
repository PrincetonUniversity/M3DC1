module rmp

  use coils
  use read_schaffer_field

  real :: tf_tilt, tf_tilt_angle
  real :: tf_shift, tf_shift_angle
  real, dimension(maxcoils) :: pf_tilt, pf_tilt_angle
  real, dimension(maxcoils) :: pf_shift, pf_shift_angle

  real, dimension(maxfilaments), private :: xc_na, zc_na
  complex, dimension(maxfilaments), private :: ic_na
  integer, private :: nc_na
  logical, private :: read_p

  type(schaffer_field), allocatable, private :: sf(:)

contains

!=========================================================================
subroutine calculate_external_fields()
  use basic
  use math
  use mesh_mod
  use sparse
  use arrays
  use coils
  use m3dc1_nint
  use newvar_mod
  use boundary_conditions

  implicit none

  include 'mpif.h'

  type(matrix_type) :: br_mat
  type(vector_type) :: psi_vec, bz_vec, p_vec
  integer :: i, j, itri, nelms, ier

  vectype, dimension(dofs_per_element,dofs_per_element) :: temp
  vectype, dimension(dofs_per_element) :: temp2, temp3, temp4

  type(field_type) :: psi_f, bz_f, p_f

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

     call rmp_field(npoints, npoints_tor, npoints_pol, &
          x_79, phi_79, z_79, &
          temp79a, temp79b, temp79c, temp79d)

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
  call calculate_external_fields()

  ! unload data
  if(iread_ext_field.ge.1) then
     do l=1, iread_ext_field
        call unload_schaffer_field(sf(l))
     end do
     deallocate(sf)
  end if
end subroutine rmp_per

subroutine rmp_field(n, nt, np, x, phi, z, br, bphi, bz, p)
  use math
  use basic
  use coils
  use gradshafranov

  implicit none

  integer, intent(in) :: n, nt, np
  real, intent(in), dimension(n) :: x, phi, z
  vectype, intent(out), dimension(n) :: br, bphi, bz
  vectype, intent(out), dimension(n), optional :: p

  integer :: i, j, ier

  complex, dimension(int_pts_main) :: fr, fphi, fz
  complex, dimension(MAX_PTS) :: brv, bthetav, bzv, phase
  real, dimension(MAX_PTS) :: r, theta, arg
#ifdef USE3D
  real, dimension(MAX_PTS) :: gr, gphi, gz, q
#endif

  real, dimension(MAX_PTS, 1, 6) :: g
  vectype, dimension(MAX_PTS) :: bbr, bbz, drbr, dzbr, drbz, dzbz

#ifdef USECOMPLEX
  complex :: sfac
  complex :: tilt_co, tilt_sn
  complex :: shift_co, shift_sn
#else
  real, dimension(int_pts_main) :: co, sn
  real, dimension(int_pts_main) :: tilt_co, tilt_sn
  real, dimension(int_pts_main) :: shift_co, shift_sn
#endif

  br = 0.
  bphi = 0.
  bz = 0.
  p = 0.

  select case(irmp)
  case(1)
     fr   = 0.    ! B_R
     fphi = 0.    ! B_phi
     fz   = 0.    ! B_Z
     
     do i=1, nc_na, 2
        call pane(ic_na(i),xc_na(i),xc_na(i+1),zc_na(i),zc_na(i+1), &
             np,x,z,ntor,fr,fphi,fz)
     end do
     
#ifdef USECOMPLEX
     br = fr
     bphi = fphi
     bz = fz
     p = 0.
#else
     do i=1, nt
        co(1:np) = &
             cos(ntor*phi((i-1)*np+1:i*np))
        sn(1:np) = &
             sin(ntor*phi((i-1)*np+1:i*np))
        br((i-1)*np+1:i*np) = &
             real(fr(1:np))*co(1:np) + &
             aimag(fr(1:np))*sn(1:np)
        bphi((i-1)*np+1:i*np) = &
             real(fphi(1:np))*co(1:np) + &
             aimag(fphi(1:np))*sn(1:np)
        bz((i-1)*np+1:i*np) = &
             real(fz(1:np))*co(1:np) + &
             aimag(fz(1:np))*sn(1:np)
     end do
#endif

     br = -twopi*br
     bphi = -twopi*bphi
     bz = -twopi*bz
     if(present(p)) p = 0.
          
     ! m/n Vacuum fields
  case(2)
     
     if(itor.eq.1) then 
        br = 0.
        bphi = 0.
        bz = 0.
     else
        r = sqrt((x - xzero)**2 + (z - zzero)**2 + regular**2)
        theta = -atan2(z - zzero, x - xzero)
        arg = ntor*r/rzero
        phase = exp((0,1)*(ntor*phi/rzero - mpol*theta))
        do i=1, n
           brv(i)     = 0.5*(Bessel_I(mpol-1, arg(i)) + Bessel_I(mpol+1, arg(i)))
           bthetav(i) = Bessel_I(mpol, arg(i)) / r(i)
           bzv(i)     = Bessel_I(mpol, arg(i))
        end do
        brv     =        (ntor/rzero)*phase*brv     *eps
        bthetav = -(0,1)*mpol        *phase*bthetav *eps
        bzv     =  (0,1)*(ntor/rzero)*phase*bzv     *eps
#ifdef USECOMPLEX
        br =  brv*cos(theta) - bthetav*sin(theta)
        bphi =  bzv
        bz = -brv*sin(theta) - bthetav*cos(theta)
#else
        br =  real(brv)*cos(theta) - real(bthetav)*sin(theta)
        bphi =  real(bzv)
        bz = -real(brv)*sin(theta) - real(bthetav)*cos(theta)
#endif
     end if
     p = 0.
     
  end select
  
  if(tf_tilt.ne.0. .or. tf_shift.ne.0.) then
#ifdef USECOMPLEX
     tilt_co  = tf_tilt*deg2rad*exp(-(0,1)*tf_tilt_angle*deg2rad )
     tilt_sn  = tf_tilt*deg2rad*exp(-(0,1)*tf_tilt_angle*deg2rad )*(0,-1)
     shift_co = tf_shift*exp(-(0,1)*tf_shift_angle*deg2rad)
     shift_sn = tf_shift*exp(-(0,1)*tf_shift_angle*deg2rad)*(0,-1)
#else
     tilt_co  = tf_tilt*deg2rad*cos(phi - tf_tilt_angle*deg2rad)
     tilt_sn  = tf_tilt*deg2rad*sin(phi - tf_tilt_angle*deg2rad)
     shift_co = tf_shift*cos(phi - tf_shift_angle*deg2rad)
     shift_sn = tf_shift*sin(phi - tf_shift_angle*deg2rad)
#endif
     br = br + bzero*rzero * ( &
          - (z/x**2)*tilt_co    &
          - (  1./x**2)*shift_sn)
     bphi = bphi + bzero*rzero * ( &
          - (z/x**2)*tilt_sn    &
          + (  1./x**2)*shift_co)
     bz = bz + bzero*rzero * ( &
          + (  1./x   )*tilt_co)
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
             cos(phi - pf_tilt_angle(j)*deg2rad)
        tilt_sn  = pf_tilt(j)*deg2rad* &
             sin(phi - pf_tilt_angle(j)*deg2rad)
        shift_co = pf_shift(j)* &
             cos(phi - pf_shift_angle(j)*deg2rad)
        shift_sn = pf_shift(j)* &
             sin(phi - pf_shift_angle(j)*deg2rad)
#endif
        call gvect(x, z, n, xc_vac(i), zc_vac(i), 1, g, 0, ier)
        bbr =   ic_vac(i)* g(:,1,3)/x
        bbz =  -ic_vac(i)* g(:,1,2)/x
        drbr =  ic_vac(i)*(g(:,1,5)/x - g(:,1,3)/x**2)
        dzbr =  ic_vac(i)* g(:,1,6)/x
        drbz = -ic_vac(i)*(g(:,1,4)/x - g(:,1,2)/x**2)
        dzbz = -ic_vac(i)* g(:,1,5)/x
        
        br = br + &
             (-shift_co*drbr &
             + tilt_sn*(-bbz - x*dzbr + z*drbr))
        bphi = bphi + &
             (shift_sn*bbr/x &
             +tilt_co*(z/x*bbr - bbz))
        bz = bz + &
             (-shift_co*drbz &
             + tilt_sn*(bbr - x*dzbz + z*drbz))
     end if
  end do


  if(iread_ext_field.ne.0) then
     do i=1, iread_ext_field 
#if defined(USECOMPLEX)
        sfac = exp(-cmplx(0.,1.)*ntor*shift_ext_field(i)*pi/180.)
        call get_external_field_ft(sf(i),x,z,fr,fphi,fz,n)
        br = br + (1e4/b0_norm)*fr  *scale_ext_field*sfac
        bphi = bphi + (1e4/b0_norm)*fphi*scale_ext_field*sfac
        bz = bz + (1e4/b0_norm)*fz  *scale_ext_field*sfac
#elif defined(USE3D)
        call get_external_field(sf(i),&
             x,phi-shift_ext_field(i)*pi/180.,z,&
             gr,gphi,gz,q,n)
        br = br + (1e4/b0_norm)*gr  *scale_ext_field
        bphi = bphi + (1e4/b0_norm)*gphi*scale_ext_field
        bz = bz + (1e4/b0_norm)*gz  *scale_ext_field
        if(sf(i)%vmec .and. present(p)) then
           p = p + q/p0_norm + pedge
        end if
#endif
     end do
  end if

end subroutine rmp_field

end module rmp
