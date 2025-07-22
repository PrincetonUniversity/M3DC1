module rmp

  use coils
  use read_schaffer_field

  implicit none

  real :: tf_tilt, tf_tilt_angle
  real :: tf_shift, tf_shift_angle
  real, dimension(maxcoils) :: pf_tilt, pf_tilt_angle
  real, dimension(maxcoils) :: pf_shift, pf_shift_angle
  real :: rmp_atten

  real, dimension(maxfilaments), private :: xc_na, zc_na
  complex, dimension(maxfilaments), private :: ic_na
  integer, private :: nc_na
  logical, private :: read_p

  type(schaffer_field), allocatable, private :: sf(:)

contains

!==============================================================================
subroutine rmp_per(ilin)
  use basic
  use arrays
  use coils
  use boundary_conditions
  use read_schaffer_field

  implicit none

  integer :: l, ierr
  character(len=13) :: ext_field_name
  integer, intent(in), optional :: ilin

  if(type_ext_field.le.0) then ! RMP field
     if(iread_ext_field.ge.1) allocate(sf(iread_ext_field))
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

     ! calculate external fields from non-axisymmetric coils and external field
     call calculate_external_fields(ilin)

     ! unload data
     call deallocate_sf

  else if(type_ext_field.eq.2) then ! Stellarator/3D vacuum field 
    allocate(sf(iread_ext_field))
    call read_stellarator_field(file_ext_field)

     ! calculate external fields from non-axisymmetric coils and external field
     call calculate_external_fields(ilin)
     ! unload data
     call deallocate_sf

  else if(type_ext_field.eq.3) then
     ! calculate external fields from currents
     call calculate_field_from_j(1)

  else
    if(myrank.eq.1) print *, &
      'ERROR: RMP fields require type_ext_field < 0 or =2.'
    call safestop(58)
  end if


end subroutine rmp_per

!==============================================================================
! For free boundary stellarator and 3D fields 

subroutine load_stellarator_field
  use basic
  use arrays
  use coils
  use init_common 
  use boundary_conditions
  use read_schaffer_field

  implicit none
  integer :: ierr

  allocate(sf(iread_ext_field))

  if(type_ext_field.eq.1) then ! Free boundary stellarator (no field subtraction)

    call read_stellarator_field(file_total_field)
    call calculate_external_fields(0)
    call deallocate_sf

  else if(type_ext_field.eq.2 .and. extsubtract.eq.1) then ! With field substraction 

    ! First load the total field while temporarily setting extsubtract=0
    extsubtract = 0
    call read_stellarator_field(file_total_field)
    call calculate_external_fields(0)
    call deallocate_sf
    ! Then reset extsubtract=1
    extsubtract = 1

  else if(type_ext_field.eq.2 .and. extsubtract.eq.0) then
    if(myrank.eq.0) print *, &
      'ERROR: type_ext_field=2 requires extsubtract=1 when itaylor=41'
    call safestop(57)
  else
    if(myrank.eq.0) print *, &
      'ERORR: Invalid external field options for stellarator.' 
    call safestop(57)

  end if
  call init_perturbations

end subroutine load_stellarator_field

! Select stellarator data to read
subroutine read_stellarator_field(field_name)
  use basic
  use arrays
  use coils
  use boundary_conditions
  use read_schaffer_field

  implicit none

  character(len=256), intent(in) :: field_name
  integer :: ierr

  if(field_name(1:10).eq.'fieldlines') then
    call load_fieldlines_field(sf(iread_ext_field), field_name,ierr)
  else if(field_name(1:4).eq.'mips') then
    call load_mips_field(sf(iread_ext_field), field_name, ierr)
#ifdef USEST
  else if(field_name(1:4).eq.'hint') then
    call load_hint_field(sf(iread_ext_field), field_name, ierr)
  else if(field_name(1:5).eq.'mgrid') then
    call load_mgrid_field(sf(iread_ext_field), field_name, vmec_filename, ierr)
#endif
  else
! Use Schaffer field as fall back
    call load_schaffer_field(sf(iread_ext_field),field_name,isample_ext_field, &
            isample_ext_field_pol,ierr)
  end if
  if(ierr.lt.0) then 
     if(myrank.eq.0) then 
        print *, "Error: could not open file: ", field_name
        print *, "FIELDLINES filename should begin with fieldlines"
        print *, "HINT filename should begin with hint"
        print *, "MGRID filename should begin with mgrid"
        print *, "MIPS filename should begin with mips"
        print *, "Other filenames will be treated as Schaffer format"
     end if
     call safestop(50)
  end if

end subroutine read_stellarator_field


!==============================================================================
! Deallocate sf field
subroutine deallocate_sf
  use basic

  implicit none

  integer :: l

  if(iread_ext_field.ge.1) then
     do l=1, iread_ext_field
        call unload_schaffer_field(sf(l))
     end do
     deallocate(sf)
  end if

end subroutine deallocate_sf

!==============================================================================
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

  complex, dimension(n) :: fr, fphi, fz
  complex, dimension(n) :: brv, bthetav, bzv, phase
  real, dimension(n) :: r, theta, arg, fac, atten
#ifdef USE3D
  real, dimension(n) :: gr, gphi, gz, q
#endif

  real, dimension(n, 1, 6) :: g
  vectype, dimension(n) :: bbr, bbz, drbr, dzbr, drbz, dzbz

#ifdef USECOMPLEX
  complex :: sfac
  complex :: tilt_co, tilt_sn
  complex :: shift_co, shift_sn
#else
  real, dimension(n) :: co, sn
  real, dimension(n) :: tilt_co, tilt_sn
  real, dimension(n) :: shift_co, shift_sn
#endif

  br = 0.
  bphi = 0.
  bz = 0.
  if(present(p)) p = 0.

  select case(irmp)
  case(0)

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
#else
     do i=1, nt
        co(1:np) = &
             cos(ntor*phi((i-1)*np+1:i*np))
        sn(1:np) = &
             sin(ntor*phi((i-1)*np+1:i*np))
        br((i-1)*np+1:i*np) = &
             real(fr(1:np))*co(1:np) - &
             aimag(fr(1:np))*sn(1:np)
        bphi((i-1)*np+1:i*np) = &
             real(fphi(1:np))*co(1:np) - &
             aimag(fphi(1:np))*sn(1:np)
        bz((i-1)*np+1:i*np) = &
             real(fz(1:np))*co(1:np) - &
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
        if(myrank.eq.0) print *, 'ERROR: irmp==2 only implemented for itor=0'
        call safestop(1)
        br = 0.
        bphi = 0.
        bz = 0.
     else
        r = sqrt((x - xzero)**2 + (z - zzero)**2 + regular**2)
        theta = -atan2(z - zzero, x - xzero)
        arg = ntor*r/rzero
        phase = exp((0,1)*(ntor*phi/rzero - mpol*theta))
        !phase = exp((0,1)*(ntor*phi/rzero - mpol*theta-pi/2))
        do i=1, n
           brv(i)     = 0.5*(Bessel_I(mpol-1, arg(i)) + Bessel_I(mpol+1, arg(i)))
           bthetav(i) = Bessel_I(mpol, arg(i)) / r(i)
           bzv(i)     = Bessel_I(mpol, arg(i))
        end do
        fac = (ntor/rzero)*0.5*(Bessel_I(mpol-1, ntor/rzero) + Bessel_I(mpol+1, ntor/rzero))
        !fac = (ntor/rzero)/(mpol*bzero)
        !fac = eps 
        if(rmp_atten.ne.0) then
           atten = exp((r-1.)/rmp_atten)
        else
           atten = 1.
        end if
        !if(present(p)) p = real(-phase*bzv     *eps / fac * atten)
        !if(present(p)) p = pedge 
        brv     =        (ntor/rzero)*phase*brv     *eps / fac * atten
        bthetav = -(0,1)*mpol        *phase*bthetav *eps / fac * atten
        bzv     =  (0,1)*(ntor/rzero)*phase*bzv     *eps / fac * atten
#ifdef USECOMPLEX
        br =  brv*cos(theta) - bthetav*sin(theta)
        bphi =  bzv
        bz = -brv*sin(theta) - bthetav*cos(theta)
#else
        br =  real(brv)*cos(theta) - real(bthetav)*sin(theta)
        bphi =  real(bzv) !+ bzero
        bz = -real(brv)*sin(theta) - real(bthetav)*cos(theta)
#endif
     end if

  case default
     if(myrank.eq.0) print *, 'Error: Option not recognized: irmp = ', irmp
     call safestop(1)
     
  end select

!!$  if(tf_tilt.ne.0. .or. tf_shift.ne.0.) then
!!$#ifdef USECOMPLEX
!!$     tilt_co  = tf_tilt*deg2rad*exp(-(0,1)*tf_tilt_angle*deg2rad )
!!$     tilt_sn  = tf_tilt*deg2rad*exp(-(0,1)*tf_tilt_angle*deg2rad )*(0,-1)
!!$     shift_co = tf_shift*exp(-(0,1)*tf_shift_angle*deg2rad)
!!$     shift_sn = tf_shift*exp(-(0,1)*tf_shift_angle*deg2rad)*(0,-1)
!!$#else
!!$     tilt_co  = tf_tilt*deg2rad*cos(phi - tf_tilt_angle*deg2rad)
!!$     tilt_sn  = tf_tilt*deg2rad*sin(phi - tf_tilt_angle*deg2rad)
!!$     shift_co = tf_shift*cos(phi - tf_shift_angle*deg2rad)
!!$     shift_sn = tf_shift*sin(phi - tf_shift_angle*deg2rad)
!!$#endif
!!$     br = br + bzero*rzero * ( &
!!$          - (z/x**2)*tilt_co    &
!!$          - (  1./x**2)*shift_sn)
!!$     bphi = bphi + bzero*rzero * ( &
!!$          - (z/x**2)*tilt_sn    &
!!$          + (  1./x**2)*shift_co)
!!$     bz = bz + bzero*rzero * ( &
!!$          + (  1./x   )*tilt_co)
!!$  end if
  
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
           p = p + q*10/p0_norm + pedge
        end if
#endif
     end do
  end if
end subroutine rmp_field


subroutine calculate_external_fields(ilin)
  use basic
  use math
  use mesh_mod
  use sparse
  use arrays
  use coils
  use m3dc1_nint
  use newvar_mod
  use boundary_conditions
  use gradshafranov

  implicit none

  include 'mpif.h'

  type(matrix_type) :: br_mat, bf_mat
  type(vector_type) :: psi_vec, bz_vec, p_vec, bf_vec
  integer :: i, itri, nelms, ier, ibound, ipsibound
  integer, intent(in), optional :: ilin
  integer :: il

  vectype, dimension(dofs_per_element,dofs_per_element,2,2) :: temp
  vectype, dimension(dofs_per_element,2) :: temp2
  vectype, dimension(dofs_per_element) :: temp3, temp4
  vectype, dimension(dofs_per_element,dofs_per_element) :: temp5

  type(field_type) :: psi_f, bz_f, p_f, bf_f, bfp_f

  if(myrank.eq.0 .and. iprint.ge.2) print *, "Calculating error fields"

  if(.not.present(ilin)) then
     il = 1
  else 
     il = ilin
  end if

  if(irmp.eq.1) then
     call load_coils(xc_na, zc_na, ic_na, nc_na, &
          'rmp_coil.dat', 'rmp_current.dat')
  end if
  if((any(pf_tilt.ne.0.) .or. any(pf_shift.ne.0.)) &
       .and. numcoils_vac.eq.0) then
     call load_coils(xc_vac,zc_vac,ic_vac,numcoils_vac, &
          'coil.dat','current.dat',coil_mask,filaments)
  end if

  call create_vector(p_vec,1)
  call associate_field(p_f,p_vec,1)

  call create_vector(psi_vec,2)
  call associate_field(psi_f,psi_vec,1)
  call associate_field(bfp_f,psi_vec,2)

  call create_vector(bz_vec,1)
  call associate_field(bz_f,bz_vec,1)

  call create_vector(bf_vec,1)
  call associate_field(bf_f,bf_vec,1)

  call set_matrix_index(br_mat, br_mat_index)
  call create_mat(br_mat, 2, 2, icomplex, 1)

  call set_matrix_index(bf_mat, bf_mat_index)
  call create_mat(bf_mat, 1, 1, icomplex, 1)

  read_p = .false.
  if(iread_ext_field.ne.0) then
     do i=1, iread_ext_field 
        if(sf(i)%vmec) read_p = .true.
     end do
  end if

  ! boundary condition on psi
  ipsibound = BOUNDARY_NONE
  !ipsibound = BOUNDARY_DIRICHLET

  ! Boundary condition on f
  !ibound = BOUNDARY_NEUMANN
  ibound = BOUNDARY_DIRICHLET
  !ibound = BOUNDARY_NONE

  if(myrank.eq.0 .and. iprint.ge.2) print *, 'calculating field values...'
  nelms = local_elements()
!!$OMP PARALLEL DO &
!!$OMP& PRIVATE(temp,temp2,temp3,temp4)
  do itri=1,nelms
        
     call define_element_quadrature(itri,int_pts_main,5)
     call define_fields(itri,0,1,0)

     call rmp_field(npoints, npoints_tor, npoints_pol, &
          x_79, phi_79, z_79, &
          temp79a, temp79b, temp79c, temp79d)

     ! psi_equation
     ! ~~~~~~~~~~~~
     ! Minimize BR, BZ
     temp(:,:,1,1) = intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DR),ri2_79) &
                   + intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),ri2_79) &
                   + regular* intxx3(mu79(:,:,OP_1),nu79(:,:,OP_1),ri4_79)
#if defined(USECOMPLEX) || defined(USE3D)
     temp(:,:,1,2) = intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DR),ri_79) &
          -          intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DZ),ri_79)
     !temp(:,:,1,2) = 0.
#else
     temp(:,:,1,2) = 0.
#endif

     !temp2(:,1) = 0.
 
     temp2(:,1) = intx3(mu79(:,:,OP_DR),temp79c,ri_79) &
          -       intx3(mu79(:,:,OP_DZ),temp79a,ri_79)

           ! Jphi
!!$      temp(:,:,1,1) = -intxx3(mu79(:,:,OP_1),nu79(:,:,OP_GS),ri2_79)
!!$      temp(:,:,1,2) = 0.
!!$      temp2(:,1) = &
!!$               + intx3(mu79(:,:,OP_DR),temp79c,ri_79) &
!!$               - intx3(mu79(:,:,OP_DZ),temp79a,ri_79)
!!$      temp2(:,1) = 0.


     ! f equation
     ! ~~~~~~~~~~

     !temp(:,:,2,1) = intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DR),ri_79) &
     !     -          intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DZ),ri_79)
     temp(:,:,2,1) = 0.

     temp(:,:,2,2) = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_LP),r2_79) !& 
     !     + regular*intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
     !temp(:,:,2,2) = &
     !     -intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DR),r2_79) &
     !     -intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),r2_79) &
     !     -2.*intxx3(mu79(:,:,OP_1),nu79(:,:,OP_DR),r_79) 
     !     + regular*intxx3(mu79(:,:,OP_1),nu79(:,:,OP_1),ri2_79)
     !temp(:,:,2,2) = &
     !     -intxx2(mu79(:,:,OP_DR),nu79(:,:,OP_DR)) &
     !     -intxx2(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ)) &
     !     + regular*intxx3(mu79(:,:,OP_1),nu79(:,:,OP_1),ri2_79)

#if defined(USE3D)
     temp2(:,2) = -intx3(mu79(:,:,OP_DP),r_79,temp79b)
#elif defined(USECOMPLEX) 
     temp2(:,2) = intx3(mu79(:,:,OP_1),r_79,temp79b)*rfac
#else
     temp2(:,2) = 0.
#endif

     !temp2(:,2) = intx2(mu79(:,:,OP_DR),temp79a) &
     !           + intx2(mu79(:,:,OP_DZ),temp79c) 
        
     !temp3 = 0. 
     temp3 = intx3(mu79(:,:,OP_1),r_79,temp79b)
        
     if(read_p) temp4 = intx2(mu79(:,:,OP_1),temp79d)

     temp5 = temp(:,:,2,2) ! &
     !      + intxx3(mu79(:,:,OP_1),nu79(:,:,OP_DPP),ri4_79)  

     call apply_boundary_mask(itri, ipsibound, temp(:,:,1,1), &
          tags=BOUND_DOMAIN)
     call apply_boundary_mask(itri, ipsibound, temp(:,:,1,2), &
          tags=BOUND_DOMAIN)
     call apply_boundary_mask(itri, ibound, temp(:,:,2,1), &
          tags=BOUND_DOMAIN)
     call apply_boundary_mask(itri, ibound, temp(:,:,2,2), &
          tags=BOUND_DOMAIN)

     call apply_boundary_mask(itri, ibound, temp5, &
          tags=BOUND_DOMAIN)

!!$OMP CRITICAL
     call insert_block(br_mat, itri, 1, 1, temp(:,:,1,1), MAT_ADD)
     call insert_block(br_mat, itri, 1, 2, temp(:,:,1,2), MAT_ADD)
     call insert_block(br_mat, itri, 2, 1, temp(:,:,2,1), MAT_ADD)
     call insert_block(br_mat, itri, 2, 2, temp(:,:,2,2), MAT_ADD)

     call insert_block(bf_mat, itri, 1, 1, temp5, MAT_ADD)

     call vector_insert_block(psi_vec, itri, 1, temp2(:,1), MAT_ADD)
     call vector_insert_block(psi_vec, itri, 2, temp2(:,2), MAT_ADD)

     call vector_insert_block(bz_vec, itri, 1, temp3, MAT_ADD)
     call vector_insert_block(bf_vec, itri, 1, temp3, MAT_ADD)

     if(read_p) call vector_insert_block(p_vec, itri, 1, temp4, MAT_ADD)
!!$OMP END CRITICAL
  end do
!OMP END PARALLEL DO

  if(myrank.eq.0 .and. iprint.ge.2) print *, 'Finalizing...'
  call sum_shared(bz_vec)
  if(read_p) call sum_shared(p_vec)
  call sum_shared(bf_vec)
 
  !solve p
  if(read_p) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving p..."

     call newsolve(mass_mat_lhs%mat,p_vec,ier)
     p_field(il) = p_f
     pe_field(il) = p_f
     call mult(pe_field(il), pefac) 
  end if

  if(il.eq.0) then
     call add(p_field(il),pedge)
     call add(pe_field(il),pedge*pefac)
  end if

  call boundary_dc(bf_vec,mat=bf_mat)
  !call boundary_dc(bf_vec,p_vec,bf_mat)
  call finalize(bf_mat)
  if(numvar.ge.2) then
     ! Solve F = R B_phi
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving bz..."
     call newsolve(mass_mat_lhs%mat,bz_vec,ier)
     ! Solve for f
     if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving bf..."
     call newsolve(bf_mat,bf_vec,ier)
     if(extsubtract.eq.1) then
        bz_ext = bz_f  ! For all cases: RMP, error fields, ST
        bf_ext = bf_f
        if (type_ext_field.ge.1 .and. (itaylor.eq.40 .or. itaylor.eq.41)) then 
        ! subtract vacuum from total field
          call mult(bz_f, -1.)
          call mult(bf_f, -1.)
          call add(bz_field(0), bz_f)
          call add(bf_field(0), bf_f)
        end if
     else
        bz_field(il) = bz_f
        bf_field(il) = bf_f
     end if
  end if

  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving psi..."
  call sum_shared(psi_vec)
  call boundary_rmp(psi_vec,br_mat)
  call finalize(br_mat)
  call newsolve(br_mat,psi_vec,ier)
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving psi: ier = ", ier
  if(extsubtract.eq.1) then
     psi_ext = psi_f ! For all cases: RMP, error fields, ST
     bfp_ext = bfp_f
     if (type_ext_field.ge.1 .and. (itaylor.eq.40 .or. itaylor.eq.41)) then 
     ! subtract vacuum from total field
       call mult(psi_f, -1.)
       call mult(bfp_f, -1.)
       call add(psi_field(0), psi_f)
       call add(bfp_field(0), bfp_f) 
     end if
  else
     psi_field(il) = psi_f
     bfp_field(il) = bfp_f
  end if

  if(iflip_b.eq.1) then 
     call mult(bz_field(il), -1.)
  end if
  if(iflip_j.eq.1) then
     call mult(psi_field(il), -1.)
     call mult(bfp_field(il), -1.)
  end if

  call destroy_vector(psi_vec)
  call destroy_vector(bz_vec)
  call destroy_vector(bf_vec)
  call destroy_vector(p_vec)
  call destroy_mat(br_mat)
  call destroy_mat(bf_mat)

  ! Add fields from TF shift / tilt
  ! This is separated because its solution has psi=0.
  if(tf_tilt.ne.0. .or. tf_shift.ne.0.) call tf_shift_tilt

  if(myrank.eq.0 .and. iprint.ge.2) print *, "Done calculating error fields"
end subroutine calculate_external_fields

subroutine tf_shift_tilt
  use basic
  use arrays
  use math
  use field
  use mesh_mod
  use m3dc1_nint
  use newvar_mod

  implicit none

  type(matrix_type) :: bz_mat, fp_mat
  type(field_type) :: ff, bzf, fpf
  integer :: nelms, itri, ier
  vectype, dimension(dofs_per_element) :: temp1
  vectype, dimension(dofs_per_element, dofs_per_element) :: temp2, temp3

#ifdef USECOMPLEX
  complex :: shift_co, tilt_sn
#else
  real, dimension(MAX_PTS) :: shift_co, tilt_sn
#endif

  call set_matrix_index(bz_mat, 100)
  call create_mat(bz_mat, 1, 1, icomplex, 0)
  call set_matrix_index(fp_mat, 101)
  call create_mat(fp_mat, 1, 1, icomplex, 0)

  call create_field(ff)
  call create_field(fpf)
  call create_field(bzf)
  ff = 0.
  fpf = 0.
  bzf = 0.

  nelms = local_elements()
  do itri=1, nelms
     call define_element_quadrature(itri,int_pts_main,5)
     call define_fields(itri,0,1,0)

#ifdef USECOMPLEX
     tilt_sn  = tf_tilt*deg2rad*exp(-(0,1)*tf_tilt_angle*deg2rad )*(0,-1)
     shift_co = tf_shift*exp(-(0,1)*tf_shift_angle*deg2rad)
#else
     tilt_sn  = tf_tilt*deg2rad*sin(phi_79 - tf_tilt_angle*deg2rad)
     shift_co = tf_shift*cos(phi_79 - tf_shift_angle*deg2rad)
#endif

     temp79a = bzero*rzero/x_79*(shift_co - z_79*tilt_sn)
     temp1 = intx2(mu79(:,:,OP_1),temp79a)

     temp2 = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_LP),r2_79)
#if defined(USECOMPLEX) || defined(USE3D)
     temp3 = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_DP))
#else
     temp3 = 0.
#endif 

     call vector_insert_block(ff%vec, itri, 1, temp1, VEC_ADD)
     call insert_block(bz_mat, itri, 1, 1, temp2, MAT_ADD)
     call insert_block(fp_mat, itri, 1, 1, temp3, MAT_ADD)
  end do
  call finalize(bz_mat)
  call finalize(fp_mat)

  call newvar_solve(ff%vec, mass_mat_lhs)
  if(extsubtract.eq.1) then 
     call add(bf_ext, ff)
  else
     call add(bf_field(1), ff)
  end if
  
  call matvecmult(fp_mat, ff%vec, fpf%vec)
  call newsolve(mass_mat_lhs%mat, fpf%vec, ier)
  if(extsubtract.eq.1) then
     call add(bfp_ext, fpf)
  else
     call add(bfp_field(1), fpf)
  end if

  call matvecmult(bz_mat, ff%vec, bzf%vec)
  call newsolve(mass_mat_lhs%mat, bzf%vec, ier)
  if(extsubtract.eq.1) then
     call add(bz_ext, bzf)
  else 
     call add(bz_field(1), bzf)
  end if

  call destroy_field(ff)
  call destroy_field(fpf)
  call destroy_field(bzf)
  call destroy_mat(bz_mat)
  call destroy_mat(fp_mat)
  
end subroutine tf_shift_tilt

!=======================================================
! boundary_rmp
! ~~~~~~~~~~~~
!
!=======================================================
subroutine boundary_rmp(rhs, mat)
  use basic
  use vector_mod
  use matrix_mod
  use boundary_conditions

  implicit none

  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer, parameter :: numvarsm = 2
  integer :: i, izone, izonedim, i_psi, i_f, numnodes, icounter_t
  real :: normal(2), curv(3)
  real :: x, z, phi
  logical :: is_boundary
!!$  real, dimension(1) :: xv, phiv, zv
!!$  vectype, dimension(1) :: br, bphi, bz
!!$  vectype :: bn
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_rmp called"

  temp = 0.

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z)
     if(.not.is_boundary) cycle

     i_psi = node_index(rhs, i, 1)
     i_f = node_index(rhs, i, 2)

     !call set_dirichlet_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)

!!$     if(ifbound.eq.1) then
     !call set_normal_bc(i_f, rhs,temp, normal,curv,izonedim,mat)
     call set_dirichlet_bc(i_f, rhs,temp, normal,curv,izonedim,mat)
!!$     else if(ifbound.eq.2) then
!!$        call set_normal_bc(i_f,rhs,temp,normal,curv,izonedim,mat)
!!$     else if(ifbound.eq.3) then
!!$        xv(1) = x
!!$        phiv(1) = phi
!!$        zv(1) = z
!!$!        call rmp_field(1, 1, 1, xv, phiv, zv, br, bphi, bz)
!!$
!!$        ! set df/dn = -B.n
!!$!        br(1) = 0.
!!$!        bz(1) = 0.
!!$        bn = -(br(1)*normal(1) + bz(1)*normal(2))
!!$        call set_np_bc(i_f,rhs,bn,izonedim,mat)
!!$     end if
  end do

end subroutine boundary_rmp

  subroutine calculate_field_from_j(ilin)
    use basic
    use math
    use mesh_mod
    use sparse
    use arrays
    use matrix_mod
    use m3dc1_nint
    use newvar_mod
    use boundary_conditions

    implicit none

    include 'mpif.h'

    type(matrix_type) :: jx_mat
    type(vector_type) :: psi_vec
#if defined(USECOMPLEX) || defined(USE3D)
    type(matrix_type) :: bf_mat
    type(vector_type) :: bf_vec, temp_vec, temp_vec2
#endif
    integer :: itri, nelms, ier, ibound, ipsibound, fbound
    integer, intent(in), optional :: ilin
    integer :: il

    vectype, dimension(dofs_per_element,dofs_per_element,2,2) :: temp
    vectype, dimension(dofs_per_element,2) :: temp2
    type(field_type) :: psi_f, fx_f
#if defined(USECOMPLEX) || defined(USE3D)
    type(field_type) :: f_f
    vectype, dimension(dofs_per_element,dofs_per_element) :: temp5
#endif

    if(myrank.eq.0 .and. iprint.ge.2) print *, "Calculating field from J"
    if(myrank.eq.0 .and. iprint.ge.2) print *, "xmag, zmag = ", xmag, zmag

    if(.not.present(ilin)) then
       il = 1
    else
       il = ilin
    end if

    call create_vector(psi_vec,2)
    call associate_field(psi_f,psi_vec,1)
    call associate_field(fx_f, psi_vec,2)

    call set_matrix_index(jx_mat, br_mat_index)
    call create_mat(jx_mat, 2, 2, icomplex, 1)

#if defined(USECOMPLEX) || defined(USE3D)
    call set_matrix_index(bf_mat, bf_mat_index)
    call create_mat(bf_mat, 1, 1, icomplex, 1)
#endif

    psi_f = 0.
    fx_f = 0.

    ! boundary condition on psi
!    ipsibound = BOUNDARY_NONE
    ipsibound = BOUNDARY_NEUMANN
!    ipsibound = BOUNDARY_DIRICHLET

    ! Boundary condition on F*
    !ibound = BOUNDARY_NEUMANN
    ibound = BOUNDARY_DIRICHLET
    !ibound = BOUNDARY_NONE

    fbound = BOUNDARY_NEUMANN

    if(myrank.eq.0 .and. iprint.ge.2) print *, 'calculating field values...'
    nelms = local_elements()

    do itri=1,nelms

     call define_element_quadrature(itri,int_pts_main,5)
     call define_fields(itri,0,1,0)

     call read_j(npoints, &
          x_79, phi_79, z_79, &
          temp79a, temp79b, temp79c)

     ! J = grad(F*) x grad(phi) - del*(psi) grad(phi) + grad_perp(psi') / R^2
     ! F* = R^2 del^2(f)

     ! ( R Jphi )                               ( -del*(nu)               0              )
     ! ( JR     ) = Jx(nu) . ( psi )   Jx(nu) = ( (1/R^2) d^2/dRdPhi(nu) -(1/R) d/dZ(nu) )
     ! ( JZ     )            ( F*  )            ( (1/R^2) d^2/dZdPhi(nu)  (1/R) d/dR(nu) )

     ! Find least squares solution to above equation:
     !            ( R Jphi )
     ! Jx^T(mu) . ( JR     ) = Jx^T(mu) . Jx(nu) . ( psi )
     !            ( JZ     )                       ( F*  )


     ! Populate Jx^T(mu) . Jx(nu)
     temp(:,:,1,1) = intxx2(mu79(:,:,OP_LP),nu79(:,:,OP_LP))
     temp(:,:,1,2) = 0.
     temp(:,:,2,1) = 0.
     temp(:,:,2,2) = intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),ri2_79) &
          +          intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DR),ri2_79)
#if defined(USE3D)
     temp(:,:,1,1) = temp(:,:,1,1) &
                   + intxx3(mu79(:,:,OP_DRP),nu79(:,:,OP_DRP),ri4_79) &
                   + intxx3(mu79(:,:,OP_DZP),nu79(:,:,OP_DZP),ri4_79)
     temp(:,:,1,2) = temp(:,:,1,2) &
          + intxx3(mu79(:,:,OP_DZP),nu79(:,:,OP_DR ), ri3_79) &
          - intxx3(mu79(:,:,OP_DRP),nu79(:,:,OP_DZ ), ri3_79)
     temp(:,:,2,1) = temp(:,:,2,1) &
          + intxx3(mu79(:,:,OP_DR ),nu79(:,:,OP_DZP), ri3_79) &
          - intxx3(mu79(:,:,OP_DZ ),nu79(:,:,OP_DRP), ri3_79)
#elif defined(USECOMPLEX)
     temp(:,:,1,1) = temp(:,:,1,1) &
                   + (intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DR),ri4_79) &
                   +  intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),ri4_79))*(rfac*conjg(rfac))
     temp(:,:,1,2) = temp(:,:,1,2) &
          + (intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DR), ri3_79) &
          -  intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DZ), ri3_79))*conjg(rfac)
     temp(:,:,2,1) = temp(:,:,2,1) &
          + (intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DZ), ri3_79) &
          -  intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DR), ri3_79))*rfac

#endif

     ! Populate Jx^T(mu) . ( R Jphi, JR, JZ )
     temp2(:,1) = -intx3(mu79(:,:,OP_LP),temp79b,r_79)
     temp2(:,2) = &
          -       intx3(mu79(:,:,OP_DZ ),temp79a,ri_79) &
          +       intx3(mu79(:,:,OP_DR ),temp79c,ri_79)
#if defined(USE3D)
     temp2(:,1) = temp2(:,1) &
          +       intx3(mu79(:,:,OP_DRP),temp79a,ri2_79) &
          +       intx3(mu79(:,:,OP_DZP),temp79c,ri2_79)
#elif defined(USECOMPLEX)
     temp2(:,1) = temp2(:,1) &
          +       (intx3(mu79(:,:,OP_DR),temp79a,ri2_79) &
          +        intx3(mu79(:,:,OP_DZ),temp79c,ri2_79))*conjg(rfac)
#endif

!!$     ! Solve for f' from F*'
!!$     ! (F*)' = R^2 Del^2 f'
!!$#if defined(USE3D)
!!$     temp5(:,:) = -(intxx3(mu79(:,:,OP_DP),nu79(:,:,OP_LP),r2_79) &
!!$          +         intxx2(mu79(:,:,OP_DP),nu79(:,:,OP_DPP)))
!!$#elif defined(USECOMPLEX)
!!$     temp5(:,:) = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_LP),r2_79)*rfac &
!!$          +       intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))*rfac**3
!!$#endif

     ! Solve for f from F*
     ! F* = R^2 Del^2 f
#if defined(USE3D)
     temp5(:,:) = (intxx3(mu79(:,:,OP_1),nu79(:,:,OP_LP),r2_79) &
          +        intxx2(mu79(:,:,OP_1),nu79(:,:,OP_DPP)))
#elif defined(USECOMPLEX)
     temp5(:,:) = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_LP),r2_79) &
          +       intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))*rfac**2
#endif

     call apply_boundary_mask(itri, ipsibound, temp(:,:,1,1), &
          tags=BOUND_DOMAIN)
     call apply_boundary_mask(itri, ipsibound, temp(:,:,1,2), &
          tags=BOUND_DOMAIN)
     call apply_boundary_mask(itri, ibound, temp(:,:,2,1), &
          tags=BOUND_DOMAIN)
     call apply_boundary_mask(itri, ibound, temp(:,:,2,2), &
          tags=BOUND_DOMAIN)

     call insert_block(jx_mat, itri, 1, 1, temp(:,:,1,1), MAT_ADD)
     call insert_block(jx_mat, itri, 1, 2, temp(:,:,1,2), MAT_ADD)
     call insert_block(jx_mat, itri, 2, 1, temp(:,:,2,1), MAT_ADD)
     call insert_block(jx_mat, itri, 2, 2, temp(:,:,2,2), MAT_ADD)

     call vector_insert_block(psi_vec, itri, 1, temp2(:,1), MAT_ADD)
     call vector_insert_block(psi_vec, itri, 2, temp2(:,2), MAT_ADD)

#if defined(USECOMPLEX) || defined(USE3D)
     call apply_boundary_mask(itri, fbound, temp5(:,:), &
          tags=BOUND_DOMAIN)

     call insert_block(bf_mat, itri, 1, 1, temp5(:,:), MAT_ADD)
#endif
  end do

  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving for (psi, F*)..."
  call sum_shared(psi_vec)
  call boundary_j(psi_vec,jx_mat)
  call finalize(jx_mat)
  call newsolve(jx_mat,psi_vec,ier)
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving psi: ier = ", ier

  if(extsubtract.eq.1) then
     psi_ext = psi_f
     bz_ext = fx_f
  else
     psi_field(il) = psi_f
     bz_field(il) = fx_f
  end if
  
#if defined(USECOMPLEX) || defined(USE3D)
  ! Solve (F*)' = R^2 Del^2(f) for f

  !  store F* in temp_vec2
  call create_vector(temp_vec2, 1)
  call associate_field(f_f, temp_vec2, 1)
  f_f = fx_f

  !  have f_f point to solution vector (which will contain f)
  call create_vector(temp_vec, 1)
  call associate_field(f_f, temp_vec, 1)
  f_f = 0.

  ! Solve F* = R^2 del^2(f) for f
  call matvecmult(mass_mat_rhs%mat, temp_vec2, temp_vec)
  call boundary_fstar(temp_vec,bf_mat)
  call finalize(bf_mat)
  call newsolve(bf_mat,temp_vec,ier)

  if(extsubtract.eq.1) then
     bf_ext = f_f
  else
     bf_field(il) = f_f
  end if

  ! Solve fp = df/dp
#ifdef USECOMPLEX
  if(extsubtract.eq.1) then
     bfp_ext = bf_ext
     call mult(bfp_ext, rfac)
  else
     bfp_field(il) = bf_field(il)
     call mult(bfp_field(il), rfac)
  end if
#elif defined(USE3D)
  call safestop(88)
#endif

  ! Solve F = F* - f''
#ifdef USECOMPLEX
  if(extsubtract.eq.1) then
     f_f = bfp_ext
     call mult(f_f, -rfac)   
     call add(bz_ext, f_f)
  else
     f_f = bfp_ext
     call mult(f_f, -rfac)   
     call add(bz_field(il), f_f)
  end if
#elif defined(USE3D)
    call safestop(89)
#endif
  
  call destroy_vector(temp_vec)
  call destroy_vector(temp_vec2)
#endif

  call destroy_vector(psi_vec)
  call destroy_mat(jx_mat)

#if defined(USECOMPLEX) || defined(USE3D)
  call destroy_vector(bf_vec)
  call destroy_mat(bf_mat)
#endif

  if(myrank.eq.0 .and. iprint.ge.2) print *, "Done calculate_field_from_j"
end subroutine calculate_field_from_j


!=======================================================
! boundary_j
! ~~~~~~~~~~
!
!=======================================================
subroutine boundary_j(rhs, mat)
  use basic
  use vector_mod
  use matrix_mod
  use boundary_conditions

  implicit none

  type(vector_type) :: rhs
  type(matrix_type), optional :: mat

  integer, parameter :: numvarsm = 2
  integer :: i, izone, izonedim, i_psi, i_f, numnodes, icounter_t
  real :: normal(2), curv(3)
  real :: x, z, phi
  logical :: is_boundary
!!$  real, dimension(1) :: xv, phiv, zv
!!$  vectype, dimension(1) :: br, bphi, bz
!!$  vectype :: bn
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_j called"

  temp = 0.

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z)
     if(.not.is_boundary) cycle

     i_psi = node_index(rhs, i, 1)
     i_f = node_index(rhs, i, 2)

!     temp = 0.
!     call set_dirichlet_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)
     temp = 0.
     call set_normal_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)

     call set_dirichlet_bc(i_f, rhs,temp, normal,curv,izonedim,mat)
  end do

end subroutine boundary_j

!=======================================================
! boundary_fstar
! ~~~~~~~~~~~~~~
!
!=======================================================
subroutine boundary_fstar(rhs, mat)
  use basic
  use vector_mod
  use matrix_mod
  use boundary_conditions

  implicit none

  type(vector_type) :: rhs
  type(matrix_type), optional :: mat

  integer, parameter :: numvarsm = 2
  integer :: i, izone, izonedim, i_f, numnodes, icounter_t
  real :: normal(2), curv(3)
  real :: x, z, phi
  logical :: is_boundary
!!$  real, dimension(1) :: xv, phiv, zv
!!$  vectype, dimension(1) :: br, bphi, bz
!!$  vectype :: bn
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_j called"

  temp = 0.

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z)
     if(.not.is_boundary) cycle

     i_f = node_index(rhs, i, 1)

!!$     temp = 0.
!!$     call set_dirichlet_bc(i_f,rhs,temp,normal,curv,izonedim,mat)
     temp = 0.
     call set_normal_bc(i_f,rhs,temp,normal,curv,izonedim,mat)

  end do

end subroutine boundary_fstar


  subroutine read_j(n, x, phi, z, jr, jphi, jz)
    use basic

    implicit none

    integer, intent(in) :: n
    real, intent(in), dimension(n) :: x, phi, z
    vectype, intent(out), dimension(n) :: jr, jphi, jz

    real, dimension(n) :: r2, ex
    real :: b0, x0, z0, s2

    ! Axisymmetric test case with helical current density
    x0 = xmag
    z0 = zmag
    b0 = 0.5
    s2 = ln**2

    !    r2(:) = (x-xh)**2 + (z-zh)**2
    r2 = (x-x0)**2 + (z-z0)**2
    ex = b0*exp(-r2/(2.*s2))

#if defined(USECOMPLEX)
    jr(:) = ex * &
         (x-x0) * ( 2.*x**2*((z-z0)**2 - s2) - (ntor*s2)**2) &
         / (x**2 * s2)
    jz(:) = ex * &
         (z-z0) * (-2.*x**2*((x-x0)**2 - s2) + (ntor*s2)**2 + 2.*x*(x-x0)*s2) &
         / (x**2 * s2)
    jphi(:) = cmplx(0.,ntor)*ex * &
         (x*((x-x0)**2 - (z-z0)**2) + (x - x0)*s2) &
         / x**2

    ! This corresponds to
    !    BR   = i*ntor*s2*(z-z0)/x  * ex * exp[i ntor phi]
    !    BZ   = i*ntor*s2*(r-r0)/x  * ex * exp[i ntor phi]
    !    BPhi = 2*b0*(x - x0)*(z0-z0)*ex * exp[i ntor phi]
#else
    jr(:)   = ex*( (z-z0)/s2       )
    jz(:)   = ex*(-(x-x0)/s2 + 1./x)
    jphi(:) = ex*(r2/s2 - 1. - x0/x)/x

    ! This corresponds to
    !    BR   = -(z-z0)/x * ex
    !    BZ   =  (r-r0)/x * ex
    !    BPhi =  ex
#endif

  end subroutine read_j



end module rmp


