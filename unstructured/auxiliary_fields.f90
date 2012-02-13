#ifdef USECOMPLEX
#define CONJUGATE(x) (conjg(x))
#else
#define CONJUGATE(x) (x)
#endif

module auxiliary_fields
  use field

  implicit none

  type(field_type) :: bdotgradp
  type(field_type) :: bdotgradt
  type(field_type) :: torque_density_em
  type(field_type) :: torque_density_ntv
  type(field_type) :: chord_mask

  logical, private :: initialized = .false.

contains

subroutine create_auxiliary_fields
  implicit none

  call create_field(bdotgradp)
  call create_field(bdotgradt)
  call create_field(torque_density_em)
  call create_field(torque_density_ntv)
  call create_field(chord_mask)
  initialized = .true.
end subroutine create_auxiliary_fields

subroutine destroy_auxiliary_fields
  implicit none

  if(.not.initialized) return
  call destroy_field(bdotgradp)
  call destroy_field(bdotgradt)
  call destroy_field(torque_density_em)
  call destroy_field(torque_density_ntv)
  call destroy_field(chord_mask)
end subroutine destroy_auxiliary_fields
  
subroutine calculate_temperatures(ilin, te, ti)
  use math
  use basic
  use m3dc1_nint
  use newvar_mod
  use diagnostics
  use metricterms_new

  implicit none

  type(field_type) :: te, ti
  integer, intent(in) :: ilin

  integer :: def_fields
  integer :: numelms
  integer :: i, itri

  vectype, dimension(dofs_per_element) :: dofs

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating temperatures'

  te = 0.
  ti = 0.

  ! specify which primitive fields are to be evalulated
  def_fields = FIELD_N + FIELD_NI + FIELD_P + FIELD_PE

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, def_fields, 1, 0)

     ! electron temperature
     if(ilin.eq.1) then
        do i=1, dofs_per_element
           dofs(i) = int3(mu79(:,OP_1,i),pe179(:,OP_1),nei79(:,OP_1)) &
                - int5(mu79(:,OP_1,i),ne179(:,OP_1),pe079(:,OP_1),nei79(:,OP_1),nei79(:,OP_1))
        end do
     else
        do i=1, dofs_per_element
           dofs(i) = int3(mu79(:,OP_1,i),pet79(:,OP_1),nei79(:,OP_1))
        end do
     end if
     call vector_insert_block(te%vec,itri,1,dofs,VEC_ADD)

     ! ion temperature
     if(ilin.eq.1) then 
        temp79a = p179(:,OP_1) - pe179(:,OP_1)
        temp79b = p079(:,OP_1) - pe079(:,OP_1)
        do i=1, dofs_per_element
           dofs(i) = int3(mu79(:,OP_1,i),temp79a,ni79(:,OP_1)) &
                - int5(mu79(:,OP_1,i),n179(:,OP_1),temp79b,ni79(:,OP_1),ni79(:,OP_1))
        end do
     else
        temp79a = pt79(:,OP_1) - pet79(:,OP_1)
        do i=1, dofs_per_element
           dofs(i) = int3(mu79(:,OP_1,i),temp79a,nei79(:,OP_1))
        end do
     end if
     call vector_insert_block(ti%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(te%vec, mass_mat_lhs)
  call newvar_solve(ti%vec, mass_mat_lhs)

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, ' Done calculating temperatures'
  
  end subroutine calculate_temperatures

subroutine calculate_auxiliary_fields(ilin)
  use math
  use basic
  use m3dc1_nint
  use newvar_mod
  use diagnostics
  use metricterms_new

  implicit none

  integer, intent(in) :: ilin

  integer :: def_fields
  integer :: numelms
  integer :: i, itri

  vectype, dimension(dofs_per_element) :: dofs

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating diagnostic fields'

  bdotgradp = 0.
  bdotgradt = 0.
  torque_density_em = 0.
  torque_density_ntv = 0.
  chord_mask = 0.

  ! specify which primitive fields are to be evalulated
  def_fields = FIELD_N + FIELD_NI + FIELD_P + FIELD_PSI + FIELD_I
  if(amupar.ne.0) then
     def_fields = def_fields + FIELD_MU + FIELD_B2I
     def_fields = def_fields + FIELD_PHI + FIELD_V + FIELD_CHI
  end if

  numelms = local_elements()

  ! EM Torque density
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, def_fields, 1, 0)

     ! magnetic torque_density (ignoring toroidal magnetic pressure gradient)
#ifdef USECOMPLEX
     do i=1, dofs_per_element
        dofs(i) = &
             + int4(ri_79,mu79(:,OP_1,i),ps179(:,OP_DR),conjg(bz179(:,OP_DZ)))&
             - int4(ri_79,mu79(:,OP_1,i),ps179(:,OP_DZ),conjg(bz179(:,OP_DR)))&
             - int3(mu79(:,OP_1,i),bf179(:,OP_DRP),conjg(bz179(:,OP_DR))) &
             - int3(mu79(:,OP_1,i),bf179(:,OP_DZP),conjg(bz179(:,OP_DZ))) &
             + int4(ri_79,mu79(:,OP_1,i),conjg(ps179(:,OP_DR)),bz179(:,OP_DZ))&
             - int4(ri_79,mu79(:,OP_1,i),conjg(ps179(:,OP_DZ)),bz179(:,OP_DR))&
             - int3(mu79(:,OP_1,i),conjg(bf179(:,OP_DRP)),bz179(:,OP_DR)) &
             - int3(mu79(:,OP_1,i),conjg(bf179(:,OP_DZP)),bz179(:,OP_DZ))
     end do
     dofs = dofs / 2.
#else
     do i=1, dofs_per_element
        dofs(i) = &
             + int4(ri_79,mu79(:,OP_1,i),pst79(:,OP_DR),bzt79(:,OP_DZ)) &
             - int4(ri_79,mu79(:,OP_1,i),pst79(:,OP_DZ),bzt79(:,OP_DR))
#ifdef USE3D
        dofs(i) = dofs(i) &
             - int3(mu79(:,OP_1,i),bft79(:,OP_DRP),bzt79(:,OP_DR)) &
             - int3(mu79(:,OP_1,i),bft79(:,OP_DZP),bzt79(:,OP_DZ))
#endif
     end do
#endif
     call vector_insert_block(torque_density_em%vec,itri,1,dofs,VEC_ADD)

     ! NTV torque_density
     if(amupar.ne.0.) then
        do i=1, dofs_per_element
           call PVV2(mu79(:,:,i),temp79f)
!           dofs(i) = int1(temp79f)
           dofs(i) = 0.

           if(ilin.eq.1) then
              call PVS1psipsi(CONJUGATE(ph179),ps179,ps079,temp79c)
              call PVS1psib  (CONJUGATE(ph179),ps179,bz079,temp79d)
              call PVS1bb    (CONJUGATE(ph179),bz179,bz079,temp79e)
              temp79a = temp79b + temp79c + temp79d + temp79e
              dofs(i) = int3(vip79(:,OP_1),temp79a,temp79f)
              call PVS1psipsi(CONJUGATE(ph179),ps079,ps179,temp79c)
              call PVS1psib  (CONJUGATE(ph179),ps079,bz179,temp79d)
              call PVS1bb    (CONJUGATE(ph179),bz079,bz179,temp79e)
              temp79a = temp79b + temp79c + temp79d + temp79e
              dofs(i) = int3(vip79(:,OP_1),temp79a,temp79f)
              call PVS1psipsi(ph079,CONJUGATE(ps179),ps179,temp79c)
              call PVS1psib  (ph079,CONJUGATE(ps179),bz179,temp79d)
              call PVS1bb    (ph079,CONJUGATE(bz179),bz179,temp79e)
              temp79a = temp79b + temp79c + temp79d + temp79e
              dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)

              if(numvar.ge.2) then
                 call PVS2psipsi(CONJUGATE(vz179),ps179,ps079,temp79c)
                 call PVS2psib  (CONJUGATE(vz179),ps179,bz079,temp79d)
                 call PVS2bb    (CONJUGATE(vz179),bz179,bz079,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
                 call PVS2psipsi(CONJUGATE(vz179),ps079,ps179,temp79c)
                 call PVS2psib  (CONJUGATE(vz179),ps079,bz179,temp79d)
                 call PVS2bb    (CONJUGATE(vz179),bz079,bz179,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
                 call PVS2psipsi(vz079,CONJUGATE(ps179),ps179,temp79c)
                 call PVS2psib  (vz079,CONJUGATE(ps179),bz179,temp79d)
                 call PVS2bb    (vz079,CONJUGATE(bz179),bz179,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
              endif
              
              if(numvar.ge.3) then
                 call PVS3psipsi(CONJUGATE(ch179),ps179,ps079,temp79c)
                 call PVS3psib  (CONJUGATE(ch179),ps179,bz079,temp79d)
                 call PVS3bb    (CONJUGATE(ch179),bz179,bz079,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
                 call PVS3psipsi(CONJUGATE(ch179),ps079,ps179,temp79c)
                 call PVS3psib  (CONJUGATE(ch179),ps079,bz179,temp79d)
                 call PVS3bb    (CONJUGATE(ch179),bz079,bz179,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
                 call PVS3psipsi(ch079,CONJUGATE(ps179),ps179,temp79c)
                 call PVS3psib  (ch079,CONJUGATE(ps179),bz179,temp79d)
                 call PVS3bb    (ch079,CONJUGATE(bz179),bz179,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
              endif
           else
              call PVS1      (pht79,temp79b)
              call PVS1psipsi(pht79,pst79,pst79,temp79c)
              call PVS1psib  (pht79,pst79,bzt79,temp79d)
              call PVS1bb    (pht79,bzt79,bzt79,temp79e)
              temp79a = temp79b + temp79c + temp79d + temp79e
              dofs(i) = int3(vip79(:,OP_1),temp79a,temp79f)

              if(numvar.ge.2) then
                 call PVS2      (vzt79,temp79b)
                 call PVS2psipsi(vzt79,pst79,pst79,temp79c)
                 call PVS2psib  (vzt79,pst79,bzt79,temp79d)
                 call PVS2bb    (vzt79,bzt79,bzt79,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
              endif
              
              if(numvar.ge.3) then
                 call PVS3      (cht79,temp79b)
                 call PVS3psipsi(cht79,pst79,pst79,temp79c)
                 call PVS3psib  (cht79,pst79,bzt79,temp79d)
                 call PVS3bb    (cht79,bzt79,bzt79,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
              endif
           endif
        end do
     else
        dofs = 0.
     endif
     call vector_insert_block(torque_density_ntv%vec,itri,1,dofs,VEC_ADD)


     ! b dot grad p
     do i=1, dofs_per_element


        if(ilin.eq.1) then 
           dofs(i) = int4(ri_79,mu79(:,OP_1,i),p179(:,OP_DZ),ps079(:,OP_DR)) &
                -    int4(ri_79,mu79(:,OP_1,i),p179(:,OP_DR),ps079(:,OP_DZ)) &
                +    int4(ri_79,mu79(:,OP_1,i),p079(:,OP_DZ),ps179(:,OP_DR)) &
                -    int4(ri_79,mu79(:,OP_1,i),p079(:,OP_DR),ps179(:,OP_DZ))
#if defined(USECOMPLEX) || defined(USE3D)
           dofs(i) = dofs(i) &
                + int4(ri2_79,mu79(:,OP_1,i),p179(:,OP_DP),bz079(:,OP_1)) &
                - int3(mu79(:,OP_1,i),p079(:,OP_DZ),bf179(:,OP_DZP)) &
                - int3(mu79(:,OP_1,i),p079(:,OP_DR),bf179(:,OP_DRP))
#endif
        else
           dofs(i) = int4(ri_79,mu79(:,OP_1,i),pt79(:,OP_DZ),pst79(:,OP_DR)) &
                -    int4(ri_79,mu79(:,OP_1,i),pt79(:,OP_DR),pst79(:,OP_DZ))
#if defined(USECOMPLEX) || defined(USE3D)
           dofs(i) = dofs(i) &
                + int4(ri2_79,mu79(:,OP_1,i),pt79(:,OP_DP),bzt79(:,OP_1)) &
                - int3(mu79(:,OP_1,i),pt79(:,OP_DZ),bft79(:,OP_DZP)) &
                - int3(mu79(:,OP_1,i),pt79(:,OP_DR),bft79(:,OP_DRP))
#endif
        end if
     end do
     call vector_insert_block(bdotgradp%vec,itri,1,dofs,VEC_ADD)

     ! b dot grad T
     do i=1, dofs_per_element
        temp79a = mu79(:,OP_1,i)*ni79(:,OP_1)
        
        if(ilin.eq.1) then 
           dofs(i) = &
                +int4(ri_79,mu79(:,OP_1,i),p179(:,OP_DZ),ps079(:,OP_DR)) &
                -int4(ri_79,mu79(:,OP_1,i),p179(:,OP_DR),ps079(:,OP_DZ)) &
                +int4(ri_79,mu79(:,OP_1,i),p079(:,OP_DZ),ps179(:,OP_DR)) &
                -int4(ri_79,mu79(:,OP_1,i),p079(:,OP_DR),ps179(:,OP_DZ)) &
                -int5(ri_79,temp79a,n179(:,OP_DZ),ps079(:,OP_DR),p079(:,OP_1))&
                +int5(ri_79,temp79a,n179(:,OP_DR),ps079(:,OP_DZ),p079(:,OP_1))&
                -int5(ri_79,temp79a,n079(:,OP_DZ),ps179(:,OP_DR),p079(:,OP_1))&
                +int5(ri_79,temp79a,n079(:,OP_DR),ps179(:,OP_DZ),p079(:,OP_1))&
                -int5(ri_79,temp79a,n079(:,OP_DZ),ps079(:,OP_DR),p179(:,OP_1))&
                +int5(ri_79,temp79a,n079(:,OP_DR),ps079(:,OP_DZ),p179(:,OP_1))
#if defined(USECOMPLEX) || defined(USE3D)
           dofs(i) = dofs(i) &
                +int4(ri2_79,mu79(:,OP_1,i),p179(:,OP_DP),bz079(:,OP_1)) &
                -int3(mu79(:,OP_1,i),p079(:,OP_DZ),bf179(:,OP_DZP)) &
                -int3(mu79(:,OP_1,i),p079(:,OP_DR),bf179(:,OP_DRP)) &
                -int5(ri2_79,temp79a,n179(:,OP_DP),bz079(:,OP_1),p079(:,OP_1))&
                +int4(temp79a,n079(:,OP_DZ),bf179(:,OP_DZP),p079(:,OP_1)) &
                +int4(temp79a,n079(:,OP_DR),bf179(:,OP_DRP),p079(:,OP_1))
#endif
        else
           dofs(i) = &
                +int4(ri_79,mu79(:,OP_1,i),pt79(:,OP_DZ),pst79(:,OP_DR)) &
                -int4(ri_79,mu79(:,OP_1,i),pt79(:,OP_DR),pst79(:,OP_DZ)) &
                -int5(ri_79,temp79a,nt79(:,OP_DZ),pst79(:,OP_DR),pt79(:,OP_1))&
                +int5(ri_79,temp79a,nt79(:,OP_DR),pst79(:,OP_DZ),pt79(:,OP_1))
#if defined(USECOMPLEX) || defined(USE3D)
           dofs(i) = dofs(i) &
                +int4(ri2_79,mu79(:,OP_1,i),pt79(:,OP_DP),bzt79(:,OP_1)) &
                -int3(mu79(:,OP_1,i),pt79(:,OP_DZ),bft79(:,OP_DZP)) &
                -int3(mu79(:,OP_1,i),pt79(:,OP_DR),bft79(:,OP_DRP)) &
                -int5(ri2_79,temp79a,nt79(:,OP_DP),bzt79(:,OP_1),pt79(:,OP_1))&
                +int4(temp79a,n079(:,OP_DZ),bft79(:,OP_DZP),pt79(:,OP_1)) &
                +int4(temp79a,n079(:,OP_DR),bft79(:,OP_DRP),pt79(:,OP_1))
#endif
        end if
     end do
     call vector_insert_block(bdotgradt%vec,itri,1,dofs,VEC_ADD)


     ! x-ray detector signal
     if(xray_detector_enabled.eq.1) then
        do i=1, dofs_per_element
           call get_chord_mask(xray_r0, xray_phi0*pi/180., xray_z0, &
                x_79, phi_79, z_79, npoints, &
                xray_theta*pi/180., xray_sigma*pi/180., temp79a)
           dofs(i) = int2(mu79(:,OP_1,i),temp79a)
        end do
        call vector_insert_block(chord_mask%vec,itri,1,dofs,VEC_ADD)
     end if

  end do

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before bdotgradp solve'
  call newvar_solve(bdotgradp%vec, mass_mat_lhs)
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before bdotgradt solve'
  call newvar_solve(bdotgradt%vec, mass_mat_lhs)
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before torque_density solve'
  call newvar_solve(torque_density_em%vec, mass_mat_lhs)
  call newvar_solve(torque_density_ntv%vec, mass_mat_lhs)

  if(xray_detector_enabled.eq.1) then
     call newvar_solve(chord_mask%vec, mass_mat_lhs)
  end if

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Done calculating diagnostic fields'
  
  end subroutine calculate_auxiliary_fields


end module auxiliary_fields
