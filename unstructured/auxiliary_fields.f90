module auxiliary_fields
  use field

  implicit none

  type(field_type) :: bdotgradp
  type(field_type) :: bdotgradt

contains

subroutine create_auxiliary_fields
  implicit none

  call create_field(bdotgradp)
  call create_field(bdotgradt)
end subroutine create_auxiliary_fields

subroutine destroy_auxiliary_fields
  implicit none

  call destroy_field(bdotgradp)
  call destroy_field(bdotgradt)
end subroutine destroy_auxiliary_fields

  
subroutine calculate_auxiliary_fields(ilin)
  use basic
  use m3dc1_nint
  use newvar_mod

  implicit none

  integer, intent(in) :: ilin

  integer :: def_fields
  integer :: numelms
  integer :: i, itri

  vectype, dimension(dofs_per_element) :: dofs

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating auxiliary fields'
#if defined(USE3D)
!
!....scj added return because this was failing on hopper and STIX    (06/21/11)
      return
#endif

  bdotgradp = 0.
  bdotgradt = 0.

  ! specify which primitive fields are to be evalulated
  def_fields = FIELD_N + FIELD_NI + FIELD_P + FIELD_PSI + FIELD_I

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, def_fields, 1, 0)

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
  end do

  call newvar_solve(bdotgradp%vec, mass_mat_lhs)
  call newvar_solve(bdotgradt%vec, mass_mat_lhs)
  
  end subroutine calculate_auxiliary_fields

end module auxiliary_fields
