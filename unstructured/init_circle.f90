!==============================================================================
! Circular magnetic field equilibrium (for parallel conduction tests)
!==============================================================================
module circular_field

  real, private :: j0

contains

subroutine circular_field_init()

    use basic
    use arrays
    use m3dc1_nint
    use mesh_mod
    use newvar_mod
    use field
    use diagnostics

    implicit none

    type(field_type) :: psi_vec, bz_vec, den_vec, p_vec
    integer :: itri, numelms, k
    vectype, dimension(dofs_per_element) :: dofs

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining Circular Equilibrium'

    j0 = 2.*bzero/(q0*rzero)
    xlim = xmag + ln
    zlim = zmag

    if(myrank.eq.0) then
       print *, 'Analytic resonant frequency: ', j0*ntor*(q0 - mpol/ntor)/(2.*sqrt(den0))
    end if
    
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
             call circular_field_equ
          else
             ! calculate perturbed fields
             call circular_field_per
          end if

          ! populate vectors for solves

          ! psi
          dofs = intx2(mu79(:,:,OP_1),ps079(:,OP_1))
          call vector_insert_block(psi_vec%vec,itri,1,dofs,VEC_ADD)

          ! bz
          dofs = intx2(mu79(:,:,OP_1),bz079(:,OP_1))
          call vector_insert_block(bz_vec%vec,itri,1,dofs,VEC_ADD)

          ! p
          dofs = intx2(mu79(:,:,OP_1),p079(:,OP_1))
          call vector_insert_block(p_vec%vec,itri,1,dofs,VEC_ADD)

          ! den
          dofs = intx2(mu79(:,:,OP_1),n079(:,OP_1))
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

    call lcfs(psi_field(0))

end subroutine circular_field_init

subroutine circular_field_equ
  use basic
  use m3dc1_nint

  implicit none

  real, dimension(MAX_PTS) :: r2

  r2 = (x_79-xmag)**2 + (z_79-zmag)**2

  where(r2.lt.ln**2)
     ps079(:,OP_1) = j0*(ln**2 - r2)/4.

     p079(:,OP_1) = j0**2*(ln**2 - r2) + pedge
  elsewhere
     ps079(:,OP_1) = -j0*ln**2*log(r2/ln**2)/4.
     p079(:,OP_1) = pedge
  end where
  bz079(:,OP_1) = bzero
  n079(:,OP_1) = den0

end subroutine circular_field_equ


subroutine circular_field_per
  use basic
  use m3dc1_nint
  use init_common

  implicit none

  call init_random(x_79, phi_79, z_79, ps079(:,OP_1))
  n079(:,OP_1) = 0.
  p079(:,OP_1) = 0.
  bz079(:,OP_1) = 0.

end subroutine circular_field_per

end module circular_field
