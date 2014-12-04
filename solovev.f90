module solovev
  implicit none

contains

  subroutine solovev_init()
    use basic
    use arrays
    use m3dc1_nint
    use mesh_mod
    use newvar_mod
    use field

    implicit none

    type(field_type) :: psi_vec, bz_vec, den_vec, p_vec
    integer :: itri, numelms, i
    vectype, dimension(dofs_per_element) :: dofs

    if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining Solov'ev Equilibrium"

    xlim = sqrt(2.*(ln*rzero) + rzero**2)
    zlim = 0.
    if(myrank.eq.0 .and. iprint.ge.1) print *, "Limiter = ", xlim

    call create_field(psi_vec)
    call create_field(bz_vec)
    call create_field(p_vec)
    call create_field(den_vec)

    psi_vec = 0.
    bz_vec = 0.
    p_vec = 0.
    den_vec = 0.

    numelms = local_elements()
    do itri=1,numelms
       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       ! calculate equilibrium fields
       call solovev_equ

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
    psi_field(0) = psi_vec

    call newvar_solve(bz_vec%vec,mass_mat_lhs)
    bz_field(0) = bz_vec

    call newvar_solve(p_vec%vec,mass_mat_lhs)
    p_field(0) = p_vec
    pe_field(0) = p_vec
    call mult(pe_field(0), pefac)

    call newvar_solve(den_vec%vec,mass_mat_lhs)
    den_field(0) = den_vec

    call destroy_field(psi_vec)
    call destroy_field(bz_vec)
    call destroy_field(den_vec)
    call destroy_field(p_vec)

  end subroutine solovev_init

  subroutine solovev_equ()
    use basic
    use m3dc1_nint

    implicit none

    real, dimension(MAX_PTS) :: g

    bz079(:,OP_1) = 1.
    n079(:,OP_1) = den0

    g = (x_79*z_79/elongation)**2 + &
         0.25*(x_79**2 - rzero**2)**2 - (ln*rzero)**2

    ps079(:,OP_1) = elongation/(2.*rzero**3*q0) * g
    where(g .lt. 0)
       p079(:,OP_1) = pedge - (1. + elongation**2) &
            / (elongation * rzero**3*q0) * ps079(:,OP_1)
    elsewhere
       p079(:,OP_1) = pedge
    end where

  end subroutine solovev_equ

  subroutine solovev_per()
    implicit none

    
  end subroutine solovev_per

end module solovev
