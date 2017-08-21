module rwm
  implicit none

  real, private :: r02

contains

  subroutine rwm_init()
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

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining RWM Equilibrium'

    ! this is the radius (squared) where psi_norm = ln
    ! here psi_norm = psi(r) / psi(1)
    r02 = exp(1. - 1./ln)
    xlim = xmag + sqrt(r02)
    zlim = zmag
    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Limiter = ', xlim
    
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
             call rwm_equ
          else
             ! calculate perturbed fields
             call rwm_per
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

    vz_field(0) = p_field(0)
    call add(vz_field(0), -pedge)
    call mult(vz_field(0), alpha0)
      
    call destroy_field(psi_vec)
    call destroy_field(bz_vec)
    call destroy_field(den_vec)
    call destroy_field(p_vec)

    call lcfs(psi_field(0))
  end subroutine rwm_init

  subroutine rwm_equ()
    use basic
    use m3dc1_nint

    implicit none

    real, dimension(MAX_PTS) :: r2

    r2 = (x_79 - xmag)**2 + (z_79 - zmag)**2

    bz079(:,OP_1) = bzero
    n079(:,OP_1) = den0
    
    where(r2.lt.r02)
       ps079(:,OP_1) = bzero * r2 / (2.*rzero*q0)
       p079(:,OP_1) = pedge + (r02 - r2)*(bzero/(rzero*q0))**2
    elsewhere
       ps079(:,OP_1) = bzero * r02 / (2.*rzero*q0) &
            * (1. + alog(r2/r02))
       p079(:,OP_1) = pedge
    end where

    if(itor.eq.1) then
       ps079(:,OP_1) = ps079(:,OP_1)*rzero   ! Bphi ~ Grad(psi)/R
       bz079(:,OP_1) = bz079(:,OP_1)*rzero   ! F = R*Bphi
    end if

  end subroutine rwm_equ

  subroutine rwm_per()
    use basic
    use math
    use m3dc1_nint
    use init_common

    implicit none

    real, dimension(MAX_PTS) :: r2, theta

    r2 = (x_79 - xmag)**2 + (z_79 - zmag)**2
    theta = atan2(z_79 - zmag, x_79 - xmag)

    bz079(:,OP_1) = 0.
    n079(:,OP_1) = 0.
    p079(:,OP_1) = 0.

    call init_random(x_79-xmag, phi_79, z_79, ps079(:,OP_1))
!!$    call init_random(-(x_79-xmag), phi_79, z_79, ps179(:,OP_1))
!!$    ps079(:,OP_1) = (ps079(:,OP_1) + ps179(:,OP_1))/2.

!!$#ifdef USECOMPLEX
!!$    ps079(:,OP_1) = exp((0.,1.)*(mpol*theta+pi/2.))
!!$#elif
!!$    ps079(:,OP_1) = real(exp((0.,1.)*(mpol*theta+pi/2.)+rfac*phi))
!!$#endif
!!$
!!$    ps079(:,OP_1) = ps079(:,OP_1)*eps*r2*(r02 - r2)**2

    where(r2.ge.r02)
       ps079(:,OP_1) = 0.
    end where
    
  end subroutine rwm_per

end module rwm
