module init_kstar

contains

!!$  subroutine kstar_init()
!!$    use basic
!!$    use arrays
!!$    use m3dc1_nint
!!$    use mesh_mod
!!$    use newvar_mod
!!$    use field
!!$    use diagnostics
!!$
!!$    implicit none
!!$
!!$    type(field_type) :: psi_vec, bz_vec, den_vec, p_vec
!!$    integer :: itri, numelms, i, k
!!$    vectype, dimension(dofs_per_element) :: dofs
!!$
!!$    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining RWM Equilibrium'
!!$
!!$    ! this is the radius (squared) where psi_norm = ln
!!$    ! here psi_norm = psi(r) / psi(1)
!!$    r02 = exp(1. - 1./ln)
!!$    xlim = xmag + sqrt(r02)
!!$    zlim = zmag
!!$    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Limiter = ', xlim
!!$    
!!$    call create_field(psi_vec)
!!$    call create_field(bz_vec)
!!$    call create_field(p_vec)
!!$    call create_field(den_vec)
!!$
!!$    numelms = local_elements()
!!$
!!$    do k=0,1
!!$       psi_vec = 0.
!!$       bz_vec = 0.
!!$       p_vec = 0.
!!$       den_vec = 0.
!!$
!!$       do itri=1,numelms
!!$          call define_element_quadrature(itri,int_pts_main,int_pts_tor)
!!$          call define_fields(itri,0,1,0)
!!$          
!!$          if(k.eq.0) then 
!!$             ! calculate equilibrium fields
!!$             call rwm_equ
!!$          else
!!$             ! calculate perturbed fields
!!$             call rwm_per
!!$          end if
!!$
!!$          ! populate vectors for solves
!!$
!!$          ! psi
!!$          do i=1, dofs_per_element
!!$             dofs(i) = int2(mu79(:,OP_1,i),ps079(:,OP_1))
!!$          end do
!!$          call vector_insert_block(psi_vec%vec,itri,1,dofs,VEC_ADD)
!!$
!!$          ! bz
!!$          do i=1, dofs_per_element
!!$             dofs(i) = int2(mu79(:,OP_1,i),bz079(:,OP_1))
!!$          end do
!!$          call vector_insert_block(bz_vec%vec,itri,1,dofs,VEC_ADD)
!!$
!!$          ! p
!!$          do i=1, dofs_per_element
!!$             dofs(i) = int2(mu79(:,OP_1,i),p079(:,OP_1))
!!$          end do
!!$          call vector_insert_block(p_vec%vec,itri,1,dofs,VEC_ADD)
!!$
!!$          ! den
!!$          do i=1, dofs_per_element
!!$             dofs(i) = int2(mu79(:,OP_1,i),n079(:,OP_1))
!!$          end do
!!$          call vector_insert_block(den_vec%vec,itri,1,dofs,VEC_ADD)
!!$       end do
!!$
!!$       ! do solves
!!$       call newvar_solve(psi_vec%vec,mass_mat_lhs)
!!$       psi_field(k) = psi_vec
!!$
!!$       call newvar_solve(bz_vec%vec,mass_mat_lhs)
!!$       bz_field(k) = bz_vec
!!$
!!$       call newvar_solve(p_vec%vec,mass_mat_lhs)
!!$       p_field(k) = p_vec
!!$       pe_field(k) = p_vec
!!$       call mult(pe_field(k), pefac)
!!$       
!!$       call newvar_solve(den_vec%vec,mass_mat_lhs)
!!$       den_field(k) = den_vec
!!$       
!!$    end do
!!$      
!!$    call destroy_field(psi_vec)
!!$    call destroy_field(bz_vec)
!!$    call destroy_field(den_vec)
!!$    call destroy_field(p_vec)
!!$
!!$    call lcfs(psi_field(0))
!!$  end subroutine kstar_init
!!$
!!$  subroutine kstar_equ()
!!$    use basic
!!$    use m3dc1_nint
!!$
!!$    implicit none
!!$
!!$    real, dimension(MAX_PTS) :: r2
!!$
!!$    r2 = (x_79 - xmag)**2 + (z_79 - zmag)**2
!!$
!!$    bz079(:,OP_1) = bzero
!!$    n079(:,OP_1) = den0
!!$    
!!$    where(r2.lt.r02)
!!$       ps079(:,OP_1) = bzero * r2 / (2.*rzero*q0)
!!$       p079(:,OP_1) = pedge + (r02 - r2)*(bzero/(rzero*q0))**2
!!$    elsewhere
!!$       ps079(:,OP_1) = bzero * r02 / (2.*rzero*q0) &
!!$            * (1. + alog(r2/r02))
!!$       p079(:,OP_1) = pedge
!!$    end where
!!$
!!$    if(itor.eq.1) then
!!$       ps079(:,OP_1) = ps079(:,OP_1)*rzero   ! Bphi ~ Grad(psi)/R
!!$       bz079(:,OP_1) = bz079(:,OP_1)*rzero   ! F = R*Bphi
!!$    end if
!!$
!!$  end subroutine kstar_equ
!!$
!!$  subroutine kstar_per()
!!$    use basic
!!$    use math
!!$    use m3dc1_nint
!!$    use init_common
!!$
!!$    implicit none
!!$
!!$    real, dimension(MAX_PTS) :: r2, theta
!!$
!!$    r2 = (x_79 - xmag)**2 + (z_79 - zmag)**2
!!$    theta = atan2(z_79 - zmag, x_79 - xmag)
!!$
!!$    bz079(:,OP_1) = 0.
!!$    n079(:,OP_1) = 0.
!!$    p079(:,OP_1) = 0.
!!$
!!$    call init_random(x_79-xmag, phi_79, z_79, ps079(:,OP_1))
!!$
!!$    where(r2.ge.r02)
!!$       ps079(:,OP_1) = 0.
!!$    end where
!!$    
!!$  end subroutine kstar_per

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
              temp(i,j) = int4(mu79(:,OP_1,i),nu79(:,OP_DR,j),co,r) &
                   +      int4(mu79(:,OP_1,i),nu79(:,OP_DZ,j),sn,r)
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

end module init_kstar
