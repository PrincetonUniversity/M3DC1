module model

  use vector_mod
  use matrix_mod

  integer, allocatable :: global_dof_ids_1(:)
  integer, allocatable :: global_dof_ids_row(:), global_dof_ids_col(:)

  type(vector_type), target :: q4_vec, r4_vec, qp4_vec, qn4_vec

  ! matrices
  type(matrix_type), target :: s1_mat, d1_mat, q1_mat, r14_mat
  type(matrix_type), target :: o1_mat, p1_mat
  type(matrix_type), target :: q42_mat, r42_mat
  type(matrix_type), target :: s2_mat, d2_mat, r2_mat, q2_mat, o2_mat, o3_mat
  type(matrix_type), target :: s8_mat, d8_mat, r8_mat, q8_mat
  type(matrix_type), target :: s9_mat, d9_mat, r9_mat, q9_mat, o9_mat
  type(matrix_type), target :: s11_mat, d11_mat, s12_mat, d12_mat

  ! positions
  integer :: u_i, vz_i, chi_i
  integer :: psi_i, bz_i, pe_i
  integer :: den_i, p_i
  integer :: bf_i, e_i
  integer :: te_i, ti_i

  ! the offset (relative to the node offset) of the named field within
  ! their respective vectors
  integer :: u_off, vz_off, chi_off
  integer :: psi_off, bz_off, pe_off
  integer :: den_off, p_off
  integer :: bf_off, e_off

contains

subroutine calc_ni(ni_field, n0_field, n1_field)
!
!....not sure what this routine is supposed to do, but as far as I can see, it is not called.
!

  use field
  use mesh_mod

  implicit none

  type(field_type), intent(inout) :: ni_field
  type(field_type), intent(in) :: n0_field, n1_field

  integer :: numnodes, inode
  vectype, dimension(dofs_per_node) :: ninv, n0, n1

  numnodes = owned_nodes()

  do inode=1,numnodes
     call get_node_data(n0_field, inode, n0)
     call get_node_data(n1_field, inode, n1)

     ninv(1) = -n1(1)/n0(1)**2
     ninv(2) = &
          -(n1(2)*n0(1) - 2.*n1(1)*n0(2))/n0(1)**3
     ninv(3) = &
          -(n1(3)*n0(1) - 2.*n1(1)*n0(3))/n0(1)**3
     ninv(4) = &
          -n1(4)/n0(1)**2 &
          +2.*(2.*n1(2)*n0(2) &
          +   n1(1)*n0(4))/n0(1)**3 &
          -6.*n1(1)*n0(2)**2/n0(1)**4
     ninv(5) = &
          -n1(5)/n0(1)**2 &
          +2.*(n1(2)*n0(3) + n1(3)*n0(2) &
          +n1(1)*n0(5))/n0(1)**3 &
          -6.*n1(1)*n0(2)*n0(3)/n0(1)**4
     ninv(6) = &
          -n1(6)/n0(1)**2 &
          +2.*(2.*n1(3)*n0(3) &
          +   n1(1)*n0(6))/n0(1)**3 &
          -6.*n1(1)*n0(3)**2/n0(1)**4

#if defined(USE3D)
!   need coding for 7-12
#endif

     call set_node_data(ni_field, inode, ninv)
  end do
end subroutine calc_ni

subroutine calc_b2i(b2i_field, psi0_field, psi1_field, b0_field, b1_field)
  use basic
  use field
  implicit none

  type(field_type), intent(inout) :: b2i_field
  type(field_type), intent(in) :: psi0_field, psi1_field, b0_field, b1_field

  vectype, dimension(dofs_per_node) :: vec, psi0, psi1, b0, b1
  vectype, dimension(dofs_per_node) :: b2i
  vectype :: b20
  real :: x, phi, z
  integer :: inode, numnodes

  numnodes = owned_nodes()

  do inode=1,numnodes
     call get_node_data(psi0_field, inode, psi0)
     call get_node_data(psi1_field, inode, psi1)
     call get_node_data(b0_field, inode, b0)
     call get_node_data(b1_field, inode, b1)

     if(itor.eq.0) then
        x = 1.
     else
        call get_node_pos(inode, x, phi, z)
     endif

     b20 = psi0(2)**2 + psi0(3)**2 + b0(1)**2

     b2i(1) = x**2/b20
     b2i(2) = 2.*x/b20 &
          - (2.*x**2/b20**2)*(psi0(2)*psi0(4) + psi0(3)*psi0(5) + b0(1)*b0(2))
     b2i(3) = &
          - (2.*x**2/b20**2)*(psi0(2)*psi0(5) + psi0(3)*psi0(6) + b0(1)*b0(3))
     b2i(4:6) = 0.
     
     vec(1) = -2.*b2i(1)**2/x**2 * &
          (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3))
     vec(2) = -2.*b2i(2)**2/x**2 * &
          (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3)) &
          -2.*b2i(1)**2/x**2 * &
          (b0(2)*b1(1) + psi0(4)*psi1(2) + psi0(5)*psi1(3)  &
          +b0(1)*b1(2) + psi0(2)*psi1(4) + psi0(3)*psi1(5)) &
          +4.*b2i(1)**2/x**3 * &
          (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3))
     vec(3) = -2.*b2i(3)**2/x**2 * &
          (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3)) &
          -2.*b2i(1)**2/x**2 * &
          (b0(3)*b1(1) + psi0(5)*psi1(2) + psi0(6)*psi1(3)  &
          +b0(1)*b1(3) + psi0(2)*psi1(5) + psi0(3)*psi1(6))
     vec(4:6) = 0.

     call set_node_data(b2i_field, inode, vec)
  end do
end subroutine calc_b2i


subroutine get_den_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inograd_n.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(iconst_n.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_den_mask

subroutine get_temp_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inograd_t.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(iconst_t.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_temp_mask

subroutine get_pres_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inograd_p.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(iconst_p.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_pres_mask


subroutine get_flux_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  if(eta_wall.eq.0.) then 
     ibound = BOUNDARY_DIRICHLET
  else
     ibound = BOUNDARY_RESISTIVE_WALL
  end if
  if(inocurrent_tor.eq.1) ibound = ior(ibound, BOUNDARY_LAPLACIAN)
  if(inocurrent_norm.eq.1) then
     if(i3d.eq.1) then
        ibound = ior(ibound, BOUNDARY_NEUMANNP)
     endif
  endif

  call get_boundary_mask(itri, ibound, imask)
end subroutine get_flux_mask

subroutine get_bz_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
!  if(eta_wall.ne.0.) ibound = ior(ibound, BOUNDARY_RESISTIVE_WALL)
  if(inocurrent_pol.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(inocurrent_norm.eq.1 .or. iconst_bz.eq.1) &
       ibound = ior(ibound, BOUNDARY_DIRICHLET)
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_bz_mask

subroutine get_q_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = BOUNDARY_DIRICHLET
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_q_mask

subroutine get_bf_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = BOUNDARY_DIRICHLET
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_bf_mask

subroutine get_vor_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inonormalflow.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  if(inoslip_pol.eq.1)   ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(vor_bc.eq.1)        ibound = ior(ibound, BOUNDARY_LAPLACIAN)
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_vor_mask

subroutine get_vz_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inoslip_tor.eq.1)   ibound = ior(ibound, BOUNDARY_DIRICHLET)
  if(inostress_tor.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_vz_mask

subroutine get_chi_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inonormalflow.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(inoslip_pol.eq.1)   ibound = ior(ibound, BOUNDARY_DIRICHLET)
  if(com_bc.eq.1)        ibound = ior(ibound, BOUNDARY_LAPLACIAN)
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_chi_mask




!=======================================================
! boundary_vel
! ~~~~~~~~~~~~
!
! sets boundary conditions for velocity fields
!=======================================================
subroutine boundary_vel(rhs, u_v, vz_v, chi_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions

  implicit none
  
  type(vector_type), intent(inout) :: rhs
  type(field_type), intent(in) :: u_v, vz_v, chi_v
  type(matrix_type), optional :: mat

  vectype, dimension(dofs_per_node) :: temp 
  real :: normal(2), curv, x, z
  integer :: i, izone, izonedim, numnodes
  integer :: i_u, i_vz, i_chi
  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vel called"

  numnodes = owned_nodes()
  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_u = node_index(u_v, i)
     if(numvar.ge.2) i_vz = node_index(vz_v, i)
     if(numvar.ge.3) i_chi = node_index(chi_v, i)
       
     ! no normal flow
     if(inonormalflow.eq.1) then
        temp = 0.
        call set_dirichlet_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
        if(numvar.ge.3) then
           call set_normal_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
        endif
     end if
     
     ! no poloidal slip
     if(inoslip_pol.eq.1) then
        temp = 0.
        call set_normal_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
        if(numvar.ge.3) then
           call set_dirichlet_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
        endif
     end if

     ! toroidal velocity
     if(numvar.ge.2) then
        ! no slip
        if(inoslip_tor.eq. 1) then
           call get_node_data(vz_field(1), i, temp)
           call set_dirichlet_bc(i_vz,rhs,temp,normal,curv,izonedim,mat)
        end if
        
        ! no toroidal stress
        if(inostress_tor.eq.1) then
           temp = 0.
           call set_normal_bc(i_vz,rhs,temp,normal,curv,izonedim,mat)
        end if
     endif
       
     ! no vorticity
     if(vor_bc.eq.1) then
        temp = 0.
        call set_laplacian_bc(i_u,rhs,temp,normal,curv,izonedim,-x,mat)
     endif

     ! no compression
     if(com_bc.eq.1 .and. numvar.ge.3) then
        select case(ivform)
        case(0)
           call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,x,mat)
        case(1)
           call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,-x,mat)
        end select
     endif
  end do
end subroutine boundary_vel


!=======================================================
! boundary_vpol
! ~~~~~~~~~~~~~
!
! sets boundary conditions for poloidal velocity fields
!=======================================================
subroutine boundary_vpol(rhs, u_v, chi_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions

  implicit none
  
  type(vector_type), intent(inout) :: rhs
  type(field_type), intent(in) :: u_v, chi_v
  type(matrix_type), optional :: mat

  vectype, dimension(dofs_per_node) :: temp 
  real :: normal(2), curv, x, z
  integer :: i, izone, izonedim, numnodes
  integer :: i_u, i_chi
  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vpol called"

  numnodes = owned_nodes()
  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_u = node_index(u_v, i)
     if(numvar.ge.3) i_chi = node_index(chi_v, i)
       
     ! no normal flow
     if(inonormalflow.eq.1) then
        temp = 0.
        call set_dirichlet_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
        if(numvar.ge.3) then
           call set_normal_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
        endif
     end if
     
     ! no poloidal slip
     if(inoslip_pol.eq.1) then
        temp = 0.
        call set_normal_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
        if(numvar.ge.3) then
           call set_dirichlet_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
        endif
     end if
      
     ! no vorticity
     if(vor_bc.eq.1) then
        temp = 0.
        call set_laplacian_bc(i_u,rhs,temp,normal,curv,izonedim,-x,mat)
     endif

     ! no compression
     if(com_bc.eq.1 .and. numvar.ge.3) then
        select case(ivform)
        case(0)
           call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,x,mat)
        case(1)
           call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,-x,mat)
        end select
     endif
  end do
end subroutine boundary_vpol



!=======================================================
! boundary_mag
! ~~~~~~~~~~~~
!
! sets boundary conditions for magnetic fields
! and electron pressure 
!=======================================================
subroutine boundary_mag(rhs, psi_v, bz_v, bf_v, e_v, mat)
  use math
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions

  implicit none
  
  type(vector_type), intent(inout) :: rhs
  type(field_type) :: psi_v, bz_v, bf_v, e_v
  type(matrix_type), optional :: mat

  vectype, dimension(dofs_per_node) :: temp, temp2, temp3
  real :: normal(2), curv, x, z
  integer :: i, izone, izonedim,  numnodes
  integer :: i_psi, i_bz, i_pe, i_e, i_bf
  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_mag called"

  numnodes = owned_nodes()
  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_psi = node_index(psi_v, i)
     if(numvar.ge.2) i_bz = node_index(bz_v, i)
!     if(numvar.ge.3 .and. ipressplit.eq.0) i_pe = node_index(pe_v, i)
     if(jadv.eq.0 .and. i3d.eq.1) i_e = node_index(e_v, i)
     if(imp_bf.eq.1) i_bf = node_index(bf_v, i)

     ! constant normal field = -n.grad(psi)/R - n.grad(f')
     if(iconst_bn.eq.1) then
        call get_node_data(psi_field(1), i, temp)
        ! add loop voltage
        if(igauge.eq.0) temp(1) = temp(1) + dt*vloop/twopi
        call set_dirichlet_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)
     end if

     ! no toroidal current = -Delta*(psi)/R
     if(inocurrent_tor.eq.1) then       
        temp = 0.
        call set_laplacian_bc(i_psi,rhs,temp,normal,curv,izonedim,-x,mat)
     end if

     ! no tangential current = n.Grad(F)/R + t.Grad(psi')/R^2
     if(inocurrent_pol.eq.1 .and. numvar.ge.2) then
        call get_node_data(bz_field(1), i, temp)
        call set_normal_bc(i_bz,rhs,temp,normal,curv,izonedim,mat)
     end if

     ! no normal current = n.Grad(psi')/R^2 - t.Grad(F)/R
     if(inocurrent_norm.eq.1) then
        if(i3d.eq.1) then
!        if(numvar.ge.2) then
!           call get_node_data(bz_field(1), i, temp2)
!        else
!           temp2 = 0.
!        endif
           call get_node_data(psi_field(1), i, temp)
           call set_normalp_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)
        endif

        if(numvar.ge.2) then
           call get_node_data(bz_field(1), i, temp)
           call set_dirichlet_bc(i_bz,rhs,temp,normal,curv,izonedim,mat)
        endif
     else if(iconst_bz.eq.1 .and. numvar.ge.2) then
        call get_node_data(bz_field(1), i, temp)
        call set_dirichlet_bc(i_bz,rhs,temp,normal,curv,izonedim,mat)
     endif

!!$     if(numvar.ge.3 .and. ipressplit.eq.0) then
!!$        if(inograd_p.eq.1) then
!!$           temp = 0.
!!$           call set_normal_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
!!$        end if
!!$        if(ipres.eq.1) then
!!$           call get_node_data(pe_field(1), i, temp)
!!$        else
!!$           call get_node_data(p_field(1), i, temp)
!!$        end if
!!$        if(iconst_p.eq.1) then
!!$           call set_dirichlet_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
!!$        else if(iconst_t.eq.1) then
!!$           call get_node_data(den_v, i, temp2)
!!$           call get_node_data(den_field(1), i, temp3)
!!$           temp = temp*temp2(1)/temp3(1)
!!$           call set_dirichlet_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
!!$        end if
!!$     endif

     if(jadv.eq.0 .and. i3d.eq.1) then
        ! electrostatic potential
        temp = 0.
        call set_dirichlet_bc(i_e,rhs,temp,normal,curv,izonedim,mat)
     endif

     if(imp_bf.eq.1) then
        temp = 0.
        call set_dirichlet_bc(i_bf,rhs,temp,normal,curv,izonedim,mat)
     end if
  end do


  if(eta_wall.ne.0.) call insert_resistive_wall(rhs,psi_v,bz_v,mat)
end subroutine boundary_mag


!=======================================================
! boundary_den
! ~~~~~~~~~~~~
!
! sets boundary conditions for density
!=======================================================
subroutine boundary_den(rhs, den_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(field_type) :: den_v
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x,z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_den called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_n = node_index(den_v, i)

     if(inograd_n.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_n.eq.1) then
        call get_node_data(den_field(1), i, temp)
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_den

!=======================================================
! boundary_te
! ~~~~~~~~~~~
!
! sets boundary conditions for electron temperature
!=======================================================
subroutine boundary_te(rhs, te_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(field_type) :: te_v
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x,z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_te called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_n = node_index(te_v, i)

     if(inograd_t.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_t.eq.1) then
        call get_node_data(te_field(1), i, temp)
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_te


subroutine boundary_ti(rhs, ti_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(field_type) :: ti_v
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x,z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_ti called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_n = node_index(ti_v, i)

     if(inograd_t.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_t.eq.1) then
        call get_node_data(ti_field(1), i, temp)
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_ti


!=======================================================
! boundary_p
! ~~~~~~~~~~
!
! sets boundary conditions for total pressure
!=======================================================
subroutine boundary_p(rhs, p_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(field_type) :: p_v
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_p

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_p called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_p = node_index(p_v, i)

     if(inograd_p.eq.1) then
        temp = 0.
        call set_normal_bc(i_p,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_p.eq.1) then
        call get_node_data(p_field(1), i, temp)
        call set_dirichlet_bc(i_p,rhs,temp,normal,curv,izonedim,mat)
     end if
   end do

 end subroutine boundary_p

!=======================================================
! boundary_pe
! ~~~~~~~~~~~
!
! sets boundary conditions for electron pressure
!=======================================================
subroutine boundary_pe(rhs, pe_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(field_type) :: pe_v
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_pe

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_pe called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_pe = node_index(pe_v, i)
     if(inograd_p.eq.1) then
        temp = 0.
        call set_normal_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_p.eq.1) then
        call get_node_data(pe_field(1), i, temp)
        call set_dirichlet_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
     end if
  end do

end subroutine boundary_pe


  subroutine insert_resistive_wall(rhs,psi_v,bz_v,mat)
    use basic
    use arrays
    use sparse
    use vacuum_interface
    use vector_mod
    use matrix_mod
    use boundary_conditions

    implicit none

    include 'mpif.h'
    
    type(vector_type), intent(inout) :: rhs
    type(field_type) :: psi_v, bz_v
    type(matrix_type), optional :: mat
    
    type(vector_type) :: tempout
    type(vector_type), save :: tempout0
    integer :: i, j, ii, jj, ir, m, n, ip, jp
    integer :: irow(2,dofs_per_node), icol(3,dofs_per_node)
    integer :: jrow(2,dofs_per_node), jcol(3,dofs_per_node)
    integer :: ibegin, jbegin, ierr
    real :: fac, thimprw
    real :: xii, zii, xjj, normi(2), normj(2), curvi, curvj
    integer, parameter :: num_rows = 1
    integer :: rows(3)
    
    logical :: is_boundary
    integer :: izone, izonedim

    logical, save :: first_time = .true.
    logical :: use_resistive_bz
    integer, allocatable :: buf(:)

    vectype, dimension(dofs_per_node,dofs_per_node,2,3) :: ssi, ddi
    vectype, dimension(dofs_per_node,dofs_per_node,2,3) :: ssj, ddj, rrj
    vectype, dimension(dofs_per_node,dofs_per_node,2,3) :: ss, dd
    vectype, dimension(dofs_per_node) :: temp, temp0

    type(vector_type) :: mag
    type(field_type) :: mag_psi, mag_bz, mag_bf

    real :: mp, a, k
    mp = 2.
    a = 1.
    k = ntor/100.

    use_resistive_bz = numvar.ge.2 .and. iconst_bz.ne.1
!    use_resistive_bz = .false.

    rows(1) = 1
    rows(2) = 3

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Inserting resistive wall'

    ! create magnetic field vector
    call create_vector(mag, 3)
    call associate_field(mag_psi, mag, 1)
    mag_psi = psi_v
    call add(mag_psi, external_psi_field, -1.)
    if(numvar.ge.2) then
       call associate_field(mag_bz, mag, 2)
       mag_bz = bz_v
       call add(mag_bz, external_bz_field, -1.)
    endif
    if(i3d.eq.1 .and. numvar.ge.2) then
       call associate_field(mag_bf, mag, 3)
       mag_bf = bf_field(1)
       call add(mag_bf, external_bf_field, -1.)
    endif

    thimprw = thimp
    
    fac = dt*eta_wall/delta_wall
    
    if(first_time) then

       call set_matrix_index(rw_rhs_mat, rw_rhs_mat_index)
       call create_mat(rw_rhs_mat, 2, 3, icomplex, .false.)
       call set_matrix_index(rw_lhs_mat, rw_lhs_mat_index)
       call create_mat(rw_lhs_mat, 2, 3, icomplex, .false.)
#ifdef USERW
       call setmatrixrwb(rw_rhs_mat_index, 1)
       call setmatrixrwb(rw_lhs_mat_index, 1)
#endif
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step rw_rhs_mat", rw_rhs_mat%imatrix     
       print *, "create_mat time_step rw_lhs_mat", rw_lhs_mat%imatrix     
#endif 

       ! calculate global dof ids
       allocate(global_dof_ids_1(nodes))
       allocate(global_dof_ids_row(nodes))
       allocate(global_dof_ids_col(nodes))
       allocate(buf(nodes))

       global_dof_ids_1 = 0
       global_dof_ids_row = 0
       global_dof_ids_col = 0

       do i=1, nodes
          if(local_id(i).le.0) cycle
          global_dof_ids_1(i) = global_node_index(rhs, global_id(i), 1)
          call get_global_node_indices(rw_rhs_mat, global_id(i), irow, icol)
          global_dof_ids_row(i) = irow(1,1)
          global_dof_ids_col(i) = icol(1,1)
       end do
       call mpi_allreduce(global_dof_ids_1, buf, nodes, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ierr)
       global_dof_ids_1 = buf
       call mpi_allreduce(global_dof_ids_row, buf, nodes, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ierr)
       global_dof_ids_row = buf
       call mpi_allreduce(global_dof_ids_col, buf, nodes, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ierr)
       global_dof_ids_col = buf

       ! populate matrices
       ssi = 0.
       ddi = 0.
       ssj = 0.
       ddj = 0.
       rrj = 0.
       ss = 0.
       dd = 0.

       do i=1, nodes
          if(local_id(i).le.0) cycle

          call boundary_node(local_id(i), is_boundary, izone, izonedim, normi,&
               curvi, xii, zii)
          if(itor.eq.0) xii = 1.

          ibegin = global_dof_ids_1(i)

          if(ibegin.le.0) then
             print *, myrank, ': Error: ibegin <= 0', i, local_id(i), global_dof_ids_1(i)
             call safestop(0)
          endif

          do ii=1, 2
             do jj=1, dofs_per_node
                irow(ii,jj) = global_dof_ids_row(i) &
                     + (ii-1)*dofs_per_node + jj-1
             end do
          end do
          do ii=1, 3
             do jj=1, dofs_per_node
                icol(ii,jj) = global_dof_ids_col(i) &
                     + (ii-1)*dofs_per_node + jj-1
             end do
          end do


          do j=1, nodes
             normj(1) = nxnode(j)
             normj(2) = nznode(j)
             xjj = xnode(j)
             if(itor.eq.0) xjj = 1.
          
             jbegin = global_dof_ids_1(j)
             do ii=1, 2
                do jj=1, dofs_per_node
                   jrow(ii,jj) = global_dof_ids_row(j) &
                        + (ii-1)*dofs_per_node + jj-1
                end do
             end do
             do ii=1, 3
                do jj=1, dofs_per_node
                   jcol(ii,jj) = global_dof_ids_col(j) &
                        + (ii-1)*dofs_per_node + jj-1
                end do
             end do

             ! Four indices are:
             ! (dof, operator, equation, field) 

             ! dpsi/dt equation
             ! ~~~~~~~~~~~~~~~~
             ! coefficients of t_j.grad(psi_j)
             rrj(1,3,1,1) = -xii*zgrbth (i,j)/xjj
             rrj(3,3,1,1) = -xii*zgrbthp(i,j)/xjj
             if(itor.eq.1) then
                rrj(3,3,1,1) = rrj(3,3,1,1) + normi(2)*zgrbth(i,j)/xjj
             endif
             
             ! coefficients of n_j.grad(f'_j)
             rrj(1,2,1,3) = -xii*zgrbth (i,j)
             rrj(3,2,1,3) = -xii*zgrbthp(i,j)
             if(itor.eq.1) then
                rrj(3,2,1,3) = rrj(3,2,1,3) + normi(2)*zgrbth(i,j)
             endif

             rrj = rrj*fac
             ssj(:,:,:,1:2) =    -thimprw *rrj(:,:,:,1:2)
             ddj(:,:,:,1:2) = (1.-thimprw)*rrj(:,:,:,1:2)
             ddj(:,:,:,3)   =         rfac*rrj(:,:,:,3)

             if(i.eq.j) then              
                ! dpsi/dt
                ssi(1,1,1,1) = 1.
                ssi(3,3,1,1) = 1.
                ddi(1,1,1,1) = 1.
                ddi(3,3,1,1) = 1.

                ssi(1,2,1,1) =  fac*thimprw
                ssi(3,5,1,1) =  fac*thimprw
                ddi(1,2,1,1) = -fac*(1.-thimprw)
                ddi(3,5,1,1) = -fac*(1.-thimprw)
             endif

!             ! cylindrical approximation
!             ! ===========================
!             ssi = 0.
!             ssj = 0.
!             ddi = 0.
!             ddj = 0.
!
!             if(i.eq.j) then
!                ssi(1,1,1,1) = 1. - fac*(mp/a)**2/k*(-5.e-3)
!                
!                ssi(1,2,1,1) = fac
!                ddi(1,1,1,1) = 1.
!             endif
!             ! ===========================

             do ir=1, num_rows
                ii = rows(ir)
                ip = ibegin + ii - 1
                
                ! transform from boundary coordinates to global coordinates
                ss(ii,:,:,:) = ssj(ii,:,:,:)
                dd(ii,:,:,:) = ddj(ii,:,:,:)
                if(i.eq.j) then
                   ss(ii,:,:,:) = ss(ii,:,:,:) + ssi(ii,:,:,:)
                   dd(ii,:,:,:) = dd(ii,:,:,:) + ddi(ii,:,:,:)
                end if

                ! insert values into the appropriate matrices
                do jj=1, dofs_per_node
                   jp = jbegin + jj - 1
                   
!                   write(90+myrank,'(2I6,2f12.5)') i, j, ss(ii,jj,1,1)

                   if(present(mat)) then
                      call insert_global(mat,ss(ii,jj,1,1), &
                           ip+psi_off,jp+psi_off,MAT_SET)
                   endif
                   call insert_global(rw_lhs_mat,ss(ii,jj,1,1), &
                        irow(1,ii),jcol(1,jj),MAT_SET)
                   call insert_global(rw_rhs_mat,dd(ii,jj,1,1), &
                        irow(1,ii),jcol(1,jj),MAT_SET)
                   if(numvar.ge.2) then
                      call insert_global(rw_rhs_mat,dd(ii,jj,1,2), &
                        irow(1,ii),jcol(2,jj),MAT_SET)
                   endif
                   if(i3d.eq.1 .and. numvar.ge.2) then
                      call insert_global(rw_rhs_mat,dd(ii,jj,1,3), &
                        irow(1,ii),jcol(3,jj),MAT_SET)
                   end if
                end do
             end do
          end do
       end do

       if(use_resistive_bz) then
          call insert_resistive_wall_bz(rw_lhs_mat, rw_rhs_mat, mat)
       endif
       call finalize(rw_rhs_mat)
       call finalize(rw_lhs_mat)

!       call write_matrix(rw_rhs_mat, 'rw_rhs_mat')
!       call write_matrix(rw_lhs_mat, 'rw_lhs_mat')

       if(present(mat)) first_time = .false.

       deallocate(global_dof_ids_1, global_dof_ids_row, global_dof_ids_col)

       call create_vector(tempout0, 2)
       call matvecmult(rw_lhs_mat, external_field, tempout0)
    endif

    call create_vector(tempout, 2)
    call matvecmult(rw_rhs_mat, mag, tempout)

    ! add contributions to rhs vector
    do i=1, nodes
       if(local_id(i).le.0) cycle

       ibegin = node_index(rhs, local_id(i), psi_i)
       call get_node_data(tempout, 1, local_id(i), temp, rotate=.false.)
       call get_node_data(tempout0, 1, local_id(i), temp0, rotate=.false.)
       temp = temp + temp0
       do j=1, num_rows
          call insert(rhs, ibegin+rows(j)-1, temp(rows(j)), VEC_SET)
       end do

       if(use_resistive_bz) then
          ibegin = node_index(rhs, local_id(i), bz_i)
          call get_node_data(tempout, 2, local_id(i), temp, rotate=.false.)
          call get_node_data(tempout0, 2, local_id(i), temp0, rotate=.false.)
          temp = temp + temp0

          do j=1, dofs_per_node
             call insert(rhs, ibegin+j-1, temp(j), VEC_ADD)
          end do
       endif
    end do

    call destroy_vector(tempout)
    call destroy_vector(mag)
  end subroutine insert_resistive_wall

  subroutine insert_resistive_wall_bz(lhs_mat,rhs_mat,mat)
    use basic
    use arrays
    use sparse
    use vacuum_interface
    use vector_mod
    use matrix_mod
    use boundary_conditions
    use m3dc1_nint

    implicit none

    type(matrix_type), intent(inout) :: rhs_mat, lhs_mat
    type(matrix_type) :: mat

    integer :: i, j, k, itri, iedge, inode, jnode, idof, jdof
    vectype :: val, temp

    logical :: is_edge(3), is_boundary
    real :: norm(2,3)
    real :: curvi, xii, zii, xjj
    real, dimension(2) :: normi, normj
    integer :: idim(3), izone, izonedim
    integer :: ii, jj
    real :: thimprw

    integer, dimension(nodes_per_element) :: inodes
    integer :: iii

    integer :: ir_lh, ic_lh, ir_rh, ic_rh
    integer :: jr_lh, jc_lh, jr_rh, jc_rh

    thimprw = thimp

    ! cycle through each boundary edge
    do i=1, boundary_edges
       itri = edge_elm(i)
       iedge = edge_edge(i)
       call boundary_edge(itri, is_edge, norm, idim)
       
       call get_element_nodes(itri,inodes)

       call define_boundary_quadrature(itri, iedge, 5, 5, norm, idim)
       call define_fields(itri, 0, 1, linear)

       ! cycle through each node of the element containing the present edge
       do iii=1,nodes_per_element
          inode = global_node_id(inodes(iii))
          call boundary_node(inodes(iii), is_boundary, izone, izonedim, &
               normi, curvi, xii, zii)
          if(itor.eq.0) xjj = 1.

          ! determine the boundary-index of inode
          do j=1, nodes
             if(inodes(iii).eq.local_id(j)) exit
          end do
          if(j.gt.nodes) then
             if(is_boundary_node(inodes(iii))) then
                print *, 'error! boundary node not found'
                call safestop(9)
             endif
             cycle
          endif

          ! these are the global row and column id's of the first dof
          ! associated with boundary-node j for the two matrices
          ir_lh = global_dof_ids_1(j)
          ic_lh = global_dof_ids_1(j)
          ir_rh = global_dof_ids_row(j)
          ic_rh = global_dof_ids_col(j)

          ! cycle through each trial function associated with inode
          do ii=1,dofs_per_node

          ! idof is the element-based dof index
          ! of global node inode
          idof = (iii-1)*dofs_per_node + ii

          ! for each boundary edge, there are 2 associated nodes
          do j=1, 2
             ! we are only calculating the value (not derivatives) of
             ! the vacuum toroidal field; choose the basis function
             ! associated with the value at inode.

             ! jdof is the element-based dof index
             ! of the first dof associated with the current boundary node
             jdof = mod(iedge+j-2,nodes_per_element)*dofs_per_node+1
             
             ! plasma toroidal field term
             do jj=1, dofs_per_node
                temp = dt*eta_wall/delta_wall * &
                     int3(ri2_79,mu79(:,OP_1,idof),nu79(:,OP_1,jdof+jj-1))

                val = thimprw*temp
                call insert_global(mat,val, &
                     dof_index(ir_lh,bz_i,ii), &
                     dof_index(ic_lh,bz_i,jj),MAT_ADD)
                call insert_global(lhs_mat,val, &
                     dof_index(ir_rh,2,ii),dof_index(ic_rh,2,jj),MAT_ADD)
                
                val = (thimprw-1.)*temp
                call insert_global(rhs_mat,val, &
                     dof_index(ir_rh,2,ii),dof_index(ic_rh,2,jj),MAT_ADD)
             end do
            
             ! The coefficient of the chosen basis function is simply
             ! the value of the vacuum toroidal field at inode,
             ! which is a linear combination of the normal field
             ! at every boundary node.
             do k=1, nodes
                ! int[ mu V(x) dS]
                temp = dt*eta_wall/delta_wall * &
                     (zgrbph(edge_nodes(j,i),k) &
                     *int3(ri_79,mu79(:,OP_1,idof),nu79(:,OP_1,jdof)) )! &
!                     +zgrbphp(edge_nodes(j,i),k) &
!                     *(normi(1) &
!                     *int3(ri_79,mu79(:,OP_1,idof),nu79(:,OP_1,jdof+2)) &
!                     -normi(2) &
!                     *int3(ri_79,mu79(:,OP_1,idof),nu79(:,OP_1,jdof+1))))

                jnode = global_id(k)
                
                normj(1) = nxnode(k)
                normj(2) = nznode(k)
                xjj = xnode(k)
                if(itor.eq.0) xjj = 1.

                ! these are the global row and column id's of the first dof
                ! associated with boundary-node k for the two matrices
                jr_lh = global_dof_ids_1(k)
                jc_lh = global_dof_ids_1(k)
                jr_rh = global_dof_ids_row(k)
                jc_rh = global_dof_ids_col(k)
                
                ! -(1/R_j) * (dpsi/dt)
                val = thimprw*temp/xjj
                call insert_global(mat, val, &
                     dof_index(ir_lh,bz_i,ii),dof_index(jc_lh,psi_i,3),MAT_ADD)
                call insert_global(lhs_mat, val, &
                     dof_index(ir_rh,2,ii),dof_index(jc_rh,1,3),MAT_ADD)

                val = (thimprw-1.)*temp/xjj
                call insert_global(rhs_mat, val, &
                     dof_index(ir_rh,2,ii),dof_index(jc_rh,1,3),MAT_ADD)

                ! -df'/dn
                val = -temp*rfac
                call insert_global(rhs_mat, val, &
                     dof_index(ir_rh,2,ii),dof_index(jc_rh,3,2),MAT_ADD)
             end do
          end do ! on j
       end do    ! on ii
       end do    ! on iii
    end do
  end subroutine insert_resistive_wall_bz


 subroutine get_temperatures(den_v, p_v, pe_v, te_v, ti_v)

  use basic
  use arrays
  use field
  use mesh_mod
  implicit none
  integer :: i, numnodes
  type(field_type) :: den_v, p_v, pe_v
  type(field_type), optional :: te_v, ti_v
   if(numvar.lt.3 .and. ipres.eq.0) return

   numnodes = owned_nodes()
   do i=1,numnodes

     if(idens.eq.1) then
        call get_node_data(den_v,i,den1_l)
     else
        if(eqsubtract.eq.1) then
           den1_l = 0.
        else
           den1_l(1) = 1.
           den1_l(2:dofs_per_node) = 0.
        endif
     endif

     if(numvar.eq.3) then
        if(ipres.eq.1) then
           call get_node_data(p_v,i,p1_l)
           call get_node_data(pe_v,i,pe1_l)
        else
           if(ipressplit.eq.0) then
              call get_node_data(pe_v,i,p1_l)
           else
              call get_node_data(p_v,i,p1_l)
           endif
           pe1_l = pefac*p1_l
        endif
     else
        if(ipres.eq.1) then
           call get_node_data(p_v,i,p1_l)
           pe1_l = pefac*p1_l
           call get_node_data(p_field(0),i,p0_l)
           pe0_l = pefac*p0_l
        endif
     endif

     if(eqsubtract.eq.1) then
        call get_node_data(p_field(0),i,p0_l)
        call get_node_data(pe_field(0),i,pe0_l)
        call get_node_data(den_field(0),i,den0_l)
     end if

     if(eqsubtract.eq.1) then
        call calc_lin_electron_temperature(te1_l, pe0_l, den0_l, pe1_l, den1_l)
        call calc_lin_ion_temperature(ti1_l, p0_l, pe0_l, den0_l, p1_l, pe1_l, den1_l)
     else
        call calc_electron_temperature(te1_l,pe1_l,den1_l)
        call calc_ion_temperature(ti1_l, p1_l, pe1_l, den1_l)
     endif

     if(present(te_v)) call set_node_data(te_v,i,te1_l)
     if(present(ti_v)) then
        if(ipressplit.eq.1 .and. ipres.eq.1) call set_node_data(ti_v,i,ti1_l)
     end if
   enddo
 end subroutine get_temperatures

 subroutine get_pressures(den_v, te_v, ti_v, p_v, pe_v)

  use basic
  use arrays
  use field
  use mesh_mod
  implicit none
  type(field_type) :: den_v, te_v, ti_v, p_v, pe_v

  integer :: i, numnodes
  
   if(numvar.lt.3 .and. ipres.eq.0) return

   numnodes = owned_nodes()
   do i=1,numnodes

     if(idens.eq.1) then
        call get_node_data(den_v,i,den1_l)
        call get_node_data(den_field(0),i,den0_l)
     else
        if(eqsubtract.eq.1) then
           den1_l = 0.
           den0_l(1) = 1.
           den0_l(2:dofs_per_node) = 0.
        else
           den0_l = 0.
           den1_l(1) = 1.
           den1_l(2:dofs_per_node) = 0.
        endif
     endif

     call get_node_data(te_v,i,te1_l)
     call get_node_data(te_field(0),i,te0_l)
     call get_node_data(ti_field(0),i,ti0_l)
     if(ipressplit.eq.1 .and. ipres.eq.1) then
        call get_node_data(ti_v,i,ti1_l)
     else
        ti1_l = te1_l*(1.-pefac)/pefac
     endif

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calc_lin_electron_pressure(pe1_l, te0_l, den0_l, te1_l, den1_l)
        call calc_lin_pressure(p1_l, te0_l, ti0_l, den0_l, te1_l, ti1_l, den1_l)
     else
        call calc_tot_electron_pressure(pe1_l,te1_l,den1_l)
        call calc_tot_pressure(p1_l, te1_l, ti1_l, den1_l)
     endif

     call set_node_data(p_v,i,p1_l)
     if(ipres.eq.1) call set_node_data(pe_v,i,pe1_l)
   enddo


   return
 end subroutine get_pressures
! calc_electron_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the electron temperature
!======================================================================
subroutine calc_electron_temperature(te, pe0, n0)
  use basic

  implicit none

  vectype, intent(out), dimension(dofs_per_node) :: te
  vectype, intent(in), dimension(dofs_per_node) :: pe0, n0

     te(1) = pe0(1)/n0(1)
     te(2) = pe0(2)/n0(1) - pe0(1)*n0(2)/n0(1)**2
     te(3) = pe0(3)/n0(1) - pe0(1)*n0(3)/n0(1)**2
     te(4) = pe0(4)/n0(1) - 2.*pe0(2)*n0(2)/n0(1)**2                       &
           - pe0(1)*n0(4)/n0(1)**2 + 2*pe0(1)*n0(2)*n0(2)/n0(1)**3
     te(5) = pe0(5)/n0(1) - pe0(3)*n0(2)/n0(1)**2  - pe0(2)*n0(3)/n0(1)**2  &
           - pe0(1)*n0(5)/n0(1)**2 + 2*pe0(1)*n0(3)*n0(2)/n0(1)**3
     te(6) = pe0(6)/n0(1) - 2.*pe0(3)*n0(3)/n0(1)**2                       &
           - pe0(1)*n0(6)/n0(1)**2 + 2*pe0(1)*n0(3)*n0(3)/n0(1)**3
#if defined(USE3D)
!    added 08/22/2012
     te(7) = pe0(7)/n0(1) - pe0(1)*n0(7)/n0(1)**2
     te(8) = pe0(8)/n0(1)         - pe0(7)*n0(2)/n0(1)**2   &
           - pe0(2)*n0(7)/n0(1)**2 - pe0(1)*n0(8)/n0(1)**2   &
                                + 2*pe0(1)*n0(2)*n0(7)/n0(1)**3
     te(9) = pe0(9)/n0(1)         - pe0(7)*n0(3)/n0(1)**2   &
           - pe0(2)*n0(7)/n0(1)**2 - pe0(1)*n0(9)/n0(1)**2   &
                                + 2*pe0(1)*n0(3)*n0(7)/n0(1)**3
     te(10) = pe0(10)/n0(1)          - pe0(8)*n0(2)/n0(1)**2   &
            - pe0(8)*n0(2)/n0(1)**2  - pe0(7)*n0(4)/n0(1)**2   &
                                     +2.*pe0(7)*n0(2)*n0(2)/n0(1)**3   &
            - pe0(4)*n0(7)/n0(1)**2         - pe0(2)*n0(8)/n0(1)**2   &
            - pe0(2)*n0(8)/n0(1)**2         - pe0(1)*n0(10)/n0(1)**2   &
            +2.*pe0(2)*n0(7)*n0(2)/n0(1)**3 +2.*pe0(1)*n0(8)*n0(2)/n0(1)**2   &
                                + 2*pe0(2)*n0(2)*n0(7)/n0(1)**3  &
                                + 2*pe0(1)*n0(4)*n0(7)/n0(1)**3  &
                                + 2*pe0(1)*n0(2)*n0(8)/n0(1)**3  &
                                - 6*pe0(1)*n0(2)*n0(7)*n0(2)/n0(1)**4
     te(11) = pe0(11)/n0(1)          - pe0(9)*n0(2)/n0(1)**2   &
            - pe0(8)*n0(3)/n0(1)**2  - pe0(7)*n0(5)/n0(1)**2   &
                                     +2.*pe0(7)*n0(3)*n0(2)/n0(1)**3   &
            - pe0(5)*n0(7)/n0(1)**2         - pe0(2)*n0(9)/n0(1)**2   &
            - pe0(3)*n0(8)/n0(1)**2         - pe0(1)*n0(11)/n0(1)**2   &
            +2.*pe0(2)*n0(7)*n0(3)/n0(1)**3 +2.*pe0(1)*n0(9)*n0(2)/n0(1)**2   &
                                + 2*pe0(3)*n0(2)*n0(7)/n0(1)**3  &
                                + 2*pe0(1)*n0(5)*n0(7)/n0(1)**3  &
                                + 2*pe0(1)*n0(3)*n0(8)/n0(1)**3  &
                                - 6*pe0(1)*n0(3)*n0(7)*n0(2)/n0(1)**4
     te(12) = pe0(12)/n0(1)          - pe0(9)*n0(3)/n0(1)**2   &
            - pe0(9)*n0(3)/n0(1)**2  - pe0(7)*n0(6)/n0(1)**2   &
                                     +2.*pe0(7)*n0(3)*n0(3)/n0(1)**3   &
            - pe0(6)*n0(7)/n0(1)**2         - pe0(3)*n0(9)/n0(1)**2   &
            - pe0(2)*n0(9)/n0(1)**2         - pe0(1)*n0(12)/n0(1)**2   &
            +2.*pe0(3)*n0(7)*n0(3)/n0(1)**3 +2.*pe0(1)*n0(9)*n0(3)/n0(1)**2   &
                                + 2*pe0(3)*n0(3)*n0(7)/n0(1)**3  &
                                + 2*pe0(1)*n0(6)*n0(7)/n0(1)**3  &
                                + 2*pe0(1)*n0(3)*n0(9)/n0(1)**3  &
                                - 6*pe0(1)*n0(3)*n0(7)*n0(3)/n0(1)**4
#endif

     ! account for fact that n0 is ion (not electron) density
     te = te/zeff
end subroutine calc_electron_temperature
! calc_ion_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the ion temperature
!======================================================================
subroutine calc_ion_temperature(ti, pres0, pe0, n0)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: pres0, pe0, n0
  vectype, intent(out), dimension(dofs_per_node) :: ti
  vectype, dimension(dofs_per_node) :: pion0

     pion0 = pres0 - pe0

     ti(1) = pion0(1)/n0(1)
     ti(2) = pion0(2)/n0(1) - pion0(1)*n0(2)/n0(1)**2
     ti(3) = pion0(3)/n0(1) - pion0(1)*n0(3)/n0(1)**2
     ti(4) = pion0(4)/n0(1) - 2.*pion0(2)*n0(2)/n0(1)**2                       &
           - pion0(1)*n0(4)/n0(1)**2 + 2*pion0(1)*n0(2)*n0(2)/n0(1)**3
     ti(5) = pion0(5)/n0(1) - pion0(3)*n0(2)/n0(1)**2  - pion0(2)*n0(3)/n0(1)**2  &
           - pion0(1)*n0(5)/n0(1)**2 + 2*pion0(1)*n0(3)*n0(2)/n0(1)**3
     ti(6) = pion0(6)/n0(1) - 2.*pion0(3)*n0(3)/n0(1)**2                       &
           - pion0(1)*n0(6)/n0(1)**2 + 2*pion0(1)*n0(3)*n0(3)/n0(1)**3
#if defined(USE3D)
!    added 08/22/2012
     ti(7) = pion0(7)/n0(1) - pion0(1)*n0(7)/n0(1)**2
     ti(8) = pion0(8)/n0(1)         - pion0(7)*n0(2)/n0(1)**2   &
           - pion0(2)*n0(7)/n0(1)**2 - pion0(1)*n0(8)/n0(1)**2   &
                                + 2*pion0(1)*n0(2)*n0(7)/n0(1)**3
     ti(9) = pion0(9)/n0(1)         - pion0(7)*n0(3)/n0(1)**2   &
           - pion0(2)*n0(7)/n0(1)**2 - pion0(1)*n0(9)/n0(1)**2   &
                                + 2*pion0(1)*n0(3)*n0(7)/n0(1)**3
     ti(10) = pion0(10)/n0(1)          - pion0(8)*n0(2)/n0(1)**2   &
            - pion0(8)*n0(2)/n0(1)**2  - pion0(7)*n0(4)/n0(1)**2   &
                                     +2.*pion0(7)*n0(2)*n0(2)/n0(1)**3   &
            - pion0(4)*n0(7)/n0(1)**2         - pion0(2)*n0(8)/n0(1)**2   &
            - pion0(2)*n0(8)/n0(1)**2         - pion0(1)*n0(10)/n0(1)**2   &
            +2.*pion0(2)*n0(7)*n0(2)/n0(1)**3 +2.*pion0(1)*n0(8)*n0(2)/n0(1)**3   &
                                + 2*pion0(2)*n0(2)*n0(7)/n0(1)**3  &
                                + 2*pion0(1)*n0(4)*n0(7)/n0(1)**3  &
                                + 2*pion0(1)*n0(2)*n0(8)/n0(1)**3  &
                                - 6*pion0(1)*n0(2)*n0(7)*n0(2)/n0(1)**4
     ti(11) = pion0(11)/n0(1)          - pion0(9)*n0(2)/n0(1)**2   &
            - pion0(8)*n0(3)/n0(1)**2  - pion0(7)*n0(5)/n0(1)**2   &
                                     +2.*pion0(7)*n0(3)*n0(2)/n0(1)**3   &
            - pion0(5)*n0(7)/n0(1)**2         - pion0(2)*n0(9)/n0(1)**2   &
            - pion0(3)*n0(8)/n0(1)**2         - pion0(1)*n0(11)/n0(1)**2   &
            +2.*pion0(2)*n0(7)*n0(3)/n0(1)**3 +2.*pion0(1)*n0(9)*n0(2)/n0(1)**3   &
                                + 2*pion0(3)*n0(2)*n0(7)/n0(1)**3  &
                                + 2*pion0(1)*n0(5)*n0(7)/n0(1)**3  &
                                + 2*pion0(1)*n0(3)*n0(8)/n0(1)**3  &
                                - 6*pion0(1)*n0(3)*n0(7)*n0(2)/n0(1)**4
     ti(12) = pion0(12)/n0(1)          - pion0(9)*n0(3)/n0(1)**2   &
            - pion0(9)*n0(3)/n0(1)**2  - pion0(7)*n0(6)/n0(1)**2   &
                                     +2.*pion0(7)*n0(3)*n0(3)/n0(1)**3   &
            - pion0(6)*n0(7)/n0(1)**2         - pion0(3)*n0(9)/n0(1)**2   &
            - pion0(3)*n0(9)/n0(1)**2         - pion0(1)*n0(12)/n0(1)**2   &
            +2.*pion0(3)*n0(7)*n0(3)/n0(1)**3 +2.*pion0(1)*n0(9)*n0(3)/n0(1)**3   &
                                + 2*pion0(3)*n0(3)*n0(7)/n0(1)**3  &
                                + 2*pion0(1)*n0(6)*n0(7)/n0(1)**3  &
                                + 2*pion0(1)*n0(3)*n0(9)/n0(1)**3  &
                                - 6*pion0(1)*n0(3)*n0(7)*n0(3)/n0(1)**4
#endif

end subroutine calc_ion_temperature
! calc_lin_electron_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the linearized electron temperature
!======================================================================
subroutine calc_lin_electron_temperature(te, pe0, n0, pe1, n1)
  use basic

  implicit none

  vectype, intent(out), dimension(dofs_per_node) :: te
  vectype, intent(in), dimension(dofs_per_node) :: pe0, n0, pe1, n1
  vectype, dimension(dofs_per_node) :: petot, ntot

  petot = pe0 + pe1
  ntot = n0 + n1
  te(1) = petot(1)/ntot(1)  &
       -pe0(1)/n0(1)
  te(2) = petot(2)/ntot(1) &
       - petot(1)*ntot(2)/ntot(1)**2  &
       -pe0(2)/n0(1) + pe0(1)*n0(2)/n0(1)**2
  te(3) = petot(3)/ntot(1) &
       - petot(1)*ntot(3)/ntot(1)**2  &
       -pe0(3)/n0(1) + pe0(1)*n0(3)/n0(1)**2
  te(4) = petot(4)/ntot(1) &
       - 2.*petot(2)*ntot(2)/ntot(1)**2    &
       - petot(1)*ntot(4)/ntot(1)**2  &
       + 2*petot(1)*ntot(2)*ntot(2)/ntot(1)**3   &
       - pe0(4)/n0(1) + 2.*pe0(2)*n0(2)/n0(1)**2                         &
       +  pe0(1)*n0(4)/n0(1)**2 - 2*pe0(1)*n0(2)*n0(2)/n0(1)**3
  te(5) = petot(5)/ntot(1) &
       - petot(3)*ntot(2)/ntot(1)**2  &
       - petot(2)*ntot(3)/ntot(1)**2  &
       - petot(1)*ntot(5)/ntot(1)**2 &
       + 2*petot(1)*ntot(3)*ntot(2)/ntot(1)**3  &
       - pe0(5)/n0(1) + pe0(3)*n0(2)/n0(1)**2  + pe0(2)*n0(3)/n0(1)**2  &
       +  pe0(1)*n0(5)/n0(1)**2 - 2*pe0(1)*n0(3)*n0(2)/n0(1)**3
  te(6) = petot(6)/ntot(1) &
       - 2.*petot(3)*ntot(3)/ntot(1)**2  &
       - petot(1)*ntot(6)/ntot(1)**2   &
       + 2*petot(1)*ntot(3)*ntot(3)/ntot(1)**3   &
       -  pe0(6)/n0(1) + 2.*pe0(3)*n0(3)/n0(1)**2                         &
       +  pe0(1)*n0(6)/n0(1)**2 - 2*pe0(1)*n0(3)*n0(3)/n0(1)**3
#if defined(USE3D)
!    added 08/22/2012
     te(7) = petot(7)/ntot(1) - petot(1)*ntot(7)/ntot(1)**2
     te(8) = petot(8)/ntot(1)         - petot(7)*ntot(2)/ntot(1)**2   &
           - petot(2)*ntot(7)/ntot(1)**2 - petot(1)*ntot(8)/ntot(1)**2   &
                                + 2*petot(1)*ntot(2)*ntot(7)/ntot(1)**3
     te(9) = petot(9)/ntot(1)         - petot(7)*ntot(3)/ntot(1)**2   &
           - petot(2)*ntot(7)/ntot(1)**2 - petot(1)*ntot(9)/ntot(1)**2   &
                                + 2*petot(1)*ntot(3)*ntot(7)/ntot(1)**3
     te(10) = petot(10)/ntot(1)          - petot(8)*ntot(2)/ntot(1)**2   &
            - petot(8)*ntot(2)/ntot(1)**2  - petot(7)*ntot(4)/ntot(1)**2   &
                                     +2.*petot(7)*ntot(2)*ntot(2)/ntot(1)**3   &
            - petot(4)*ntot(7)/ntot(1)**2         - petot(2)*ntot(8)/ntot(1)**2   &
            - petot(2)*ntot(8)/ntot(1)**2         - petot(1)*ntot(10)/ntot(1)**2   &
          +2.*petot(2)*ntot(7)*ntot(2)/ntot(1)**3 +2.*petot(1)*ntot(8)*ntot(2)/ntot(1)**3 &
                                + 2*petot(2)*ntot(2)*ntot(7)/ntot(1)**3  &
                                + 2*petot(1)*ntot(4)*ntot(7)/ntot(1)**3  &
                                + 2*petot(1)*ntot(2)*ntot(8)/ntot(1)**3  &
                                - 6*petot(1)*ntot(2)*ntot(7)*ntot(2)/ntot(1)**4
     te(11) = petot(11)/ntot(1)          - petot(9)*ntot(2)/ntot(1)**2   &
            - petot(8)*ntot(3)/ntot(1)**2  - petot(7)*ntot(5)/ntot(1)**2   &
                                     +2.*petot(7)*ntot(3)*ntot(2)/ntot(1)**3   &
            - petot(5)*ntot(7)/ntot(1)**2         - petot(2)*ntot(9)/ntot(1)**2   &
            - petot(3)*ntot(8)/ntot(1)**2         - petot(1)*ntot(11)/ntot(1)**2   &
          +2.*petot(2)*ntot(7)*ntot(3)/ntot(1)**3 +2.*petot(1)*ntot(9)*ntot(2)/ntot(1)**3 &
                                + 2*petot(3)*ntot(2)*ntot(7)/ntot(1)**3  &
                                + 2*petot(1)*ntot(5)*ntot(7)/ntot(1)**3  &
                                + 2*petot(1)*ntot(3)*ntot(8)/ntot(1)**3  &
                                - 6*petot(1)*ntot(3)*ntot(7)*ntot(2)/ntot(1)**4
     te(12) = petot(12)/ntot(1)          - petot(9)*ntot(3)/ntot(1)**2   &
            - petot(9)*ntot(3)/ntot(1)**2  - petot(7)*ntot(6)/ntot(1)**2   &
                                     +2.*petot(7)*ntot(3)*ntot(3)/ntot(1)**3   &
            - petot(6)*ntot(7)/ntot(1)**2         - petot(3)*ntot(9)/ntot(1)**2   &
            - petot(3)*ntot(9)/ntot(1)**2         - petot(1)*ntot(12)/ntot(1)**2   &
          +2.*petot(3)*ntot(7)*ntot(3)/ntot(1)**3 +2.*petot(1)*ntot(9)*ntot(3)/ntot(1)**3 &
                                + 2*petot(3)*ntot(3)*ntot(7)/ntot(1)**3  &
                                + 2*petot(1)*ntot(6)*ntot(7)/ntot(1)**3  &
                                + 2*petot(1)*ntot(3)*ntot(9)/ntot(1)**3  &
                                - 6*petot(1)*ntot(3)*ntot(7)*ntot(3)/ntot(1)**4
#endif

  ! account for fact that n0 is ion (not electron) density
  te = te/zeff
end subroutine calc_lin_electron_temperature
! calc_lin_ion_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the linearized ion temperature
!======================================================================
subroutine calc_lin_ion_temperature(ti, pres0, pe0, n0, pres1, pe1, n1)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: pres0, pe0, n0, pres1, pe1, n1
  vectype, intent(out), dimension(dofs_per_node) :: ti
  vectype, dimension(dofs_per_node) :: pion0, pion1, ptot, ntot

     pion0 = pres0 - pe0
     pion1 = pres1 - pe1
     ptot = pion0 + pion1
     ntot = n0 + n1

     ti(1) = ptot(1)/ntot(1)  &
           -pion0(1)/n0(1)
     ti(2) = ptot(2)/ntot(1) &
           - ptot(1)*ntot(2)/ntot(1)**2   &
            -pion0(2)/n0(1) + pion0(1)*n0(2)/n0(1)**2
     ti(3) = ptot(3)/ntot(1) &
           - ptot(1)*ntot(3)/ntot(1)**2   &
             -pion0(3)/n0(1) + pion0(1)*n0(3)/n0(1)**2
     ti(4) = ptot(4)/ntot(1) &
           - 2.*ptot(2)*ntot(2)/ntot(1)**2   &
           - ptot(1)*ntot(4)/ntot(1)**2 &
           + 2*ptot(1)*ntot(2)*ntot(2)/ntot(1)**3   &
           - pion0(4)/n0(1) + 2.*pion0(2)*n0(2)/n0(1)**2       &              
           + pion0(1)*n0(4)/n0(1)**2 - 2*pion0(1)*n0(2)*n0(2)/n0(1)**3

     ti(5) = ptot(5)/ntot(1) &
           - ptot(3)*ntot(2)/ntot(1)**2  &
           - ptot(2)*ntot(3)/ntot(1)**2  &
           - ptot(1)*ntot(5)/ntot(1)**2  &
         + 2*ptot(1)*ntot(3)*ntot(2)/ntot(1)**3   &
            - pion0(5)/n0(1) + pion0(3)*n0(2)/n0(1)**2  &
            + pion0(2)*n0(3)/n0(1)**2  &
            + pion0(1)*n0(5)/n0(1)**2 - 2*pion0(1)*n0(3)*n0(2)/n0(1)**3

     ti(6) = ptot(6)/ntot(1) &
           - 2.*ptot(3)*ntot(3)/ntot(1)**2   &
           - ptot(1)*ntot(6)/ntot(1)**2 &
           + 2*ptot(1)*ntot(3)*ntot(3)/ntot(1)**3  &
           -  pion0(6)/n0(1) + 2.*pion0(3)*n0(3)/n0(1)**2                       &
           +  pion0(1)*n0(6)/n0(1)**2 - 2*pion0(1)*n0(3)*n0(3)/n0(1)**3
#if defined(USE3D)
!    added 08/22/2012
     ti(7) = ptot(7)/ntot(1) - ptot(1)*ntot(7)/ntot(1)**2
     ti(8) = ptot(8)/ntot(1)         - ptot(7)*ntot(2)/ntot(1)**2   &
           - ptot(2)*ntot(7)/ntot(1)**2 - ptot(1)*ntot(8)/ntot(1)**2   &
                                + 2*ptot(1)*ntot(2)*ntot(7)/ntot(1)**3
     ti(9) = ptot(9)/ntot(1)         - ptot(7)*ntot(3)/ntot(1)**2   &
           - ptot(2)*ntot(7)/ntot(1)**2 - ptot(1)*ntot(9)/ntot(1)**2   &
                                + 2*ptot(1)*ntot(3)*ntot(7)/ntot(1)**3
     ti(10) = ptot(10)/ntot(1)          - ptot(8)*ntot(2)/ntot(1)**2   &
            - ptot(8)*ntot(2)/ntot(1)**2  - ptot(7)*ntot(4)/ntot(1)**2   &
                                     +2.*ptot(7)*ntot(2)*ntot(2)/ntot(1)**3   &
            - ptot(4)*ntot(7)/ntot(1)**2         - ptot(2)*ntot(8)/ntot(1)**2   &
            - ptot(2)*ntot(8)/ntot(1)**2         - ptot(1)*ntot(10)/ntot(1)**2   &
          +2.*ptot(2)*ntot(7)*ntot(2)/ntot(1)**3 +2.*ptot(1)*ntot(8)*ntot(2)/ntot(1)**3 &
                                + 2*ptot(2)*ntot(2)*ntot(7)/ntot(1)**3  &
                                + 2*ptot(1)*ntot(4)*ntot(7)/ntot(1)**3  &
                                + 2*ptot(1)*ntot(2)*ntot(8)/ntot(1)**3  &
                                - 6*ptot(1)*ntot(2)*ntot(7)*ntot(2)/ntot(1)**4
     ti(11) = ptot(11)/ntot(1)          - ptot(9)*ntot(2)/ntot(1)**2   &
            - ptot(8)*ntot(3)/ntot(1)**2  - ptot(7)*ntot(5)/ntot(1)**2   &
                                     +2.*ptot(7)*ntot(3)*ntot(2)/ntot(1)**3   &
            - ptot(5)*ntot(7)/ntot(1)**2         - ptot(2)*ntot(9)/ntot(1)**2   &
            - ptot(3)*ntot(8)/ntot(1)**2         - ptot(1)*ntot(11)/ntot(1)**2   &
          +2.*ptot(2)*ntot(7)*ntot(3)/ntot(1)**3 +2.*ptot(1)*ntot(9)*ntot(2)/ntot(1)**3 &
                                + 2*ptot(3)*ntot(2)*ntot(7)/ntot(1)**3  &
                                + 2*ptot(1)*ntot(5)*ntot(7)/ntot(1)**3  &
                                + 2*ptot(1)*ntot(3)*ntot(8)/ntot(1)**3  &
                                - 6*ptot(1)*ntot(3)*ntot(7)*ntot(2)/ntot(1)**4
     ti(12) = ptot(12)/ntot(1)          - ptot(9)*ntot(3)/ntot(1)**2   &
            - ptot(9)*ntot(3)/ntot(1)**2  - ptot(7)*ntot(6)/ntot(1)**2   &
                                     +2.*ptot(7)*ntot(3)*ntot(3)/ntot(1)**3   &
            - ptot(6)*ntot(7)/ntot(1)**2         - ptot(3)*ntot(9)/ntot(1)**2   &
            - ptot(3)*ntot(9)/ntot(1)**2         - ptot(1)*ntot(12)/ntot(1)**2   &
          +2.*ptot(3)*ntot(7)*ntot(3)/ntot(1)**3 +2.*ptot(1)*ntot(9)*ntot(3)/ntot(1)**3 &
                                + 2*ptot(3)*ntot(3)*ntot(7)/ntot(1)**3  &
                                + 2*ptot(1)*ntot(6)*ntot(7)/ntot(1)**3  &
                                + 2*ptot(1)*ntot(3)*ntot(9)/ntot(1)**3  &
                                - 6*ptot(1)*ntot(3)*ntot(7)*ntot(3)/ntot(1)**4
#endif
    
     return

end subroutine calc_lin_ion_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
!======================================================================
subroutine calc_tot_electron_pressure(pe, te0, n0)
  use basic

  implicit none

  vectype, intent(out), dimension(dofs_per_node) :: pe
  vectype, intent(in), dimension(dofs_per_node) :: te0, n0

     pe(1) = te0(1)*n0(1)
     pe(2) = te0(2)*n0(1) + te0(1)*n0(2)
     pe(3) = te0(3)*n0(1) + te0(1)*n0(3)
     pe(4) = te0(4)*n0(1) + 2.*te0(2)*n0(2) + te0(1)*n0(4) 
     pe(5) = te0(5)*n0(1) + te0(3)*n0(2)  &
           + te0(2)*n0(3) + te0(1)*n0(5)
     pe(6) = te0(6)*n0(1) + 2.*te0(3)*n0(3) + te0(1)*n0(6)

#if defined(USE3D)
!   added 8/21/2012
     pe(7) = te0(7)*n0(1) &
          +  te0(1)*n0(7)
     pe(8) = te0(8)*n0(1) + te0(7)*n0(2) &
           + te0(1)*n0(8) +  te0(2)*n0(7) 
     pe(9) = te0(9)*n0(1) + te0(3)*n0(7) &
           + te0(1)*n0(9) + te0(7)*n0(3)
     pe(10)= te0(10)*n0(1) + 2.*te0(8)*n0(2) + te0(7)*n0(4) &
           + te0(1)*n0(10) + 2.*te0(2)*n0(8) + te0(4)*n0(7) 
     pe(11)= te0(11)*n0(1) + te0(9)*n0(2)  + te0(8)*n0(3) + te0(7)*n0(5) &
           + te0(1)*n0(11) + te0(2)*n0(9)  + te0(3)*n0(8) + te0(5)*n0(7) 
         
     pe(12)= te0(12)*n0(1) + 2.*te0(9)*n0(3) + te0(7)*n0(6) &
           + te0(1)*n0(12) + 2.*te0(3)*n0(9) + te0(6)*n0(7)
#endif

     return

end subroutine calc_tot_electron_pressure
! calc_tot_pressure
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the total pressure
!======================================================================
subroutine calc_tot_pressure(pres0, te0, ti0, n0)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: te0, ti0, n0
  vectype, intent(out), dimension(dofs_per_node) :: pres0
  vectype, dimension(dofs_per_node) :: tepti

     tepti = te0 + ti0

     pres0(1) = tepti(1)*n0(1)
     pres0(2) = tepti(2)*n0(1) + tepti(1)*n0(2)
     pres0(3) = tepti(3)*n0(1) + tepti(1)*n0(3)
     pres0(4) = tepti(4)*n0(1) + 2.*tepti(2)*n0(2) + tepti(1)*n0(4) 
     pres0(5) = tepti(5)*n0(1) + tepti(3)*n0(2)  &
              + tepti(2)*n0(3) + tepti(1)*n0(5)
     pres0(6) = tepti(6)*n0(1) + 2.*tepti(3)*n0(3) + tepti(1)*n0(6)

#if defined(USE3D)
!   added 8/21/2012
     pres0(7) = tepti(7)*n0(1) &
             +  tepti(1)*n0(7)
     pres0(8) = tepti(8)*n0(1) + tepti(7)*n0(2) &
              + tepti(1)*n0(8) +  tepti(2)*n0(7) 
     pres0(9) = tepti(9)*n0(1) + tepti(3)*n0(7) &
              + tepti(1)*n0(9) + tepti(7)*n0(3)
     pres0(10)= tepti(10)*n0(1) + 2.*tepti(8)*n0(2) + tepti(7)*n0(4) &
              + tepti(1)*n0(10) + 2.*tepti(2)*n0(8) + tepti(4)*n0(7) 
     pres0(11)= tepti(11)*n0(1) + tepti(9)*n0(2)  + tepti(8)*n0(3) + tepti(7)*n0(5) &
              + tepti(1)*n0(11) + tepti(2)*n0(9)  + tepti(3)*n0(8) + tepti(5)*n0(7) 
         
     pres0(12)= tepti(12)*n0(1) + 2.*tepti(9)*n0(3) + tepti(7)*n0(6) &
              + tepti(1)*n0(12) + 2.*tepti(3)*n0(9) + tepti(6)*n0(7)
#endif

     return

end subroutine calc_tot_pressure
subroutine calc_lin_electron_pressure(pe, te0, n0, te1, n1)
  use basic

  implicit none

  vectype, intent(out), dimension(dofs_per_node) :: pe
  vectype, intent(in), dimension(dofs_per_node) :: te0, n0, te1, n1
  vectype, dimension(dofs_per_node) :: ttot, ntot

     ttot = te0 + te1
     ntot = n0 + n1
     pe(1) = ttot(1)*ntot(1) &
           -  te0(1)*n0(1)

     pe(2) = ttot(2)*ntot(1) &
           + ttot(1)*ntot(2) &
           -  te0(2)*n0(1) - te0(1)*n0(2)
 

     pe(3) = ttot(3)*ntot(1) &
           + ttot(1)*ntot(3) &
           -  te0(3)*n0(1) - te0(1)*n0(3)

     pe(4) = ttot(4)*ntot(1) &
           + 2.*ttot(2)*ntot(2) &
           + ttot(1)*ntot(4)    &
           - te0(4)*n0(1) - 2.*te0(2)*n0(2) - te0(1)*n0(4)

     pe(5) = ttot(5)*ntot(1) &
           + ttot(3)*ntot(2) &
           + ttot(2)*ntot(3) &
           + ttot(1)*ntot(5) &
           - te0(5)*n0(1) - te0(3)*n0(2) -  te0(2)*n0(3) - te0(1)*n0(5)

     pe(6) = ttot(6)*ntot(1) &
           + 2.*ttot(3)*ntot(3) &
           + ttot(1)*ntot(6) &
           - te0(6)*n0(1) - 2.*te0(3)*n0(3) - te0(1)*n0(6)
#if defined(USE3D)
!   added 8/21/2012
     pe(7) = te1(7)*ntot(1) &
          +  ttot(1)*n1(7)
     pe(8) = te1(8)*ntot(1) + te1(7)*ntot(2) &
           + ttot(1)*n1(8) +  ttot(2)*n1(7)
     pe(9) = te1(9)*ntot(1) + ttot(3)*n1(7) &
           + ttot(1)*n1(9) + te1(7)*ntot(3)
     pe(10)= te1(10)*ntot(1) + 2.*te1(8)*ntot(2) + te1(7)*ntot(4) &
           + ttot(1)*n1(10) + 2.*ttot(2)*n1(8) + ttot(4)*n1(7)
     pe(11)= te1(11)*ntot(1) + te1(9)*ntot(2)  + te1(8)*ntot(3) + te1(7)*ntot(5) &
           + ttot(1)*n1(11) + ttot(2)*n1(9)  + ttot(3)*n1(8) + ttot(5)*n1(7)
     pe(12)= te1(12)*ntot(1) + 2.*te1(9)*ntot(3) + te1(7)*ntot(6) &
           + ttot(1)*n1(12) + 2.*ttot(3)*n1(9) + ttot(6)*n1(7)
#endif

     return

end subroutine calc_lin_electron_pressure
! calc_linearized pressure
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the linearized pressure
!======================================================================
subroutine calc_lin_pressure(pres1, te0, ti0, n0, te1, ti1, n1)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: te0, ti0, n0, te1, ti1, n1
  vectype, intent(out), dimension(dofs_per_node) :: pres1
  vectype, dimension(dofs_per_node) :: tepti, tepti0, nt

     tepti0 = te0 + ti0
     tepti  = te0 + te1 + ti0 + ti1
     nt = n0 + n1

     pres1(1) = tepti(1)*nt(1) - tepti0(1)*n0(1)
     pres1(2) = tepti(2)*nt(1) + tepti(1)*nt(2) &
             -  tepti0(2)*n0(1)- tepti0(1)*n0(2)
     pres1(3) = tepti(3)*nt(1) + tepti(1)*nt(3) &
              - tepti0(3)*n0(1)- tepti0(1)*n0(3)
     pres1(4) = tepti(4)*nt(1) + 2.*tepti(2)*nt(2) + tepti(1)*nt(4) &
              - tepti0(4)*n0(1)- 2.*tepti0(2)*n0(2)- tepti0(1)*n0(4)
     pres1(5) = tepti(5)*nt(1) + tepti(3)*nt(2)  &
              + tepti(2)*nt(3) + tepti(1)*nt(5)  &
              - tepti0(5)*n0(1)- tepti0(3)*n0(2) &
              - tepti0(2)*n0(3)- tepti0(1)*n0(5)
     pres1(6) = tepti(6)*nt(1) + 2.*tepti(3)*nt(3) + tepti(1)*nt(6) &
              - tepti0(6)*n0(1)- 2.*tepti0(3)*n0(3)- tepti0(1)*n0(6)

#if defined(USE3D)
!   added 8/21/2012
     pres1(7) = tepti(7)*nt(1) &
             +  tepti(1)*nt(7)
     pres1(8) = tepti(8)*nt(1) + tepti(7)*nt(2) &
              + tepti(1)*nt(8) +  tepti(2)*nt(7) 
     pres1(9) = tepti(9)*nt(1) + tepti(3)*nt(7) &
              + tepti(1)*nt(9) + tepti(7)*nt(3)
     pres1(10)= tepti(10)*nt(1) + 2.*tepti(8)*nt(2) + tepti(7)*nt(4) &
              + tepti(1)*nt(10) + 2.*tepti(2)*nt(8) + tepti(4)*nt(7) 
     pres1(11)= tepti(11)*nt(1) + tepti(9)*nt(2)  + tepti(8)*nt(3) + tepti(7)*nt(5) &
              + tepti(1)*nt(11) + tepti(2)*nt(9)  + tepti(3)*nt(8) + tepti(5)*nt(7) 
         
     pres1(12)= tepti(12)*nt(1) + 2.*tepti(9)*nt(3) + tepti(7)*nt(6) &
              + tepti(1)*nt(12) + 2.*tepti(3)*nt(9) + tepti(6)*nt(7)
#endif




     return

end subroutine calc_lin_pressure

subroutine calcnorm(temp, nsize, l2norm)
  use basic
  use arrays
  use field
  use mesh_mod
  implicit none
  integer :: nsize, numnodes, i,ii
  real :: l2norm, sum
  type(vector_type), intent(in) :: temp

  numnodes = owned_nodes()

  sum = 0.
  do i = 1,numnodes*nsize
    ii = 1 + (i-1)*dofs_per_node
#ifdef USECOMPLEX
    sum = sum + temp%data(ii)*conjg(temp%data(ii))
#else
    sum = sum + temp%data(ii)**2
#endif
  enddo
  l2norm = sqrt(sum)
  return
  

end subroutine calcnorm



end module model
