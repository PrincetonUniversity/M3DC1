module model

  use vector_mod
  use matrix_mod

  integer, allocatable :: global_dof_ids_1(:)
  integer, allocatable :: global_dof_ids_row(:), global_dof_ids_col(:)

  type(vector_type), target :: q4_vec, r4_vec, qp4_vec, qn4_vec, qn5_vec

  ! matrices
  type(matrix_type), target :: s1_mat, d1_mat, q1_mat, r14_mat
  type(matrix_type), target :: o1_mat, p1_mat
  type(matrix_type), target :: q42_mat, r42_mat, q43_mat, r43_mat
  type(matrix_type), target :: q44_mat, r44_mat
  type(matrix_type), target :: s2_mat, d2_mat, r2_mat, q2_mat, o2_mat, o3_mat
  type(matrix_type), target :: s8_mat, d8_mat, r8_mat, q8_mat
  type(matrix_type), target :: s9_mat, d9_mat, r9_mat, q9_mat, o9_mat
  type(matrix_type), target :: qp42_mat, rp42_mat
  type(matrix_type), target :: s11_mat, d11_mat, s12_mat, d12_mat
  type(matrix_type), target :: s15_mat, d15_mat, r15_mat, q15_mat, k15_mat

  ! positions
  integer :: u_i, vz_i, chi_i
  integer :: psi_i, bz_i, pe_i
  integer :: den_i, p_i
  integer :: bf_i, e_i
  integer :: te_i, ti_i
  integer :: nre_i

contains

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
  
  call get_boundary_mask(itri, ibound, imask, BOUND_ANY)
end subroutine get_den_mask

subroutine get_nre_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(0.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(1.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)

  call get_boundary_mask(itri, ibound, imask, BOUND_ANY)
end subroutine get_nre_mask

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
  if(iconst_t.eq.1 .or. iconst_p.ge.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  
  call get_boundary_mask(itri, ibound, imask, BOUND_ANY)
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
  if(iconst_p.ge.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  
  call get_boundary_mask(itri, ibound, imask, BOUND_ANY)
end subroutine get_pres_mask


subroutine get_flux_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inocurrent_tor.eq.1) ibound = ior(ibound, BOUNDARY_LAPLACIAN)
  if(inocurrent_norm.eq.1) then
     if(i3d.eq.1) then
        ibound = ior(ibound, BOUNDARY_NEUMANNP)
     endif
  endif
  if(iconst_bn.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)

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
  if(inocurrent_pol.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(iconst_bz.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
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

  if(ifbound.eq.1) then
     ibound = BOUNDARY_DIRICHLET
  else if(ifbound.eq.2) then 
     ibound = BOUNDARY_NEUMANN
  else
     ibound = 0
  end if
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_bf_mask

subroutine get_esp_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = BOUNDARY_DIRICHLET
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_esp_mask

subroutine get_j_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = BOUNDARY_DIRICHLET
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_j_mask

subroutine get_vor_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0

  if(inonormalflow.eq.1) then
     ibound = ior(ibound, BOUNDARY_DIRICHLET)
  elseif(inonormalflow.eq.2) then
     ibound = ior(ibound, BOUNDARY_MULTI_DT)
  end if

  if(inoslip_pol.eq.1) then
     ibound = ior(ibound, BOUNDARY_NEUMANN)
  elseif(inoslip_pol.eq.2) then
     ibound = ior(ibound, BOUNDARY_MULTI_DN)
  end if

  if(vor_bc.eq.1)        ibound = ior(ibound, BOUNDARY_LAPLACIAN)

  call get_boundary_mask(itri, ibound, imask, BOUND_ANY)

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
  call get_boundary_mask(itri, ibound, imask, BOUND_ANY)
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

  ! inonormalflow=2 and inoslip_pol=2 conditions are included in U equation

  if(com_bc.eq.1)        ibound = ior(ibound, BOUNDARY_LAPLACIAN)
  call get_boundary_mask(itri, ibound, imask, BOUND_ANY)
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
  real :: normal(2), curv(3), x, z, phi
  integer :: i, izone, izonedim, numnodes, icounter_t
  integer :: i_u, i_vz, i_chi
  logical :: is_boundary
  integer :: ibegin(2), ibc(2)
  vectype :: coeff(2), xp(2)

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) then
     print *, "boundary_vel called"
     if(inonormalflow.eq.2) print *, "  Using multibc for normal flow"
     if(inoslip_pol.eq.2)   print *, "  Using multibc for poloidal flow"
  end if

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_ANY)
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
     elseif(inonormalflow.eq.2) then

        temp = 0.

        ! U
        ibegin(1) = i_u
        ibc(1) = BOUND_DT
        coeff(1) = -1.0
        xp(1) = itor

        ! chi
        ibegin(2) = i_chi
        ibc(2) = BOUND_DN
        coeff(2) = 1.0
        xp(2) = -2*itor

        call set_multi_bc(2,ibegin,ibc,coeff,xp,rhs,temp,normal,curv,izonedim,x,mat)

     end if
     
     ! no poloidal slip
     if(inoslip_pol.eq.1) then
        temp = 0.
        call set_normal_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
        if(numvar.ge.3) then
           call set_dirichlet_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
        endif
     elseif(inoslip_pol.eq.2) then

        temp = 0.

        ! U
        ibegin(1) = i_u
        ibc(1) = BOUND_DN
        coeff(1) = 1.0
        xp(1) = itor

        ! chi
        ibegin(2) = i_chi
        ibc(2) = BOUND_DT
        coeff(2) = 1.0
        xp(2) = -2*itor

        call set_multi_bc(2,ibegin,ibc,coeff,xp,rhs,temp,normal,curv,izonedim,x,mat)

     end if

     ! toroidal velocity
     if(numvar.ge.2) then
        ! no slip
        if(inoslip_tor.eq. 1) then
#ifdef USEST
           temp = 0.
#else
           call get_node_data(vz_field(1), i, temp)
#endif
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
        call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,-x,mat)
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
  real :: normal(2), curv(3), x, z, phi
  integer :: i, izone, izonedim, numnodes, icounter_t
  integer :: i_u, i_chi
  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vpol called"

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_ANY)
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
        call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,-x,mat)
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
subroutine boundary_mag(rhs, psi_v, bz_v, bfp_v, e_v, mat)
  use math
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions

  implicit none
  
  type(vector_type), intent(inout) :: rhs
  type(field_type) :: psi_v, bz_v, bfp_v, e_v
  type(matrix_type), optional :: mat

  vectype, dimension(dofs_per_node) :: temp !, temp2, temp3
  real :: normal(2), curv(3), x, z, phi
  integer :: i, izone, izonedim,  numnodes, icounter_t
  integer :: i_psi, i_bz, i_e, i_bf !, i_pe
  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_mag called"

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z)
     if(.not.is_boundary) cycle

     i_psi = node_index(psi_v, i)
     if(numvar.ge.2) i_bz = node_index(bz_v, i)
!     if(numvar.ge.3 .and. ipressplit.eq.0) i_pe = node_index(pe_v, i)
     if((jadv.eq.0).or.(jadv.eq.1 .and. imp_hyper.ge.1)) i_e = node_index(e_v, i)
     if(imp_bf.eq.1) i_bf = node_index(bfp_v, i)

     ! constant normal field = -t.grad(psi)/R - n.grad(f')
     if(iconst_bn.eq.1) then
        call get_node_data(psi_field(1), i, temp)
        ! add loop voltage
        if(idiff .gt. 0) temp = 0.
        if(igauge.eq.0) temp(1) = temp(1) + dt*vloop/toroidal_period
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
        if(idiff .gt. 0) temp = 0
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
           if(idiff .gt. 0) temp = 0
           call set_normalp_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)
        endif

        if(numvar.ge.2) then
           call get_node_data(bz_field(1), i, temp)
           if(idiff .gt. 0) temp = 0
           call set_dirichlet_bc(i_bz,rhs,temp,normal,curv,izonedim,mat)
        endif
     else if(iconst_bz.eq.1 .and. numvar.ge.2) then
        call get_node_data(bz_field(1), i, temp)
        if(idiff .gt. 0) temp = 0
        temp(1) = temp(1) + dt*vloopRZ
        call set_dirichlet_bc(i_bz,rhs,temp,normal,curv,izonedim,mat)
     endif

     if((jadv.eq.0).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
        ! electrostatic potential or del_star_psi
        temp = 0.
        call set_dirichlet_bc(i_e,rhs,temp,normal,curv,izonedim,mat)
     endif

     if(imp_bf.eq.1) then
#ifdef USEST
        call get_node_data(bfp_field(1), i, temp)
#else
        temp = 0.
#endif
        if(ifbound.eq.1) then 
           call set_dirichlet_bc(i_bf,rhs,temp,normal,curv,izonedim,mat)
        else if(ifbound.eq.2) then 
           call set_normal_bc(i_bf,rhs,temp,normal,curv,izonedim,mat)
        end if
     end if
  end do

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
  
  integer :: i, izone, izonedim, numnodes, icounter_t
  real :: normal(2), curv(3), x,z, phi
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_den called"

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_ANY)
     if(.not.is_boundary) cycle

     i_n = node_index(den_v, i)

     if(inograd_n.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_n.eq.1) then
        call get_node_data(den_field(1), i, temp)
        if(idiff .gt. 0) temp = 0.   ! this is for change in density from n to n+1
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_den

!=======================================================
! boundary_nre
! ~~~~~~~~~~~~
!
! sets boundary conditions for RE
!=======================================================
subroutine boundary_nre(rhs, nre_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none

  type(vector_type) :: rhs
  type(field_type) :: nre_v
  type(matrix_type), optional :: mat

  integer :: i, izone, izonedim, numnodes, icounter_t
  real :: normal(2), curv(3), x,z, phi
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_nre called"

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_ANY)
     if(.not.is_boundary) cycle

     i_n = node_index(nre_v, i)

     if(0.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(1.eq.1) then
        call get_node_data(nre_field(1), i, temp)
        temp = 0.
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_nre


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
  
  integer :: i, izone, izonedim, numnodes, icounter_t
  real :: normal(2), curv(3), x,z, phi
  logical :: is_boundary, is_inner
  vectype, dimension(dofs_per_node) :: temp, temp2, temp3

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_te called"

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_ANY)
     if(.not.is_boundary) cycle
     call boundary_node(i,is_inner,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_FIRSTWALL)
     
     i_n = node_index(te_v, i)

     if(inograd_t.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_t.eq.1) then
        if((tebound.gt.0.) .and. is_inner) then
           temp = 0.
           temp(1) = tebound
        else
           call get_node_data(te_field(1), i, temp)
        end if

        if(idiff .gt. 0) temp = 0.

        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     else if(iconst_t.eq.2) then
        temp = 0.
        if(eqsubtract.eq.0 .and. idiff.eq.0) temp(1) = tedge
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     else if(iconst_p.ge.1) then
        call get_node_data(pe_field(1), i, temp)
        call get_node_data(ne_field(1), i, temp2)
        if(eqsubtract.eq.1) then
           call get_node_data(ne_field(0), i, temp3)
           temp2 = temp2 + temp3
        end if

        temp3 = 0.
        temp3(1) = temp(1) / temp2(1)
        temp3(2) = temp(2) / temp2(1) - temp2(2) * temp(1) / temp2(1)**2
        temp3(3) = temp(3) / temp2(1) - temp2(3) * temp(1) / temp2(1)**2

        if(eqsubtract.eq.1 .and. idiff.eq.0) then
           call get_node_data(te_field(0),i,temp2)
           temp3 = temp3 - temp2
        end if

        call set_dirichlet_bc(i_n,rhs,temp3,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_te


subroutine boundary_ti(rhs, ti_v, mat)
  use basic
  use field
  use arrays
  use matrix_mod
  use boundary_conditions
  use kprad_m3dc1
  implicit none
  
  type(vector_type) :: rhs
  type(field_type) :: ti_v
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes, icounter_t
  real :: normal(2), curv(3), x,z, phi
  logical :: is_boundary, is_inner
  vectype, dimension(dofs_per_node) :: temp, temp2, temp3

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_ti called"

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_ANY)
     if(.not.is_boundary) cycle
     call boundary_node(i,is_inner,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_FIRSTWALL)
     
     i_n = node_index(ti_v, i)

     if(inograd_t.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_t.eq.1 .or. iconst_p.ge.1) then
        if((tibound.gt.0.) .and. is_inner) then
           temp = 0.
           temp(1) = tibound
        else
           call get_node_data(ti_field(1), i, temp)
        end if

        if(idiff .gt. 0) temp = 0.

        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
!!$     else if(iconst_p.eq.1) then
!!$        call get_node_data(p_field(1), i, temp)
!!$        call get_node_data(pe_field(1), i, temp2)
!!$        temp = temp - temp2
!!$
!!$        call get_node_data(den_field(1), i, temp2)
!!$        if(eqsubtract.eq.1) then
!!$           call get_node_data(den_field(0), i, temp3)
!!$           temp2 = temp2 + temp3
!!$        end if
!!$        if(ikprad.ne.0) then
!!$           do i=1, kprad_z
!!$              call get_node_data(kprad_n(i), i, temp3)
!!$              temp2 = temp2 + temp3
!!$           end do
!!$        end if
!!$
!!$        temp3 = 0.
!!$        temp3(1) = temp(1) / temp2(1)
!!$        temp3(2) = temp(2) / temp2(1) - temp2(2) * temp(1) / temp2(1)**2 
!!$        temp3(3) = temp(3) / temp2(1) - temp2(3) * temp(1) / temp2(1)**2 
!!$
!!$        if(eqsubtract.eq.1 .and. idiff.eq.0) then
!!$           call get_node_data(ti_field(0),i,temp2)
!!$           temp3 = temp3 - temp2
!!$        end if
!!$
!!$        call set_dirichlet_bc(i_n,rhs,temp3,normal,curv,izonedim,mat)
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
  
  integer :: i, izone, izonedim, numnodes, icounter_t
  real :: normal(2), curv(3), x, z, phi
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_p

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_p called"

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_ANY)
     if(.not.is_boundary) cycle

     i_p = node_index(p_v, i)

     if(inograd_p.eq.1) then
        temp = 0.
        call set_normal_bc(i_p,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_p.eq.1) then
        call get_node_data(p_field(1), i, temp)
       
        if(idiff .gt. 0) temp = 0.

        call set_dirichlet_bc(i_p,rhs,temp,normal,curv,izonedim,mat)
     else if(iconst_p.eq.2) then
        temp = 0.
        if(eqsubtract.eq.0 .and. idiff.eq.0) temp(1) = pedge

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
  
  integer :: i, izone, izonedim, numnodes, icounter_t
  real :: normal(2), curv(3), x, phi, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_pe

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_pe called"

  numnodes = owned_nodes()
  do icounter_t=1,numnodes
     i = nodes_owned(icounter_t)
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
          BOUND_ANY)
     if(.not.is_boundary) cycle

     i_pe = node_index(pe_v, i)
     if(inograd_p.eq.1) then
        temp = 0.
        call set_normal_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_p.eq.1) then
        call get_node_data(pe_field(1), i, temp)

        if(idiff .gt. 0) temp = 0.

        call set_dirichlet_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
     else if(iconst_p.eq.2) then
        temp = 0.
        if(eqsubtract.eq.0 .and. idiff.eq.0) temp(1) = pedge*pefac
     end if
  end do

end subroutine boundary_pe


end module model
