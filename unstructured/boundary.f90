module boundary_conditions

  integer, parameter :: BOUNDARY_NONE           =  0
  integer, parameter :: BOUNDARY_DIRICHLET      =  1
  integer, parameter :: BOUNDARY_NEUMANN        =  2
  integer, parameter :: BOUNDARY_LAPLACIAN      =  4
  integer, parameter :: BOUNDARY_RESISTIVE_WALL =  8
  integer, parameter :: BOUNDARY_AXISYMMETRIC   = 16

contains


!======================================================================
! get_boundary_mask
! ~~~~~~~~~~~~~~~~~
! returns a vector imask, such that 
! imask(i) = 0 if the row associated with dof i of element itri
!              will be overwritten by boundary condition ibound
! imask(i) = 1 otherwise
!======================================================================
subroutine get_boundary_mask(itri, ibound, imask)
  use element
  use mesh_mod
  implicit none

  integer, intent(in) :: itri, ibound
  integer, intent(out), dimension(dofs_per_element) :: imask

  integer :: inode(nodes_per_element)
  real :: norm(2), curv, x, z
  integer :: izone, izonedim
  logical :: is_boundary
  integer :: i, k
#ifdef USE3D
  integer :: iplane
#endif
  
  imask = 1
  if(ibound.eq.BOUNDARY_NONE) return
  
  call get_element_nodes(itri, inode)
  do i=1, nodes_per_element
     k = (i-1)*dofs_per_node + 1

#ifdef USE3D
     if(iand(ibound, BOUNDARY_AXISYMMETRIC).eq.BOUNDARY_AXISYMMETRIC) then
        call getplaneid(iplane)
           
        ! in the axisymmetric case, the only elements not set by bcs
        ! are those on associated with nodes on the first plane
        ! that are not also associated with toroidal derivatives
        if(iplane.eq.1 .and. i.le.pol_nodes_per_element) then
           imask(k+6:k+11) = 0
        else
           imask(k:k+11) = 0
        endif
     endif
#endif

     call boundary_node(inode(i),is_boundary,izone,izonedim,norm,curv,x,z)
     if(.not.is_boundary) cycle

     if(iand(ibound, BOUNDARY_RESISTIVE_WALL).eq.BOUNDARY_RESISTIVE_WALL) then
        imask(k) = 0
     endif
     if(iand(ibound, BOUNDARY_DIRICHLET).eq.BOUNDARY_DIRICHLET) then
        imask(k  ) = 0
        imask(k+2) = 0
        imask(k+5) = 0
        if(izonedim.eq.0) then
           imask(k+1) = 0
           imask(k+3) = 0
        end if
#ifdef USE3D
        imask(k+6 ) = 0
        imask(k+8 ) = 0
        imask(k+11) = 0
        if(izonedim.eq.0) then
           imask(k+7 ) = 0
           imask(k+9 ) = 0
        end if
#endif
     endif
     if(iand(ibound, BOUNDARY_NEUMANN).eq.BOUNDARY_NEUMANN) then
        imask(k+1) = 0
        imask(k+4) = 0
        if(izonedim.eq.0) then
           imask(k+2) = 0
        endif
#ifdef USE3D
        imask(k+7) = 0
        imask(k+10) = 0
        if(izonedim.eq.0) then
           imask(k+8) = 0
        endif
#endif
     endif
     if(iand(ibound, BOUNDARY_LAPLACIAN).eq.BOUNDARY_LAPLACIAN) then
        imask(k+3) = 0
     endif
  end do
end subroutine get_boundary_mask


subroutine apply_boundary_mask(itri, ibound, vals, imaskin)
  use element
  integer, intent(in) :: itri, ibound
  vectype, intent(inout), dimension(dofs_per_element,dofs_per_element) :: vals
  integer, dimension(dofs_per_element), optional :: imaskin

  integer, dimension(dofs_per_element) :: imask
  integer :: i

  if(present(imaskin)) then
     imask = imaskin
  else
     call get_boundary_mask(itri, ibound, imask)
  end if
  do i=1, dofs_per_element
     vals(i,:) = vals(i,:)*imask(i)
  end do
end subroutine apply_boundary_mask


!======================================================================
! set_total_bc
!======================================================================
subroutine set_total_bc(ibegin,rhs,bv,normal,curv,izonedim,mat)
  use basic
  use vector_mod
  use matrix_mod
  implicit none

  type(matrix_type), optional :: mat          ! matrix handle
  integer, intent(in) :: ibegin               ! first dof of field
  type(vector_type) :: rhs                    ! right-hand-side of equation
  vectype, intent(in), dimension(dofs_per_node) :: bv     ! boundary values
  real, intent(in) :: normal(2), curv
  integer, intent(in) :: izonedim             ! dimension of boundary
  
  ! clamp value
  if(present(mat)) then
     call identity_row(mat, ibegin)
     call identity_row(mat, ibegin+1)
     call identity_row(mat, ibegin+2)
     call identity_row(mat, ibegin+3)
     call identity_row(mat, ibegin+4)
     call identity_row(mat, ibegin+5)
  end if
  call insert(rhs, ibegin  , bv(1), VEC_SET)
  call insert(rhs, ibegin+1, bv(2), VEC_SET)
  call insert(rhs, ibegin+2, bv(3), VEC_SET)
  call insert(rhs, ibegin+3, bv(4), VEC_SET)
  call insert(rhs, ibegin+4, bv(5), VEC_SET)
  call insert(rhs, ibegin+5, bv(6), VEC_SET)
end subroutine set_total_bc




!======================================================================
! set_dirichlet_bc
!======================================================================
subroutine set_dirichlet_bc(ibegin,rhs,bv,normal,curv,izonedim,mat)
  use basic
  use vector_mod
  use matrix_mod
  implicit none

  integer, intent(in) :: ibegin               ! first dof of field
  type(vector_type) :: rhs                    ! right-hand-side of equation
  vectype, intent(in), dimension(dofs_per_node) :: bv  ! boundary values
  real, intent(in) :: normal(2), curv
  integer, intent(in) :: izonedim             ! dimension of boundary
  type(matrix_type), optional :: mat
  
  ! clamp value
  if(present(mat)) call identity_row(mat, ibegin)
  call insert(rhs, ibegin, bv(1), VEC_SET)
#ifdef USE3D
  if(present(mat)) call identity_row(mat, ibegin+6)
  call insert(rhs, ibegin+6, bv(7), VEC_SET)
#endif
     
  ! clamp tangential derivative
  call set_tangent_bc(ibegin,rhs,bv,normal,curv,izonedim,mat)

end subroutine set_dirichlet_bc


!======================================================================
! set_tangent_bc
!======================================================================
subroutine set_tangent_bc(ibegin,rhs,bv,normal,curv,izonedim,mat)
  use basic
  use vector_mod
  use matrix_mod

  implicit none
  
  integer, intent(in) :: ibegin               ! first dof of field
  real, intent(in) :: normal(2), curv
  type(vector_type), intent(inout) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(dofs_per_node) :: bv     ! boundary values
  integer, intent(in) :: izonedim             ! dimension of boundary
  type(matrix_type), optional :: mat

  vectype, dimension(dofs_per_node) :: bv_rotated

  call rotate_dofs(bv, bv_rotated, normal, curv, 1)

  ! t
  if(present(mat)) call identity_row(mat, ibegin+2)
  call insert(rhs, ibegin+2, bv_rotated(3), VEC_SET)
#ifdef USE3D
  if(present(mat)) call identity_row(mat, ibegin+8)
  call insert(rhs, ibegin+8, bv_rotated(9), VEC_SET)
#endif
  

  ! tt
  if(present(mat)) call identity_row(mat, ibegin+5)
  call insert(rhs, ibegin+5, bv_rotated(6), VEC_SET)
#ifdef USE3D
  if(present(mat)) call identity_row(mat, ibegin+11)
  call insert(rhs, ibegin+11, bv_rotated(12), VEC_SET)
#endif

  if(izonedim.eq.0) then
     ! n
     if(present(mat)) call identity_row(mat, ibegin+1)
     call insert(rhs, ibegin+1, bv_rotated(2), VEC_SET)
#ifdef USE3D
     if(present(mat)) call identity_row(mat, ibegin+7)
     call insert(rhs, ibegin+7, bv_rotated(8), VEC_SET)
#endif

     ! nn
     if(present(mat)) call identity_row(mat, ibegin+3)
     call insert(rhs, ibegin+3, bv_rotated(4), VEC_SET)
#ifdef USE3D
     if(present(mat)) call identity_row(mat, ibegin+9)
     call insert(rhs, ibegin+9, bv_rotated(10), VEC_SET)
#endif
  endif

end subroutine set_tangent_bc


!======================================================================
! set_normal_bc
!======================================================================
subroutine set_normal_bc(ibegin,rhs,bv,normal,curv,izonedim,mat)
  use basic
  use vector_mod
  use matrix_mod

  implicit none
  
  integer, intent(in) :: ibegin               ! first dof of field
  real, intent(in) :: normal(2), curv
  type(vector_type), intent(inout) :: rhs     ! right-hand-side of equation
  vectype, intent(in), dimension(dofs_per_node) :: bv     ! boundary values
  integer, intent(in) :: izonedim             ! dimension of boundary
  type(matrix_type), optional :: mat

  vectype, dimension(dofs_per_node) :: bv_rotated

  call rotate_dofs(bv, bv_rotated, normal, curv, 1)
     
  ! n
  if(present(mat)) call identity_row(mat, ibegin+1)
  call insert(rhs, ibegin+1, bv_rotated(2), VEC_SET)
#ifdef USE3D
     if(present(mat)) call identity_row(mat, ibegin+7)
     call insert(rhs, ibegin+7, bv_rotated(8), VEC_SET)
#endif

  ! nt
  if(present(mat)) call identity_row(mat, ibegin+4)
  call insert(rhs, ibegin+4, bv_rotated(5), VEC_SET)
#ifdef USE3D
     if(present(mat)) call identity_row(mat, ibegin+10)
     call insert(rhs, ibegin+10, bv_rotated(11), VEC_SET)
#endif

  if(izonedim.eq.0) then
     ! t
     if(present(mat)) call identity_row(mat, ibegin+2)
     call insert(rhs, ibegin+2, bv_rotated(3), VEC_SET)
#ifdef USE3D
     if(present(mat)) call identity_row(mat, ibegin+8)
     call insert(rhs, ibegin+8, bv_rotated(9), VEC_SET)
#endif
  endif
     
end subroutine set_normal_bc
   
   
!======================================================================
! set_laplacian_bc
!======================================================================
subroutine set_laplacian_bc(ibegin,rhs,bv,normal,curv,izonedim,radius,mat)
  use basic
  use vector_mod
  use matrix_mod

  implicit none

  integer, intent(in) :: ibegin      ! first dof of field
  real, intent(in) :: normal(2), curv
  type(vector_type), intent(inout) :: rhs ! right-hand-side of equation
  vectype, intent(in), dimension(dofs_per_node) :: bv
  integer, intent(in) :: izonedim    ! dimension of boundary
  real, intent(in) :: radius         ! radial coordinate of node
                                     ! for gs operator, let radius -> -radius
  type(matrix_type), optional :: mat
  
  integer :: numvals, irow
  integer, dimension(4) :: cols
  vectype, dimension(4) :: vals

  if(itor.eq.1) then
     numvals = 4
     cols(1) = ibegin + 1
     cols(2) = ibegin + 2
     cols(3) = ibegin + 3
     cols(4) = ibegin + 5
     vals(1) =  normal(1)/radius
     vals(2) = -normal(2)/radius
     vals(3) =  1.
     vals(4) =  1.
  else
     numvals = 2
     cols(1) = ibegin + 3
     cols(2) = ibegin + 5
     vals(1) = 1.
     vals(2) = 1.
  endif

  if(izonedim.eq.1) then
     ! edge points
     ! ~~~~~~~~~~~
     irow = ibegin + 3
     if(present(mat)) &
          call set_row_vals(mat, irow, numvals, cols, vals)
     call insert(rhs, irow, 0., VEC_SET)
  end if
  
end subroutine set_laplacian_bc


!=======================================================
! boundary_dc
! ~~~~~~~~~~~
!
! sets homogeneous dirichlet boundary condition
!=======================================================
subroutine boundary_dc(rhs, bvec, mat)
  use basic
  use vector_mod
  use matrix_mod

  implicit none
  
  type(vector_type) :: rhs
  type(vector_type), intent(in) :: bvec
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim
  integer :: ibegin, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_dc called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     ibegin = node_index(rhs, i, 1)

     call get_node_data(bvec, 1, i, temp)
     call set_dirichlet_bc(ibegin,rhs,temp,normal,curv,izonedim,mat)
  end do

  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_dc done"

end subroutine boundary_dc


!=======================================================
! boundary_nm
! ~~~~~~~~~~~
!
! sets homogeneous Neumann boundary condition
!=======================================================
subroutine boundary_nm(rhs, bvec, mat)
  use basic
  use vector_mod
  use matrix_mod

  implicit none

  type(vector_type) :: rhs
  type(vector_type), intent(in) :: bvec
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim
  integer :: ibegin, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_nm called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     ibegin = node_index(rhs, i, 1)

     call get_node_data(bvec, 1, i, temp)
     call set_normal_bc(ibegin,rhs,temp,normal,curv,izonedim,mat)
  end do

end subroutine boundary_nm


!=======================================================
! boundary_cy
! ~~~~~~~~~~~
!
! sets cauchy boundary condition
!=======================================================
subroutine boundary_cy(rhs, mat)
  use basic
  use vector_mod
  use matrix_mod

  implicit none
  
  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, ind, numnodes
  real :: normal(2), curv, x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_cy called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     ind = node_index(rhs, i, 1)

     temp = 0.
     call set_dirichlet_bc(ind,rhs,temp,normal,curv,izonedim,mat)
     temp = 0.
     call set_normal_bc(ind,rhs,temp,normal,curv,izonedim,mat)
  end do

end subroutine boundary_cy



!=======================================================
! boundary_vor
! ~~~~~~~~~~~~
!
! sets boundary conditions on Delta*(phi) 
! in the smoother
!=======================================================
subroutine boundary_vor(rhs, mat)
  use basic
  use vector_mod
  use matrix_mod

  implicit none

  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer, parameter :: numvarsm = 2
  integer :: i, izone, izonedim, i_u, i_vor, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vor called"

  temp = 0.

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_vor = node_index(rhs, i, 1)
     i_u = node_index(rhs, i, 2)

     if(inonormalflow.eq.1) then
        call set_dirichlet_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
     end if
     
     if(inoslip_pol.eq.1) then
        call set_normal_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
     end if

     ! no vorticity
     call set_dirichlet_bc(i_vor,rhs,temp,normal,curv,izonedim,mat)

     if(vor_bc.eq.1) then
        call set_laplacian_bc(i_u,rhs,temp,normal,curv,izonedim,-x,mat)
     endif
  end do

end subroutine boundary_vor



!=======================================================
! boundary_jphi
! ~~~~~~~~~~~~~
!
! sets boundary conditions on Delta*(phi) 
! in the smoother
!=======================================================
subroutine boundary_jphi(rhs, mat)
  use basic
  use arrays
  use vector_mod
  use matrix_mod

  implicit none

  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer, parameter :: numvarsm = 2
  integer :: i, izone, izonedim, i_psi, i_jphi, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_jphi called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_jphi = node_index(rhs, i, 1)
     i_psi = node_index(rhs, i, 2)

     call get_node_data(psi_field(1), i, temp)
     call set_dirichlet_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)

     ! no vorticity
     temp = 0.
     call set_dirichlet_bc(i_jphi,rhs,temp,normal,curv,izonedim,mat)
  end do

end subroutine boundary_jphi


!=======================================================
! boundary_com
! ~~~~~~~~~~~~
!
! sets boundary conditions on Del^2(chi) in the smoother
!=======================================================
subroutine boundary_com(rhs, mat)
  use basic
  use vector_mod
  use matrix_mod

  implicit none

  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer, parameter :: numvarsm = 2
  integer :: i, izone, izonedim, i_com, i_chi, numnodes
  real :: normal(2), curv
  real :: x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_com called"

  temp = 0.

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_com = node_index(rhs, i, 1)
     i_chi = node_index(rhs, i, 2)

     ! clamp compression
     call set_dirichlet_bc(i_com,rhs,temp,normal,curv,izonedim,mat)

     if(inonormalflow.eq.1) then
        call set_normal_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
     end if

     if(inoslip_pol.eq.1) then
        call set_dirichlet_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
     end if

     ! no compression
     if(com_bc.eq.1) then
        call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,x,mat)
     endif
  end do

end subroutine boundary_com

end module boundary_conditions
