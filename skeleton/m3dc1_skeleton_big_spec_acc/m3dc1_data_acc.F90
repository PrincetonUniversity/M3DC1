module m3dc1_data

  use field

  implicit none

  type(field_type) :: psi_field, eta_field

  real,parameter :: dt=0.1 

  integer :: myrank
  integer :: maxrank

  real :: current, volume
!acc declare create (psi_field, eta_field)
!acc declare copyin(dt) create(myrank,maxrank,current,volume)
!$acc declare copyin(dt) 
  
  contains

    subroutine allocate_data()
      use element
      use mesh
      implicit none

      integer :: elms

      elms = local_elements()

      allocate(gtri(coeffs_per_tri,dofs_per_tri,elms))
      allocate(gtri_old(coeffs_per_tri,dofs_per_tri,elms))
      allocate(htri(coeffs_per_dphi,dofs_per_dphi,elms))

!$acc enter data create(psi_field,eta_field)
      call create_field(psi_field)
      call create_field(eta_field)
    end subroutine allocate_data


    subroutine deallocate_data()
      use element
      implicit none
      
      deallocate(gtri, gtri_old, htri)

      call destroy_field(psi_field)
      call destroy_field(eta_field)
!$acc exit data delete(psi_field,eta_field)
    end subroutine deallocate_data

    subroutine initialize_fields()
      use element
      use mesh
      implicit none

      integer :: inode, numnodes
      real :: R, Phi, Z
      real, dimension(dofs_per_node) :: psi

      current = 0.
      volume = 0.

      numnodes = local_nodes()
      do inode=1, numnodes
         call get_node_pos(inode, R, Phi, Z)

         psi = 0.

         psi(1) = R**2 + Z**2
         psi(2) = 2.*R
         psi(3) = 2.*Z
         psi(4) = 2.
         psi(5) = 0.
         psi(6) = 2.

         call set_node_data(psi_field, inode, psi)
      end do
    end subroutine initialize_fields

end module m3dc1_data
