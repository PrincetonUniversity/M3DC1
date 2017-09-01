!====================
! M3D-C1 KPRAD module
!====================
module kprad_m3dc1

  use kprad
  use field

  implicit none

  type(field_type), allocatable :: kprad_n(:)
  type(field_type) :: kprad_rad

  integer :: ikprad     ! 1 = use kprad model
  integer :: kprad_z    ! Z of impurity species

  ! variables for setting initial conditions
  real :: kprad_fz
  real :: kprad_nz

contains

  !==================================
  ! kprad_init
  ! ~~~~~~~~~~
  !==================================
  subroutine kprad_init(z, ierr)
    implicit none

    integer, intent(in) :: z
    integer, intent(out) :: ierr

    integer :: i
    
    kprad_z = z

    call kprad_atomic_data_sub(kprad_z, ierr)
    if(ierr.ne.0) return
    
    allocate(kprad_n(0:kprad_z))
    
    do i=0, kprad_z
       call create_field(kprad_n(i))
    end do
    call create_field(kprad_rad)

  end subroutine kprad_init
    
  !==================================
  ! kprad_destroy
  ! ~~~~~~~~~~~~~
  !==================================
  subroutine kprad_destroy
    implicit none

    integer :: i

    if(allocated(kprad_n)) then
       do i=0, kprad_z
          call destroy_field(kprad_n(i))
       end do
       deallocate(kprad_n)
       call destroy_field(kprad_rad)
    end if

    call kprad_deallocate()

  end subroutine kprad_destroy


  !===========================
  ! kprad_onestep
  ! ~~~~~~~~~~~~~
  !===========================
  subroutine kprad_onestep()
    use math
    use basic
    use newvar_mod
    use m3dc1_nint

    implicit none
    
    real, dimension(MAX_PTS) :: ne, te
    real, dimension(MAX_PTS,0:kprad_z) :: nz
    real, dimension(MAX_PTS) :: dw_brem
    real, dimension(MAX_PTS,0:kprad_z) :: dw_rad

    integer :: i, itri, nelms
    vectype, dimension(dofs_per_element) :: dofs

    ! convert ne, te to cm^-3, 
    ne = net79(:,OP_1)*n0_norm
    te = tet79(:,OP_1)*p0_norm/n0_norm * 1.6022e-12

    call kprad_advance_densities(dt, MAX_PTS, kprad_z, &
         ne, te, nz, dw_rad, dw_brem)

    ! convert nz, dw_rad, dw_brem to normalized units
    nz = nz / n0_norm
    dw_rad = dw_rad / (b0_norm**2 / (4.*pi))
    dw_brem = dw_brem / (b0_norm**2 / (4.*pi))

    ! Do ionization / recombination step
    nelms = local_elements()
    do i=0, kprad_z
       kprad_n(i) = 0.

       do itri=i, nelms
          temp79a = nz(:,i)
          dofs = intx2(mu79(:,:,OP_1), temp79a)
          call vector_insert_block(kprad_n(i)%vec,itri,1,dofs,VEC_ADD)
       end do
       
       call newvar_solve(kprad_n(i)%vec, mass_mat_lhs)
    end do
  end subroutine kprad_onestep

end module kprad_m3dc1
