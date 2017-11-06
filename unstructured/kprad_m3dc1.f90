!====================
! M3D-C1 KPRAD module
!====================
module kprad_m3dc1

  use kprad
  use field

  implicit none

  type(field_type), allocatable :: kprad_n(:)
  type(field_type), allocatable, private :: kprad_temp(:)
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
  subroutine kprad_init(ierr)
    implicit none

    integer, intent(out) :: ierr

    integer :: i

    ierr = 0
    if(ikprad.eq.0) return
    
    call kprad_atomic_data_sub(kprad_z, ierr)
    if(ierr.ne.0) return
    
    allocate(kprad_n(0:kprad_z))
    allocate(kprad_temp(0:kprad_z))
    
    do i=0, kprad_z
       call create_field(kprad_n(i))
       call create_field(kprad_temp(i))
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

    if(ikprad.eq.0) return

    if(allocated(kprad_n)) then
       do i=0, kprad_z
          call destroy_field(kprad_n(i))
          call destroy_field(kprad_temp(i))
       end do
       deallocate(kprad_n, kprad_temp)
       call destroy_field(kprad_rad)
    end if

    call kprad_deallocate()

  end subroutine kprad_destroy


  !==================================
  ! kprad_init_conds
  ! ~~~~~~~~~~~~~~~~
  !==================================
  subroutine kprad_init_conds()
    use arrays

    implicit none

    if(ikprad.eq.0) return

    ! initial neutral impurity density is
    ! nz_0 = kprad_fz*ni + kprad_nz
    kprad_n(0) = den_field(0)
    call mult(kprad_n(0), kprad_fz)
    call add(kprad_n(0), kprad_nz)

  end subroutine kprad_init_conds


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
    
    real :: dt_s
    real, dimension(MAX_PTS) :: ne, te
    real, dimension(MAX_PTS,0:kprad_z) :: nz
    real, dimension(MAX_PTS) :: dw_brem
    real, dimension(MAX_PTS,0:kprad_z) :: dw_rad

    integer :: i, itri, nelms, def_fields
    vectype, dimension(dofs_per_element) :: dofs

    if(ikprad.ne.1) return

    if(myrank.eq.0) print *, 'Advancing KPRAD.'

    do i=0, kprad_z
       kprad_temp(i) = 0.
    end do
    kprad_rad = 0.

    def_fields = FIELD_N + FIELD_TE

    ! Do ionization / recombination step
    nelms = local_elements()
    do itri=1, nelms
       call define_element_quadrature(itri,int_pts_main,5)
       call define_fields(itri,def_fields,1,0)

       ! evaluate impurity density
       do i=0, kprad_z
          call eval_ops(itri, kprad_n(i), ph079, rfac)
          nz(:,i) = ph079(:,OP_1)
       end do
       where(nz.lt.0.) nz = 0.

       ! convert nz, ne, te, dt to cgs / eV
       nz = nz*n0_norm
       ne = net79(:,OP_1)*n0_norm
       te = tet79(:,OP_1)*p0_norm/n0_norm / 1.6022e-12
       dt_s = dt*t0_norm

       ! advance densities at each integration point
       ! for one MHD timestep (dt_s)
       call kprad_advance_densities(dt_s, MAX_PTS, kprad_z, &
            ne, te, nz, dw_rad, dw_brem)

       ! convert nz, dw_rad, dw_brem to normalized units
       nz = nz / n0_norm
       dw_rad = dw_rad / p0_norm
       dw_brem = dw_brem / p0_norm

       do i=0, kprad_z
          temp79a = nz(:,i)
          dofs = intx2(mu79(:,:,OP_1), temp79a)
          call vector_insert_block(kprad_temp(i)%vec,itri,1,dofs,VEC_ADD)
       end do
       temp79b = (dw_rad(:,kprad_z) + dw_brem) / dt
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_rad%vec, itri,1,dofs,VEC_ADD)
    end do

    do i=0, kprad_z
       call newvar_solve(kprad_temp(i)%vec, mass_mat_lhs)
       kprad_n(i) = kprad_temp(i)
    end do
    call newvar_solve(kprad_rad%vec, mass_mat_lhs)

    if(myrank.eq.0) print *, ' Done advancing KPRAD'
  end subroutine kprad_onestep

end module kprad_m3dc1
