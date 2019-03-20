!====================
! M3D-C1 KPRAD module
!====================
module kprad_m3dc1

  use kprad
  use field

  implicit none

  type(field_type), allocatable :: kprad_n(:)
  type(field_type), allocatable, private :: kprad_temp(:)
  type(field_type) :: kprad_rad      ! power lost to line radiation
  type(field_type) :: kprad_brem     ! power lost to bremsstrahlung
  type(field_type) :: kprad_ion      ! power lost to ionization
  type(field_type) :: kprad_reck     ! power lost to recombination (kinetic)
  type(field_type) :: kprad_recp     ! power lost to recombination (potential)
  type(field_type) :: kprad_sigma_e  ! electron source / sink due to ionization / recomb
  type(field_type) :: kprad_sigma_i  ! total ion source / sink due to ionization / recomb

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
    call create_field(kprad_brem)
    call create_field(kprad_ion)
    call create_field(kprad_reck)
    call create_field(kprad_recp)
    call create_field(kprad_sigma_e)
    call create_field(kprad_sigma_i)

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
       call destroy_field(kprad_brem)
       call destroy_field(kprad_ion)
       call destroy_field(kprad_reck)
       call destroy_field(kprad_recp)
       call destroy_field(kprad_sigma_e)
       call destroy_field(kprad_sigma_i)
    end if

    call kprad_deallocate()

  end subroutine kprad_destroy


  !==================================
  ! kprad_init_conds
  ! ~~~~~~~~~~~~~~~~
  !==================================
  subroutine kprad_init_conds()
    use basic
    use arrays
    use pellet
    use m3dc1_nint
    use newvar_mod

    implicit none

    integer :: itri, nelms, def_fields
    vectype, dimension(dofs_per_element) :: dofs
    real, dimension(MAX_PTS) :: p
    integer :: ip

    if(ikprad.eq.0) return

    kprad_n(0) = 0.

    def_fields = FIELD_N + FIELD_TE
    if(ipellet.lt.0. .and. ipellet_z.eq.kprad_z) &
         def_fields = def_fields + FIELD_P

    nelms = local_elements()
    do itri=1, nelms
       call define_element_quadrature(itri,int_pts_main,5)
       call define_fields(itri,def_fields,1,1,1)

       temp79a = kprad_nz +  kprad_fz*nt79(:,OP_1)

       if(ipellet.lt.0. .and. ipellet_z.eq.kprad_z) then
          p = pt79(:,OP_1)
          do ip=1,npellets
             temp79a = temp79a + &
                  pellet_rate(ip)*pellet_distribution(ip, x_79, phi_79, z_79, p, 1)
          end do
       end if

       dofs = intx2(mu79(:,:,OP_1),temp79a)
       call vector_insert_block(kprad_n(0)%vec,itri,1,dofs,VEC_ADD)
    end do

    call newvar_solve(kprad_n(0)%vec, mass_mat_lhs)

  end subroutine kprad_init_conds


  !=======================================================
  ! boundary_kprad
  ! ~~~~~~~~~~~~~~
  !
  ! sets boundary conditions for density
  !=======================================================
  subroutine boundary_kprad(rhs, den_v, mat)
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
    real :: normal(2), curv, x,z, phi
    logical :: is_boundary
    vectype, dimension(dofs_per_node) :: temp
    
    integer :: i_n
    
    if(iper.eq.1 .and. jper.eq.1) return
    if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_kprad called"
    
    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       i = nodes_owned(icounter_t)
       call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
            all_boundaries)
       if(.not.is_boundary) cycle
       
       i_n = node_index(den_v, i)
       
       if(inograd_n.eq.1) then
          temp = 0.
          call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
       end if
       if(iconst_n.eq.1) then
          if(idiff .gt. 0) then
             temp = 0.
          else
             call get_node_data(den_v, i, temp)
          end if
          call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
       end if
    end do
    
  end subroutine boundary_kprad
  

  subroutine kprad_advect(dti)
    use basic
    use matrix_mod
    use m3dc1_nint
    use boundary_conditions
    use metricterms_new
    use sparse

    implicit none

    real, intent(in) :: dti
    type(matrix_type) :: nmat_lhs, nmat_rhs
    type(field_type) :: rhs
    integer :: itri, j, numelms, ierr, def_fields, izone
!    integer, dimension(dofs_per_element) :: imask
    vectype, dimension(dofs_per_element) :: tempx
    vectype, dimension(dofs_per_element,dofs_per_element) :: tempxx
    vectype, dimension(dofs_per_element,dofs_per_element) :: ssterm, ddterm

    if(ikprad.eq.0) return

    if(myrank.eq.0 .and. iprint.ge.1) print *, ' In kprad_advect'

    call set_matrix_index(nmat_lhs, kprad_lhs_index)
    call set_matrix_index(nmat_rhs, kprad_rhs_index)
    call create_mat(nmat_lhs, 1, 1, icomplex, 1)
    call create_mat(nmat_rhs, 1, 1, icomplex, 0)
    call clear_mat(nmat_lhs)
    call clear_mat(nmat_rhs)
    call create_field(rhs)


    def_fields = FIELD_PHI + FIELD_V + FIELD_CHI

    if(myrank.eq.0 .and. iprint.ge.2) print *, '  populating matrix'

    numelms = local_elements()
    do itri=1, numelms

       call get_zone(itri, izone)

       call define_element_quadrature(itri, int_pts_main, 5)
       call define_fields(itri, def_fields, 1, linear)
       
       tempxx = n1n(mu79,nu79)
       ssterm = tempxx
       ddterm = tempxx

       if(izone.ne.1) goto 100

       do j=1, dofs_per_element
          ! NUMVAR = 1
          ! ~~~~~~~~~~      
          tempx = n1ndenm(mu79,nu79(j,:,:),denm,vzt79) &
               +  n1nu   (mu79,nu79(j,:,:),pht79)
          ssterm(:,j) = ssterm(:,j) -     thimp *dti*tempx
          ddterm(:,j) = ddterm(:,j) + (1.-thimp)*dti*tempx
          
#if defined(USECOMPLEX) || defined(USE3D)
          ! NUMVAR = 2
          ! ~~~~~~~~~~
          if(numvar.ge.2) then
             tempx = n1nv(mu79,nu79(j,:,:),vzt79)
             ssterm(:,j) = ssterm(:,j) -     thimp *dti*tempx
             ddterm(:,j) = ddterm(:,j) + (1.-thimp)*dti*tempx
          endif
#endif
          
          ! NUMVAR = 3
          ! ~~~~~~~~~~
          if(numvar.ge.3) then
             tempx = n1nchi(mu79,nu79(j,:,:),cht79)
             ssterm(:,j) = ssterm(:,j) -     thimp *dti*tempx
             ddterm(:,j) = ddterm(:,j) + (1.-thimp)*dti*tempx
          endif
       end do

100    continue

!       call get_den_mask(itri, imask)
!       call apply_boundary_mask(itri, 0, ssterm, imask)
!       call apply_boundary_mask(itri, 0, ddterm, imask)

       call insert_block(nmat_lhs,itri,1,1,ssterm,MAT_ADD)
       call insert_block(nmat_rhs,itri,1,1,ddterm,MAT_ADD)
    end do

    if(myrank.eq.0 .and. iprint.ge.2) print *, '  finalizing'

!    call boundary_kprad(rhs%vec, kprad_n(0), nmat_lhs)
    call finalize(nmat_rhs)
    call finalize(nmat_lhs)
    
    if(myrank.eq.0 .and. iprint.ge.2) print *, '  solving'

    do j=1, kprad_z
       rhs = 0.
       call matvecmult(nmat_rhs, kprad_n(j)%vec, rhs%vec)
!       call boundary_kprad(rhs%vec, kprad_n(j))
       ierr = 0
       call newsolve(nmat_lhs, rhs%vec, ierr)
       kprad_n(j) = rhs
    end do

    if(myrank.eq.0 .and. iprint.ge.2) print *, '  destroying'

    call destroy_field(rhs)
    call destroy_mat(nmat_rhs)
    call destroy_mat(nmat_lhs)

    if(myrank.eq.0 .and. iprint.ge.2) print *, '  done'

  end subroutine kprad_advect


  !===========================
  ! kprad_ionize
  ! ~~~~~~~~~~~~
  !===========================
  subroutine kprad_ionize(dti)
    use math
    use basic
    use newvar_mod
    use m3dc1_nint
    use pellet

    implicit none

    real, intent(in) :: dti
    
    real :: dt_s
    real, dimension(MAX_PTS) :: ne, te, n0_old, p
    real, dimension(MAX_PTS,0:kprad_z) :: nz
    real, dimension(MAX_PTS) :: dw_brem
    real, dimension(MAX_PTS,0:kprad_z) :: dw_rad, dw_ion, dw_reck, dw_recp
    real, dimension(MAX_PTS) :: source    ! neutral particle source

    integer :: i, itri, nelms, def_fields, izone
    real, parameter :: min_te = .01
    real, parameter :: min_ne = 1e8
    vectype, dimension(dofs_per_element) :: dofs
    integer :: ip

    if(ikprad.ne.1) return

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Advancing KPRAD.'

    do i=0, kprad_z
       kprad_temp(i) = 0.
    end do
    kprad_rad = 0.
    kprad_brem = 0.
    kprad_ion = 0.
    kprad_reck = 0.
    kprad_recp = 0.
    kprad_sigma_e = 0.
    kprad_sigma_i = 0.
    source = 0.

    def_fields = FIELD_N + FIELD_TE
    if(ipellet.ge.1 .and. ipellet_z.eq.kprad_z) &
         def_fields = def_fields + FIELD_P

    if(myrank.eq.0 .and. iprint.ge.2) print *, ' populating matrix'

    ! Do ionization / recombination step
    nelms = local_elements()
    do itri=1, nelms
       call get_zone(itri, izone)

       call define_element_quadrature(itri,int_pts_main,5)
       call define_fields(itri,def_fields,1,0)

       ! evaluate impurity density
       do i=0, kprad_z
          call eval_ops(itri, kprad_n(i), ph079, rfac)
          nz(:,i) = ph079(:,OP_1)
       end do
       ne = net79(:,OP_1)

       if(ipellet.ge.1 .and. ipellet_z.eq.kprad_z) then
          p = pt79(:,OP_1)
          do ip=1,npellets
             source = source + pellet_rate(ip)*pellet_distribution(ip, x_79, phi_79, z_79, p, 1)
          end do
       end if

       n0_old = sum(nz(:,1:kprad_z),2)

       ! convert nz, ne, te, dt to cgs / eV
       nz = nz*n0_norm
       ne = ne*n0_norm
       source = source*n0_norm/t0_norm
       te = tet79(:,OP_1)*p0_norm/n0_norm / 1.6022e-12
       dt_s = dti*t0_norm

       where(te.lt.min_te .or. te.ne.te) te = min_te
       where(nz.lt.0.) nz = 0.
       where(ne.lt.min_ne) ne = min_ne

       ! advance densities at each integration point
       ! for one MHD timestep (dt_s)
       if(izone.eq.1) then
          call kprad_advance_densities(dt_s, MAX_PTS, kprad_z, &
               ne, te, nz, dw_rad, dw_brem, dw_ion, dw_reck, dw_recp, source)
       else
          dw_rad = 0.
          dw_brem = 0.
          dw_ion = 0.
          dw_reck = 0.
          dw_recp = 0.
       end if
       
       ! convert nz, dw_rad, dw_brem to normalized units
       ! nz is given in /cm^3
       nz = nz / n0_norm
       ne = ne / n0_norm

       ! dw is given in J/cm^3
       ! factor of 1e7 needed to convert J to erg
       dw_rad = dw_rad * 1.e7 / p0_norm
       dw_brem = dw_brem * 1.e7 / p0_norm
       dw_ion  = dw_ion * 1.e7 / p0_norm
       dw_reck = dw_reck * 1.e7 / p0_norm
       dw_recp = dw_recp * 1.e7 / p0_norm
       
       ! New charge state densities
       do i=0, kprad_z
          temp79a = nz(:,i)
          dofs = intx2(mu79(:,:,OP_1), temp79a)
          call vector_insert_block(kprad_temp(i)%vec,itri,1,dofs,VEC_ADD)
       end do

       ! Line Radiation
       temp79b = dw_rad(:,kprad_z) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_rad%vec, itri,1,dofs,VEC_ADD)

       ! Bremsstrahlung
       temp79b = dw_brem / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_brem%vec, itri,1,dofs,VEC_ADD)

       ! Ionization
       temp79b = dw_ion(:,kprad_z) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_ion%vec, itri,1,dofs,VEC_ADD)

       ! Recombination (kinetic)
       temp79b = dw_reck(:,kprad_z) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_reck%vec, itri,1,dofs,VEC_ADD)

       ! Recombination (potential)
       temp79b = dw_recp(:,kprad_z) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_recp%vec, itri,1,dofs,VEC_ADD)

       ! Electron source
       temp79c = ne - net79(:,OP_1)
       dofs = intx2(mu79(:,:,OP_1),temp79c) / dti
       call vector_insert_block(kprad_sigma_e%vec, itri,1,dofs,VEC_ADD)

       ! Total ion source
       temp79d = sum(nz(:,1:kprad_z),2) - n0_old
       dofs = intx2(mu79(:,:,OP_1),temp79d) / dti
       call vector_insert_block(kprad_sigma_i%vec, itri,1,dofs,VEC_ADD)
    end do

    if(myrank.eq.0 .and. iprint.ge.2) print *, ' solving'

    do i=0, kprad_z
       call newvar_solve(kprad_temp(i)%vec, mass_mat_lhs)
       kprad_n(i) = kprad_temp(i)
    end do
    call newvar_solve(kprad_rad%vec, mass_mat_lhs)
    call newvar_solve(kprad_brem%vec, mass_mat_lhs)
    call newvar_solve(kprad_ion%vec, mass_mat_lhs)
    call newvar_solve(kprad_reck%vec, mass_mat_lhs)
    call newvar_solve(kprad_recp%vec, mass_mat_lhs)
    call newvar_solve(kprad_sigma_e%vec, mass_mat_lhs)
    call newvar_solve(kprad_sigma_i%vec, mass_mat_lhs)

    call kprad_rebase_dt

    if(myrank.eq.0) print *, ' Done advancing KPRAD'
  end subroutine kprad_ionize

end module kprad_m3dc1
