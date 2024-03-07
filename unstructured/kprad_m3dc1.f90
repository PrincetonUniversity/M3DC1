!====================
! M3D-C1 KPRAD module
!====================
module kprad_m3dc1

  use kprad
  use field

  implicit none

  type(field_type), allocatable :: kprad_n(:), kprad_n_prev(:)
  type(field_type), allocatable, private :: kprad_temp(:)
  type(field_type) :: kprad_rad      ! power lost to line radiation
  type(field_type) :: kprad_brem     ! power lost to bremsstrahlung
  type(field_type) :: kprad_ion      ! power lost to ionization
  type(field_type) :: kprad_reck     ! power lost to recombination (kinetic)
  type(field_type) :: kprad_recp     ! power lost to recombination (potential)
  type(field_type) :: kprad_sigma_e  ! electron source / sink due to ionization / recomb
  type(field_type) :: kprad_sigma_i  ! total ion source / sink due to ionization / recomb

  integer :: kprad_z    ! Z of impurity species

  ! Model for advection/diffusion of neutrals
  ! 0 = no adv/diff, 1 = advect & diff, 2 = diffuse only
  integer :: ikprad_evolve_neutrals

  ! scaling factor for neutral particle diffusion
  real :: kprad_n0_denm_fac

  ! variables for setting initial conditions
  real :: kprad_fz
  real :: kprad_nz

  integer :: iread_lp_source
  real :: lp_source_dt
  real :: lp_source_mass
  real, allocatable :: lp_source_rate(:)
  type(field_type), allocatable :: kprad_particle_source(:)

contains

  !==================================
  ! kprad_init
  ! ~~~~~~~~~~
  !==================================
  subroutine kprad_init(ierr)
    use basic
#ifdef USEADAS
    use adas_m3dc1
#endif
    implicit none
    integer, intent(out) :: ierr
    integer :: i
    character(len=32) :: fname

    ierr = 0
    if (ikprad.eq.0) return

    call kprad_allocate(kprad_z)
    
    if (ikprad.eq.1) then
       call kprad_atomic_data_sub(kprad_z, ierr)
    elseif (ikprad.eq.-1) then
#ifdef USEADAS
       call load_adf11(kprad_z)
#endif
    end if
    if(ierr.ne.0) return
    
    allocate(kprad_n(0:kprad_z))
    if(isolve_with_guess==1) allocate(kprad_n_prev(0:kprad_z))
    allocate(kprad_temp(0:kprad_z))
    allocate(kprad_particle_source(0:kprad_z))
    allocate(lp_source_rate(0:kprad_z))

    if (ispradapt .eq. 1) then
       do i=0, kprad_z
          write(fname,"(A5,I2.2,A)")  "kprn", i, 0
          call create_field(kprad_n(i), trim(fname))
          write(fname,"(A5,I2.2,A)")  "kprt", i, 0
          call create_field(kprad_temp(i), trim(fname))
          write(fname,"(A5,I2.2,A)")  "kprp", i, 0
          call create_field(kprad_particle_source(i), trim(fname))
          kprad_particle_source(i) = 0.
       end do
       call create_field(kprad_rad, "kprad_rad")
       call create_field(kprad_brem, "kprad_brem")
       call create_field(kprad_ion, "kprad_ion")
       call create_field(kprad_reck, "kprad_reck")
       call create_field(kprad_recp, "kprad_recp")
       call create_field(kprad_sigma_e,"kprad_sigma_e")
       call create_field(kprad_sigma_i, "kprad_sigma_i")
    else
       do i=0, kprad_z
          call create_field(kprad_n(i))
          call create_field(kprad_temp(i))
          call create_field(kprad_particle_source(i))
          kprad_particle_source(i) = 0.
       end do

       if(isolve_with_guess==1) then
          do i=0, kprad_z
             call create_field(kprad_n_prev(i))
          enddo
       endif

       call create_field(kprad_rad)
       call create_field(kprad_brem)
       call create_field(kprad_ion)
       call create_field(kprad_reck)
       call create_field(kprad_recp)
       call create_field(kprad_sigma_e)
       call create_field(kprad_sigma_i)
    endif
    if(ikprad_min_option.eq.2 .or. ikprad_min_option.eq.3) then
       kprad_nemin = kprad_nemin*n0_norm
       kprad_temin = kprad_temin*p0_norm/n0_norm / 1.6022e-12
    end if

  end subroutine kprad_init
    
  !==================================
  ! kprad_destroy
  ! ~~~~~~~~~~~~~
  !==================================
  subroutine kprad_destroy
    use basic, only:isolve_with_guess
    implicit none

    integer :: i

    if(ikprad.eq.0) return

    if(allocated(kprad_n)) then
       do i=0, kprad_z
          call destroy_field(kprad_n(i))
          call destroy_field(kprad_temp(i))
          call destroy_field(kprad_particle_source(i))
       end do
       deallocate(kprad_n, kprad_temp)
       if(isolve_with_guess==1) then
          do i=0, kprad_z
             call destroy_field(kprad_n_prev(i))
          enddo
          deallocate(kprad_n_prev)
       endif       
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
    integer :: ip, izone

    if(ikprad.eq.0) return

    kprad_n(0) = 0.

    def_fields = FIELD_N + FIELD_TE
    if(ipellet.lt.0. .and. ipellet_z.eq.kprad_z) &
         def_fields = def_fields + FIELD_P

    nelms = local_elements()
    do itri=1, nelms
       call define_element_quadrature(itri,int_pts_main,5)
       call define_fields(itri,def_fields,1,1,1)
       call get_zone(itri, izone)

       temp79a = kprad_nz +  kprad_fz*nt79(:,OP_1)

       if(ipellet.lt.0. .and. ipellet_z.eq.kprad_z) then
          p = pt79(:,OP_1)
          do ip=1,npellets
             temp79a = temp79a + &
                  pellet_rate(ip)*pellet_distribution(ip, x_79, phi_79, z_79, p, 1, izone)
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
    real :: normal(2), curv(3), x,z, phi
    logical :: is_boundary
    vectype, dimension(dofs_per_node) :: temp
    
    integer :: i_n
    
    if(iper.eq.1 .and. jper.eq.1) return
    if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_kprad called"
    
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

    include 'mpif.h'

    real, intent(in) :: dti
    type(matrix_type) :: nmat_lhs, nmat_rhs
    type(field_type) :: rhs
    integer :: itri, j, numelms, ierr, def_fields, izone, itmp, ier
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

    def_fields = FIELD_PHI + FIELD_V + FIELD_CHI + FIELD_DENM

    if(myrank.eq.0 .and. iprint.ge.2) print *, '  populating matrix'

    ! Evolve charged impurities
    ! =========================
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
          tempx = n1ndenm(mu79,nu79(j,:,:),denm79,vzt79)
          ssterm(:,j) = ssterm(:,j) -     thimp *dti*tempx
          ddterm(:,j) = ddterm(:,j) + (1.-thimp)*dti*tempx
          
          tempx = n1nu(mu79,nu79(j,:,:),pht79)
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
       ierr = 0
       rhs = 0.
       call matvecmult(nmat_rhs, kprad_n(j)%vec, rhs%vec)
!       call boundary_kprad(rhs%vec, kprad_n(j))
       if(isolve_with_guess==1) then
         call newsolve_with_guess(nmat_lhs, rhs%vec, kprad_n_prev(j)%vec, ierr)
       else
         call newsolve(nmat_lhs, rhs%vec, ierr)
       endif 
       if(is_nan(rhs%vec)) ierr = 1
       call mpi_allreduce(ierr, itmp, 1, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ier)
       ierr = itmp

       if(ierr.ne.0) then
          if(myrank.eq.0) &
               print *, 'Error in impurity ion solve ', j
       else
          kprad_n(j) = rhs
       end if
       if(isolve_with_guess==1) kprad_n_prev(j) = kprad_n(j)
    end do


    ! Evolve neutrals
    ! ===============
    if(ikprad_evolve_neutrals.ge.1) then
       call clear_mat(nmat_lhs)
       call clear_mat(nmat_rhs)
    
       if(myrank.eq.0 .and. iprint.ge.2) print *, '  populating neutral matrix'
    
       do itri=1, numelms
       
          call get_zone(itri, izone)
       
          call define_element_quadrature(itri, int_pts_main, 5)
          call define_fields(itri, def_fields, 1, linear)
          
          tempxx = n1n(mu79,nu79)
          ssterm = tempxx
          ddterm = tempxx
          
          if(izone.ne.1) goto 200
    
          do j=1, dofs_per_element
             tempx = n1ndenm(mu79,nu79(j,:,:),denm79,vzt79) &
                  * kprad_n0_denm_fac
             ssterm(:,j) = ssterm(:,j) -     thimp *dti*tempx
             ddterm(:,j) = ddterm(:,j) + (1.-thimp)*dti*tempx
          end do

          if(ikprad_evolve_neutrals.eq.2) goto 200

          do j=1, dofs_per_element
             tempx = n1nu(mu79,nu79(j,:,:),pht79)
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
          
200       continue

          call insert_block(nmat_lhs,itri,1,1,ssterm,MAT_ADD)
          call insert_block(nmat_rhs,itri,1,1,ddterm,MAT_ADD)
       end do
    
       if(myrank.eq.0 .and. iprint.ge.2) print *, '  finalizing neutral'
       
       call finalize(nmat_rhs)
       call finalize(nmat_lhs)
       
       if(myrank.eq.0 .and. iprint.ge.2) print *, '  solving neutral'
       
       rhs = 0.
       call matvecmult(nmat_rhs, kprad_n(0)%vec, rhs%vec)
       ierr = 0
       if(isolve_with_guess==1) then
         call newsolve_with_guess(nmat_lhs, rhs%vec, kprad_n_prev(0)%vec, ierr)
       else
         call newsolve(nmat_lhs, rhs%vec, ierr)
       endif 
       if(is_nan(rhs%vec)) ierr = 1
       call mpi_allreduce(ierr, itmp, 1, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ier)
       ierr = itmp

       if(ierr.ne.0) then
          if(myrank.eq.0) &
               print *, 'Error in impurity neutral solve'
       else
          kprad_n(0) = rhs          
       end if
       if(isolve_with_guess==1) kprad_n_prev(0) = kprad_n(0)
    end if


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
    real, dimension(MAX_PTS) :: ne, te, den, ti, n0_old, p
    real, dimension(MAX_PTS,0:kprad_z) :: nz, nz_nokprad
    real, dimension(MAX_PTS) :: dw_brem
    real, dimension(MAX_PTS,0:kprad_z) :: dw_rad, dw_ion, dw_reck, dw_recp
    real, dimension(MAX_PTS,0:kprad_z) :: source    ! particle source
    logical, dimension(MAX_PTS) :: advance_kprad

    integer :: i, itri, nelms, def_fields, izone
    vectype, dimension(dofs_per_element) :: dofs
    integer :: ip
    character(len=5):: imp_model

    if(ikprad.eq.0) return

    if(ikprad.eq.1) then
       imp_model = 'KPRAD'
    elseif(ikprad.eq.-1) then
       imp_model = 'ADAS'
    end if

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Advancing ', trim(imp_model)

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

    def_fields = FIELD_N + FIELD_P + FIELD_TE + FIELD_TI + FIELD_DENM

    if(myrank.eq.0 .and. iprint.ge.2) print *, ' populating matrix'

    ! Do ionization / recombination step
    nelms = local_elements()
    do itri=1, nelms
       source = 0.

       call get_zone(itri, izone)

       call define_element_quadrature(itri,int_pts_main,5)
       call define_fields(itri,def_fields,1,0)

       ! evaluate impurity density
       do i=0, kprad_z
          call eval_ops(itri, kprad_n(i), ph079, rfac)
          nz(:,i) = ph079(:,OP_1)
       end do

       ne = net79(:,OP_1)
       te = tet79(:,OP_1)
       den = nt79(:,OP_1)
       ti = tit79(:,OP_1)
       if(ikprad_te_offset .gt. 0) te = te - eta_te_offset
       p = pt79(:,OP_1)

       if(ipellet.ge.1 .and. ipellet_z.eq.kprad_z .and. iread_lp_source.eq.0) then
          p = pt79(:,OP_1)
          source = 0.
          do ip=1,npellets
             source(:,0) = source(:,0) + pellet_rate(ip)*pellet_distribution(ip, x_79, phi_79, z_79, p, 1, izone)
          end do
       end if

       if(iread_lp_source.eq.1) then
          do i=0, kprad_z
             call eval_ops(itri, kprad_particle_source(i), ch079, rfac)
             source(:,i) = source(:,i) + ch079(:,OP_1)
          end do
       else if (iread_lp_source.eq.2) then
          p = pt79(:,OP_1)
          do i=0, kprad_z
             ! Deposit over distribution of pellet #1
             source(:,i) = source(:,i) + lp_source_rate(i)*pellet_distribution(1, x_79, phi_79, z_79, p, 1, izone)
          end do
       end if

       n0_old = sum(nz(:,1:kprad_z),2)

       where(nz.lt.0.) nz = 0.

       ! determine where KPRAD advance will be used
       ! and impurity densities if charge states don't advance
       !  old nz (with source added)
       if (ikprad_min_option.eq.1) then
          advance_kprad = .not.(te.lt.kprad_temin .or. te.ne.te .or. &
                                ne.lt.kprad_nemin .or. ne.ne.ne)
       else
          advance_kprad = .not.(te.ne.te .or. ne.ne.ne)
       end if
       nz_nokprad = nz
       nz_nokprad = nz_nokprad + dti*source

       if(any(advance_kprad)) then ! skip if no KPRAD advance in this triangle

          ! convert nz, ne, te, dt to cgs / eV
          p = p*p0_norm
          nz = nz*n0_norm
          ne = ne*n0_norm
          den = den*n0_norm
          source = source*n0_norm/t0_norm
          te = te*p0_norm/n0_norm / 1.6022e-12
          ti = ti*p0_norm/n0_norm / 1.6022e-12
          dt_s = dti*t0_norm
       
          ! advance densities at each integration point
          ! for one MHD timestep (dt_s)
          if(izone.eq.1) then
             call kprad_advance_densities(dt_s, MAX_PTS, kprad_z, p, ne, &
                  te, den, ti, nz, dw_rad, dw_brem, dw_ion, dw_reck, dw_recp, source)
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

       end if

       ! New charge state densities
       ! Old nz or source added to neutrals if not advancing KPRAD at that point
       do i=0, kprad_z
          temp79a = merge(nz(:,i), nz_nokprad(:,i), advance_kprad)
          dofs = intx2(mu79(:,:,OP_1), temp79a)
          call vector_insert_block(kprad_temp(i)%vec,itri,1,dofs,VEC_ADD)
       end do

       ! Line Radiation (0 if not advancing KPRAD at that point)
       temp79b = merge(dw_rad(:,kprad_z), 0., advance_kprad) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
!Check for and delete spurius values   (scj 6/21/20)
       do i=1,MAX_PTS
          if(abs(temp79b(i)) .gt. 1.e0) then
            temp79b(i) = 0.
          endif
       enddo
!
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_rad%vec, itri,1,dofs,VEC_ADD)

       ! Bremsstrahlung  (0 if not advancing KPRAD at that point)
       temp79b = merge(dw_brem, 0., advance_kprad) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_brem%vec, itri,1,dofs,VEC_ADD)

       ! Ionization (0 if not advancing KPRAD at that point)
       temp79b = merge(dw_ion(:,kprad_z), 0., advance_kprad) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_ion%vec, itri,1,dofs,VEC_ADD)

       ! Recombination (kinetic) (0 if not advancing KPRAD at that point)
       temp79b = merge(dw_reck(:,kprad_z), 0., advance_kprad) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_reck%vec, itri,1,dofs,VEC_ADD)

       ! Recombination (potential) (0 if not advancing KPRAD at that point)
       temp79b = merge(dw_recp(:,kprad_z), 0., advance_kprad) / dti
       where(temp79b.ne.temp79b) temp79b = 0.
       dofs = intx2(mu79(:,:,OP_1),temp79b)
       call vector_insert_block(kprad_recp%vec, itri,1,dofs,VEC_ADD)

       ! Electron source (0 if not advancing KPRAD at that point)
       temp79c = merge(ne - real(net79(:,OP_1)), 0., advance_kprad)
       dofs = intx2(mu79(:,:,OP_1),temp79c) / dti
       call vector_insert_block(kprad_sigma_e%vec, itri,1,dofs,VEC_ADD)

       ! Total ion source (0 if not advancing KPRAD at that point)
       temp79d = merge(sum(nz(:,1:kprad_z),2) - n0_old, 0., advance_kprad)
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

    if(myrank.eq.0) print *, ' Done advancing ', trim(imp_model)
  end subroutine kprad_ionize

  subroutine read_lp_source(filename, ierr)
    use basic
    use read_ascii
    use pellet
    use newvar_mod
    use math

    implicit none

    include 'mpif.h'

    integer :: mpierr
    integer, intent(out) :: ierr
    character(len=*), intent(in) :: filename
    integer, parameter :: ifile = 113
    integer :: n, i, j, m
    real, allocatable :: x_vals(:), y_vals(:), z_vals(:), phi_vals(:), temp_vals(:)
    real, allocatable :: n_vals(:,:)
    real :: N_per_LP, ntot
    

    if(iprint.ge.1 .and. myrank.eq.0) print *, 'Reading LP source...'
    ierr = 0

    if(ikprad.eq.0 .or. kprad_z.le.0) then
       if(myrank.eq.0) print *, 'Error: kprad must be enabled to read LP source'
       call safestop(8)
    end if
       
    ! Read the data
    if(iprint.ge.2 .and. myrank.eq.0) print *, ' reading file ', filename

    ! Read lp_source_dt and lp_source_mass
    if(myrank.eq.0) then
       open(unit=ifile, file=filename, status='old', action='read')
       read(ifile, *) lp_source_dt, lp_source_mass
       close(ifile)
    end if
    call MPI_bcast(lp_source_dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
    call MPI_bcast(lp_source_mass,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)
    
    n = 0
    call read_ascii_column(filename, x_vals, n, skip=1, icol=1)
    call read_ascii_column(filename, y_vals, n, skip=1, icol=2)
    call read_ascii_column(filename, z_vals, n, skip=1, icol=3)
    allocate(n_vals(n,0:kprad_z))
    do i=0, kprad_z
       if(iprint.ge.2 .and. myrank.eq.0) print *, ' reading column', 12+i
       call read_ascii_column(filename, temp_vals, n, skip=1, icol=12+i)
       n_vals(:,i) = temp_vals
       deallocate(temp_vals)
    end do
    if(n.eq.0) then
       if(myrank.eq.0) print *, 'Error: no data in LP source file ', filename
       ierr = 1
       return
    end if
    allocate(phi_vals(n))

    ! convert from cgs to normalized units
    x_vals = x_vals / l0_norm
    y_vals = y_vals / l0_norm
    z_vals = z_vals / l0_norm
    n_vals = n_vals / n0_norm
    lp_source_dt = lp_source_dt / t0_norm

    ! Each LP has same impurity mass, defined by lp_source_mass
    N_per_LP = (N_Avo/(n0_norm*l0_norm**3))*lp_source_mass/M_table(kprad_z)
    if(iprint.ge.1 .and. myrank.eq.0) then
       print *, 'PELLET N_per_LP: ', N_per_LP
    end if
    do j=1, n
       ntot = sum(n_vals(j,:))
       n_vals(j,:) = n_vals(j,:)*N_per_LP/ntot
    end do

    if(iread_lp_source.eq.1) then

       ! Read LP distribution directly
       ! Assumes a single pellet

       ! convert from local to (R,phi,Z)
       x_vals = x_vals + pellet_r(1)
       y_vals = y_vals + pellet_z(1)

       ! convert z from length to angle
       phi_vals = z_vals / x_vals + pellet_phi(1)
       where(phi_vals.lt.0.) phi_vals = phi_vals + toroidal_period
       where(phi_vals.gt.toroidal_period) phi_vals = phi_vals - toroidal_period

       if(iprint.ge.1 .and. myrank.eq.0) then
          print *, 'PELLET sum(n_vals): ', sum(n_vals)
          print *, 'PELLET lp_source_dt: ', lp_source_dt
       end if
    
       ! convert density to rate
       n_vals = n_vals / lp_source_dt
#ifndef USE3D
       ! need density rate per radian in 2D
       n_vals = n_vals / toroidal_period
#endif

       ! construct fields using data
       if(iprint.ge.2 .and. myrank.eq.0) print *, ' constructing fields'
       do i=0, kprad_z
          kprad_particle_source(i) = 0.
       end do
       m = kprad_z+1
       call deltafuns(n,x_vals,phi_vals,y_vals,m,n_vals,kprad_particle_source,ierr)
       if(iprint.ge.2 .and. myrank.eq.0) print *, ' solving fields'
       do i=0, kprad_z
          call newvar_solve(kprad_particle_source(i)%vec,mass_mat_lhs)
       end do

       ! free temporary arrays
       if(iprint.ge.2 .and. myrank.eq.0) print *, ' freeing data'
       deallocate(x_vals, y_vals, z_vals, n_vals, phi_vals)

    else if(iread_lp_source.eq.2) then
       ! Sum LP density to get impurity rate
       lp_source_rate = 0.
       do j=1, n
          lp_source_rate = lp_source_rate + n_vals(j,:)
       end do

       ! convert particle # to rate
       lp_source_rate = lp_source_rate / lp_source_dt

    end if
    
    if(iprint.ge.1 .and. myrank.eq.0) print *, 'Done reading LP source'

  end subroutine read_lp_source

  

! ===========================================================
! deltafuns
! ~~~~~~~~~
! sets jout_i =  <mu_i | -val*delta(R-x)delta(Z-z)> 
! ===========================================================
subroutine deltafuns(n,x,phi,z,m,val,jout, ier)

  use mesh_mod
  use basic
  use arrays
  use field
  use m3dc1_nint
  use math

  implicit none

  include 'mpif.h'

  integer, intent(in) :: n, m
  real, intent(in), dimension(n) :: x, phi, z
  real, intent(in), dimension(n,m) :: val
  type(field_type), intent(inout), dimension(m) :: jout
  integer, intent(out) :: ier

  type(element_data) :: d
  integer, dimension(n) :: itri, in_domain, in_domains
  integer :: i, j, k, it
  real, dimension(n) :: x1, z1
  real :: si, zi, eta
  vectype, dimension(dofs_per_element) :: temp, temp2
  real, dimension(dofs_per_element,coeffs_per_element) :: c

  ier = 0
  itri = 0
  it = 0
  do j=1, n
     call whattri(x(j), phi(j), z(j), it, x1(j), z1(j))
     itri(j) = it
  end do
  where(itri.gt.0)
     in_domain = 1
  elsewhere 
     in_domain = 0
  end where

  call mpi_allreduce(in_domain,in_domains,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
  if(ier.ne.0) return

  ! Check to make sure all points are in domain
  do j=1, n
     if(in_domains(j).eq.0) then
        if(myrank.eq.0) &
             print *, 'Error: point not found in domain', j, x(j), phi(j), z(j)
        ier = 1
        return
     end if
  end do

  do j=1, n
     if(itri(j).le.0) cycle

     ! calculate local coordinates
     call get_element_data(itri(j), d)
     call global_to_local(d, x(j), phi(j), z(j), si, zi, eta)
        
     ! calculate temp_i = val*mu_i(si,zi,eta)
     if(iprecompute_metric.eq.1) then
        c = ctri(:,:,itri(j))
     else
        call local_coeff_vector(itri(j), c)
     endif 
     temp = 0.
     do k=1, coeffs_per_tri
#ifdef USE3D
        temp(:) = temp(:) + c(:,k)*si**mi(k)*eta**ni(k)*zi**li(k)
#else
        temp(:) = temp(:) + c(:,k)*si**mi(k)*eta**ni(k)
#endif
     end do
     if(equilibrate.ne.0) temp(:) = temp(:)*equil_fac(:, itri(j))

     do i=1, m
        temp2 = temp*val(j,i)/in_domains(j)
        call vector_insert_block(jout(i)%vec, itri(j), jout(i)%index, temp2, VEC_ADD)
     end do
  end do
end subroutine deltafuns


end module kprad_m3dc1
