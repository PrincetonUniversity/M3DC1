! Kinetic energetic ion module, J. Breslau, 2015
module particles
   use mesh_mod
   use field
   use matrix_mod
   use newvar_mod
   use gradshafranov
   !use mpi_f08
   use ieee_arithmetic
   implicit none
   private
#ifdef USEPARTICLES
   public :: particle_test, particle_step, particle_scaleback, finalize_particles, hdf5_read_particles, hdf5_write_particles
   public :: set_parallel_velocity
#endif

   real, parameter :: e_mks = 1.6022e-19      !Elementary charge, in Coulombs
   real, parameter :: m_proton = 1.6726e-27   !Proton mass in kg
   real, parameter :: A_deuteron = 1.999      !Deuteron/proton mass ratio
   real, parameter :: A_triton = 2.9942615    !Triton/proton mass ratio
   real, parameter :: A_alpha = 3.972599689   !alpha/proton mass ratio
   real, parameter :: qm_proton = e_mks/m_proton  !Proton charge/mass ratio in C/kg
   real, parameter :: vp1eV = 1.3841e+04      !Thermal speed of a 1 eV proton in m/s
   real, parameter :: twopi = 6.283185307179586476925286766559
   integer, parameter :: vspdims = 2          !Dimensions in velocity space (2=gk; 3=full)
#ifdef USE3D
   integer, parameter :: nneighbors = 5       !Max # of nearest neighbors of a prism
#else
   integer, parameter :: nneighbors = 3       !Max # of nearest neighbors of a triangle
#endif
   integer, parameter :: nglayers = 2 !# of layers of ghost zones around each mesh domain

   type elfield
      sequence
      vectype, dimension(coeffs_per_element) :: psiv0, psiv1, Bzv0, Bzv1
      vectype, dimension(coeffs_per_element) :: Bfpv0, Bfpv1
      vectype, dimension(coeffs_per_element) :: er, ephi, ez
      vectype, dimension(coeffs_per_element) :: rho, nf, tf
      vectype, dimension(coeffs_per_element) :: nfi, tfi
      vectype, dimension(coeffs_per_element) :: B0, B1, B1_last
      vectype, dimension(coeffs_per_element) :: pe, ne, pe0, ne0, te0
      vectype, dimension(coeffs_per_element) :: rst, zst
   end type elfield

   type xgeomterms
      sequence
      real, dimension(coeffs_per_element) :: g, dr, dz
      real, dimension(coeffs_per_element) :: drr, drz, dzz
#ifdef USE3D
      real, dimension(coeffs_per_element) :: dphi, drphi, dzphi
      !real, dimension(coeffs_per_element) :: drrphi, drzphi, dzzphi
      !real, dimension(coeffs_per_element) :: drphiphi, dzphiphi
#endif
   end type xgeomterms

   type particle
      sequence
      real, dimension(3)       :: x           !Position in cylindrical coords
      real, dimension(vspdims) :: v           !Velocity
      real                     :: wt    !Particle weighting in delta-f scheme
      real                     :: dB      !Time
      real                     :: B0
      real                     :: dE
      integer                  :: sps
      integer                  :: gid         !Unique global particle index
      integer                  :: jel         !Predicted element of residence
      real, dimension(3, 4)     :: kx
      integer, dimension(4)     :: kel
      real                     :: f0
      real, dimension(3)       :: x0           !Position in cylindrical coords
      real, dimension(vspdims) :: v0           !Velocity
      logical                  :: deleted
   end type particle

   type(particle), dimension(:), pointer :: pdata  !Particle arrays
   integer :: win_pdata
   !type(elfield), dimension(:), allocatable, target :: elfieldcoefs        !Receive buffer for jumping particles
   type(elfield), dimension(:), pointer :: elfieldcoefs       !Receive buffer for jumping particles
!$acc declare link(elfieldcoefs)
   integer :: win_elfieldcoefs
   real, dimension(2) :: m_ion, q_ion, qm_ion
!$acc declare create(m_ion,q_ion,qm_ion)
   real, dimension(2) :: nrmfac
   integer :: particle_linear_particle, psubsteps_particle, kinetic_thermal_ion_particle
!$acc declare create(particle_linear_particle, psubsteps_particle, kinetic_thermal_ion_particle)
   integer :: iconst_f0_particle
!$acc declare create(iconst_f0_particle)
   real :: dt_particle, t0_norm_particle, v0_norm_particle, b0_norm_particle
!$acc declare create(dt_particle, t0_norm_particle, v0_norm_particle,b0_norm_particle)
   complex :: rfac_particle
!$acc declare create(rfac_particle)
   real, dimension(:, :, :), pointer :: mesh_coord !Neighbor tracking arrays
!$acc declare link(mesh_coord)
   integer :: win_mesh_coord
   integer, dimension(:, :), pointer :: mesh_nodes !Neighbor tracking arrays
   integer :: win_mesh_nodes
   integer, dimension(:, :), allocatable :: neighborlist !Neighbor tracking arrays
!$acc declare link(neighborlist)
   integer, dimension(:), allocatable :: localmeshid
   integer :: nparticles, locparts
   integer :: mpi_particle !User-defined MPI datatype for particle communication
   integer :: mpi_elfield
   type(matrix_type) :: diff2_mat, diff3_mat
   type(matrix_type) :: s1_0_mat
   !type(newvar_matrix) :: diff_mat
   integer :: hostcomm, rowcomm
   integer :: hostrank, rowrank, ncols, nrows
   integer, dimension(:), allocatable :: nelm_row, displs_elm
   vectype, dimension(:, :), allocatable :: coeffspaf, coeffspef, coeffspaf_local, coeffspef_local
   vectype, dimension(:, :), allocatable :: coeffspai, coeffspei, coeffspai_local, coeffspei_local
   vectype, dimension(:, :), allocatable :: coeffsdei0, coeffsdei0_local
   vectype, dimension(:, :), allocatable :: coeffsdef0, coeffsdef0_local
   vectype, dimension(:, :), allocatable :: coeffsdef, coeffsdef_local
   vectype, dimension(:, :), allocatable :: coeffsdei, coeffsdei_local
   vectype, dimension(:, :), allocatable :: coeffsvpi, coeffsvpi_local
   vectype, dimension(:, :), allocatable :: coeffsvpf, coeffsvpf_local
   integer :: ielm_min, ielm_max, ipart_begin, ipart_end
   integer :: nelms, nelms_global, nnodes_global
   real :: psi_axis, nf_axis, nfi_axis, toroidal_period_particle
!$acc declare create(psi_axis, nf_axis, nfi_axis, toroidal_period_particle)
   integer :: gyroaverage_particle, fast_ion_dist_particle
!$acc declare create(gyroaverage_particle,fast_ion_dist_particle)
   real :: fast_ion_max_energy_particle, kinetic_rhomax_particle
!$acc declare create(fast_ion_max_energy_particle,kinetic_rhomax_particle)
   integer :: num_energy, num_pitch, num_r
!$acc declare create(num_energy,num_pitch,num_r)
   real, dimension(:), allocatable :: energy_array, pitch_array, r_array
!$acc declare link(energy_array, pitch_array, r_array)
   real, dimension(:, :, :), allocatable :: f_array, f_array2
!$acc declare link(f_array, f_array2)
contains

#ifdef USEPARTICLES

!Define MPI datatype for particle communication
! Note: any changes to the "particle" user-defined datatype must be reflected
!       in the definitions of pnvars, pblklen, ptyps, and pdspls below.
subroutine define_mpi_particle(ierr)
   implicit none

   include 'mpif.h'

   integer, intent(out) :: ierr
   integer, parameter :: pnvars = 15
   integer, dimension(pnvars), parameter :: pblklen = (/3, vspdims, 1, 1, 1, 1, 1, 1, 1, 12, 4, 1, 3, vspdims, 1/)
   integer(kind=MPI_ADDRESS_KIND), dimension(pnvars) :: pdspls
   integer, dimension(pnvars), parameter :: ptyps = (/MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
 MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
 MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, &
 MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_LOGICAL/)

   type(particle) :: dum_par

   !Set up component displacements array
   call mpi_get_address(dum_par%x, pdspls(1), ierr)
   call mpi_get_address(dum_par%v, pdspls(2), ierr)
   pdspls(2) = pdspls(2) - pdspls(1)
   call mpi_get_address(dum_par%wt, pdspls(3), ierr)
   pdspls(3) = pdspls(3) - pdspls(1)
   call mpi_get_address(dum_par%dB, pdspls(4), ierr)
   pdspls(4) = pdspls(4) - pdspls(1)
   call mpi_get_address(dum_par%B0, pdspls(5), ierr)
   pdspls(5) = pdspls(5) - pdspls(1)
   call mpi_get_address(dum_par%dE, pdspls(6), ierr)
   pdspls(6) = pdspls(6) - pdspls(1)
   call mpi_get_address(dum_par%sps, pdspls(7), ierr)
   pdspls(7) = pdspls(7) - pdspls(1)
   call mpi_get_address(dum_par%gid, pdspls(8), ierr)
   pdspls(8) = pdspls(8) - pdspls(1)
   call mpi_get_address(dum_par%jel, pdspls(9), ierr)
   pdspls(9) = pdspls(9) - pdspls(1)
   call mpi_get_address(dum_par%kx, pdspls(10), ierr)
   pdspls(10) = pdspls(10) - pdspls(1)
   call mpi_get_address(dum_par%kel, pdspls(11), ierr)
   pdspls(11) = pdspls(11) - pdspls(1)
   call mpi_get_address(dum_par%f0, pdspls(12), ierr)
   pdspls(12) = pdspls(12) - pdspls(1)
   call mpi_get_address(dum_par%x0, pdspls(13), ierr)
   pdspls(13) = pdspls(13) - pdspls(1)
   call mpi_get_address(dum_par%v0, pdspls(14), ierr)
   pdspls(14) = pdspls(14) - pdspls(1)
   call mpi_get_address(dum_par%deleted, pdspls(15), ierr)
   pdspls(15) = pdspls(15) - pdspls(1)
   pdspls(1) = 0

   call mpi_type_create_struct(pnvars, pblklen, pdspls, ptyps, mpi_particle, ierr)
   if (ierr .ne. 0) return
   call mpi_type_commit(mpi_particle, ierr)
end subroutine define_mpi_particle

subroutine define_mpi_elfield(ierr)
   implicit none

   include 'mpif.h'

   integer, intent(out) :: ierr
   integer, parameter :: pnvars = 24
   integer :: i
   integer, dimension(pnvars), parameter :: pblklen = [(coeffs_per_element, i=1, pnvars)]
   integer(kind=MPI_ADDRESS_KIND), dimension(pnvars) :: pdspls
   integer, dimension(pnvars) :: ptyps

   type(elfield) :: dum_elfield

#ifdef USECOMPLEX
   ptyps(:) = MPI_DOUBLE_complex
#else
   ptyps(:) = MPI_DOUBLE_precision
#endif
   !Set up component displacements array
   call mpi_get_address(dum_elfield%psiv0, pdspls(1), ierr)
   call mpi_get_address(dum_elfield%psiv1, pdspls(2), ierr)
   pdspls(2) = pdspls(2) - pdspls(1)
   call mpi_get_address(dum_elfield%Bzv0, pdspls(3), ierr)
   pdspls(3) = pdspls(3) - pdspls(1)
   call mpi_get_address(dum_elfield%Bzv1, pdspls(4), ierr)
   pdspls(4) = pdspls(4) - pdspls(1)
   call mpi_get_address(dum_elfield%Bfpv0, pdspls(5), ierr)
   pdspls(5) = pdspls(5) - pdspls(1)
   call mpi_get_address(dum_elfield%Bfpv1, pdspls(6), ierr)
   pdspls(6) = pdspls(6) - pdspls(1)
   call mpi_get_address(dum_elfield%er, pdspls(7), ierr)
   pdspls(7) = pdspls(7) - pdspls(1)
   call mpi_get_address(dum_elfield%ephi, pdspls(8), ierr)
   pdspls(8) = pdspls(8) - pdspls(1)
   call mpi_get_address(dum_elfield%ez, pdspls(9), ierr)
   pdspls(9) = pdspls(9) - pdspls(1)
   call mpi_get_address(dum_elfield%rho, pdspls(10), ierr)
   pdspls(10) = pdspls(10) - pdspls(1)
   call mpi_get_address(dum_elfield%nf, pdspls(11), ierr)
   pdspls(11) = pdspls(11) - pdspls(1)
   call mpi_get_address(dum_elfield%tf, pdspls(12), ierr)
   pdspls(12) = pdspls(12) - pdspls(1)
   call mpi_get_address(dum_elfield%nfi, pdspls(13), ierr)
   pdspls(13) = pdspls(13) - pdspls(1)
   call mpi_get_address(dum_elfield%tfi, pdspls(14), ierr)
   pdspls(14) = pdspls(14) - pdspls(1)
   call mpi_get_address(dum_elfield%B0, pdspls(15), ierr)
   pdspls(15) = pdspls(15) - pdspls(1)
   call mpi_get_address(dum_elfield%B1, pdspls(16), ierr)
   pdspls(16) = pdspls(16) - pdspls(1)
   call mpi_get_address(dum_elfield%B1_last, pdspls(17), ierr)
   pdspls(17) = pdspls(17) - pdspls(1)
   call mpi_get_address(dum_elfield%pe, pdspls(18), ierr)
   pdspls(18) = pdspls(18) - pdspls(1)
   call mpi_get_address(dum_elfield%ne, pdspls(19), ierr)
   pdspls(19) = pdspls(19) - pdspls(1)
   call mpi_get_address(dum_elfield%pe0, pdspls(20), ierr)
   pdspls(20) = pdspls(20) - pdspls(1)
   call mpi_get_address(dum_elfield%ne0, pdspls(21), ierr)
   pdspls(21) = pdspls(21) - pdspls(1)
   call mpi_get_address(dum_elfield%te0, pdspls(22), ierr)
   pdspls(22) = pdspls(22) - pdspls(1)
   call mpi_get_address(dum_elfield%rst, pdspls(23), ierr)
   pdspls(23) = pdspls(23) - pdspls(1)
   call mpi_get_address(dum_elfield%zst, pdspls(24), ierr)
   pdspls(24) = pdspls(24) - pdspls(1)
   pdspls(1) = 0

   call mpi_type_create_struct(pnvars, pblklen, pdspls, ptyps, mpi_elfield, ierr)
   if (ierr .ne. 0) return
   call mpi_type_commit(mpi_elfield, ierr)
end subroutine define_mpi_elfield

!---------------------------------------------------------------------------
subroutine particle_test
   use basic
   use diagnostics
   use auxiliary_fields
   use m3dc1_nint
   use arrays
   implicit none
   include 'mpif.h'

   integer, parameter :: trunit = 120
   real :: tstart, tend
   integer :: ierr
   vectype, dimension(coeffs_per_element) :: denfi0
   real, dimension(2) :: nrmfac_temp
   type(particle) :: dpar
   type(xgeomterms) :: geomterms
   integer :: itri, ielm, izone



   if (myrank .eq. 0) then
      print *, 'xlim2 = ', xlim2
      print *, 'nonrect = ', nonrect
      print *, 'iwall_is_limiter = ', iwall_is_limiter
      print *, 'ifixedb = ', ifixedb
      print *, 'GS magnetic axis at ', xmag, ', ', zmag
      print *, 'jadv = ', jadv
      print *, 'imp_hyper = ', imp_hyper
      !print *,'xzero, zzero = ',xzero,', ',zzero
      print *, 'rfac = ', rfac
   end if

   !Precompute electric field components (do it this way for testing only!)
   call calculate_auxiliary_fields(eqsubtract)

   !Initialize particle population
   call second(tstart)
   call init_particles(irestart .gt. 0, ierr)
   if (ierr .ne. 0) return
   if (myrank .eq. 0) then
      call second(tend)
      write (0, '(I12,A,f9.2,A)') nparticles, ' particles initialized in', &
         tend - tstart, ' seconds.'
   end if
   !call particle_step(dt*t0_norm)

   if (irestart.eq.0) then
      nrmfac(:)=1.0
      call update_particle_pressure
      nrmfac_temp(:)=0.
      dpar%x(1) = xmag
      dpar%x(3) = zmag
      dpar%x(2) = 0.
      itri = 0
      call get_geom_terms(dpar%x, itri, geomterms, .false., ierr)
      if (localmeshid(itri)>0) then
         write(0,*) itri
         ielm=localmeshid(itri)
         call get_zone(ielm, izone)
         call define_element_quadrature(ielm, int_pts_main, int_pts_tor)
         call define_fields(ielm, FIELD_KIN, 1, 0)

         if (kinetic_fast_ion.eq.1) then
            call eval_ops(ielm, den_f_0, nf79)
            nrmfac_temp(2)=nrmfac(2)/sum(nf79(:,OP_1))*sum(nf079(:,OP_1))
         endif
         if (kinetic_thermal_ion.eq.1) then
            call eval_ops(ielm, den_i_0, nfi79)
            nrmfac_temp(1)=nrmfac(1)/sum(nfi79(:,OP_1))*sum(nfi079(:,OP_1))
         endif
      endif
      call mpi_allreduce(nrmfac_temp, nrmfac, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   endif
   nrmfac = nrmfac*kinetic_nrmfac_scale
   call update_particle_pressure

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
end subroutine particle_test

!---------------------------------------------------------------------------
! Particle physics parameters
subroutine get_ion_physics_params(speed)
   use basic
   implicit none

   real, intent(out) :: speed
   real, parameter :: c_mks = 2.9979e+8
   real EmaxeV

   EmaxeV = 10000.0                            !Peak ion kinetic energy in eV
   m_ion(1) = ion_mass * m_proton                 !Ion mass in kg
   m_ion(2) = fast_ion_mass * m_proton                 !Ion mass in kg
   !m_ion = (0.111)*m_proton                   !fishbone
   !q_ion = Z_ion * e_mks                       !Ion charge in C
   q_ion(1) = z_ion*e_mks                       !fishbone
   q_ion(2) = fast_ion_z*e_mks                       !fishbone
   qm_ion = q_ion/m_ion                      !Ion charge/mass ratio
   speed = sqrt(EmaxeV/ion_mass)*vp1eV     !Peak ion speed in m/s
   speed = (v0_norm/100.0)*4.0
   !speed = v0_norm/100.*2.58
   !if (myrank.eq.0) print *,'peak ion speed = ',speed,' m/s = ',speed/c_mks,' c.'
end subroutine get_ion_physics_params
!---------------------------------------------------------------------------
real function tri_area(x1, z1, x2, z2, x3, z3)
!$acc routine seq
   implicit none

   real :: x1, z1, x2, z2, x3, z3

   tri_area = 0.5*(x1*(z2 - z3) + x2*(z3 - z1) + x3*(z1 - z2))
end function tri_area
!---------------------------------------------------------------------------
subroutine init_particles(lrestart, ierr)
   use basic
   use arrays
   use hdf5_output
   use m3dc1_nint
   use gradshafranov
   use read_ascii
   use model
   use boundary_conditions

   implicit none
   include 'mpif.h'

   logical, intent(in) :: lrestart
   integer, intent(out) :: ierr

   type(particle) :: dpar, dparold  !Dummy particle
   type(elfield), dimension(nneighbors + 1) :: elcoefs
   type(xgeomterms) :: geomterms
   real, dimension(3) :: Bcyl, deltaB
   integer, dimension(nodes_per_element) :: nodeids
#ifndef USE3D
   real, dimension(2) :: mmsa
#endif
   real    :: x1, x2, x3, z1, z2, z3, phi1, phi2, phi3, pdx, pdphi, pdz, pdl, vpar, vperp
   real    :: x_min, x_max, phi_min, phi_max, z_min, z_max, area
   real    :: maxspeed, lambda_min, lambda_max, B0, Bmax, pphi, lambda
   real    :: gyroperiod, gfrac = 5.0e-3, gkfrac = 800.0, dtp, ldtmin
   real    :: par_r, par_theta, x_temp, z_temp, ran_temp, ran_temp2, y1, y2
   integer :: npr, npphi, npz, npe, npmu, nphiv, ir, iphi, iz, ip, ipar, nptheta
   integer(kind=8) :: npar_local
   integer :: ielm, ielm_global, inode, inode_global, lc, noc, tridex, itri, ielmold = 0
   integer :: gfx, isghost, xid, nle = 0, nge = 0
   integer :: ie, imu, pperrow
   character(len=32) :: part_file_name
   type(particle), dimension(:), allocatable :: pdata_local  !Particle arrays
   type(element_data) :: eldat
   real :: xi, zi, eta, xi2, eta2
   vectype, dimension(dofs_per_element, dofs_per_element) :: tempxx
   integer :: blocky, starty, endy
   TYPE(C_PTR) :: baseptr
   INTEGER(KIND=MPI_ADDRESS_KIND) :: arraysize
   TYPE(elfield) :: elfieldtest
   integer :: disp_unit
   integer :: icol, irow
   integer, dimension(:), allocatable :: recvcounts, displs
   integer :: sendcount
   character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
   integer :: hostname_len
   integer :: info
   integer :: child_pid
   integer(kind=MPI_ADDRESS_KIND) :: pppp
   real, dimension(3) :: gradpsi, gradrho, gradf, gradT
   real :: gradcoef
   real :: rho, f0, fx, T00, f_mesh, df0de
   real, allocatable, dimension(:) :: rho_vals, psi_vals, nf_vals, tf_vals
   integer :: ipoint, dotloc

   real, dimension(3) :: bhat, Bstar, Jcyl, svec, dBdR, dBdphi, dBdz, gradB0, gradB1, BxgrdB, dvdt
   real :: gradcoef0, df0de0, df0dxi0
   real :: Rinv, Bss, dB1, di, R_axis, di_axis
   integer :: num_f
   real, dimension(:), allocatable :: f_array_temp
   real :: radi, pitch, energy
   integer :: radi_i, pitch_i, energy_i
   !real, dimension(:, :, :), allocatable :: f_array2
   real :: f1, f2, f3, f4, f5, f6
   integer :: ran_size
   integer, allocatable :: seed(:)
   integer, dimension(dofs_per_element) :: imask
   type(vector_type) :: b1_vel
   real, dimension(:, :, :), allocatable :: mesh_coord_temp
         type(field_type) ::   u_v
   type(field_type) ::  vz_v
   type(field_type) :: chi_v
  type(vector_type) :: vel_vec

   integer, dimension(:, :), allocatable :: mesh_nodes_temp
   integer :: nelm_row_temp
   integer :: ielm_min_temp, ielm_max_temp
   integer :: sps
   real :: npar_ratio, npar_ratio_temp, npar_ratio_local, npar_fac

   !Allocate particle pressure tensor components

   call get_ion_physics_params(maxspeed)

   !nf_field=te_field(0)
   !tf_field=te_field(0)
   !call mult(tf_field,496.69)

   call m3dc1_mesh_getnumglobalent(0, nnodes_global)
#ifdef USE3D
   call m3dc1_mesh_getnumglobalent(3, nelms_global)
#else
   call m3dc1_mesh_getnumglobalent(2, nelms_global)
#endif
   call define_mpi_particle(ierr)
   call define_mpi_elfield(ierr)

   CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                            MPI_INFO_NULL, hostcomm, ierr)
   CALL MPI_Comm_rank(hostcomm, hostrank, ierr)
   CALL MPI_Comm_size(hostcomm, ncols, ierr)
   call MPI_Comm_split(MPI_COMM_WORLD, hostrank, myrank, rowcomm, ierr)
   CALL MPI_Comm_size(rowcomm, nrows, ierr)
   CALL MPI_Comm_rank(rowcomm, rowrank, ierr)
   disp_unit = 1
   arraysize = 0
   if (hostrank == 0) arraysize = nelms_global*nodes_per_element*3*sizeof(1.d0)
   CALL MPI_Win_allocate_shared(arraysize, disp_unit, MPI_INFO_NULL, &
                                hostcomm, baseptr, win_mesh_coord, ierr)
   if (hostrank /= 0) CALL MPI_Win_shared_query(win_mesh_coord, 0, arraysize, disp_unit, baseptr, ierr)
   !CALL C_F_POINTER(baseptr, mesh_coord_shared, [2,3,nelms])
   CALL C_F_POINTER(baseptr, mesh_coord, [3, nodes_per_element, nelms_global])
   !allocate(mesh_coord(2,3,nelms))

   disp_unit = 1
   arraysize = 0
   if (hostrank == 0) arraysize = nelms_global*nodes_per_element*sizeof(1)
   CALL MPI_Win_allocate_shared(arraysize, disp_unit, MPI_INFO_NULL, &
                                hostcomm, baseptr, win_mesh_nodes, ierr)
   if (hostrank /= 0) CALL MPI_Win_shared_query(win_mesh_nodes, 0, arraysize, disp_unit, baseptr, ierr)
   CALL C_F_POINTER(baseptr, mesh_nodes, [nodes_per_element, nelms_global])
   !allocate(mesh_nodes(3,nelms))

   allocate (localmeshid(nelms_global))
   localmeshid = 0
   nelms = local_elements()
   ielm_min = 1e9
   ielm_max = -1
   do ielm = 1, nelms
#ifdef USE3D
      call m3dc1_ent_getglobalid(3, ielm - 1, ielm_global)
#else
      call m3dc1_ent_getglobalid(2, ielm - 1, ielm_global)
#endif
      ielm_global = ielm_global + 1
      if (ielm_min > ielm_global) ielm_min = ielm_global
      if (ielm_max < ielm_global) ielm_max = ielm_global
      localmeshid(ielm_global) = ielm
      call get_element_nodes(ielm, nodeids)
      do inode = 1, nodes_per_element
         call m3dc1_ent_getglobalid(0, nodeids(inode) - 1, inode_global)
         mesh_nodes(inode, ielm_global) = inode_global + 1
         call get_node_pos(nodeids(inode), mesh_coord(1,inode,ielm_global), mesh_coord(2,inode,ielm_global),&
            mesh_coord(3,inode,ielm_global))
      end do
   end do
   call mpi_barrier(mpi_comm_world, ierr)
   ielm_min_temp=ielm_min
   ielm_max_temp=ielm_max
   call mpi_allreduce(ielm_min_temp, ielm_min, 1, MPI_INTEGER, MPI_MIN, hostcomm, ierr)
   call mpi_allreduce(ielm_max_temp, ielm_max, 1, MPI_INTEGER, MPI_MAX, hostcomm, ierr)
   allocate (nelm_row(nrows), displs_elm(nrows))
   nelm_row(rowrank + 1) = ielm_max - ielm_min + 1
   nelm_row_temp=nelm_row(rowrank + 1)
   call mpi_allgather(nelm_row_temp, 1, MPI_INTEGER, nelm_row, 1, MPI_INTEGER, rowcomm, ierr)
   displs_elm(1) = 0
   do irow = 2, nrows
      displs_elm(irow) = displs_elm(irow - 1) + nelm_row(irow - 1)
   end do
    if (hostrank == 1) then
      allocate (recvcounts(nrows))
      allocate (displs(nrows))
      sendcount = nodes_per_element*nelm_row(rowrank + 1)
      recvcounts = nodes_per_element*nelm_row
      displs = nodes_per_element*displs_elm
      allocate (mesh_nodes_temp(nodes_per_element, nelm_row(rowrank+1)))
      mesh_nodes_temp=mesh_nodes(:,ielm_min:ielm_max)
      call MPI_ALLGATHERV(mesh_nodes_temp,sendcount,MPI_INTEGER, mesh_nodes,recvcounts,displs,&
         MPI_INTEGER, rowcomm,ierr)
      deallocate (mesh_nodes_temp)
      sendcount = nodes_per_element*3*nelm_row(rowrank + 1)
      recvcounts = nodes_per_element*3*nelm_row
      displs = nodes_per_element*3*displs_elm
      allocate (mesh_coord_temp(3, nodes_per_element, nelm_row(rowrank+1)))
      mesh_coord_temp=mesh_coord(:,:,ielm_min:ielm_max)
      call MPI_ALLGATHERV(mesh_coord_temp,sendcount,MPI_DOUBLE_PRECISION, mesh_coord,recvcounts,&
         displs,MPI_DOUBLE_PRECISION, rowcomm,ierr) 
      deallocate (mesh_coord_temp)
      deallocate (recvcounts)
      deallocate (displs)
   end if
   call mpi_barrier(mpi_comm_world, ierr)
   disp_unit = 1
   arraysize = 0
   if (hostrank == 0) arraysize = nelms_global*sizeof(elfieldtest)
   CALL MPI_Win_allocate_shared(arraysize, disp_unit, MPI_INFO_NULL, &
                                hostcomm, baseptr, win_elfieldcoefs, ierr)
   if (hostrank /= 0) CALL MPI_Win_shared_query(win_elfieldcoefs, 0, arraysize, disp_unit, baseptr, ierr)
   CALL C_F_POINTER(baseptr, elfieldcoefs, [nelms_global])
   !allocate(elfieldcoefs(nelms))
   call get_field_coefs(1)
   !call MPI_Win_fence(0, win_elfieldcoefs)

   !npar = nplanes
   !if (npar < 4) npar = 4
   !npar = npar*2e6

   !Allocate local storage for particle data
   disp_unit = 1
   arraysize = 0
   if (hostrank == 0) then
      arraysize = num_par_max
      arraysize = arraysize*sizeof(dpar)
   end if
   CALL MPI_Win_allocate_shared(arraysize, disp_unit, MPI_INFO_NULL, &
                                hostcomm, baseptr, win_pdata, ierr)
   if (hostrank /= 0) CALL MPI_Win_shared_query(win_pdata, 0, arraysize, disp_unit, baseptr, ierr)
   CALL C_F_POINTER(baseptr, pdata, [num_par_max])

   !Normalize integrals
   !nrmfac(1) = 3207919.*4.*0.0876/twopi/npar*4.e6*143 !nstx
!#ifdef USEST
   !nrmfac(2) = 3207919.*4.*0.1349267*36.36*0.00025/npar*4.e6
!#else
   !nrmfac(2) = 3207919.*4.*0.1349267*36.36*0.0297/npar*4.e6
!#endif

   !kinetic_nrmfac = kinetic_nrmfac/num_par_fac
!#ifdef USE3D
   !kinetic_nrmfac = kinetic_nrmfac*toroidal_period
!#endif

   !Set up 'neighborlist' table of element neighbors for ensemble tracking
   call find_element_neighbors

   dt_particle = dt
   t0_norm_particle = t0_norm
   v0_norm_particle = v0_norm
   b0_norm_particle = b0_norm
   rfac_particle = rfac
   particle_linear_particle = particle_linear
   toroidal_period_particle = toroidal_period
   gyroaverage_particle = igyroaverage
   iconst_f0_particle = iconst_f0
   psubsteps_particle = particle_substeps
   kinetic_thermal_ion_particle=kinetic_thermal_ion
   fast_ion_dist_particle = fast_ion_dist
   fast_ion_max_energy_particle = fast_ion_max_energy
   kinetic_rhomax_particle = kinetic_rhomax
   dpar%x(1) = xmag
   dpar%x(3) = zmag
   dpar%x(2) = 0.
   itri = 0
   call get_geom_terms(dpar%x, itri, geomterms, .false., ierr)
   psi_axis = dot_product(elfieldcoefs(itri)%psiv0, geomterms%g)
   nf_axis = dot_product(elfieldcoefs(itri)%nf, geomterms%g)
   nfi_axis = dot_product(elfieldcoefs(itri)%nfi, geomterms%g)
#ifdef USEST
   R_axis=dot_product(elfieldcoefs(itri)%rst, geomterms%g)
   di_axis = 1./(dot_product(elfieldcoefs(itri)%rst,geomterms%dr)*dot_product(elfieldcoefs(itri)%zst,geomterms%dz)&
      - dot_product(elfieldcoefs(itri)%zst,geomterms%dr)*dot_product(elfieldcoefs(itri)%rst,geomterms%dz))
#endif

   if ((kinetic_fast_ion.eq.1).and.(fast_ion_dist.eq.0)) then
      num_energy=0
      call read_ascii_column('ep_energy_array', energy_array, num_energy, icol=1)
        !energy_array=energy_array*1.7
      num_pitch=0
      call read_ascii_column('ep_pitch_array', pitch_array, num_pitch, icol=1)
      num_r=0
      call read_ascii_column('ep_r_array', r_array, num_r, icol=1)
      num_f=0
      call read_ascii_column('ep_f_array', f_array_temp, num_f, icol=1)
      f_array_temp=f_array_temp/maxval(f_array_temp)
      allocate(f_array(num_energy,num_pitch,num_r))
      f_array=reshape(f_array_temp,[num_energy,num_pitch,num_r])
      deallocate(f_array_temp)
      call read_ascii_column('ep_f_array2', f_array_temp, num_f, icol=1)
      f_array_temp=f_array_temp/maxval(f_array_temp)
      allocate(f_array2(num_energy,num_pitch,num_r))
      f_array2=reshape(f_array_temp,[num_energy,num_pitch,num_r])
       !allocate(f_array2(num_energy,num_pitch,num_r))
      !do pitch_i=1,num_pitch
      !   f_array2(:,pitch_i,:)=f_array(:,1+num_pitch-pitch_i,:)
      !enddo
      !f_array=f_array2
   endif

   if (lrestart) then
      write (part_file_name, '("ions_",I4.4,".h5")') times_output
      write (0, *) part_file_name
      call hdf5_read_particles(part_file_name, ierr)
      if (myrank .eq. 0) print *, 'read_parts returned with ierr=', ierr
   else
      p_f_par%vec = 0.; p_f_perp%vec = 0.
      p_i_par%vec = 0.; p_i_perp%vec = 0.
      den_i_0%vec = 0.
      den_f_0%vec = 0.
      allocate (pdata_local(num_par_max/maxrank*10))
      !First pass: assign particles to processors, elements
      locparts = 0
      !allocate(ran(npar*5))
      !call random_number(ran)
      !call random_seed(size=ran_size)
      !allocate (seed(ran_size))
      !call random_seed(get=seed)
      !seed = seed + myrank
      !call random_seed(put=seed)

      do sps = 1,2
         if ((sps==1).and.(kinetic_thermal_ion==0)) cycle
         if ((sps==2).and.(kinetic_fast_ion==0)) cycle
         npar_ratio=0.

      do ielm = 1, nelms
#ifdef USE3D
         call m3dc1_ent_getglobalid(3, ielm - 1, itri)
#else
         call m3dc1_ent_getglobalid(2, ielm - 1, itri)
#endif
         itri = itri + 1
         x_min = minval(mesh_coord(1, :, itri))
         x_max = maxval(mesh_coord(1, :, itri))
         phi_min = mesh_coord(2, 1, itri)
         phi_max = mesh_coord(2, 4, itri)
         if (phi_max == 0) phi_max = toroidal_period
         z_min = minval(mesh_coord(3, :, itri))
         z_max = maxval(mesh_coord(3, :, itri))

         area=tri_area(mesh_coord(1,1,itri),mesh_coord(3,1,itri),mesh_coord(1,2,itri),&
            mesh_coord(3,2,itri),mesh_coord(1,3,itri),mesh_coord(3,3,itri))
#ifdef USE3D
         dpar%x = (mesh_coord(:, 1, itri) + mesh_coord(:, 2, itri) + mesh_coord(:, 3, itri) + &
            & mesh_coord(:, 4, itri) + mesh_coord(:, 5, itri) + mesh_coord(:, 6, itri))/6.
#else
         dpar%x = (mesh_coord(:, 1, itri) + mesh_coord(:, 2, itri) + mesh_coord(:, 3, itri))/3.
#endif
         call get_geom_terms(dpar%x, itri, geomterms, .false., ierr)
#ifdef USEST
         di = 1./(dot_product(elfieldcoefs(itri)%rst,geomterms%dr)*dot_product(elfieldcoefs(itri)%zst,geomterms%dz)&
            - dot_product(elfieldcoefs(itri)%zst,geomterms%dr)*dot_product(elfieldcoefs(itri)%rst,geomterms%dz))
         call update_geom_terms_st(geomterms, elfieldcoefs(itri), .false.)
#endif
         if (sps==1) then
            f_mesh = dot_product(elfieldcoefs(itri)%nfi,geomterms%g)/nfi_axis
         elseif (fast_ion_dist.eq.0) then
            f_mesh = 1.
         else
            f_mesh = dot_product(elfieldcoefs(itri)%nf,geomterms%g)/nf_axis
         endif
#ifdef USEST
         !npar_ratio_local=(x_max-x_min)*(dot_product(elfieldcoefs(itri)%rst,geomterms%g)/R_axis)/(di/di_axis)*(z_max-z_min)*f_mesh*2/nplanes
         npar_ratio_local=area*(dot_product(elfieldcoefs(itri)%rst,geomterms%g)/R_axis)/(di/di_axis)*f_mesh*2/nplanes
#else
         npar_ratio_local = area*(x_max + x_min)/2.*f_mesh*2/nplanes
#endif
         npar_ratio=npar_ratio+npar_ratio_local
      end do
      npar_ratio_temp=npar_ratio
      call mpi_allreduce(npar_ratio_temp, npar_ratio, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      npar_fac=0.95/npar_ratio*num_par_scale(sps)
      if ((kinetic_fast_ion.eq.1).and.(kinetic_thermal_ion.eq.1)) then
         npar_fac=npar_fac*0.5
      endif
 
      do ielm = 1, nelms
#ifdef USE3D
         call m3dc1_ent_getglobalid(3, ielm - 1, itri)
#else
         call m3dc1_ent_getglobalid(2, ielm - 1, itri)
#endif
         itri = itri + 1
         x_min = minval(mesh_coord(1, :, itri))
         x_max = maxval(mesh_coord(1, :, itri))
         phi_min = mesh_coord(2, 1, itri)
         phi_max = mesh_coord(2, 4, itri)
         if (phi_max == 0) phi_max = toroidal_period
         z_min = minval(mesh_coord(3, :, itri))
         z_max = maxval(mesh_coord(3, :, itri))

         !area=x1*(z2-z3)+x2*(z3-z1)+x3*(z1-z2)
         area=tri_area(mesh_coord(1,1,itri),mesh_coord(3,1,itri),mesh_coord(1,2,itri),&
            mesh_coord(3,2,itri),mesh_coord(1,3,itri),mesh_coord(3,3,itri))
         !hou!npar_local=int(npar*area*mesh_coord(1,1,itri)*0.03/nplanes)
#ifdef USE3D
         dpar%x = (mesh_coord(:, 1, itri) + mesh_coord(:, 2, itri) + mesh_coord(:, 3, itri) + &
            & mesh_coord(:, 4, itri) + mesh_coord(:, 5, itri) + mesh_coord(:, 6, itri))/6.
#else
         dpar%x = (mesh_coord(:, 1, itri) + mesh_coord(:, 2, itri) + mesh_coord(:, 3, itri))/3.
#endif
         dpar%v = 0.
         B0 = 1.
         call get_geom_terms(dpar%x, itri, geomterms, .false., ierr)
#ifdef USEST
         di = 1./(dot_product(elfieldcoefs(itri)%rst,geomterms%dr)*dot_product(elfieldcoefs(itri)%zst,geomterms%dz)&
            - dot_product(elfieldcoefs(itri)%zst,geomterms%dr)*dot_product(elfieldcoefs(itri)%rst,geomterms%dz))
         call update_geom_terms_st(geomterms, elfieldcoefs(itri), .false.)
#endif
         if (sps==1) then
            f_mesh = dot_product(elfieldcoefs(itri)%nfi,geomterms%g)/nfi_axis
         elseif (fast_ion_dist.eq.0) then
            f_mesh = 1.
         else
            f_mesh = dot_product(elfieldcoefs(itri)%nf,geomterms%g)/nf_axis
         endif
#ifdef USEST
         npar_local=int(num_par_max*area*(dot_product(elfieldcoefs(itri)%rst,geomterms%g)/R_axis)/(di/di_axis)&
            *f_mesh*2/nplanes*npar_fac)!fullf
#else
         npar_local = int(num_par_max*area*(x_max + x_min)/2.*f_mesh*2/nplanes*npar_fac)!fullf
#endif
         ipar = 1
         !do while (ipar<=npar_local)
         do ipar=1,npar_local
            call random_number(ran_temp)
#ifdef USEST
            x_temp = ran_temp*(x_max-x_min)+x_min
#else
            x_temp = sqrt(ran_temp*(x_max**2 - x_min**2) + x_min**2)
#endif
            !dpar%x(1) = x_temp*(x2-x1)+x1
            call random_number(ran_temp)
            z_temp = z_min + ran_temp*(z_max - z_min)
            if ((tri_area(mesh_coord(1, 1, itri), mesh_coord(3, 1, itri), mesh_coord(1, 2, itri),&
               mesh_coord(3, 2, itri), x_temp, z_temp) < 0) &
            .or.(tri_area(x_temp,z_temp,mesh_coord(1,2,itri),&
            mesh_coord(3,2,itri),mesh_coord(1,3,itri),mesh_coord(3,3,itri))<0) &
            .or.(tri_area(mesh_coord(1,1,itri),mesh_coord(3,1,itri), x_temp,&
            z_temp,mesh_coord(1,3,itri),mesh_coord(3,3,itri))<0)) cycle
            dpar%x(1) = x_temp
            dpar%x(3) = z_temp
            call random_number(ran_temp)
#ifdef USE3D
            dpar%x(2) = phi_min + (phi_max - phi_min)*ran_temp
#else
            dpar%x(2) = twopi*ran_temp
#endif
            call get_geom_terms(dpar%x, itri, geomterms, vspdims .eq. 5, ierr)
#ifdef USEST
            di = 1./(dot_product(elfieldcoefs(itri)%rst,geomterms%dr)*dot_product(elfieldcoefs(itri)%zst,geomterms%dz)&
               - dot_product(elfieldcoefs(itri)%zst,geomterms%dr)*dot_product(elfieldcoefs(itri)%rst,geomterms%dz))
            call update_geom_terms_st(geomterms, elfieldcoefs(itri), vspdims.eq.5)
#endif
            if (ierr .ne. 0) then
               print *, myrank, ': get_geom_terms call failed for particle', ip, &
                  ' of element', itri
               cycle
            end if !ierr
            if (sps==1) then
               T00 = dot_product(elfieldcoefs(itri)%tfi,geomterms%g)
            else
               T00 = dot_product(elfieldcoefs(itri)%tf,geomterms%g)
            endif
            dpar%jel = itri
            dpar%sps = sps
              !Rinv = 1.0/dpar%x(1)
               !Rinv = 1.0/dot_product(elfieldcoefs(itri)%rst,geomterms%g)
               !dpar%df0dpsi=gradcoef
               !!monoenergy
               !y2 = acos(ran_temp*2-1)
               !vpar = maxspeed * cos(y2)
               !vperp = maxspeed * sin(y2)
               !Pphi = dot_product(elfieldcoefs(itri)%psiv0, geomterms%g) + (1./qm_ion) * vpar *  dot_product(elfieldcoefs(itri)%Bzv0, geomterms%g)/B0
               !dpar%f0=Pphi
               !dpar%f0=maxspeed**2
            if ((sps.eq.1).or.(fast_ion_dist.eq.1)) then
               !for Maxwellian
               call random_number(ran_temp)
               y1 = sqrt( - 2*log(ran_temp) )
               vperp=y1*sqrt(T00*1.6e-19/m_ion(sps))
               call random_number(ran_temp)
               call random_number(ran_temp2)
               y1 = sqrt( - 2*log(ran_temp) )* cos( twopi*ran_temp2)
               vpar=y1*sqrt(T00*1.6e-19/m_ion(sps))
               !dpar%f0=dpar%f0*T00**(-1.5)*exp(-m_ion*(vpar**2+vperp**2)/(2*T00*1.6e-19))
            elseif (fast_ion_dist.eq.2) then
               !for slowingdown
               call random_number(ran_temp)
               y1 = (exp(ran_temp*log((fast_ion_max_energy/T00)**(1.5)+1)) - 1.)**(1./3.)
               call random_number(ran_temp)
               y2 = acos(ran_temp*2 - 1)
               vpar = y1*sqrt(2*T00*1.6e-19/m_ion(sps))*cos(y2)
               vperp = y1*sqrt(2*T00*1.6e-19/m_ion(sps))*sin(y2)
            elseif (fast_ion_dist.eq.0) then
               radi = dot_product(elfieldcoefs(itri)%rho, geomterms%g)
               !radi_i=int((radi-r_array(1))/(r_array(2)-r_array(1)))+1
               call random_number(ran_temp)
               pitch = ran_temp*2-1
               call random_number(ran_temp)
               !energy = (ran_temp*energy_array(num_energy)**1.5)**(1/1.5)
               energy = (ran_temp*100000**1.5)**(1/1.5)
               y1 = sqrt(energy*2*1.6e-19/m_ion(sps))
               vpar = y1 * pitch
               vperp = y1 * sqrt(1-pitch**2)

               call getBcyl(dpar%x, elfieldcoefs(itri), geomterms, Bcyl, deltaB, gradB0, gradB1, dB1)
               B0 = sqrt(dot_product(Bcyl, Bcyl))
               pphi = (dot_product(elfieldcoefs(itri)%psiv0, geomterms%g)+(1./qm_ion(sps)) * vpar *  dot_product(elfieldcoefs(itri)%Bzv0, geomterms%g)/B0)-psimin
               radi_i=int((pphi-r_array(1))/(r_array(2)-r_array(1)))+1
               lambda = vperp**2/B0*1.99/y1**2
               pitch_i=int((lambda-pitch_array(1))/(pitch_array(2)-pitch_array(1)))+1
               energy_i=int((energy-energy_array(1))/(energy_array(2)-energy_array(1)))+1

               if (vpar<0) then
               f1=f_array(energy_i,pitch_i,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
               f1=f1+f_array(energy_i+1,pitch_i,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
               f2=f_array(energy_i,pitch_i+1,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
               f2=f2+f_array(energy_i+1,pitch_i+1,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
               f3=f1*(pitch_array(pitch_i+1)-lambda)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
               f3=f3+f2*(lambda-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
               f4=f_array(energy_i,pitch_i,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
               f4=f4+f_array(energy_i+1,pitch_i,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
               f5=f_array(energy_i,pitch_i+1,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
               f5=f5+f_array(energy_i+1,pitch_i+1,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
               f6=f4*(pitch_array(pitch_i+1)-lambda)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
               f6=f6+f5*(lambda-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
               f0=f3*(r_array(radi_i+1)-pphi)/(r_array(radi_i+1)-r_array(radi_i))
               f0=f0+f6*(pphi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
!#ifdef USEST
!               f0=f0*dot_product(geomterms%g,elfieldcoefs(itri)%rst)/di
!               f0=f0/(f_mesh*R_axis/di_axis*2)
!#else
!               f0 = f0/(f_mesh*2)
!#endif
!               !if (f0<0) write(0,*) f0
               else
               f1=f_array2(energy_i,pitch_i,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
               f1=f1+f_array2(energy_i+1,pitch_i,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
               f2=f_array2(energy_i,pitch_i+1,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
               f2=f2+f_array2(energy_i+1,pitch_i+1,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
               f3=f1*(pitch_array(pitch_i+1)-lambda)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
               f3=f3+f2*(lambda-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
               f4=f_array2(energy_i,pitch_i,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
               f4=f4+f_array2(energy_i+1,pitch_i,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
               f5=f_array2(energy_i,pitch_i+1,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
               f5=f5+f_array2(energy_i+1,pitch_i+1,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
               f6=f4*(pitch_array(pitch_i+1)-lambda)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
               f6=f6+f5*(lambda-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
               f0=f3*(r_array(radi_i+1)-pphi)/(r_array(radi_i+1)-r_array(radi_i))
               f0=f0+f6*(pphi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
!#ifdef USEST
!               f0=f0*dot_product(geomterms%g,elfieldcoefs(itri)%rst)/di
!               f0=f0/(f_mesh*R_axis/di_axis*2)
!#else
!               f0 = f0/(f_mesh*2)
!#endif
!               !if (f0<0) write(0,*) f0
               endif
               call random_number(ran_temp)
               !if (f0>0.8) write(0,*) pphi,lambda,energy
               if (ran_temp > f0) cycle
               y1 = sqrt(energy*2*1.6e-19/m_ion(sps))
               vpar = y1 * pitch
               vperp = y1 * sqrt(1-pitch**2)
            endif
            !call evalf0(dpar%x, vpar, vperp, elfieldcoefs(itri), geomterms, sps, f0, gradcoef0, df0de0, df0dxi0)
            f0 = 1.0
            dpar%f0=f0
            if (vspdims .eq. 2) then
               call getBcyl(dpar%x, elfieldcoefs(itri), geomterms, Bcyl, deltaB, gradB0, gradB1, dB1)
            else
               !call getBcylprime(dpar%x, elfieldcoefs(itri), geomterms, Bcyl, deltaB, dBdR, dBdphi, dBdz, .false.)
            end if
            !write(0,*) Bcyl
            B0 = sqrt(dot_product(Bcyl, Bcyl))
            bhat = Bcyl/B0                         !Unit vector in b direction
            !dpar%f0=dpar%f0*1./((vpar**2+vperp**2)**(1.5)+T00**3)

            !if (vspdims.eq.5) then !full orbit
            !dpar%v(4) = vpar
            !dpar%v(5) = (0.5*vperp**2)/(qm_ion*B0)
            !gradB0(1) = dot_product(bhat, dBdR)
            !gradB0(2) = Rinv*dot_product(bhat, dBdphi)
            !gradB0(3) = dot_product(bhat, dBdz)

            !BxgrdB(1) = Bcyl(2)*gradB0(3) - Bcyl(3)*gradB0(2)
            !BxgrdB(2) = Bcyl(3)*gradB0(1) - Bcyl(1)*gradB0(3)
            !BxgrdB(3) = Bcyl(1)*gradB0(2) - Bcyl(2)*gradB0(1)

            !Jcyl(1) = Rinv*dBdphi(3) - dBdz(2)
            !Jcyl(2) = dBdz(1) - dBdR(3)
            !Jcyl(3) = dBdR(2) - Rinv*dBdphi(1)
    !!if (itor.eq.1) Jcyl(3) = Jcyl(3) + Rinv*B_cyl(2)
            !Jcyl(3) = Jcyl(3) + Rinv*Bcyl(2)

            !svec = (Jcyl + BxgrdB/B0)/B0
            !Bstar = Bcyl + (vpar/qm_ion)*svec
            !Bss = dot_product(Bstar, bhat)
            !svec = dpar%v(5)*gradB0 ! - g_mks/qm_ion
            !dpar%v(1) = (vpar*Bstar(1) + bhat(2)*svec(3) - bhat(3)*svec(2))/Bss
            !dpar%v(2) = (vpar*Bstar(2) + bhat(3)*svec(1) - bhat(1)*svec(3))/Bss
            !dpar%v(3) = (vpar*Bstar(3) + bhat(1)*svec(2) - bhat(2)*svec(1))/Bss
            !call advancex(dpar%x, dpar%v, dt*t0_norm/psubsteps/2.)
            !else !gyro- or drift kinetic
            dpar%v(1) = vpar                        !v_parallel
            !pdata(ielm)%ion(ip)%v(1) = 100000.                        !v_parallel
            dpar%v(2) = (0.5*vperp**2)/(qm_ion(sps)*B0)  !mu/q
            !pdata(ielm)%ion(ip)%wt = 1.
            dpar%x0=dpar%x
            dpar%v0=dpar%v
            dpar%dB = 0.
            dpar%dE = 0.
            !endif !vspdims

            dpar%wt = 0.
            dpar%B0 = B0
            dpar%kel(:) = itri
            do ipoint = 1, 4
               dpar%kx(:, ipoint) = dpar%x
            end do
            dpar%deleted = .false.
            locparts = locparts + 1
            dpar%gid = itri*10000 + ipar
            pdata_local(locparts) = dpar
            !ipar = ipar + 1
         end do !iz
      end do
      end do

      allocate (recvcounts(ncols))
      allocate (displs(ncols))
      sendcount = locparts
      recvcounts(hostrank + 1) = sendcount
      call MPI_ALLGATHER(sendcount, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, hostcomm, ierr)
      displs(1) = 0
      do icol = 2, ncols
         displs(icol) = displs(icol - 1) + recvcounts(icol - 1)
      end do
      call MPI_GATHERV(pdata_local(1:locparts), sendcount, MPI_particle, pdata,&
         recvcounts, displs, MPI_particle, 0, hostcomm, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
      ipart_begin = 1
      ipart_end = sum(recvcounts)
      deallocate (recvcounts)
      deallocate (displs)
      deallocate (pdata_local)

   end if !!lrestart

   !if ((kinetic_fast_ion.eq.1).and.(fast_ion_dist.eq.0)) then
   !   if (energy_array(1).eq.0) energy_array=energy_array+0.1*(energy_array(2)-energy_array(1))
   !   do energy_i=1,num_energy
   !      f_array(energy_i,:,:)=f_array(energy_i,:,:)/sqrt(energy_array(energy_i)/energy_array(num_energy))
   !   enddo
   !endif

   if (hostrank == 0) then
      call delete_particle(.true.)
   end if
   call MPI_Bcast(ipart_begin, 1, MPI_INTEGER, 0, hostcomm, ierr)
   call MPI_Bcast(ipart_end, 1, MPI_INTEGER, 0, hostcomm, ierr)
   call MPI_Bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
#ifdef _OPENACC
   if (hostrank < num_devices) then
      write (0, *) 'num_device', num_devices
!$acc enter data create(mesh_coord,neighborlist)
!!$acc enter data create(pdata(starty:endy)) async(blocky)
!$acc update device(mesh_coord,neighborlist)
!!$acc update device(pdata(starty:endy)) async(blocky)
!$acc update device(m_ion,q_ion,qm_ion)
!$acc update device(dt_particle,t0_norm_particle,v0_norm_particle,b0_norm_particle,rfac_particle)
!$acc update device(particle_linear_particle,psi_axis,nf_axis,nfi_axis,toroidal_period_particle)
!$acc update device(gyroaverage_particle,psubsteps_particle,iconst_f0_particle,kinetic_rhomax_particle)
!$acc update device(kinetic_thermal_ion_particle,fast_ion_dist_particle,fast_ion_max_energy_particle)
      if ((kinetic_fast_ion.eq.1).and.(fast_ion_dist.eq.0)) then
!$acc update device(num_energy,num_pitch,num_r) async(blocky)
!$acc enter data create(energy_array,pitch_array,r_array,f_array) async(blocky)
!$acc update device(energy_array,pitch_array,r_array,f_array) async(blocky)
      endif
   end if !hostrank
#endif

   allocate (coeffspaf(coeffs_per_element, nelms_global))
   allocate (coeffspef(coeffs_per_element, nelms_global))
   allocate (coeffspaf_local(coeffs_per_element, nelms_global))
   allocate (coeffspef_local(coeffs_per_element, nelms_global))
   allocate (coeffspai(coeffs_per_element, nelms_global))
   allocate (coeffspei(coeffs_per_element, nelms_global))
   allocate (coeffspai_local(coeffs_per_element, nelms_global))
   allocate (coeffspei_local(coeffs_per_element, nelms_global))
   allocate (coeffsdef0(coeffs_per_element, nelms_global))
   allocate (coeffsdef0_local(coeffs_per_element, nelms_global))
   allocate (coeffsdei0(coeffs_per_element, nelms_global))
   allocate (coeffsdei0_local(coeffs_per_element, nelms_global))
   allocate (coeffsdef(coeffs_per_element, nelms_global))
   allocate (coeffsdef_local(coeffs_per_element, nelms_global))
   allocate (coeffsdei(coeffs_per_element, nelms_global))
   allocate (coeffsdei_local(coeffs_per_element, nelms_global))
   allocate (coeffsvpf(coeffs_per_element, nelms_global))
   allocate (coeffsvpf_local(coeffs_per_element, nelms_global))
   allocate (coeffsvpi(coeffs_per_element, nelms_global))
   allocate (coeffsvpi_local(coeffs_per_element, nelms_global))
   call set_matrix_index(diff2_mat, 170)
   call create_mat(diff2_mat, 1, 1, icomplex, 1)
   do itri = 1, nelms
      call define_element_quadrature(itri, int_pts_main, int_pts_tor)
      call define_fields(itri, 0, 1, 0)
      tempxx = intxx2(mu79(:, :, OP_1), nu79(:, :, OP_1))
      tempxx = tempxx + smooth_par*(intxx2(mu79(:, :, OP_DZZ), nu79(:, :, OP_DZZ)) + intxx2(mu79(:, :, OP_DRR), nu79(:, :, OP_DRR)))
#ifdef USE3D
      tempxx = tempxx + smooth_par*intxx3(mu79(:, :, OP_DPP), nu79(:, :, OP_DPP), ri4_79)
#endif
      call insert_block(diff2_mat, itri, 1, 1, tempxx, MAT_ADD)
   end do
   call finalize(diff2_mat)
   call set_matrix_index(diff3_mat, 171)
   call create_mat(diff3_mat, 1, 1, icomplex, 1)
   do itri = 1, nelms
      call define_element_quadrature(itri, int_pts_main, int_pts_tor)
      call define_fields(itri, FIELD_PSI+FIELD_I+FIELD_B2I, 1, 0)
      tempxx = intxx2(mu79(:, :, OP_1), nu79(:, :, OP_1))
      !tempxx = tempxx + smooth_pres*(intxx2(mu79(:, :, OP_DZZ), nu79(:, :, OP_DZZ)) + intxx2(mu79(:, :, OP_DRR), nu79(:, :, OP_DRR)))
!#ifdef USE3D
      !tempxx = tempxx + smooth_pres*intxx3(mu79(:, :, OP_DPP), nu79(:, :, OP_DPP), ri4_79)
!#endif
      tempxx = tempxx + smooth_dens_parallel*(&
            + intxx5(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),ri_79*pstx79(:,OP_DR)-bfptx79(:,OP_DZ),ri_79*pstx79(:,OP_DR)-bfptx79(:,OP_DZ),b2i79(:,OP_1)) &
            - intxx5(mu79(:,:,OP_DZ),nu79(:,:,OP_DR),ri_79*pstx79(:,OP_DR)-bfptx79(:,OP_DZ),ri_79*pstx79(:,OP_DZ)+bfptx79(:,OP_DR),b2i79(:,OP_1)) &
            + intxx5(mu79(:,:,OP_DZ),nu79(:,:,OP_DP),ri_79*pstx79(:,OP_DR)-bfptx79(:,OP_DZ),ri2_79*bztx79(:,OP_1),b2i79(:,OP_1)) &
            - intxx5(mu79(:,:,OP_DR),nu79(:,:,OP_DZ),ri_79*pstx79(:,OP_DZ)+bfptx79(:,OP_DR),ri_79*pstx79(:,OP_DR)-bfptx79(:,OP_DZ),b2i79(:,OP_1)) &
            + intxx5(mu79(:,:,OP_DR),nu79(:,:,OP_DR),ri_79*pstx79(:,OP_DZ)+bfptx79(:,OP_DR),ri_79*pstx79(:,OP_DZ)+bfptx79(:,OP_DR),b2i79(:,OP_1)) &
            - intxx5(mu79(:,:,OP_DR),nu79(:,:,OP_DP),ri_79*pstx79(:,OP_DZ)+bfptx79(:,OP_DR),ri2_79*bztx79(:,OP_1),b2i79(:,OP_1)) &
            + intxx5(mu79(:,:,OP_DP),nu79(:,:,OP_DZ),ri2_79*bztx79(:,OP_1),ri_79*pstx79(:,OP_DR)-bfptx79(:,OP_DZ),b2i79(:,OP_1)) &
            - intxx5(mu79(:,:,OP_DP),nu79(:,:,OP_DR),ri2_79*bztx79(:,OP_1),ri_79*pstx79(:,OP_DZ)+bfptx79(:,OP_DR),b2i79(:,OP_1)) &
            + intxx5(mu79(:,:,OP_DP),nu79(:,:,OP_DP),ri2_79*bztx79(:,OP_1),ri2_79*bztx79(:,OP_1),b2i79(:,OP_1)) &
            )

   call insert_block(diff3_mat, itri, 1, 1, tempxx, MAT_ADD)
   end do
   call finalize(diff3_mat)

end subroutine init_particles

subroutine advance_particles(tinc)
   use basic  !For MPI variables
   use omp_lib
   implicit none

   real, intent(in) :: tinc  !Time increment for particle advance
   type(elfield), dimension(nneighbors + 1) :: elcoefs
   real, parameter :: twopi = 6.283185307179586476925286766559
   real    :: dtp, trem
   integer :: ielm, itri, ipart, ierr, istep
   integer :: nlost, ipe, lunf, gunf, nstep
   integer :: nhop, thop, nreas, treas  !Stats on ptcle movement w/in local domain
   integer :: isghost, totin
   integer :: ndevice, blocky, starty, endy
   integer :: npart_local

   if (tinc .le. 0.0) return

   npart_local = ipart_end - ipart_begin + 1
#ifdef _OPENACC
   blocky = hostrank + 1
   starty = npart_local/num_devices*(blocky - 1) + ipart_begin
   endy = npart_local/num_devices*blocky + ipart_begin - 1
   if (blocky == num_devices) endy = ipart_end
!$acc enter data create(elfieldcoefs)
!$acc update device(elfieldcoefs)
!$acc enter data create(pdata(starty:endy))
!$acc update device(pdata(starty:endy))
#endif
   do istep = 1, psubsteps_particle
#ifdef _OPENACC
      starty = npart_local/num_devices*(blocky - 1) + ipart_begin
      endy = npart_local/num_devices*blocky + ipart_begin - 1
      if (blocky == num_devices) endy = ipart_end
#else
      starty = npart_local/ncols*hostrank + ipart_begin
      endy = npart_local/ncols*(hostrank+1) + ipart_begin - 1
      if (hostrank == ncols-1) endy = ipart_end
#endif
#ifndef _OPENACC
!$omp parallel do default(none) shared(pdata,tinc,nparticles,starty,endy,psubsteps_particle,istep) PRIVATE(ierr)
#endif
!$acc parallel loop present(pdata(starty:endy))
      do ipart = starty, endy
         if (pdata(ipart)%deleted) cycle  !Skip if element is empty
         !call rk4(pdata(ielm)%ion(ipart), tinc, itri, ierr)
         call rk4(pdata(ipart), tinc, istep .eq. psubsteps_particle, ierr)
         if (ierr .eq. 1) then ! Particle exited local+ghost domain -> lost
            pdata(ipart)%deleted = .true.
            !pdata(ipart)%x = pdata(ipart)%x0
            !pdata(ipart)%v = pdata(ipart)%v0
            !pdata(ipart)%wt = 0.
            !call mesh_search(pdata(ipart)%jel, pdata(ipart)%x, itri)
            !pdata(ipart)%jel = itri
            !pdata(ipart)%kel(:) = itri
            !cycle !Break out of tinc loop, go to next particle.
         end if
      end do!ielm
#ifndef _OPENACC
!$OMP END PARALLEL do
#endif
   end do
#ifdef _OPENACC
   starty = npart_local/num_devices*(blocky - 1) + ipart_begin
   endy = npart_local/num_devices*blocky + ipart_begin - 1
   if (blocky == num_devices) endy = ipart_end
!$acc exit data copyout(pdata(starty:endy))
!$acc exit data detach(elfieldcoefs)
#endif
end subroutine advance_particles

!---------------------------------------------------------------------------
! 4th-order Runge-Kutta integrator, no adaptive time step control.
!  Four derivative evaluations per step.
! (Adapted from Numerical Recipes in C, 2nd edition, 1994, pp. 712-713).
!
subroutine rk4(part, dt, last_step, ierr)
!$acc routine seq
   implicit none

   type(particle), intent(inout) :: part
   real, intent(in) :: dt
   logical, intent(in) :: last_step
   !type(elfield), dimension(nneighbors+1), intent(in) :: fh
   integer :: itri, itri2
   integer, dimension(4) :: kel
   integer, intent(out) :: ierr

   real, parameter :: onethird = 1.0/3.0
   real, dimension(3) :: k1, k2, k3, k4, y1, x2, lr, lr2
   real, dimension(vspdims) :: l1, l2, l3, l4, z1
   real :: hh, m1, m2, m3, m4, w1, n1, n2, n3, n4
   real :: B0, B0inv
   real, dimension(3) :: B_cyl, deltaB, bhat, gradB0, gradB1, xtemp
   real, dimension(vspdims) :: vtemp
   type(xgeomterms)   :: geomterms
   integer ktri, ipoint
   real :: ran_temp, dB1
   real :: x, y, dphi, vR, vphi
   real :: wtt, wt2, wt3
   real :: gradcoef, df0de, df0dxi, f0, f00
   real :: energy, spsq, xi

   !ierr = 0
   hh = 0.5*dt
   !B0=part%B0
   kel = part%kel
   itri = part%jel
   !if (ierr .eq. 1) return
   !1st step
   call fdot(part%x, part%v, part%wt, k1, l1, m1, n1, itri, kel, part%f0, ierr, part%sps, part%B0)
   if (ierr .eq. 1) return
   y1 = part%x + hh*k1; z1 = part%v + hh*l1; w1 = part%wt + hh*m1
   !write(0,*) k1(1),k1(2),k1(3)

   !2nd step
   call fdot(y1, z1, w1, k2, l2, m2, n2, itri, kel, part%f0, ierr, part%sps, part%B0)
   if (ierr .eq. 1) return
   y1 = part%x + hh*k2; z1 = part%v + hh*l2; w1 = part%wt + hh*m2

   !3rd step
   call fdot(y1, z1, w1, k3, l3, m3, n3, itri, kel, part%f0, ierr, part%sps, part%B0)
   if (ierr .eq. 1) return
   y1 = part%x + dt*k3; z1 = part%v + dt*l3; w1 = part%wt + dt*m3

   !4th step
   call fdot(y1, z1, w1, k4, l4, m4, n4, itri, kel, part%f0, ierr, part%sps, part%B0)
   if (ierr .eq. 1) return
   part%x = part%x + onethird*dt*(k2 + k3 + 0.5*(k1 + k4))
   part%v = part%v + onethird*dt*(l2 + l3 + 0.5*(l1 + l4))
   part%wt = part%wt + onethird*dt*(m2 + m3 + 0.5*(m1 + m4))
   if ((.not. particle_linear_particle .eq. 1) .and. (part%wt < -10.)) then
      part%wt = 0.
      !ierr=1
      !return
   end if
   !if ((abs(part%wt) > 0.05)) then
   !   part%wt = 0.
   !   !ierr=1
   !   !return
   !end if

   if (last_step) then
  !Determine final particle element location
   xtemp = part%x
   vtemp = part%v
   call get_geom_terms(xtemp, itri, geomterms, .false., ierr)
   if (ierr .eq. 1) return
#ifdef USEST
   call update_geom_terms_st(geomterms, elfieldcoefs(itri), .false.)
#endif
   call getBcyl(xtemp, elfieldcoefs(itri), geomterms, B_cyl, deltaB, gradB0, gradB1, dB1)
   B0inv = 1.0/sqrt(dot_product(B_cyl, B_cyl))  !1/magnitude of B
   bhat = B_cyl*B0inv                         !Unit vector in b direction
   !write(0,*) part%gid
   !if (part%gid==2117425261) then
   !if (abs(dot_product(elfieldcoefs(itri)%rho, geomterms%g))>kinetic_rhomax_particle) part%wt=0.
   ! write(0,*) (abs(dot_product(elfieldcoefs(itri)%rho, geomterms%g)))
   !write(0,*) part%v(1)**2+2.*qm_ion*part%v(2)/B0inv, 1./B0inv, part%v(1), itri
   !if (part%x(2)<0) write(0,*) 'cc',dot_product(geomterms%g,elfieldcoefs(itri)%rst),dot_product(geomterms%g,elfieldcoefs(itri)%zst)
   !write(0,*) xtemp(1),xtemp(2),xtemp(3),dot_product(geomterms%g,elfieldcoefs(itri)%Bzv0),itri
   !endif
   x2 = xtemp
   part%dB = 0
   do ipoint = 1, 4
      if (gyroaverage_particle.eq.1) then
         select case (ipoint)
         case (1)
            ran_temp = mod(xtemp(1)*1.e6, 10.)/10.*twopi
            lr(1) = cos(ran_temp)
            lr(3) = sin(ran_temp)
            lr(2) = -(lr(1)*bhat(1) + lr(3)*bhat(3))/bhat(2)
            if (vspdims .eq. 2) then
               lr = lr/sqrt(dot_product(lr, lr))*sqrt(2.0*qm_ion(part%sps)*part%v(2)/B0inv)/qm_ion(part%sps)*B0inv
            else
               lr = lr/sqrt(dot_product(lr, lr))*sqrt(2.0*qm_ion(part%sps)*part%v(5)/B0inv)/qm_ion(part%sps)*B0inv
            end if
            lr(2) = lr(2)/x2(1)
            x2 = x2 + lr
         case (2)
            x2 = x2 - 2*lr
         case (3)
            x2 = x2 + lr
            lr(2) = lr(2)*x2(1)
            lr2(1) = lr(2)*bhat(3) - lr(3)*bhat(2)
            lr2(2) = (lr(3)*bhat(1) - lr(1)*bhat(3))/x2(1)
            lr2(3) = lr(1)*bhat(2) - lr(2)*bhat(1)
            x2 = x2 + lr2
         case (4)
            x2 = x2 - 2*lr2
         end select
         call mesh_search(kel(ipoint), x2, itri2)
         if (itri2 < 0) then
            ierr = 1
            return
         end if
         part%kel(ipoint) = itri2
         part%kx(:, ipoint) = x2
         call get_geom_terms(x2, itri2, geomterms, .false., ierr)
         if (ierr .eq. 1) return
         call getBcyl(x2, elfieldcoefs(itri2), geomterms, B_cyl, deltaB, gradB0, gradB1, dB1)
         part%dB = part%dB + dot_product(deltaB, bhat)*B0inv
      else
         part%kel(ipoint) = itri
         part%kx(:, ipoint) = x2
      end if
   end do
   if (gyroaverage_particle.eq.1) then
      part%dB = part%dB/4.
   else
      part%dB = dot_product(deltaB, bhat)*B0inv
   end if
   part%B0 = 1./B0inv ! fluid particle
   part%jel = itri

      !spsq = part%v(1)*part%v(1) + 2.0*qm_ion(part%sps)*part%v(2)*part%B0
      !xi = part%v(1)/sqrt(spsq)
      !energy = m_ion(part%sps)*spsq/2./1.6e-19
      !part%x0(1)=energy
      !part%x0(2)=xi
      !write(0,*) xi
   ! end if



    !call evalf0(part%x, part%v(1), sqrt(2.0*qm_ion(part%sps)*part%v(2)/B0inv), elfieldcoefs(itri), geomterms, part%sps, f0, gradcoef, df0de, df0dxi)
    !!call evalf0(part%x, part%v(1), sqrt(2.0*qm_ion(part%sps)*part%v(2)*part%B0), elfieldcoefs(itri), geomterms, part%sps, f0, gradcoef, df0de, df0dxi) ! fluid particle
    !if (part%f0/f0>10) then
    !   !if (floor(mod(part%x(1)*100000,500.0))==0) then
    !      write(0,*) "33333333333333333333",part%f0/f0
    !      part%x=part%x0
    !      part%v=part%v0
    !      part%wt=0.
    !      call mesh_search(part%jel, part%x, itri)
    !      part%jel=itri
    !      part%kel(:)=itri
    !      !endif
    !endif
    !if (f0/part%f0>10) then
    !   part%x=part%x0
    !   part%v=part%v0
    !   part%wt=0.
    !   call mesh_search(part%jel, part%x, itri)
    !   part%jel=itri
    !   part%kel(:)=itri
    !   write(0,*) "55555555555555",part%f0/f0
    !endif
   endif
end subroutine rk4

subroutine fdot(x, v, w, dxdt, dvdt, dwdt, dEpdt, itri, kel, f00, ierr, sps, B00)
!$acc routine seq
   use basic
   implicit none

   real, dimension(3), intent(inout)                          :: x
   !real, dimension(3)                             :: x
   real, dimension(3), intent(out)                            :: dxdt
   real, dimension(vspdims), intent(in)                       :: v
   real, dimension(vspdims), intent(out)                      :: dvdt
   real, intent(in)                                           :: w
   real, intent(in)                                           :: f00
   real, intent(out)                                          :: dwdt
   real, intent(out)                                          :: dEpdt
   !type(elfield), dimension(nneighbors+1), target, intent(in) :: fh
   integer, intent(inout)                                     :: itri
   integer, dimension(4), intent(inout)                       :: kel
   integer, intent(out)                                       :: ierr
   !integer, intent(in)                                        :: gid
   !real, intent(out)                                        :: deltaB0
   integer, intent(in)                                       :: sps
   real, intent(in)                                       :: B00
   !real, intent(in)                                        :: df0de, df0dpsi
   real, parameter :: g_mks = 9.8067 ! earth avg surf grav accel in m/s/s
   type(elfield), target  :: fh_hop
   real, dimension(3) :: lr, lr2
   integer :: itri2
   type(xgeomterms)   :: geomterms, geomterms2
   real, dimension(3) :: x2, B_cyl, B0_cyl, B_cyl2, Jcyl, BxgrdB, deltaB, deltaB2, E_cyl, E_cyl2
   real, dimension(3) :: bhat, bhat0, svec, svec0, Bstar, Bstar0
   real, dimension(vspdims)                                  :: v2, vs, vu
   real, dimension(3) :: dBdR, dBdphi, dBdz, dB0dR, dB0dphi, dB0dz, dB1dR, dB1dphi, dB1dz
   real, dimension(3) :: dB0dR2, dB0dphi2, dB0dz2, dB1dR2, dB1dphi2, dB1dz2
   real, dimension(3) :: gradB0, gradB, gradB02, gradB1, gradB12, dEdR, dEdphi, dEdz
   real, dimension(3) :: weqv0, weqv1, weqvD, weqvD1, gradpsi, gradf, gradpe, gradrho
   real, dimension(3) :: gradpsi0, gradpsi1, gradBz
   vectype, dimension(3) :: temp
   real f0, T0, tmp1, tmp2, tmp3, tmp4, tmp5, df0de, df0dxi, spd, gradcoef, dB1, dB12, j0xb, ne0, te0, dBdt, dEdt, dxidt
   real dpsidt, df0dpsi
   real :: Rinv, Rinv2, B0inv, Binv, Bss, Bss0
   real :: dRdphi, dZdphi, di, dxdR, dxdZ, dydR, dydZ
   real, dimension(3)                            :: dxdt2, dxdt0
   real, dimension(vspdims)                      :: dvdt0, weqa1
!!$OMP THREADPRIVATE(Rinv)
   integer :: tridex
   integer :: ipoint
   real :: Pphi, ran_temp, hh

   !x=x_in
   !ierr = 0

   !Need terms to compute fields to calculate acceleration
   call get_geom_terms(x, itri, geomterms, vspdims .eq. 2, ierr)
   if (ierr .ne. 0) then
      !write(0,*) 'outoutout'
      return
   end if

#ifdef USEST
   !Get electric field components
   Rinv = 1.0/dot_product(elfieldcoefs(itri)%rst,geomterms%g)
   !Rinv = 1.0
   dRdphi = dot_product(elfieldcoefs(itri)%rst,geomterms%dphi)
   dZdphi = dot_product(elfieldcoefs(itri)%zst,geomterms%dphi)
   di = 1./(dot_product(elfieldcoefs(itri)%rst,geomterms%dr)*dot_product(elfieldcoefs(itri)%zst,geomterms%dz) -&
      dot_product(elfieldcoefs(itri)%zst,geomterms%dr)*dot_product(elfieldcoefs(itri)%rst,geomterms%dz))
   dxdR = di*dot_product(elfieldcoefs(itri)%zst,geomterms%dz)
   dxdZ = -di*dot_product(elfieldcoefs(itri)%rst,geomterms%dz)
   dydR = -di*dot_product(elfieldcoefs(itri)%zst,geomterms%dr)
   dydZ = di*dot_product(elfieldcoefs(itri)%rst,geomterms%dr)
   call update_geom_terms_st(geomterms, elfieldcoefs(itri), vspdims.eq.2)
#else
   !if (itor.eq.1) Rinv = 1.0/x(1)
   Rinv = 1.0/x(1)
#endif

   !Calculate time derivatives
   !call getBcyl(x, fhptr, geomterms, B_cyl, deltaB, gradB0)
   call getBcyl(x, elfieldcoefs(itri), geomterms, B0_cyl, deltaB, gradB0, gradB1, dB1)
   B0inv = 1.0/sqrt(dot_product(B0_cyl, B0_cyl))  !1/magnitude of B
   bhat0 = B0_cyl*B0inv                         !Unit vector in b direction

   if (gyroaverage_particle.eq.1) then
      x2 = x
      deltaB = 0.
      E_cyl = 0.
      B0_cyl = 0.
      dB0dR = 0.
      dB0dphi = 0.
      dB0dz = 0.
      dB1dR = 0.
      dB1dphi = 0.
      dB1dz = 0.
      gradpe = 0.
      j0xb = 0.
      do ipoint = 1, 4
         select case (ipoint)
         case (1)
            ran_temp = mod(x(1)*1.e6, 10.)/10.*twopi
            lr(1) = cos(ran_temp)
            lr(3) = sin(ran_temp)
            lr(2) = -(lr(1)*bhat0(1) + lr(3)*bhat0(3))/bhat0(2)
            lr = lr/sqrt(dot_product(lr, lr))*sqrt(2.0*qm_ion(sps)*v(2)/B0inv)/qm_ion(sps)*B0inv
            lr(2) = lr(2)/x2(1)
            x2 = x2 + lr
         case (2)
            x2 = x2 - 2*lr
         case (3)
            x2 = x2 + lr
            lr(2) = lr(2)*x2(1)
            lr2(1) = lr(2)*bhat0(3) - lr(3)*bhat0(2)
            lr2(2) = (lr(3)*bhat0(1) - lr(1)*bhat0(3))/x2(1)
            lr2(3) = lr(1)*bhat0(2) - lr(2)*bhat0(1)
            x2 = x2 + lr2
         case (4)
            x2 = x2 - 2*lr2
         end select
         !call mesh_search(kel(ipoint), x2, itri2)
         !kel(ipoint)=itri2
         call get_geom_terms(x2, kel(ipoint), geomterms2, .false., ierr)
         if (ierr .ne. 0) then
            return
         end if
#ifdef USEST
         !Get electric field components
         Rinv2 = 1.0/dot_product(elfieldcoefs(ipoint)%rst,geomterms2%g)
         !Rinv = 1.0
         !dRdphi = dot_product(elfieldcoefs(itri)%rst,geomterms%dphi)
         !dZdphi = dot_product(elfieldcoefs(itri)%zst,geomterms%dphi)
         !di = 1./(dot_product(elfieldcoefs(itri)%rst,geomterms%dr)*dot_product(elfieldcoefs(itri)%zst,geomterms%dz) -&
         !   dot_product(elfieldcoefs(itri)%zst,geomterms%dr)*dot_product(elfieldcoefs(itri)%rst,geomterms%dz))
         !dxdR = di*dot_product(elfieldcoefs(itri)%zst,geomterms%dz)
         !dxdZ = -di*dot_product(elfieldcoefs(itri)%rst,geomterms%dz)
         !dydR = -di*dot_product(elfieldcoefs(itri)%zst,geomterms%dr)
         !dydZ = di*dot_product(elfieldcoefs(itri)%rst,geomterms%dr)
         !call update_geom_terms_st(geomterms, elfieldcoefs(itri), vspdims.eq.2)
#else
         Rinv2 = 1.0/x2(1)
#endif
         call getEcyl(x2, elfieldcoefs(kel(ipoint)), geomterms2, E_cyl2)
         !call getBcyl(x2, elfieldcoefs(kel(ipoint)), geomterms2, B_cyl2, deltaB2, gradB02, gradB12, dB12)
         call getBcylprime(x2, elfieldcoefs(kel(ipoint)), geomterms2, B_cyl2, deltaB2, dB0dR2, dB0dphi2, dB0dz2, dB1dR2, dB1dphi2, dB1dz2)
         !call getBcyl_last(x2, elfieldcoefs(kel(ipoint)), geomterms2, B_cyl2, deltaB_last2)
         E_cyl = E_cyl + E_cyl2
         B0_cyl = B0_cyl + B_cyl2
         deltaB = deltaB + deltaB2
         dB0dR = dB0dR + dB0dR2
         dB0dphi = dB0dphi + dB0dphi2
         dB0dz = dB0dz + dB0dz2
         dB1dR = dB1dR + dB1dR2
         dB1dphi = dB1dphi + dB1dphi2
         dB1dz = dB1dz + dB1dz2

         if (kinetic_thermal_ion_particle.eq.1) then
            temp(1) = dot_product(geomterms2%dr, elfieldcoefs(kel(ipoint))%pe)
            temp(3) = dot_product(geomterms2%dz, elfieldcoefs(kel(ipoint))%pe)
#ifdef USECOMPLEX
            temp(2) = dot_product(geomterms2%g, elfieldcoefs(kel(ipoint))%pe)*rfac_particle/x2(1)
#elif defined(USE3D)
            temp(2) = Rinv2*dot_product(geomterms2%dphi, elfieldcoefs(kel(ipoint))%pe)
#else
            temp(2) = 0.
#endif
#ifdef USECOMPLEX
            gradpe = gradpe + real(temp*exp(rfac_particle*x2(2)))
#else
            gradpe = gradpe + temp
#endif
            temp(1) = dot_product(geomterms2%dr, elfieldcoefs(kel(ipoint))%pe0)
            temp(3) = dot_product(geomterms2%dz, elfieldcoefs(kel(ipoint))%pe0)
#ifdef USEST
            temp(2) = Rinv2*dot_product(geomterms2%dphi, elfieldcoefs(kel(ipoint))%pe0)
#else
            temp(2) = 0.
#endif
            j0xb = j0xb - dot_product(temp, deltaB2)*B0inv
         endif
      end do
      E_cyl = E_cyl/4.
      B0_cyl = B0_cyl/4.
      deltaB = deltaB/4.
      dB0dR = dB0dR/4.
      dB0dphi = dB0dphi/4.
      dB0dz = dB0dz/4.
      dB1dR = dB1dR/4.
      dB1dphi = dB1dphi/4.
      dB1dz = dB1dz/4.
      if (kinetic_thermal_ion_particle.eq.1) then
         gradpe = gradpe/4.
         j0xb = j0xb/4.
      endif
   else
      kel(:) = itri
      !call getEcylprime(x, fhptr, geomterms, E_cyl, dEdR, dEdphi, dEdz)
      call getEcyl(x, elfieldcoefs(itri), geomterms, E_cyl)
      call getBcylprime(x, elfieldcoefs(itri), geomterms, B0_cyl, deltaB, dB0dR, dB0dphi, dB0dz, dB1dR, dB1dphi, dB1dz)
         !call getBcyl_last(x, fhptr, geomterms, B_cyl2, deltaB_last)
      if (kinetic_thermal_ion_particle.eq.1) then
         !temp(1) = dot_product(geomterms%dr, elfieldcoefs(itri)%pe)
         !temp(3) = dot_product(geomterms%dz, elfieldcoefs(itri)%pe)
         temp(1) = dot_product(geomterms%dr, elfieldcoefs(itri)%ne)
         temp(3) = dot_product(geomterms%dz, elfieldcoefs(itri)%ne)
#ifdef USECOMPLEX
         !temp(2) = dot_product(geomterms%g, elfieldcoefs(itri)%pe)*rfac_particle/x(1)
         temp(2) = dot_product(geomterms%g, elfieldcoefs(itri)%ne)*rfac_particle/x(1)
#elif defined(USE3D)
         !temp(2) = Rinv*dot_product(geomterms%dphi, elfieldcoefs(itri)%pe)
         temp(2) = Rinv*dot_product(geomterms%dphi, elfieldcoefs(itri)%ne)
#else
         temp(2) = 0.
#endif
#ifdef USECOMPLEX
         gradpe=real(temp * exp(rfac_particle*x(2)))
#else
         gradpe=temp
#endif
         gradpe=gradpe*dot_product(geomterms%g, elfieldcoefs(itri)%te0)

         !temp(1) = dot_product(geomterms%dr, elfieldcoefs(itri)%pe0)
         !temp(3) = dot_product(geomterms%dz, elfieldcoefs(itri)%pe0)
         temp(1) = dot_product(geomterms%dr, elfieldcoefs(itri)%ne0)
         temp(3) = dot_product(geomterms%dz, elfieldcoefs(itri)%ne0)
#ifdef USEST
         temp(2) = Rinv*dot_product(geomterms%dphi, elfieldcoefs(itri)%ne0)
#else
         temp(2) = 0.
#endif
         !j0xb=-dot_product(temp,deltaB)*B0inv
         j0xb=-dot_product(temp,deltaB)*B0inv*dot_product(geomterms%g, elfieldcoefs(itri)%te0)
         !if (real(dot_product(geomterms%g, fhptr%psiv0))<0.21) then
         !   gradpe=0.
            !j0xb=0.
         !endif
      endif
   end if

   B_cyl = B0_cyl + deltaB
   dBdR = dB0dR + dB1dR
   dBdphi = dB0dphi + dB1dphi
   dBdz = dB0dz + dB1dz
   !call getBcyl(x, elfieldcoefs(itri), geomterms, B_cyl, deltaB, gradB0, gradB1, dB1)
   !write(0,*) B_cyl(1), B_cyl(1), B_cyl(2)
   B0inv = 1.0/sqrt(dot_product(B0_cyl, B0_cyl))  !1/magnitude of B
   bhat0 = B0_cyl*B0inv                         !Unit vector in b direction
   Binv = 1.0/sqrt(dot_product(B_cyl, B_cyl))  !1/magnitude of B
   bhat = B_cyl*Binv                         !Unit vector in b direction

   !if (particle_linear_particle .eq. 1) then
      E_cyl=E_cyl-dot_product(E_cyl,B0_cyl)*B0inv
   !else
      !E_cyl=E_cyl-dot_product(E_cyl,B_cyl)*Binv
   !endif
   ! E_cyl=0.

   ! Gradient of B0 = grad(B.B)/(2 B0) = (B . grad B)/B0
   gradB0(1) = dot_product(bhat0, dB0dR)
   gradB0(2) = Rinv*dot_product(bhat0, dB0dphi)
   gradB0(3) = dot_product(bhat0, dB0dz)

   ! Curl of bhat = curl(B/B0) = curl(B)/B0 - (grad B0 x B)/(B0**2)
   BxgrdB(1) = B0_cyl(2)*gradB0(3) - B0_cyl(3)*gradB0(2)
   BxgrdB(2) = B0_cyl(3)*gradB0(1) - B0_cyl(1)*gradB0(3)
   BxgrdB(3) = B0_cyl(1)*gradB0(2) - B0_cyl(2)*gradB0(1)

   Jcyl(1) = Rinv*dB0dphi(3) - dB0dz(2)
   Jcyl(2) = dB0dz(1) - dB0dR(3)
   Jcyl(3) = dB0dR(2) - Rinv*dB0dphi(1)
   !if (itor.eq.1) Jcyl(3) = Jcyl(3) + Rinv*B_cyl(2)
   Jcyl(3) = Jcyl(3) + Rinv*B0_cyl(2)

   !tmp1 = (v(1)*v(1)) * (B0inv*B0inv)/qm_ion
   !tmp2 = tmp1*B0inv + v(2)*(B0inv*B0inv)
   !weqvD = tmp2*BxgrdB + tmp1*(Jcyl - dot_product(bhat, Jcyl)*bhat)
   svec = (Jcyl + BxgrdB*B0inv)*B0inv

   Bstar0 = B0_cyl + (v(1)/qm_ion(sps))*svec
   !Bstar0 = B0_cyl ! fluid particle
   Bss0 = dot_product(Bstar0, bhat0)

   svec0 = v(2)*gradB0 ! - g_mks/qm_ion
   !svec = svec - E_cyl  ! - g_mks/qm_ion
   !svec = v(2)*gradB0  ! - g_mks/qm_ion
   !svec0 = 0. ! fluid particle

   dxdt0(1) = (v(1)*Bstar0(1) + bhat0(2)*svec0(3) - bhat0(3)*svec0(2))/Bss0
   dxdt0(2) = (v(1)*Bstar0(2) + bhat0(3)*svec0(1) - bhat0(1)*svec0(3))/Bss0
   dxdt0(3) = (v(1)*Bstar0(3) + bhat0(1)*svec0(2) - bhat0(2)*svec0(1))/Bss0

   dvdt0(1) = -qm_ion(sps)*dot_product(Bstar0, svec0)/Bss0
   !dvdt0(1) = -qm_ion*dot_product(bhat0, svec)
   !dvdt0(1) = 0 ! fluid particle
   dvdt0(2) = 0. !magnetic moment is conserved.

   !call getBcylprime(x, elfieldcoefs(itri), geomterms, B_cyl, deltaB, dBdR, dBdphi, dBdz, .false.)
   !call getBcyl(x, fhptr, geomterms, B_cyl, deltaB, gradB0, gradB1,1)

   ! Gradient of B0 = grad(B.B)/(2 B0) = (B . grad B)/B0
   if (particle_linear_particle .eq. 1) then
      gradB(1) = dot_product(bhat0, dBdR)
      gradB(2) = Rinv*dot_product(bhat0, dBdphi)
      gradB(3) = dot_product(bhat0, dBdz)
   else
      gradB(1) = dot_product(bhat, dBdR)
      gradB(2) = Rinv*dot_product(bhat, dBdphi)
      gradB(3) = dot_product(bhat, dBdz)
   end if

   ! Curl of bhat = curl(B/B0) = curl(B)/B0 - (grad B0 x B)/(B0**2)
   if (particle_linear_particle .eq. 1) then
      BxgrdB(1) = B0_cyl(2)*gradB(3) - B0_cyl(3)*gradB(2)
      BxgrdB(2) = B0_cyl(3)*gradB(1) - B0_cyl(1)*gradB(3)
      BxgrdB(3) = B0_cyl(1)*gradB(2) - B0_cyl(2)*gradB(1)
   else
      BxgrdB(1) = B_cyl(2)*gradB(3) - B_cyl(3)*gradB(2)
      BxgrdB(2) = B_cyl(3)*gradB(1) - B_cyl(1)*gradB(3)
      BxgrdB(3) = B_cyl(1)*gradB(2) - B_cyl(2)*gradB(1)
   end if

   Jcyl(1) = Rinv*dBdphi(3) - dBdz(2)
   Jcyl(2) = dBdz(1) - dBdR(3)
   Jcyl(3) = dBdR(2) - Rinv*dBdphi(1)
   !if (itor.eq.1) Jcyl(3) = Jcyl(3) + Rinv*B_cyl(2)
   Jcyl(3) = Jcyl(3) + Rinv*B_cyl(2)

   !tmp1 = (v(1)*v(1)) * (B0inv*B0inv)/qm_ion
   !tmp2 = tmp1*B0inv + v(2)*(B0inv*B0inv)
   !weqvD = tmp2*BxgrdB + tmp1*(Jcyl - dot_product(bhat, Jcyl)*bhat)
   if (particle_linear_particle .eq. 1) then
      svec = (Jcyl + BxgrdB*B0inv)*B0inv
   else
      svec = (Jcyl + BxgrdB*Binv)*Binv
   end if

   Bstar = B_cyl + (v(1)/qm_ion(sps))*svec
   !Bstar = B_cyl ! fluid particle
   Bss = dot_product(Bstar, bhat)

   svec = v(2)*gradB ! - g_mks/qm_ion
   svec = svec - E_cyl  ! - g_mks/qm_ion
   !svec = -E_cyl ! fluid particle
   !svec = 0.
   !svec = v(2)*gradB0  ! - g_mks/qm_ion

   if (particle_linear_particle .eq. 1) then
      dxdt(1) = (v(1)*Bstar(1) + bhat0(2)*svec(3) - bhat0(3)*svec(2))/Bss0
      dxdt(2) = (v(1)*Bstar(2) + bhat0(3)*svec(1) - bhat0(1)*svec(3))/Bss0
      dxdt(3) = (v(1)*Bstar(3) + bhat0(1)*svec(2) - bhat0(2)*svec(1))/Bss0

      dvdt(1) = -qm_ion(sps)*(dot_product(Bstar0, svec)+dot_product(Bstar, svec0)-dot_product(Bstar0,svec0))/Bss0
      !dvdt(1) = -qm_ion*dot_product(bhat0, svec)
      !dvdt(1) = 0
      dvdt(2) = 0. !magnetic moment is conserved.
   else
      dxdt(1) = (v(1)*Bstar(1) + bhat(2)*svec(3) - bhat(3)*svec(2))/Bss
      dxdt(2) = (v(1)*Bstar(2) + bhat(3)*svec(1) - bhat(1)*svec(3))/Bss
      dxdt(3) = (v(1)*Bstar(3) + bhat(1)*svec(2) - bhat(2)*svec(1))/Bss

      dvdt(1) = -qm_ion(sps)*dot_product(Bstar, svec)/Bss
      !dvdt(1) = -qm_ion*dot_product(bhat0, svec)
      !dvdt(1) = 0 ! fluid particle
      dvdt(2) = 0. !magnetic moment is conserved.
   end if

   !Weights evolve in delta-f method only.
   ! V1 = (ExB)/(B**2) + U deltaB/B
   ! weqv1(1) = ((E_cyl(2)*B_cyl(3) - E_cyl(3)*B_cyl(2))*B0inv + 0*v(1)*deltaB(1))*B0inv
   ! weqv1(2) = ((E_cyl(3)*B_cyl(1) - E_cyl(1)*B_cyl(3))*B0inv + 0*v(1)*deltaB(2))*B0inv
   ! weqv1(3) = ((E_cyl(1)*B_cyl(2) - E_cyl(2)*B_cyl(1))*B0inv + 0*v(1)*deltaB(3))*B0inv
   spd = sqrt(v(1)*v(1) + 2.0*qm_ion(sps)*v(2)/B0inv)
   weqv1 = dxdt - dxdt0
   weqa1 = dvdt - dvdt0
   gradpsi0 = 0.0
   gradpsi0(1) = dot_product(elfieldcoefs(itri)%psiv0, geomterms%dr)
   gradpsi0(3) = dot_product(elfieldcoefs(itri)%psiv0, geomterms%dz)
   !gradpsi1 = 0.0
   !gradpsi1(1) = dot_product(geomterms%dr, elfieldcoefs(itri)%psiv1)
   !gradpsi1(3) = dot_product(geomterms%dz, elfieldcoefs(itri)%psiv1)
!#ifdef USECOMPLEX
   !gradpsi1(2) = dot_product(geomterms%g, elfieldcoefs(itri)%psiv1)*rfac_particle
!#else
   !gradpsi1(2) = dot_product(geomterms%dphi, elfieldcoefs(itri)%psiv1)
!#endif

   dpsidt = dot_product(weqv1,gradpsi0)!+dot_product(dxdt0, gradpsi1)
   dpsidt = dpsidt + (1./qm_ion(sps)) * weqa1(1) *  dot_product(elfieldcoefs(itri)%Bzv0, geomterms%g)*B0inv
   gradbz = 0.0
   gradbz(1) = dot_product(elfieldcoefs(itri)%Bzv0, geomterms%dr)
   gradbz(3) = dot_product(elfieldcoefs(itri)%Bzv0, geomterms%dz)
   dpsidt = dpsidt + (1./qm_ion(sps)) * v(1) * dot_product(weqv1,gradbz)*B0inv
   dBdt = dot_product(weqv1, gradB0)+0*dot_product(dxdt0, gradB-gradB0)
   dpsidt = dpsidt - (1./qm_ion(sps)) * v(1) * dot_product(elfieldcoefs(itri)%Bzv0, geomterms%g)*B0inv**2*dBdt
   dEdt = m_ion(sps)*v(1)*weqa1(1) + q_ion(sps)*v(2)*dBdt
   dxidt = -q_ion(sps)*v(2)*1.99/(m_ion(sps)*v(1)**2*0.5+q_ion(sps)*v(2)/B0inv)**2*dEdt

   !dEdt = 0. ! fluid particle
   !dxidt = 0. ! fluid particle

   ! vD = (1/(e B**3))(M_i U**2 + mu B)(B x grad B) + ((M_i U**2)/(eB**2))*J_perp
   tmp1 = (v(1)*v(1))*(B0inv*B0inv)/qm_ion(sps)
   tmp2 = tmp1*B0inv + v(2)*(B0inv*B0inv)
   weqvD = tmp2*BxgrdB + tmp1*(Jcyl - dot_product(bhat, Jcyl)*bhat)
   weqv0 = v(1)*bhat + weqvD
   tmp2 = (v(4)*v(4))*(B0inv*B0inv)/qm_ion(sps)*B0inv
   tmp2 = tmp1*B0inv
   weqvD1 = tmp2*BxgrdB + tmp1*(Jcyl - dot_product(bhat, Jcyl)*bhat)

   ne0 = dot_product(geomterms%g, elfieldcoefs(itri)%ne0)

   if (kinetic_thermal_ion_particle.eq.1) dEdt = dEdt + 1*q_ion(sps)*v(1)*(dot_product(-gradpe,bhat0)+j0xb)/ne0

   ! tmp1=dot_product(-gradpe,bhat)+j0xb
   ! temp(1) = dot_product(geomterms%dr, fhptr%ne)
   ! temp(3) = dot_product(geomterms%dz, fhptr%ne)
   ! temp(2) = dot_product(geomterms%g, fhptr%ne)*rfac_particle/x(1)
   ! te0 = dot_product(geomterms%g, fhptr%te0)
   ! gradpe=real(temp * te0 * exp(rfac_particle*x(2)))
   ! temp(1) = dot_product(geomterms%dr, fhptr%ne0)
   ! temp(3) = dot_product(geomterms%dz, fhptr%ne0)
   ! temp(2) = 0.
   ! j0xb=-dot_product(temp*te0,deltaB)*B0inv

   ! write(0,*) tmp1,(dot_product(-gradpe,bhat)+j0xb)
   gradrho = 0.0
   gradrho(1) = dot_product(elfieldcoefs(itri)%rho, geomterms%dr)
   gradrho(3) = dot_product(elfieldcoefs(itri)%rho, geomterms%dz)
#ifdef USEST
   gradrho(2) = Rinv*dot_product(elfieldcoefs(itri)%rho, geomterms%dphi)
#endif
   call evalf0(x, v(1), sqrt(2.0*qm_ion(sps)*v(2)/B0inv), elfieldcoefs(itri), geomterms, sps, 1./B0inv, f0, df0dpsi, df0de, df0dxi)
   !call evalf0(x, v(1), sqrt(2.0*qm_ion(sps)*v(2)*B00), elfieldcoefs(itri), geomterms, sps, f0, gradcoef, df0de, df0dxi) ! fluid particle
   gradf = gradrho*gradcoef
   if (iconst_f0_particle.eq.1) then
      gradf = gradf*f0/f00
      df0de = df0de*f0/f00
      df0dxi = df0dxi*f0/f00
   endif
   ! write(0,*) v(1), sqrt(2.0*qm_ion(sps)*v(2)/B0inv)
   ! gradcoef, df0de
   ! call evalf0_advance(x, v, 1.0/B0inv, B_cyl(2), elfieldcoefs(itri), geomterms, sps, f0, gradf, df0de, df0dxi)

   ! vD = (1/(e B**3))(M_i U**2 + mu B)(B x grad B) + ((M_i U**2)/(eB**2))*J_perp
   tmp3 = tmp2*dot_product(gradB0, deltaB)
   !tmp4 = tmp3 - v(1)*dot_product(Jcyl,E_cyl)/qm_ion*(B0inv*B0inv*B0inv)
   tmp4 = tmp3
   tmp5 = -tmp4*dot_product(elfieldcoefs(itri)%Bzv0, geomterms%g)*gradcoef
  
   !dBdt(1) = Rinv*dEdphi(3) - dEdz(2)
   !dBdt(1) = 0.
   !dBdt(2) = dEdz(1) - dEdR(3)
   !dBdt(3) = dEdR(2) - Rinv*dEdphi(1)
   !dBdt(3) = 0.
   !dBdt(3) = dBdt(3) + Rinv*E_cyl(2)
   !dB0dt=dot_product(bhat, deltaB-deltaB_last)/dt_particle/t0_norm_particle
   !dB0dt=dB1/dt_particle/t0_norm_particle
   !dB0dt=dot_product(bhat, -dBdt)
   !dwdt = (w - 1.0)*(dot_product(weqv1, gradf0) - &
   !     q_ion*dot_product(weqv0, E_cyl)*df0de)/f0
   !if (sps == 1) then
      !dwdt = (-1.0)*(1.*dot_product(weqv1, gradpsi*gradcoef) + &
      !    !0*dot_product(real(fhptr%Bzv0), geomterms%g)*B0inv*dot_product(bhat,E_cyl-v(2)*gradB1)*gradcoef+&
 !!dwdt = (-1.0)*(deltaB(1) + &
 !!dwdt = (-1.0)*(((E_cyl(1)*B_cyl(2)-E_cyl(2)*B_cyl(1))*B0inv+v(1)*deltaB(3))*B0inv*1e6 + &
      !     !1.*q_ion*dot_product(weqvD, E_cyl)*df0de + 1.*tmp5)*f0/f00
      !     !1.*q_ion*(dot_product(weqvD, E_cyl)+1*v(5)*dB0dt)*df0de + 0.*tmp5)
      !     !1.*q_ion*(dot_product(weqv0, E_cyl)+0*v(5)*dB0dt)*df0de + 0.*tmp5)
      !     !1.*q_ion*(1*dot_product(weqvD, E_cyl)+1*dot_product(weqv0, -gradpe)+0*v(5)*dB0dt)*df0de + 0.*tmp5)
      !     1.*q_ion*(1*dot_product(weqvD, E_cyl)+1.*v(1)*dot_product(E_cyl,bhat)+1*v(1)*(dot_product(-gradpe,bhat)+j0xb)/ne0+0*v(2)*dB0dt)*df0de + 0.*tmp5)
      !     !1.*q_ion*(1*dot_product(weqvD, E_cyl)+1*v(4)*(dot_product(-gradpe,bhat)+j0xb)+1*v(5)*dB0dt)*df0de + 0.*tmp5)
      !     !1.*q_ion*(dot_product(weqv0, bhat)*dot_product(E_cyl,bhat)+0*v(5)*dB0dt)*df0de + 0.*tmp5)
      !     ! dxdt=0.
      !     ! dvdt=0.
   !else
      dwdt = (-1.0)*(dpsidt*df0dpsi + &
                     !0*dot_product(real(fhptr%Bzv0), geomterms%g)*B0inv*dot_product(bhat,E_cyl-v(2)*gradB1)*gradcoef+&
                     !dwdt = (-1.0)*(deltaB(1) + &
                     !dwdt = (-1.0)*(((E_cyl(1)*B_cyl(2)-E_cyl(2)*B_cyl(1))*B0inv+v(1)*deltaB(3))*B0inv*1e6 + &
                     !1.*q_ion*dot_product(weqvD, E_cyl)*df0de + 1.*tmp5)*f0/f00
                     !1.*q_ion*(dot_product(weqvD, E_cyl)+1*v(5)*dB0dt)*df0de + 0.*tmp5)
                     !1.*q_ion*(dot_product(weqv0, E_cyl)+0*v(5)*dB0dt)*df0de + 0.*tmp5)
                     !1.*q_ion*(1*dot_product(weqvD, E_cyl)+1*dot_product(weqv0, -gradpe)+0*v(5)*dB0dt)*df0de + 0.*tmp5)
                     !1.*q_ion*(1*dot_product(weqvD, E_cyl) + 1.*v(1)*dot_product(E_cyl, bhat))*(df0de - v(1)/spd**3/m_ion*df0dxi) + 0.*tmp5 + &
                     !1.*qm_ion*(1*dot_product(weqvD1,E_cyl)/v(1)+1.*dot_product(E_cyl,bhat)-v(2)*(dot_product(gradB0,deltaB)*B0inv+dot_product(gradB1,bhat)))*df0dxi/spd)
                     +dEdt*df0de + dxidt*df0dxi)

      !1.*q_ion*(dot_product(weqv0, bhat)*dot_product(E_cyl,bhat)+0*v(5)*dB0dt)*df0de + 0.*tmp5)
   !end if
   !if (deltaB_last(2).ne.0.) write(0,*) deltaB(2),deltaB_last(2)
   !dEpdt = q_ion*(dot_product(weqvD, E_cyl)+v(2)*dB0dt)
   !write(0,*) dt_particle, t0_norm_particlk
   !Pphi = dot_product(fhptr%psiv0, geomterms%g) + (1./qm_ion) * v(1) *  dot_product(fhptr%Bzv0, geomterms%g) *B0inv!orbit
   !Pphi = dot_product(fhptr%psiv0, geomterms%g) !orbit
   !dwdt = abs(Pphi/f00-1)/(dt_particle*t0_norm_particle)!orbit
   !dwdt = abs((v(1)**2+2.*qm_ion*v(2)/B0inv)/f00-1)/(dt_particle*t0_norm_particle)!orbit
   !dwdt = 1./(dt_particle*t0_norm_particle)
   !dwdt = f00/(dt_particle*t0_norm_particle)
   !dwdt = 0
   !if (linear_particle) then
   !dwdt = dwdt * f00
   !else
   !dwdt = dwdt * (f00-w)
   !endif
   if (particle_linear_particle .eq. 0) then
      dwdt = dwdt*(1 - w)
      !dEpdt = dEpdt *(w+(1-w)*dot_product(deltaB,bhat)*B0inv)
      !dvdt(1) = dvdt(1)+qm_ion*(dot_product(-gradpe,bhat)+j0xb)/ne0
   else
      dxdt = dxdt0
      dvdt = dvdt0
   end if
   !dxdt=0.
   !dvdt=0.
   dxdt(2) = Rinv*dxdt(2)  !phi-dot = (v_phi / R) for cylindrical case
#ifdef USEST
   dxdt(1) = dxdt(1)-dRdphi*dxdt(2)
   dxdt(3) = dxdt(3)-dZdphi*dxdt(2)
   dxdt2(1) = dxdR*dxdt(1)+dxdZ*dxdt(3)
   dxdt2(3) = dydR*dxdt(1)+dydZ*dxdt(3)
   dxdt2(2) = dxdt(2)
   dxdt=dxdt2
#endif

   !dwdt = (1.-w)/(dt_particle*t0_norm_particle)
   !B0=1./B0inv
   !deltaB0=dot_product(deltaB,bhat)
   !deltaB0=f0
end subroutine fdot

subroutine particle_scaleback(scalefac)
   use basic
   use field
   use arrays
   use mesh_mod
   implicit none

   include 'mpif.h'

   vectype, intent(in) :: scalefac
   integer :: ipart, ielm, ierr

   if (hostrank == 0) then
      do ipart = ipart_begin, ipart_end
         pdata(ipart)%wt = pdata(ipart)%wt*scalefac
         pdata(ipart)%dB = pdata(ipart)%dB*scalefac
      end do
      do ielm = ielm_min, ielm_max
         !elfieldcoefs(ielm)%psiv1=elfieldcoefs(ielm)%psiv1*scalefac
         !elfieldcoefs(ielm)%Bzv1=elfieldcoefs(ielm)%Bzv1*scalefac
         !elfieldcoefs(ielm)%Bfv=elfieldcoefs(ielm)%Bfv*scalefac
         elfieldcoefs(ielm)%B1 = elfieldcoefs(ielm)%B1*scalefac
      end do
      call delete_particle(.true.)
   end if
   call MPI_Bcast(ipart_begin, 1, MPI_INTEGER, 0, hostcomm, ierr)
   call MPI_Bcast(ipart_end, 1, MPI_INTEGER, 0, hostcomm, ierr)
   call MPI_Bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
   call mult(p_f_par, scalefac)
   call mult(p_f_perp, scalefac)
   call mult(den_f_1, scalefac)
   call mult(p_i_par, scalefac)
   call mult(p_i_perp, scalefac)
   call mult(v_i_par, scalefac)
   call mult(den_i_1, scalefac)
end subroutine particle_scaleback
!---------------------------------------------------------------------------
subroutine delete_particle(exchange)
   use basic
   use diagnostics
   implicit none
   include 'mpif.h'

   logical, intent(in) :: exchange
   !integer, intent(in) :: ipart
   !type(elplist), intent(inout) :: pbuf
   !integer, intent(out) :: ierr

   integer :: ipart, npart
   integer :: ierr
   integer, dimension(:), allocatable :: recvcounts, displs
   integer :: sendcount, irow
   real :: depar
   type(particle), dimension(:), allocatable :: pdata_temp

   !Error checking
   !np = pbuf%np
   !if (ipart.lt.1.or.ipart.gt.nparticles) then
   !   ierr = 2
   !   return
   !endif

   !Replace the particle with the last one in the array
   if (particle_nodelete.eq.0) then
   npart = ipart_begin - 1
   depar = 0.
   do ipart = ipart_begin, ipart_end
      !write(0,*) pdata(ipart)%deleted
      if (pdata(ipart)%deleted) then
         cycle
      end if
      npart = npart + 1
      !spsq = pdata(ipart)%v(1)**2 + 2.0*qm_ion*pdata(ipart)%v(2)*pdata(ipart)%B0
      !depar=depar+pdata(ipart)%wt*pdata(ipart)%wt2
      depar = depar + pdata(ipart)%wt*pdata(ipart)%dE
      !depar=depar+pdata(ipart)%dE
      pdata(ipart)%dE = 0
      if (ipart .ne. npart) pdata(npart) = pdata(ipart)
   end do
   ipart_end = npart
   endif
   write (0, *) ipart_begin, ipart_end
   if (exchange) then
      allocate (recvcounts(nrows))
      allocate (displs(nrows))
      sendcount = ipart_end - ipart_begin + 1
      recvcounts(rowrank + 1) = sendcount
      call MPI_ALLGATHER(sendcount, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, rowcomm, ierr)
      displs(1) = 0
      do irow = 2, nrows
         displs(irow) = displs(irow - 1) + recvcounts(irow - 1)
      end do
      if ((ipart_begin .eq. 1) .and. (rowrank .gt. 0)) then
         pdata(displs(rowrank + 1) + 1:displs(rowrank + 1) + sendcount) = pdata(ipart_begin:ipart_end)
         ipart_begin = displs(rowrank + 1) + 1
         ipart_end = displs(rowrank + 1) + sendcount
      end if
      allocate (pdata_temp(ipart_end-ipart_begin+1))
      pdata_temp=pdata(ipart_begin:ipart_end)
      call MPI_ALLGATHERV(pdata_temp, sendcount, MPI_particle, pdata, recvcounts, displs, MPI_particle, rowcomm, ierr)
      deallocate (pdata_temp)
      !call mpi_allreduce(depar, depar, 1, MPI_DOUBLE_PRECISION, MPI_SUM, rowcomm, ierr)
      call mpi_barrier(rowcomm, ierr)
      nparticles = sum(recvcounts)
      ipart_begin = nparticles/nrows*rowrank + 1
      ipart_end = nparticles/nrows*(rowrank + 1)
      if (rowrank == nrows - 1) ipart_end = nparticles
      write (0, *) "nparticles", rowrank, ipart_begin, ipart_end
      deallocate (recvcounts)
      deallocate (displs)
   end if
   !ierr = 0
end subroutine delete_particle

!---------------------------------------------------------------------------
subroutine particle_step(pdt)
   use basic
   use diagnostics
   use auxiliary_fields
   implicit none
   include 'mpif.h'

   real, intent(in) :: pdt

   real    :: tstart, tend
   integer :: istep, ierr, ipart
   integer :: isubcycle

   call mpi_barrier(mpi_comm_world, ierr)
   if (kinetic_thermal_ion.eq.1) then
      call set_parallel_velocity
   endif
   call calculate_electric_fields(linear)
   do isubcycle=1,particle_subcycles
      if (kinetic_thermal_ion.eq.1) then
         call set_den_smooth
      endif
      !Advance particle positions
      call get_field_coefs(0)
      call mpi_barrier(mpi_comm_world, ierr)
      !call MPI_Win_fence(0, win_elfieldcoefs)
#ifdef _OPENACC
      if (hostrank < num_devices) then
#else
      if (hostrank < ncols) then
#endif
         call second_new(tstart)
         call advance_particles(pdt/particle_substeps/particle_subcycles)
         call second_new(tend)
         write (0, '(A,I7,A,f9.2,A)') 'Particle advance completed', particle_substeps, ' steps in', &
            tend - tstart, ' seconds.'
      end if
      call mpi_barrier(mpi_comm_world, ierr)
      if (hostrank < 1) then
         call second(tstart)
         if ((isubcycle==particle_subcycles).and.(mod(ntime - ntime0, ntimepr) == 0)) then
            call delete_particle(.true.)
         else
            call delete_particle(.false.)
         end if
         call second(tend)
         write (0, '(A,f9.2,A)') 'Particle delete in', tend - tstart, ' seconds.'
      end if
      call MPI_Bcast(ipart_begin, 1, MPI_INTEGER, 0, hostcomm, ierr)
      call MPI_Bcast(ipart_end, 1, MPI_INTEGER, 0, hostcomm, ierr)
      call MPI_Bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
      !Compute particle pressure tensor components
      call update_particle_pressure
      call mpi_barrier(mpi_comm_world, ierr)

      if (kinetic_thermal_ion.eq.1) then
         call set_density
      endif
      call mpi_barrier(mpi_comm_world, ierr)
   enddo

end subroutine particle_step
!---------------------------------------------------------------------------
subroutine update_particle_pressure
   use basic
   use arrays
   use diagnostics
   implicit none
   include 'mpif.h'

   real    :: tstart, tend
   integer :: ierr

   !call MPI_Bcast(ipart_begin, 1, MPI_INTEGER, 0, hostcomm, ierr)
   !call MPI_Bcast(ipart_end, 1, MPI_INTEGER, 0, hostcomm, ierr)
   !call MPI_Bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
   !call mpi_barrier(mpi_comm_world, ierr)
   !call MPI_Win_fence(0, win_pdata)
   !Deposit particles to compute RHS.
   call second(tstart)
   call particle_pressure_rhs

   !Invert mass matrix to solve for field components
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   if (myrank .eq. 0) then
      call second(tend)
      write (0, '(A,f9.2,A)') 'Pressure tensor RHS vecs calculated in', tend - tstart, ' sec.'
   end if

   call second(tstart)
   call solve_pi_tensor
   if (myrank .eq. 0) then
      call second(tend)
      write (0, '(A,f9.2,A)') 'Pressure tensor LHS vecs calculated in', tend - tstart, ' sec.'
   end if
end subroutine update_particle_pressure
!---------------------------------------------------------------------------
subroutine finalize_particles
   use arrays
   implicit none
   include 'mpif.h'

   !integer :: nelms, ielm
   integer :: ierr

   if (allocated(neighborlist)) deallocate (neighborlist)

   !if (allocated(pdata)) then
   !nelms = size(pdata)
   !do ielm=1,nelms
   !if (allocated(pdata(ielm)%ion)) deallocate(pdata(ielm)%ion)
   !enddo !ielm

   !deallocate(pdata)
   !endif

   CALL MPI_Win_free(win_pdata, ierr)
   CALL MPI_Win_free(win_elfieldcoefs, ierr)
   CALL MPI_Win_free(win_mesh_coord, ierr)
   CALL MPI_Win_free(win_mesh_nodes, ierr)
!    if (allocated(jmppar)) then
   !nelms = size(jmppar)
   !do ielm=1,nelms
   !if (allocated(jmppar(ielm)%ion)) deallocate(jmppar(ielm)%ion)
   !enddo !ielm

   !deallocate(jmppar)
   !endif

   !if (allocated(jinbuf)) deallocate(jinbuf)
   !if (allocated(dnbr)) deallocate(dnbr)
   !if (allocated(dnlist)) deallocate(dnlist)

!    call mpi_type_free(mpi_particle, ielm)
   call destroy_field(p_i_par); call destroy_field(p_i_perp)
end subroutine finalize_particles
!---------------------------------------------------------------------------
subroutine find_element_neighbors
   use basic
   implicit none

#ifdef USE3D
   integer, parameter :: maxconnect = 96 !Max # of faces converging on any mesh node
   integer, parameter :: ifverts = 4     !3 or 4 verts define a face
#else
   integer, parameter :: maxconnect = 24  !Max # of edges converging on any mesh node
   integer, parameter :: ifverts = 2     !Two verts define an edge
#endif

   type d_face
      integer :: el0, side
      integer, dimension(ifverts - 1) :: v
   end type d_face

   type face
      type(d_face), dimension(maxconnect) :: o
      integer :: n = 0
   end type face

   type(face), dimension(:), allocatable  :: facelist
   integer, dimension(nodes_per_element)  :: enode
   integer, dimension(ifverts, nneighbors) :: sidevecsub !Vector subscripts for iface IDs
   integer, dimension(ifverts) :: iface
   integer, dimension(1) :: ml
   integer :: ielm, side, v1, ivrt
   logical :: ep

   !if (nelms.lt.1) return

#ifdef USE3D
   if (nodes_per_element .eq. 6) then !prisms
      !Define prism faces in terms of nodes
      sidevecsub(:, 1) = [2, 5, 6, 3]
      sidevecsub(:, 2) = [3, 6, 4, 1]
      sidevecsub(:, 3) = [1, 4, 5, 2]
      sidevecsub(:, 4) = [1, 2, 3, 1]
      sidevecsub(:, 5) = [4, 6, 5, 4]
#else
   if (nodes_per_element .eq. 3) then !triangles
      !Define triangle edges in terms of nodes
      sidevecsub(:, 1) = [2, 3]
      sidevecsub(:, 2) = [3, 1]
      sidevecsub(:, 3) = [1, 2]
#endif
      !Allocate storage, initialize
      allocate (neighborlist(nneighbors, nelms_global), facelist(nnodes_global))
      neighborlist = -1
      !Loop over local elements
      do ielm = 1, nelms_global
         enode = mesh_nodes(:, ielm)
         !call get_element_nodes(ielm, enode)

         !Loop over faces of this element
         do side = 1, nneighbors
            iface = enode(sidevecsub(:, side))
            ml = minloc(iface)
#ifdef USE3D
            if (iface(1) == iface(4)) then
               iface(1:3) = cshift(iface(1:3), ml(1) - 1)  !Arrange so lowest-indexed node comes 1st
               iface(4) = iface(1)
               v1 = iface(1)
               if (iface(2) > iface(3)) iface = iface([1, 3, 2, 1])
            else
               iface = cshift(iface, ml(1) - 1)  !Arrange so lowest-indexed node comes 1st
               v1 = iface(1)
               if (iface(2) > iface(4)) iface = iface([1, 4, 3, 2])
            end if
#else
            iface = cshift(iface, ml(1) - 1)  !Arrange so lowest-indexed node comes 1st
            v1 = iface(1)
#endif

            !Search if the face is already present
            ep = .false.
            do ivrt = 1, facelist(v1)%n
               if (veceq(facelist(v1)%o(ivrt)%v, iface(2:))) then
                  !Yes, update neighbor table
                  neighborlist(side, ielm) = facelist(v1)%o(ivrt)%el0
                  neighborlist(facelist(v1)%o(ivrt)%side, facelist(v1)%o(ivrt)%el0) = ielm
                  !if (side==4) write(0,*) iface, neighborlist(:,ielm)
                  ep = .true.
                  exit
               end if !veceq...
            end do !ivrt

            if (.not. ep) then !Face was not present; add it.
               facelist(v1)%n = facelist(v1)%n + 1
               if (facelist(v1)%n .gt. maxconnect) then !out of range
                  print *, 'Error: too many connections in find_element_neighbors.'
                  deallocate (facelist, neighborlist)
                  return
               end if !n out-of-range`
               facelist(v1)%o(facelist(v1)%n)%v = iface(2:)
               facelist(v1)%o(facelist(v1)%n)%el0 = ielm
               facelist(v1)%o(facelist(v1)%n)%side = side
               !if (side==4) write(0,*) iface, facelist(v1)%o(facelist(v1)%n)%v
            end if
         end do !side
      end do !ielm

      deallocate (facelist)
   else
      if (myrank .eq. 0) print *, nodes_per_element, ' nodes per element; cannot find neighbors.'
   end if !nodes_per_element...
end subroutine find_element_neighbors
!---------------------------------------------------------------------------
logical function veceq(v1, v2)
   use basic
   implicit none

#ifdef USE3D
   integer, dimension(3), intent(in) :: v1, v2
   veceq = (v1(1) .eq. v2(1) .and. v1(2) .eq. v2(2) .and. v1(3) .eq. v2(3))
#else
   integer, dimension(1), intent(in) :: v1, v2
   veceq = (v1(1) .eq. v2(1))
#endif
end function veceq
!---------------------------------------------------------------------------
subroutine mesh_search(initial_simplex, final_position, final_simplex)
!$acc routine seq
   implicit none

   integer, intent(in) :: initial_simplex
   real, dimension(3), intent(inout) :: final_position
   integer, intent(out) :: final_simplex
   integer, parameter :: max_count = 1000
   real, dimension(3) :: pos
   real, dimension(3, nodes_per_element) :: x
   real, dimension(3) :: area, h
   integer :: simplex, min_edge, search_count
   logical :: located

   simplex = initial_simplex
   located = .false.
   search_count = 0
   !pos=final_position
   !#ifdef USE3D
   if (final_position(2) .lt. 0.0) then
      final_position(2) = mod(final_position(2), toroidal_period_particle) + toroidal_period_particle
   else
      final_position(2) = mod(final_position(2), toroidal_period_particle)
   end if
   !if (final_position(2).lt.0.0) then
!!if (neighborlist(4,simplex)==-1) write(0,*) '444',neighborlist(1,simplex),neighborlist(2,simplex),neighborlist(3,simplex),neighborlist(4,simplex),neighborlist(5,simplex)
   !simplex=neighborlist(4,simplex)
   !elseif (final_position(2) .ge.toroidal_period_particle) then
   !final_position(2) = final_position(2) - toroidal_period_particle
!!if (neighborlist(5,simplex)==-1) write(0,*) '555',neighborlist(1,simplex),neighborlist(2,simplex),neighborlist(3,simplex),neighborlist(4,simplex),neighborlist(5,simplex)
   !simplex=neighborlist(5,simplex)
   !endif
   !#endif
   do while (.not. located)
      !write(0,*) simplex
      !if (search_count > max_count) then
      !write(0,*) simplex, final_position(1), final_position(2), final_position(3)
      !endif
      x = mesh_coord(:, :, simplex)
#ifdef USE3D
      if (x(2, 4) == 0) x(2, 4) = toroidal_period_particle
      if (final_position(2) < x(2, 1)) then
         !if (neighborlist(4,simplex)==-1) write(0,*) '444',neighborlist(1,simplex),neighborlist(2,simplex),neighborlist(3,simplex),neighborlist(4,simplex),neighborlist(5,simplex)
         simplex = neighborlist(4, simplex)
         search_count = search_count + 1
         cycle
      elseif (final_position(2) > x(2, 4)) then
         !if (neighborlist(5,simplex)==-1) write(0,*) '555',neighborlist(1,simplex),neighborlist(2,simplex),neighborlist(3,simplex),neighborlist(4,simplex),neighborlist(5,simplex)
         simplex = neighborlist(5, simplex)
         search_count = search_count + 1
         cycle
      end if
#endif
      area(1) = tri_area(final_position(1), final_position(3), x(1, 2), x(3, 2), x(1, 3), x(3, 3))
      area(2) = tri_area(x(1, 1), x(3, 1), final_position(1), final_position(3), x(1, 3), x(3, 3))
      area(3) = tri_area(x(1, 1), x(3, 1), x(1, 2), x(3, 2), final_position(1), final_position(3))
      if (all(area >= 0.)) then
         located = .true.
         final_simplex = simplex
         exit
      end if
      h(1) = area(1)/sqrt((x(1, 2) - x(1, 3))**2 + (x(3, 2) - x(3, 3))**2)
      h(2) = area(2)/sqrt((x(1, 1) - x(1, 3))**2 + (x(3, 1) - x(3, 3))**2)
      h(3) = area(3)/sqrt((x(1, 1) - x(1, 2))**2 + (x(3, 1) - x(3, 2))**2)

      min_edge = 1
      if (h(2) < h(1)) min_edge = 2
      if (h(3) < h(min_edge)) min_edge = 3
      if (neighborlist(min_edge, simplex) == -1) then
         final_simplex = -1
         exit
      end if
      ! do while (neighborlist(min_edge,simplex)==-1)
      !min_edge=min_edge+1
      !if (min_edge>3) min_edge=1
      ! enddo
      !located = .true.
      !final_simplex = simplex
      simplex = neighborlist(min_edge, simplex)

      search_count = search_count + 1
      if (search_count > max_count) then
         !  write(0,*) simplex, final_position(1), final_position(2), final_position(3)
         final_simplex = -1
         exit
      end if
      !if (search_count == max_count) then
      !final_simplex = -1
      !exit
      !endif
   end do
end subroutine
!---------------------------------------------------------------------------
subroutine element_local(ielm, R, Phi, Z, xi, eta, zi, co, sn)
!$acc routine seq
   implicit none

   integer, intent(in) :: ielm
   real, intent(in) :: R, Z, Phi
   real, intent(out) :: xi, eta, zi, co, sn
   real, dimension(3, nodes_per_element) :: x
   real :: x2, z2, hi, b

   x = mesh_coord(:, :, ielm)
   x2 = x(1, 2) - x(1, 1)
   z2 = x(3, 2) - x(3, 1)
   hi = 1./sqrt(x2**2 + z2**2)
   co = x2*hi
   sn = z2*hi
   b = co*(x(1, 3) - x(1, 1)) + sn*(x(3, 3) - x(3, 1))
   xi = (R - x(1, 1))*co + (Z - x(3, 1))*sn - b
   eta = -(R - x(1, 1))*sn + (Z - x(3, 1))*co
   zi = Phi - x(2, 1)
end subroutine element_local
!---------------------------------------------------------------------------
! Compute terms for a function and its partial derivatives with respect to
!  R and z in the reduced quintic expansion at position x.
subroutine get_geom_terms(x, ielm, gh, ic2, ierr)
!$acc routine seq
   implicit none

   real, dimension(3), intent(inout) :: x
   integer, intent(inout) :: ielm
   !type(elfield), dimension(nneighbors+1), intent(in) :: fh
   integer, dimension(20) :: mi = [0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 3, 2, 1, 0]
   integer, dimension(20) :: ni = [0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 2, 3, 4, 5]
   integer, dimension(4) :: li = [0, 1, 2, 3]
   type(xgeomterms), intent(out) :: gh  !Geometric terms handle
   logical, intent(in)  :: ic2          !Compute 2nd derivative terms?
   integer, intent(out) :: ierr

   real :: xi, zi, eta, dxi, deta, co, sn, co2, sn2
   real :: d2xi, d2eta, dxieta
   integer pp
#ifdef USE3D
   real    :: gtmp, drtmp, dztmp, zpow
   integer :: ii, jj
#endif
   real, dimension(2) :: mmsa
   integer ktri
   real, dimension(-2:5) :: xipow, etapow

   ierr = 0
   !tridex = -1

   if (ielm .lt. 1) ielm = 1
   call mesh_search(ielm, x, ktri)
   ielm = ktri
   if (ielm .le. 0) then !The triangle is not in the local partition
      ierr = 1
      return
   end if

   !tridex = ensemble_index(fh, ielm)

   !call get_element_data(ielm, eldat)
   !call global_to_local(eldat, x(1), x(2), x(3), xi, zi, eta)
   call element_local(ielm, x(1), x(2), x(3), xi, eta, zi, co, sn)

   !Precompute powers
   xipow(-1) = 0.0; etapow(-1) = 0.0
   xipow(0) = 1.0; etapow(0) = 1.0
   do pp = 1, 5
      xipow(pp) = xi*xipow(pp - 1)
      etapow(pp) = eta*etapow(pp - 1)
   end do

   !Compute terms for function & 1st derivatives
   gh%dr = 0.; gh%dz = 0.
   do pp = 1, coeffs_per_tri
      gh%g(pp) = xipow(mi(pp))*etapow(ni(pp))

      dxi = mi(pp)*xipow(mi(pp) - 1)*etapow(ni(pp))
      deta = xipow(mi(pp))*ni(pp)*etapow(ni(pp) - 1)

      gh%dr(pp) = gh%dr(pp) + co*dxi - sn*deta
      gh%dz(pp) = gh%dz(pp) + sn*dxi + co*deta

#ifdef USE3D
      gtmp = gh%g(pp); drtmp = gh%dr(pp); dztmp = gh%dz(pp)
      do ii = 1, coeffs_per_dphi
         jj = pp + (ii - 1)*coeffs_per_tri
         zpow = zi**li(ii)

         gh%g(jj) = gtmp*zpow
         gh%dr(jj) = drtmp*zpow
         gh%dz(jj) = dztmp*zpow

         !First toroidal derivative
         if (li(ii) .gt. 0) then
            zpow = li(ii)*zi**(li(ii) - 1)
            gh%dphi(jj) = gtmp*zpow
            gh%drphi(jj) = drtmp*zpow
            gh%dzphi(jj) = dztmp*zpow

       !!Second toroidal derivative
            !if (li(ii).gt.1) then
            !   zpow = li(ii) * (li(ii) - 1) * zi**(li(ii) - 2)
            !   gh%drphiphi(jj) = drtmp * zpow
            !   gh%dzphiphi(jj) = dztmp * zpow
            !else
            !   gh%drphiphi(jj) = 0.
            !   gh%dzphiphi(jj) = 0.
            !endif
         else
            gh%dphi(jj) = 0.
            gh%drphi(jj) = 0.
            gh%dzphi(jj) = 0.
            !gh%drphiphi(jj) = 0.
            !gh%dzphiphi(jj) = 0.
         end if
      end do !ii
#endif
   end do !pp

   if (ic2) then !2nd derivative terms
      xipow(-2) = 0.0; etapow(-2) = 0.0
      gh%drr = 0.; gh%drz = 0.; gh%dzz = 0.
      co2 = co**2; sn2 = sn**2
      do pp = 1, coeffs_per_tri
         d2xi = mi(pp)*(mi(pp) - 1)*xipow(mi(pp) - 2)*etapow(ni(pp))
         dxieta = mi(pp)*ni(pp)*xipow(mi(pp) - 1)*etapow(ni(pp) - 1)
         d2eta = xipow(mi(pp))*ni(pp)*(ni(pp) - 1)*etapow(ni(pp) - 2)

         gh%drr(pp) = gh%drr(pp) + d2xi*co2 - &
                      (2.0*dxieta*co - d2eta*sn)*sn
         gh%drz(pp) = gh%drz(pp) + (d2xi - d2eta)*co*sn + dxieta*(co2 - sn2)
         gh%dzz(pp) = gh%dzz(pp) + d2xi*sn2 + &
                      (2.0*dxieta*sn + d2eta*co)*co

#ifdef USE3D
         gtmp = gh%drr(pp); drtmp = gh%drz(pp); dztmp = gh%dzz(pp)
         do ii = 1, coeffs_per_dphi
            jj = pp + (ii - 1)*coeffs_per_tri
            zpow = zi**li(ii)

            gh%drr(jj) = gtmp*zpow
            gh%drz(jj) = drtmp*zpow
            gh%dzz(jj) = dztmp*zpow

       !!First toroidal derivative
            !if (li(ii).gt.0) then
            !   zpow = li(ii) * zi**(li(ii) - 1)
            !   gh%drrphi(jj) = gh%drr(pp) * zpow
            !   gh%drzphi(jj) = gh%drz(pp) * zpow
            !   gh%dzzphi(jj) = gh%dzz(pp) * zpow
            !else
            !   gh%drrphi(jj) = 0.
            !   gh%drzphi(jj) = 0.
            !   gh%dzzphi(jj) = 0.
            !endif
         end do !ii
#endif
      end do !pp
   end if !ic2
end subroutine get_geom_terms

subroutine update_geom_terms_st(gh, fh, ic2)
!$acc routine seq
   implicit none

   type(elfield), intent(in) :: fh
   type(xgeomterms), intent(inout) :: gh  !Geometric terms handle
   type(xgeomterms) :: gh2  !Geometric terms handle
   logical, intent(in)  :: ic2          !Compute 2nd derivative terms?

   real :: rstDR, rstDZ, zstDR, zstDZ, rstDP, zstDP
   real :: rstDRR, rstDRZ, rstDZZ, zstDRR, zstDRZ, zstDZZ
   real :: tempa, tempb, tempc, tempd, tempe
   real :: di, di2
#ifdef USE3D
   real :: rstDRP, rstDZP, zstDRP, zstDZP
#endif

   rstDR = dot_product(fh%rst, gh%dr)
   rstDZ = dot_product(fh%rst, gh%dz)
   zstDR = dot_product(fh%zst, gh%dr)
   zstDZ = dot_product(fh%zst, gh%dz)
   di = 1./(rstDR*zstDZ - rstDZ*zstDR)
   di2 = di*di
#ifdef USE3D
   rstDP = dot_product(fh%rst, gh%dphi)
   zstDP = dot_product(fh%zst, gh%dphi)
#endif

   gh2%g = gh%g
   gh2%dr = di*zstDZ*gh%dr - di*zstDR*gh%dz
   gh2%dz = di*rstDR*gh%dz - di*rstDZ*gh%dr
#ifdef USE3D
   gh2%dphi = -rstDP*gh2%dr - zstDP*gh2%dz + gh%dphi
#endif
   if (ic2) then !2nd derivative terms
      rstDRR = dot_product(fh%rst, gh%drr)
      rstDRZ = dot_product(fh%rst, gh%drz)
      rstDZZ = dot_product(fh%rst, gh%dzz)
      zstDRR = dot_product(fh%zst, gh%drr)
      zstDRZ = dot_product(fh%zst, gh%drz)
      zstDZZ = dot_product(fh%zst, gh%dzz)
      tempa = rstDR*zstDRZ + rstDRR*zstDZ &
              - rstDZ*zstDRR - rstDRZ*zstDR
      tempb = rstDR*zstDZZ + rstDRZ*zstDZ &
              - rstDZ*zstDRZ - rstDZZ*zstDR
      tempc = rstDR*tempb - rstDZ*tempa
      tempd = zstDZ*tempa - zstDR*tempb

      gh2%drr = (zstDZ**2*gh%drr + zstDR**2*gh%dzz - 2*zstDR*zstDZ*gh%drz &
                 + (zstDZ*zstDRZ - zstDR*zstDZZ)*gh%dr &
                 + gh%dz*(zstDR*zstDRZ - zstDZ*zstDRR) &
                 - tempd*gh2%dr)*di2
      gh2%dzz = (rstDZ**2*gh%drr + rstDR**2*gh%dzz - 2*rstDR*rstDZ*gh%drz &
                 + (rstDZ*rstDRZ - rstDR*rstDZZ)*gh%dr &
                 + gh%dz*(rstDR*rstDRZ - rstDZ*rstDRR) &
                 - tempc*gh2%dz)*di2
      gh2%drz = ((rstDR*zstDZ + zstDR*rstDZ)*gh%drz &
                 - rstDR*zstDR*gh%dzz - rstDZ*zstDZ*gh%drr &
                 - (zstDZ*rstDRZ - zstDR*rstDZZ)*gh%dr &
                 - gh%dz*(zstDR*rstDRZ - zstDZ*rstDRR) &
                 - tempd*gh2%dz)*di2
#ifdef USE3D
      rstDRP = dot_product(fh%rst, gh%drphi)
      rstDZP = dot_product(fh%rst, gh%dzphi)
      zstDRP = dot_product(fh%zst, gh%drphi)
      zstDZP = dot_product(fh%zst, gh%dzphi)
      tempe = rstDR*zstDZP + rstDRP*zstDZ &
              - rstDZ*zstDRP - rstDZP*zstDR
      gh2%drphi = (zstDZ*gh%drphi - zstDR*gh%dzphi &
                   + zstDZP*gh%dr - zstDRP*gh%dz &
                   - tempe*gh2%dr)*di
      gh2%dzphi = (rstDR*gh%dzphi - rstDZ*gh%drphi &
                   + rstDRP*gh%dz - rstDZP*gh%dr &
                   - tempe*gh2%dz)*di
      gh2%drphi = -rstDP*gh2%drr - zstDP*gh2%drz &
                  + gh2%drphi
      gh2%dzphi = -rstDP*gh2%drz - zstDP*gh2%dzz &
                  + gh2%dzphi
#endif
   end if !ic2

   gh = gh2

end subroutine update_geom_terms_st
!---------------------------------------------------------------------------
subroutine get_field_coefs(eq)
   use arrays
   use basic
   use auxiliary_fields
   implicit none
   include 'mpif.h'

   !type(elfield), intent(out) :: fh  !Field handle
   integer, intent(in) :: eq
   integer :: ielm, ielm_global
   !logical, intent(in) :: getE
   integer :: isghost
   integer, dimension(:), allocatable :: recvcounts, displs
   integer :: sendcount
!    logical :: use_f = .false.
   integer :: ierr
   real :: factor
   type(elfield), dimension(:), allocatable :: elfieldcoefs_temp

   do ielm = 1, nelms!Always get magnetic field components

#ifdef USE3D
      call m3dc1_ent_getglobalid(3, ielm - 1, ielm_global)
#else
      call m3dc1_ent_getglobalid(2, ielm - 1, ielm_global)
#endif
      ielm_global = ielm_global + 1
      if (eq == 1) then
         call calcavector(ielm, psi_field(0), elfieldcoefs(ielm_global)%psiv0)
         call calcavector(ielm, bz_field(0), elfieldcoefs(ielm_global)%Bzv0)
         call calcavector(ielm, bfp_field(0), elfieldcoefs(ielm_global)%Bfpv0)
         call calcavector(ielm, rho_field, elfieldcoefs(ielm_global)%rho)
         call calcavector(ielm, nf_field, elfieldcoefs(ielm_global)%nf)
         call calcavector(ielm, tf_field, elfieldcoefs(ielm_global)%tf)
         call calcavector(ielm, nfi_field, elfieldcoefs(ielm_global)%nfi)
         call calcavector(ielm, tfi_field, elfieldcoefs(ielm_global)%tfi)
         elfieldcoefs(ielm_global)%Bfpv1 = 0.
         elfieldcoefs(ielm_global)%psiv1 = 0.
         elfieldcoefs(ielm_global)%Bzv1 = 0.
#ifdef USEST
         call calcavector(ielm, rst, elfieldcoefs(ielm_global)%rst)
         call calcavector(ielm, zst, elfieldcoefs(ielm_global)%zst)
#endif
      end if
      call calcavector(ielm, bfp_field(1), elfieldcoefs(ielm_global)%Bfpv1)
      call calcavector(ielm, psi_field(1), elfieldcoefs(ielm_global)%psiv1)
      call calcavector(ielm, bz_field(1), elfieldcoefs(ielm_global)%Bzv1)
      if (kinetic_thermal_ion.eq.1) then
         factor = 1*c_light/ &
                  sqrt(4.*3.14159*n0_norm*(z_ion*e_c)**2/m0_norm)/ &
                  l0_norm*(v0_norm/100.0*b0_norm/1.e4)
         call calcavector(ielm, p_field(1), elfieldcoefs(ielm_global)%pe)
         call calcavector(ielm, densmooth_field, elfieldcoefs(ielm_global)%ne)
         !call calcavector(ielm, den_field(1), elfieldcoefs(ielm_global)%ne)
         elfieldcoefs(ielm_global)%pe=elfieldcoefs(ielm_global)%pe*factor
         if (eq==1) then
            call calcavector(ielm, p_field(0), elfieldcoefs(ielm_global)%pe0)
            elfieldcoefs(ielm_global)%pe0=elfieldcoefs(ielm_global)%pe0*factor
            call calcavector(ielm, te_field(0), elfieldcoefs(ielm_global)%te0)
            elfieldcoefs(ielm_global)%te0=elfieldcoefs(ielm_global)%te0*2.0*factor
            call calcavector(ielm, den_field(0), elfieldcoefs(ielm_global)%ne0)
         endif
      endif
      !Get electric field components if needed
      call calcavector(ielm, ef_r, elfieldcoefs(ielm_global)%er)
      call calcavector(ielm, ef_phi, elfieldcoefs(ielm_global)%ephi)
      call calcavector(ielm, ef_z, elfieldcoefs(ielm_global)%ez)

   end do
   call mpi_barrier(mpi_comm_world, ierr)
   if (hostrank == 1) then
      allocate (recvcounts(nrows))
      allocate (displs(nrows))
      sendcount = nelm_row(rowrank + 1)
      recvcounts = nelm_row
      displs = displs_elm
      allocate (elfieldcoefs_temp(nelm_row(rowrank+1)))
      elfieldcoefs_temp=elfieldcoefs(ielm_min:ielm_max)
call MPI_ALLGATHERV(elfieldcoefs_temp,sendcount,MPI_elfield, elfieldcoefs,recvcounts,displs,MPI_elfield, rowcomm,ierr)
      deallocate (elfieldcoefs_temp)
      deallocate (recvcounts)
      deallocate (displs)
   end if    !call m3dc1_ent_isghost(2, ielm-1, isghost)
   call mpi_barrier(mpi_comm_world, ierr)
end subroutine get_field_coefs

!---------------------------------------------------------------------------
subroutine getBcyl(x, fh, gh, Bcyl, deltaB, gradB0, gradB1, dB1)
!$acc routine seq
   use basic
   implicit none

   real, dimension(3), intent(in) :: x       !Position
   type(elfield), intent(in) :: fh           !Field handle
   type(xgeomterms), intent(in) :: gh        !Geometric terms handle
   real, dimension(3), intent(out) :: Bcyl   !Output total magnetic field
   real, dimension(3), intent(out) :: deltaB !Output perturbed part of Bcyl
   real, dimension(3), intent(out) :: gradB0, gradB1 !Output perturbed part of Bcyl
   real, intent(out) :: dB1 !Output perturbed part of Bcyl

   vectype, dimension(3) :: temp
   real :: Rinv
!!$OMP THREADPRIVATE(Rinv)

   !if (itor.eq.1) Rinv = 1.0/x(1)
#ifdef USEST
   Rinv = 1.0/dot_product(fh%rst,gh%g)
#else
   Rinv = 1.0/x(1)
#endif
   !Rinv = 1.0

   !Total/Equilibrium part
   !B_poloidal_axisymmetric = grad psi x grad phi
   Bcyl(1) = -Rinv*dot_product(fh%psiv0, gh%dz)
   Bcyl(3) = Rinv*dot_product(fh%psiv0, gh%dr)
   !B_toroidal = B_Z / R
   Bcyl(2) = Rinv*dot_product(fh%Bzv0, gh%g)
   Bcyl(1) = Bcyl(1) - dot_product(gh%dr, fh%Bfpv0)
   Bcyl(3) = Bcyl(3) - dot_product(gh%dz, fh%Bfpv0)

   !Perturbed part
   temp(1) = -Rinv*dot_product(gh%dz, fh%psiv1)
   temp(3) = Rinv*dot_product(gh%dr, fh%psiv1)
   temp(2) = Rinv*dot_product(gh%g, fh%Bzv1)
#ifdef USECOMPLEX
   temp(1) = temp(1) - dot_product(gh%dr, fh%Bfpv1)
   temp(3) = temp(3) - dot_product(gh%dz, fh%Bfpv1)
   deltaB = real(temp*exp(rfac_particle*x(2)))
#else
   temp(1) = temp(1) - dot_product(gh%dr, fh%Bfpv1)
   temp(3) = temp(3) - dot_product(gh%dz, fh%Bfpv1)
   deltaB = temp
#endif
   if (.not. (particle_linear_particle .eq. 1)) Bcyl = Bcyl + deltaB
   gradB0(1) = dot_product(fh%B0, gh%dr)
   gradB0(3) = dot_product(fh%B0, gh%dz)
#ifdef USE3D
   gradB0(2) = Rinv*dot_product(fh%B0, gh%dphi)
#else
   gradB0(2) = 0
#endif
   temp(1) = dot_product(gh%dr, fh%B1)
   temp(3) = dot_product(gh%dz, fh%B1)
#ifdef USECOMPLEX
   temp(2) = Rinv*dot_product(gh%g, fh%B1)*rfac_particle
#elif defined(USE3D)
   temp(2) = Rinv*dot_product(fh%B1, gh%dphi)
#else
   temp(2) = 0.
#endif
#ifdef USECOMPLEX
   gradB1 = real(temp*exp(rfac_particle*x(2)))
#else
   gradB1 = temp
#endif
   temp(1) = dot_product(gh%g, fh%B1) - dot_product(gh%g, fh%B1_last)
#ifdef USECOMPLEX
   ! dB1 = real(temp(1) * exp(rfac_particle*x(2)))
   dB1 = 0.
#else
   ! dB1 = temp(1)
   dB1 = 0.
#endif
   if (.not. (particle_linear_particle .eq. 1)) gradB0 = gradB0 + gradB1

   Bcyl = Bcyl*b0_norm_particle/1.e4
   deltaB = deltaB*b0_norm_particle/1.e4
   gradB0 = gradB0*b0_norm_particle/1.e4
   gradB1 = gradB1*b0_norm_particle/1.e4
   dB1 = dB1*b0_norm_particle/1.e4
end subroutine getBcyl
!---------------------------------------------------------------------------
subroutine getBcylprime(x, fh, gh, Bcyl, deltaB, dB0dR, dB0dphi, dB0dz, dB1dR, dB1dphi, dB1dz)
!$acc routine seq
   use basic
   implicit none

   real, dimension(3), intent(in)  :: x
   type(elfield), intent(in)       :: fh
   type(xgeomterms), intent(in)    :: gh
   real, dimension(3), intent(out) :: Bcyl, deltaB, dB0dR, dB0dphi, dB0dz, dB1dR, dB1dphi, dB1dz

   vectype, dimension(3) :: temp, tempR, tempphi, tempz
   real :: Rinv
!!$OMP THREADPRIVATE(Rinv)

   !if (itor.eq.1) Rinv = 1.0/x(1)
#ifdef USEST
   Rinv = 1.0/dot_product(fh%rst,gh%g)
#else
   Rinv = 1.0/x(1)
#endif
   !Rinv = 1.0

   !Total/Equilibrium part
   !B_poloidal_axisym = grad psi x grad phi
   Bcyl(1) = -Rinv*dot_product(fh%psiv0, gh%dz)
   dB0dR(1) = -Rinv*dot_product(fh%psiv0, gh%drz)
   dB0dz(1) = -Rinv*dot_product(fh%psiv0, gh%dzz)

   Bcyl(3) = Rinv*dot_product(fh%psiv0, gh%dr)
   dB0dR(3) = Rinv*dot_product(fh%psiv0, gh%drr)
   dB0dz(3) = Rinv*dot_product(fh%psiv0, gh%drz)
   !B_toroidal = B_Z / R
   Bcyl(2) = Rinv*dot_product(fh%Bzv0, gh%g)
   dB0dR(2) = Rinv*dot_product(fh%Bzv0, gh%dr)
   dB0dz(2) = Rinv*dot_product(fh%Bzv0, gh%dz)

   !if (itor.eq.1) dBdR = dBdR - Rinv*Bcyl
   dB0dR = dB0dR - Rinv*Bcyl

#ifdef USE3D
   !Non-axisymmetric B_poloidal term: - grad f'
   dB0dphi(1) = -Rinv*dot_product(fh%psiv0, gh%dzphi)
   dB0dphi(3) = Rinv*dot_product(fh%psiv0, gh%drphi)
   dB0dphi(2) = Rinv*dot_product(fh%Bzv0, gh%dphi)

   Bcyl(1) = Bcyl(1) - dot_product(fh%Bfpv0, gh%dr)
   dB0dR(1) = dB0dR(1) - dot_product(fh%Bfpv0, gh%drr)
   dB0dphi(1) = dB0dphi(1) - dot_product(fh%Bfpv0, gh%drphi)
   dB0dz(1) = dB0dz(1) - dot_product(fh%Bfpv0, gh%drz)

   Bcyl(3) = Bcyl(3) - dot_product(fh%Bfpv0, gh%dz)
   dB0dR(3) = dB0dR(3) - dot_product(fh%Bfpv0, gh%drz)
   dB0dphi(3) = dB0dphi(3) - dot_product(fh%Bfpv0, gh%dzphi)
   dB0dz(3) = dB0dz(3) - dot_product(fh%Bfpv0, gh%dzz)

   !dBdphi = 0.
#else
   dB0dphi = 0.
#endif

   !Perturbed part
   !B_poloidal = grad psi x grad phi
   temp(1) = -Rinv*dot_product(gh%dz, fh%psiv1)
   tempR(1) = -Rinv*dot_product(gh%drz, fh%psiv1)
   tempz(1) = -Rinv*dot_product(gh%dzz, fh%psiv1)

   temp(3) = Rinv*dot_product(gh%dr, fh%psiv1)
   tempR(3) = Rinv*dot_product(gh%drr, fh%psiv1)
   tempz(3) = Rinv*dot_product(gh%drz, fh%psiv1)

   !B_toroidal = B_Z / R
   temp(2) = Rinv*dot_product(gh%g, fh%Bzv1)
   tempR(2) = Rinv*dot_product(gh%dr, fh%Bzv1)
   tempz(2) = Rinv*dot_product(gh%dz, fh%Bzv1)

   tempR = tempR - Rinv*temp

#ifdef USE3D
   tempphi(1) = -Rinv*dot_product(gh%dzphi, fh%psiv1)
   tempphi(3) = Rinv*dot_product(gh%drphi, fh%psiv1)
   tempphi(2) = Rinv*dot_product(gh%dphi, fh%Bzv1)
#else
   tempphi = temp*rfac_particle
#endif
   !if (itor.eq.1) tempR = tempR - Rinv*temp

#ifdef USECOMPLEX
   temp(1) = temp(1) - dot_product(gh%dr, fh%Bfpv1)
   temp(3) = temp(3) - dot_product(gh%dz, fh%Bfpv1)
   deltaB = real(temp*exp(rfac_particle*x(2)))

   tempR(1) = tempR(1) - dot_product(gh%drr, fh%Bfpv1)
   tempR(3) = tempR(3) - dot_product(gh%drz, fh%Bfpv1)
   dB1dR = real(tempR*exp(rfac_particle*x(2)))
   !dBdR = dBdR + real(tempR * exp(rfac*x(2)))

   tempphi(1) = tempphi(1) - dot_product(gh%dr, fh%Bfpv1)*rfac_particle
   tempphi(3) = tempphi(3) - dot_product(gh%dz, fh%Bfpv1)*rfac_particle
   dB1dphi = real(tempphi*exp(rfac_particle*x(2)))

   tempz(1) = tempz(1) - dot_product(gh%drz, fh%Bfpv1)
   tempz(3) = tempz(3) - dot_product(gh%dzz, fh%Bfpv1)
   dB1dz = real(tempz*exp(rfac_particle*x(2)))
   !dBdz = dBdz + real(tempz * exp(rfac*x(2)))

   !dBdphi = real(temp * rfac * exp(rfac*x(2)))
#else
   temp(1) = temp(1) - dot_product(gh%dr, fh%Bfpv1)
   temp(3) = temp(3) - dot_product(gh%dz, fh%Bfpv1)
   deltaB = temp
   !dBdphi(1) = -Rinv*dot_product(fh%psiv0, gh%dzphi) - dot_product(fh%Bfv, gh%drphiphi)

   tempR(1) = tempR(1) - dot_product(gh%drr, fh%Bfpv1)
   tempR(3) = tempR(3) - dot_product(gh%drz, fh%Bfpv1)
   dB1dR = tempR
   !dBdR = dBdR + tempR

#ifdef USE3D
   tempphi(1) = tempphi(1) - dot_product(fh%Bfpv1, gh%drphi)
   tempphi(3) = tempphi(3) - dot_product(fh%Bfpv1, gh%dzphi)
#endif
   dB1dphi = tempphi

   tempz(1) = tempz(1) - dot_product(gh%drz, fh%Bfpv1)
   tempz(3) = tempz(3) - dot_product(gh%dzz, fh%Bfpv1)
   dB1dz = tempz
   !dBdz = dBdz + tempz
#endif

   !if (.not. (linear_particle == 1)) then
   !   Bcyl = Bcyl + deltaB
   !   dBdR = dBdR + tempR
   !   dBdz = dBdz + tempz
   !   dBdphi = dBdphi + tempphi
   !end if
   Bcyl = Bcyl*b0_norm_particle/1.e4
   deltaB = deltaB*b0_norm_particle/1.e4
   dB0dR = dB0dR*b0_norm_particle/1.e4
   dB0dz = dB0dz*b0_norm_particle/1.e4
   dB0dphi = dB0dphi*b0_norm_particle/1.e4
   dB1dR = dB1dR*b0_norm_particle/1.e4
   dB1dz = dB1dz*b0_norm_particle/1.e4
   dB1dphi = dB1dphi*b0_norm_particle/1.e4
end subroutine getBcylprime
!---------------------------------------------------------------------------
subroutine getEcyl(x, fh, gh, Ecyl)
!$acc routine seq
   use arrays
   use basic
   use auxiliary_fields
   implicit none

   real, dimension(3), intent(in) :: x       !Position
   type(elfield), intent(in) :: fh
   type(xgeomterms), intent(in) :: gh
   real, dimension(3), intent(out) :: Ecyl

   vectype, dimension(3) :: temp
   vectype :: psitemp

   temp(1) = dot_product(gh%g, fh%er)
   temp(2) = dot_product(gh%g, fh%ephi)
   temp(3) = dot_product(gh%g, fh%ez)
   psitemp = dot_product(gh%g, fh%psiv0)

#ifdef USECOMPLEX
   !Ecyl = real(temp)
   Ecyl = real(temp*exp(rfac_particle*x(2)))
#else
   Ecyl = temp
#endif
   Ecyl = Ecyl*(v0_norm_particle/100.0*b0_norm_particle/1.e4)
   !Ecyl = Ecyl * (v0_norm_particle/100.0)
   !if (real(psitemp)<0.21) Ecyl=0.
end subroutine getEcyl

!---------------------------------------------------------------------------
subroutine evalf0(x, vpar, vperp, fh, gh, sps, B0, f0, df0dpsi, df0de, df0dxi)
!$acc routine seq
   use basic
   implicit none

   real, dimension(3), intent(in)       :: x
   real, intent(in)                     :: vpar, vperp
   type(elfield), intent(in)            :: fh
   type(xgeomterms), intent(in)         :: gh
   integer, intent(in)                  :: sps
   real, intent(in)                     :: B0
   real, intent(out)                    :: df0dpsi, df0de, df0dxi
   real, dimension(3)                   :: gradf, gradT, gradrho
   real, intent(out)      :: f0
   real :: T0, gradcoef2
   real :: lambda, gradcoef

   !Math, physics parameters for Maxwellian distribution
   real, parameter :: twopi = 6.283185307179586476925286766559
   real, parameter :: vthermal = 100.0*vp1eV  !Maxwellian temperature = 10 keV
   !real, parameter :: nrmfac = (vthermal*sqrt(twopi))**(-3)
   real, parameter :: ecoef = -0.5*(vthermal**(-2))

   !Parameters for slowing-down distribution
   real, parameter :: vDTalpha = sqrt(3.5e+6/A_alpha)*vp1eV  !Birth speed of DT alpha
   real, parameter :: vbeam = sqrt(1.0e+5)*vp1eV !Injection speed of 100 kV beam proton
   real :: vmax
   !real, parameter :: vmax = 981481 !or vbeam
   real :: vcrit
   !real :: sdc = 1.5/(twopi * log((vmax/vcrit)**3 + 1.0))

   real s, spd, spsq !|v|, v-squared
   real Ipa_d3v, Ipe_d3v
   real psi0, g0, Pphi
   real xi, radi, pitch, energy
   integer radi_i, pitch_i, energy_i
   real f1, f2, f3, f4, f5, f6
   real df0dr, Rinv

   !vmax = (v0_norm_particle/100.0)*4.0
   vmax = 1.e7
   !kim!vmax = 1.e6
   !vcrit = 0.58 * vmax
   vcrit = 0.5*vmax
   ! if (vspdims .eq. 5) then
      !spsq = dot_product(v, v)
      ! spsq = v(4)*v(4) + 2.0*qm_ion(2)*v(5)*modB
      ! xi = v(4)/sqrt(spsq)
   ! else
      spsq = vpar**2 + vperp**2
      xi = vpar/sqrt(spsq)
   ! end if

   select case (1)
   case (0) !Spatially uniform Maxwellian
      !f0 = nrmfac*exp(ecoef*spsq)
      !gradpsi = 0.0
      !gradcoef = 0.
      !df0de = (2.0*ecoef/m_ion)*f0
      !if (.true.) then
      !   Ipa_d3v = m_ion*(vthermal**2)
      !   Ipe_d3v = Ipa_d3v
      !end if
   case (1) !Slowing-down distribution (Fu 2006 formulation)
      spd = sqrt(spsq)
      if (.true.) then
         f0 = 1./(spd**3 + vcrit**3)
         !kim!df0de = (-3.0 * spd * f0) / (m_ion)
         f0 = vcrit**3*f0
         !hou!f0 = 1./(spd**3 + vcrit**3)
         df0de = -1./(400*1e3*1.6e-19)
      else
         f0 = 0.0
         df0de = 0.0
      end if
      !gradpsi = 0.0
      !gradpsi(1) = dot_product(fh%psiv0, gh%dr)
      !gradpsi(3) = dot_product(fh%psiv0, gh%dz)
      !fishbone!gradpsi(1) = gradpsi(1)+v(1)/qm_ion/modB*dot_product(fh%Bzv0, gh%dr)
      !fishbone!gradpsi(3) = gradpsi(3)+v(1)/qm_ion/modB*dot_product(fh%Bzv0, gh%dz)
      !Rinv = 1.0/dot_product(fh%rst,gh%g)
      ! psi0 = dot_product(fh%psiv0, gh%g)
      ! g0 = dot_product(fh%Bzv0, gh%g)
      !gradf0 = gradf0/(0.0555/4.)
      !Pphi = psi0 + (1./qm_ion) * v(1) * Bphi * x(1) / modB / (b0_norm_particle/1.e4)
      ! Pphi = psi0 + (1./qm_ion(2))*v(1)*g0/modB
      !hou!s=sqrt(abs(psi0/0.28036+1))
      !s=sqrt(abs(psi0/0.08+1))
      !hou!gradcoef2=-1/0.3/cosh((s-0.5)/0.2)**2/(2*s)/0.28036
      !gradcoef2=-1/0.3/cosh((s-0.5)/0.2)**2/(2*s)/0.08
      !f0=f0*exp((Pphi-0.0555)/(0.0555/4.))
      !hou!f0=exp(-m_ion*spd**2/(2*400*1e3*1.6e-19))*exp(-0.2/0.3*tanh((s-0.5)/(0.2)))
      !f0=exp(-m_ion*spd**2/(2*17*1e3*1.6e-19))*exp(-0.2/0.3*tanh((s-0.5)/(0.2)))
      !f0=exp(-0.2/0.3*tanh((s-0.5)/(0.2)))/exp(-0.2/0.3*tanh((-0.5/0.2)))
      !psi0=1-psi0/psi_axis
#ifdef USEST
      Rinv = 1.0/dot_product(fh%rst,gh%g)
#endif
 
! #ifndef USEST
      !call evaluate_spline_particle(nf_spline, psi0, f0, gradcoef)
      if (sps==1) then
         f0 = dot_product(fh%nfi, gh%g)
         gradf(1) = dot_product(fh%nfi, gh%dr)
         gradf(3) = dot_product(fh%nfi, gh%dz)
#ifdef USEST
         gradf(2) = Rinv*dot_product(fh%nfi, gh%dphi)
#else
         gradf(2) = 0.
#endif
      else
         f0 = dot_product(fh%nf, gh%g)
         gradf(1) = dot_product(fh%nf, gh%dr)
         gradf(3) = dot_product(fh%nf, gh%dz)
#ifdef USEST
         gradf(2) = Rinv*dot_product(fh%nf, gh%dphi)
#else
         gradf(2) = 0.
#endif
      endif
      gradf = gradf/f0
      !if (abs(xi)<0.02) df0dxi=0.
      ! if (abs(xi)>0.01) then
      ! df0de = df0de+0.5*(1-xi**2)/xi/(m_ion*0.5*spd**2)*df0dxi
      !df0de = df0de-gradcoef/(spd*xi)*dot_product(real(fh%Bzv0), gh%g)/modB/(q_ion)
      ! endif
      !call evaluate_spline_particle(tf_spline, psi0, T0, gradT)
      if (sps==1) then
         T0 = dot_product(fh%tfi, gh%g)
         gradT(1) = dot_product(fh%tfi, gh%dr)
         gradT(3) = dot_product(fh%tfi, gh%dz)
#ifdef USEST
         gradT(2) = Rinv*dot_product(fh%tfi, gh%dphi)
#else
         gradT(2) = 0.
#endif
      else
         T0 = dot_product(fh%tf, gh%g)
         gradT(1) = dot_product(fh%tf, gh%dr)
         gradT(3) = dot_product(fh%tf, gh%dz)
#ifdef USEST
         gradT(2) = Rinv*dot_product(fh%tf, gh%dphi)
#else
         gradT(2) = 0.
#endif
      endif
      
      if ((sps==1).or.(fast_ion_dist_particle.eq.1)) then
         !maxwellian
         gradf=gradf+(-1.5+m_ion(sps)*spd**2/(2*T0*1.6e-19))*gradT/T0
    !!if (psi0>0.64) gradcoef=0.
    !!if (psi0<0.1) gradcoef=0.
         f0=f0*T0**(-1.5)*exp(-m_ion(sps)*spd**2/(2*T0*1.6e-19))
         df0de = -1./(T0*1.6e-19)
         df0dxi = 0.
      elseif (fast_ion_dist_particle.eq.2) then
         !slowingdown
         gradf=gradf-3./(spd**3 + (2*T0*1.6e-19/m_ion(sps))**(1.5))*(2*T0*1.6e-19/m_ion(sps))**(0.5)/m_ion(2)*1.6e-19*gradT
         f0=f0/(spd**3+(T0*2*1.6e-19/m_ion(sps))**(1.5))
         df0de = -3./(spd**3 + (2*T0*1.6e-19/m_ion(sps))**(1.5))*spd/m_ion(2)
         df0dxi = 0.
         !df0dxi=-2*(xi+0.632)/0.55**2
         ! if (xi>-0.02) df0dxi=0.

      !f0=f0/nf_spline%y(1)
      !gradcoef=1.
      !T0=17
      !f0=f0*exp((psi0-0.0555)/(0.0555/4.))
      !gradf0 = gradf0/(1+4.*(1-f0temp/0.0607)**2)*8*(1-f0temp/0.0607)
      elseif (fast_ion_dist_particle.eq.0) then
       pphi = (dot_product(fh%psiv0, gh%g)+(1./qm_ion(sps)) * vpar *  dot_product(fh%Bzv0, gh%g)/B0)-psimin
       if (pphi<r_array(1)) pphi=r_array(1)
       if (pphi>r_array(num_r)) pphi=r_array(num_r)
       radi = pphi
       radi_i=int((radi-r_array(1))/(r_array(2)-r_array(1)))+1
       if (radi_i>=num_r) radi_i=num_r-1
       lambda = vperp**2/B0*1.99/spsq
       if (lambda<pitch_array(1)) lambda=pitch_array(1)
       if (lambda>pitch_array(num_pitch)) lambda=pitch_array(num_pitch)
       pitch = lambda
       pitch_i=int((lambda-pitch_array(1))/(pitch_array(2)-pitch_array(1)))+1
       if (pitch_i>=num_pitch) pitch_i=num_pitch-1
       energy = m_ion(sps)*spd**2/2./1.6e-19
       if (energy<energy_array(1)) energy=energy_array(1)
       if (energy>energy_array(num_energy)) energy=energy_array(num_energy)
       energy_i=int((energy-energy_array(1))/(energy_array(2)-energy_array(1)))+1
       if (energy_i>=num_energy) energy_i=num_energy-1
       if (vpar<0) then
       f1=f_array(energy_i,pitch_i,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f1=f1+f_array(energy_i+1,pitch_i,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f2=f_array(energy_i,pitch_i+1,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f2=f2+f_array(energy_i+1,pitch_i+1,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f3=f1*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f3=f3+f2*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f4=f_array(energy_i,pitch_i,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f4=f4+f_array(energy_i+1,pitch_i,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f5=f_array(energy_i,pitch_i+1,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f5=f5+f_array(energy_i+1,pitch_i+1,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f6=f4*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f6=f6+f5*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f0=f3*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f0=f0+f6*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       !if ((f6/f3>10.).or.(f6/f3<0.1)) then
          df0dr=(f6-f3)/(r_array(radi_i+1)-r_array(radi_i))/(f6+f3+1e-10)*2
       !else
          !df0dr=(f6-f3)/(r_array(radi_i+1)-r_array(radi_i))/(f0+1e-10)
       !endif
       if (abs(df0dr)>0.2/(r_array(2)-r_array(1))) df0dr=0
       !df0dr=df0dr*0.003
       gradf = 0.
       gradf(1) = dot_product(fh%rho, gh%dr)
       gradf(3) = dot_product(fh%rho, gh%dz)
       !gradcoef=dot_product(gradf,gradpsi)/dot_product(gradpsi,gradpsi)
       !gradcoef=-gradcoef*df0dr
       gradf=gradf*df0dr
       f1=f_array(energy_i,pitch_i,radi_i)*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f1=f1+f_array(energy_i,pitch_i+1,radi_i)*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f2=f_array(energy_i,pitch_i,radi_i+1)*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f2=f2+f_array(energy_i,pitch_i+1,radi_i+1)*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f3=f1*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f3=f3+f2*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       f4=f_array(energy_i+1,pitch_i,radi_i)*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f4=f4+f_array(energy_i+1,pitch_i+1,radi_i)*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f5=f_array(energy_i+1,pitch_i,radi_i+1)*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f5=f5+f_array(energy_i+1,pitch_i+1,radi_i+1)*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f6=f4*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f6=f6+f5*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       !if ((f6/f3>10.).or.(f6/f3<0.1)) then
          df0de=(f6-f3)/(energy_array(energy_i+1)-energy_array(energy_i))/(f6+f3+1e-10)*2
       !else
          !df0de=(f6-f3)/(energy_array(energy_i+1)-energy_array(energy_i))/(f0+1e-10)
       !endif
       !if ((energy_i==13) .or. (energy_i==20) .or. (energy_i==42)) then
       !if (energy>174000) then
       !   df0de=0.
       !endif
       if (abs(df0de)>1./(energy_array(2)-energy_array(1))) df0de=0
       !df0de=0.
       df0de=df0de/1.6e-19
       f1=f_array(energy_i,pitch_i,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f1=f1+f_array(energy_i+1,pitch_i,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f2=f_array(energy_i,pitch_i,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f2=f2+f_array(energy_i+1,pitch_i,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f3=f1*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f3=f3+f2*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       f4=f_array(energy_i,pitch_i+1,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f4=f4+f_array(energy_i+1,pitch_i+1,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f5=f_array(energy_i,pitch_i+1,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f5=f5+f_array(energy_i+1,pitch_i+1,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f6=f4*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f6=f6+f5*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       !if ((f6/f3>10.).or.(f6/f3<0.1)) then
          df0dxi=(f6-f3)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))/(f6+f3+1e-10)*2
       !else
          !df0dxi=(f6-f3)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))/(f0+1e-10)
       !endif
       if (abs(df0dxi)>1./(pitch_array(2)-pitch_array(1))) df0dxi=0
       !df0dxi=0
       !df0dxi=exp(-((xi+0.5)/0.3)**2)*(-2)*(xi+0.5)/0.3/0.3
    else
       f1=f_array2(energy_i,pitch_i,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f1=f1+f_array2(energy_i+1,pitch_i,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f2=f_array2(energy_i,pitch_i+1,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f2=f2+f_array2(energy_i+1,pitch_i+1,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f3=f1*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f3=f3+f2*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f4=f_array2(energy_i,pitch_i,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f4=f4+f_array2(energy_i+1,pitch_i,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f5=f_array2(energy_i,pitch_i+1,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f5=f5+f_array2(energy_i+1,pitch_i+1,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f6=f4*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f6=f6+f5*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f0=f3*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f0=f0+f6*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       !if ((f6/f3>10.).or.(f6/f3<0.1)) then
          df0dr=(f6-f3)/(r_array(radi_i+1)-r_array(radi_i))/(f6+f3+1e-10)*2
       !else
          !df0dr=(f6-f3)/(r_array(radi_i+1)-r_array(radi_i))/(f0+1e-10)
       !endif
       if (abs(df0dr)>0.2/(r_array(2)-r_array(1))) df0dr=0
       !df0dr=df0dr*0.003
       gradf = 0.
       gradf(1) = dot_product(fh%rho, gh%dr)
       gradf(3) = dot_product(fh%rho, gh%dz)
       !gradcoef=dot_product(gradf,gradpsi)/dot_product(gradpsi,gradpsi)
       !gradcoef=-gradcoef*df0dr
       gradf=gradf*df0dr
       f1=f_array2(energy_i,pitch_i,radi_i)*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f1=f1+f_array2(energy_i,pitch_i+1,radi_i)*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f2=f_array2(energy_i,pitch_i,radi_i+1)*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f2=f2+f_array2(energy_i,pitch_i+1,radi_i+1)*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f3=f1*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f3=f3+f2*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       f4=f_array2(energy_i+1,pitch_i,radi_i)*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f4=f4+f_array2(energy_i+1,pitch_i+1,radi_i)*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f5=f_array2(energy_i+1,pitch_i,radi_i+1)*(pitch_array(pitch_i+1)-pitch)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f5=f5+f_array2(energy_i+1,pitch_i+1,radi_i+1)*(pitch-pitch_array(pitch_i))/(pitch_array(pitch_i+1)-pitch_array(pitch_i))
       f6=f4*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f6=f6+f5*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       !if ((f6/f3>10.).or.(f6/f3<0.1)) then
          df0de=(f6-f3)/(energy_array(energy_i+1)-energy_array(energy_i))/(f6+f3+1e-10)*2
       !else
          !df0de=(f6-f3)/(energy_array(energy_i+1)-energy_array(energy_i))/(f0+1e-10)
       !endif
       !if ((energy_i==13) .or. (energy_i==20) .or. (energy_i==42)) then
       !if (energy>174000) then
       !   df0de=0.
       !endif
       if (abs(df0de)>1./(energy_array(2)-energy_array(1))) df0de=0
       !df0de=0.
       df0de=df0de/1.6e-19
       f1=f_array2(energy_i,pitch_i,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f1=f1+f_array2(energy_i+1,pitch_i,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f2=f_array2(energy_i,pitch_i,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f2=f2+f_array2(energy_i+1,pitch_i,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f3=f1*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f3=f3+f2*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       f4=f_array2(energy_i,pitch_i+1,radi_i)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f4=f4+f_array2(energy_i+1,pitch_i+1,radi_i)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f5=f_array2(energy_i,pitch_i+1,radi_i+1)*(energy_array(energy_i+1)-energy)/(energy_array(energy_i+1)-energy_array(energy_i))
       f5=f5+f_array2(energy_i+1,pitch_i+1,radi_i+1)*(energy-energy_array(energy_i))/(energy_array(energy_i+1)-energy_array(energy_i))
       f6=f4*(r_array(radi_i+1)-radi)/(r_array(radi_i+1)-r_array(radi_i))
       f6=f6+f5*(radi-r_array(radi_i))/(r_array(radi_i+1)-r_array(radi_i))
       !if ((f6/f3>10.).or.(f6/f3<0.1)) then
          df0dxi=(f6-f3)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))/(f6+f3+1e-10)*2
       !else
          !df0dxi=(f6-f3)/(pitch_array(pitch_i+1)-pitch_array(pitch_i))/(f0+1e-10)
       !endif
       if (abs(df0dxi)>1./(pitch_array(2)-pitch_array(1))) df0dxi=0
       !df0dxi=0
       !df0dxi=exp(-((xi+0.5)/0.3)**2)*(-2)*(xi+0.5)/0.3/0.3
    endif
    endif
      gradrho(1) = dot_product(fh%rho,gh%dr)
      gradrho(3) = dot_product(fh%rho,gh%dz)
#ifdef USEST
      gradrho(2) = Rinv*dot_product(fh%rho,gh%dphi)
#else
      gradrho(2) = 0.
#endif
      gradcoef = dot_product(gradf, gradrho)/dot_product(gradrho, gradrho)
      df0dpsi = df0dr
   case default
      f0 = 1.0
      !gradpsi = 0.
      !gradcoef = 0.0
      df0de = 0.0
   end select
   if (real(dot_product(fh%rho, gh%g))>kinetic_rhomax_particle) then
      gradf=0.
      df0de=0.
      df0dxi=0.
   endif

end subroutine evalf0
!---------------------------------------------------------------------------
subroutine particle_pressure_rhs
   use basic
   use arrays
   use math
   use m3dc1_nint
   implicit none
   include 'mpif.h'
   intrinsic matmul

   real, dimension(dofs_per_element, coeffs_per_element) :: cl
   real, dimension(coeffs_per_element) :: wnuhere, wnuhere2, deltaBhere
#ifdef USECOMPLEX
   complex phfac
#endif
   real, dimension(3) :: B_part, deltaB
   real, dimension(vspdims) :: vperp
   !type(elfield), dimension(nneighbors+1) :: elcoefs
   type(xgeomterms) :: geomterms, geomterms2
   real             :: B0, vpar, ppar, pperp
   integer          :: i, ierr, ielm, ielm_local, ielm_global, ipart, itri, itri2, tridex, isghost
   !integer          :: ibp, iwe, iok
   type(element_data) :: eldat
   integer :: ielm2
   integer, dimension(nodes_per_element) :: nodeids, nodeids2
   real :: x1, phi1, z1, x2, phi2, z2, x3, phi3, z3, area
   real :: x21, phi21, z21, x22, phi22, z22, x23, phi23, z23
   integer :: k
   real :: f0, df0de
   real, dimension(3) :: gradf0
   integer :: ipoint
   integer :: ipart_begin_local, ipart_end_local
   real, dimension(3) :: xtemp
   real, dimension(vspdims) :: vtemp
   real, dimension(3) :: B_cyl, gradB0, gradB1
   real :: dB1
   !nelms = size(pdata)
   !nelms = local_elements()
   !elcoefs(:)%itri = 0

   coeffspaf_local = 0.; coeffspef_local = 0.
   coeffspai_local = 0.; coeffspei_local = 0.
   coeffsdef0_local = 0.; coeffsdef0_local = 0.
   coeffsdei0_local = 0.; coeffsdei0_local = 0.
   coeffsdef_local = 0.; coeffsdef_local = 0.
   coeffsdei_local = 0.; coeffsdei_local = 0.
   coeffsvpi_local = 0.

   ipart_begin_local = (ipart_end - ipart_begin + 1)/ncols*hostrank + ipart_begin
   ipart_end_local = (ipart_end - ipart_begin + 1)/ncols*(hostrank + 1) - 1 + ipart_begin
   if (hostrank == ncols - 1) ipart_end_local = ipart_end
   do ipart = ipart_begin_local, ipart_end_local
      !if (pdata(ipart)%deleted) cycle
      do ipoint = 1, 4
         if (igyroaverage.eq.1) then
            itri = pdata(ipart)%kel(ipoint)
         else
            itri = pdata(ipart)%jel
         end if
         if (igyroaverage.eq.1) then
            call get_geom_terms(pdata(ipart)%kx(:, ipoint), itri, &
                                geomterms, .false., ierr)
         else
            xtemp = pdata(ipart)%x
            vtemp = pdata(ipart)%v
            if (vspdims .eq. 5) call advancex(xtemp, vtemp, -t0_norm*dt/particle_substeps/2)
            call get_geom_terms(xtemp, itri, &
                                geomterms, .false., ierr)
         end if
#ifdef USEST
            call update_geom_terms_st(geomterms, elfieldcoefs(itri), .false.)
#endif

         if (ierr .ne. 0) then
            print *, myrank, ': Bad particle in pressure tensor integral; skipping.'
            write (0, *) itri, pdata(ipart)%x
            !ibp = ibp + 1
            cycle !next particle
         end if
         if (itri .ne. pdata(ipart)%kel(ipoint)) then
            print *, myrank, ': Particle in wrong element in pressure tensor integral:', pdata(ipart)%gid
            print *, ipart, ': ', itri, '.ne.', pdata(ipart)%kel(ipoint)
       !!print *,myrank,': pcoord ',pdata(ielm)%ion(ipart)%gid,' = ',pdata(ielm)%ion(ipart)%x
       !!iwe = iwe + 1
            !cycle !next particle
         end if
         if (pdata(ipart)%B0.ne.0) then
            B0 = pdata(ipart)%B0
         else
            call getBcyl(pdata(ipart)%x, elfieldcoefs(itri), geomterms, B_cyl, deltaB, gradB0, gradB1, dB1)
            !itri2 = itri !fluid particle
            !call get_geom_terms(pdata(ipart)%x0, itri2, geomterms2, .false., ierr) !fluid particle
            !call getBcyl(pdata(ipart)%x0, elfieldcoefs(itri2), geomterms2, B_cyl, deltaB, gradB0, gradB1, dB1) !fluid particle
            B0 = sqrt(dot_product(B_cyl, B_cyl))  !1/magnitude of B
            pdata(ipart)%B0 = B0
         endif
         !Use B and v to get parallel and perp components of particle velocity
         if (vspdims .eq. 2) then ! drift-kinetic: v_|| = v(1),  mu = q * v(2)
            vpar = pdata(ipart)%v(1)
            pperp = q_ion(pdata(ipart)%sps)*pdata(ipart)%v(2)*B0
         else !full orbit: v_|| = v.B/|B|,  v_perp = v - v_||
            if (B0 .gt. 0.0) then !non-degenerate
               !vpar = dot_product(pdata(ipart)%v, B_part(1:vspdims))/B0
               vpar = pdata(ipart)%v(4)
               !vperp = pdata(ipart)%v - (vpar/B0)*B_part(1:vspdims)
               pperp = q_ion(pdata(ipart)%sps)*pdata(ipart)%v(5)*B0
            else !degenerate case: no B field, pressure is scalar
               !vpar = 0.0
               !vperp = pdata(ipart)%v
            end if !degenerate?
            !pperp = 0.5 * m_ion * dot_product(vperp, vperp)
         end if !full-orbit?
         ppar = m_ion(pdata(ipart)%sps)*vpar**2
         if (particle_linear == 1) then
            wnuhere = (pdata(ipart)%wt + pdata(ipart)%dB)*geomterms%g
            wnuhere2 = (pdata(ipart)%wt + pdata(ipart)%dB + pdata(ipart)%dB)*geomterms%g
         else
            wnuhere = (pdata(ipart)%wt + (1 - pdata(ipart)%wt)*pdata(ipart)%dB)*geomterms%g
            wnuhere2 = (pdata(ipart)%wt + (1 - pdata(ipart)%wt)*pdata(ipart)%dB + pdata(ipart)%dB)*geomterms%g
            !if (abs(pdata(ipart)%wt)>1) then
            !   wnuhere=0
            !   wnuhere2=0
            !endif
         end if
         wnuhere = (pdata(ipart)%wt) * geomterms%g
         wnuhere2 = (pdata(ipart)%wt) * geomterms%g
         ! wnuhere = geomterms%g
         ! wnuhere2 = geomterms%g
         wnuhere = wnuhere*nrmfac(pdata(ipart)%sps)
         wnuhere2 = wnuhere2*nrmfac(pdata(ipart)%sps)
         !if (pdata(ipart)%sps==1) then
         !wnuhere=0.
         !wnuhere2=0.
         !endif
         !wnuhere = pdata(ipart)%wt * matmul(cl, geomterms%g) * pdata(ipart)%x(1)/10.
         !deltaBhere = pdata(ielm)%ion(ipart)%f0 *dot_product(B_part,deltaB)/B0**2* matmul(cl,geomterms%g)
         if (pdata(ipart)%sps == 1) then
            coeffsdei0_local(:,itri) = coeffsdei0_local(:,itri) + geomterms%g*nrmfac(pdata(ipart)%sps)/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)
#ifndef USECOMPLEX
            coeffspai_local(:, itri) = coeffspai_local(:, itri) + ppar*wnuhere/4
            ! coeffspai_local(:,itri) = coeffspai_local(:,itri) + ppar*geomterms%g*nrmfac(pdata(ipart)%sps)/4
            coeffspei_local(:, itri) = coeffspei_local(:, itri) + pperp*wnuhere2/4
            ! coeffspei_local(:,itri) = coeffspei_local(:,itri) + pperp*geomterms%g*nrmfac(pdata(ipart)%sps)/4
            coeffsdei_local(:, itri) = coeffsdei_local(:, itri) + wnuhere/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)
            coeffsvpi_local(:, itri) = coeffsvpi_local(:, itri) + vpar/(v0_norm/100.)*wnuhere/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)
            !dofspa = intx2(mu79(:,:,OP_1),ppar79(:,OP_1))
            !dofspe = intx2(mu79(:,:,OP_1),pper79(:,OP_1))
#else
            !Extract appropriate Fourier component of particle contribution
            if (igyroaverage.eq.1) then
               phfac = exp(-rfac*pdata(ipart)%kx(2, ipoint))
            else
               phfac = exp(-rfac*xtemp(2))
            end if
            coeffspai_local(:, itri) = coeffspai_local(:, itri) + ppar*phfac*wnuhere/4*2
            !coeffspai_local(:,itri) = coeffspai_local(:,itri) + ppar*geomterms%g*nrmfac(pdata(ipart)%sps)/4
            coeffspei_local(:, itri) = coeffspei_local(:, itri) + pperp*phfac*wnuhere2/4*2
            !coeffspei_local(:,itri) = coeffspei_local(:,itri) + pperp*geomterms%g*nrmfac(pdata(ipart)%sps)/4
            coeffsdei_local(:, itri) = coeffsdei_local(:, itri) + phfac*wnuhere/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)*2
            coeffsvpi_local(:,itri) = coeffsvpi_local(:,itri) + vpar/(v0_norm/100.)*phfac*wnuhere/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)*2
            !dofspen = dofspen + pperp*phfac*deltaBhere
            !dofspan = intx2(mu79(:,:,OP_1),ppar79(:,OP_1))
            !dofspen = intx2(mu79(:,:,OP_1),pper79(:,OP_1))
#endif
         else
            coeffsdef0_local(:,itri) = coeffsdef0_local(:,itri) + geomterms%g*nrmfac(pdata(ipart)%sps)/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)
#ifndef USECOMPLEX
            !coeffspaf_local(:,itri) = coeffspaf_local(:,itri) + ppar*geomterms%g*nrmfac(pdata(ipart)%sps)/4
            coeffspaf_local(:, itri) = coeffspaf_local(:, itri) + ppar*wnuhere/4
            !coeffspef_local(:,itri) = coeffspef_local(:,itri) + pperp*geomterms%g*nrmfac(pdata(ipart)%sps)/4
            coeffspef_local(:, itri) = coeffspef_local(:, itri) + pperp*wnuhere2/4
            coeffsdef_local(:, itri) = coeffsdef_local(:, itri) + wnuhere/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)
            coeffsvpf_local(:, itri) = coeffsvpf_local(:, itri) + vpar/(v0_norm/100.)*wnuhere/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)
            !dofspa = intx2(mu79(:,:,OP_1),ppar79(:,OP_1))
            !dofspe = intx2(mu79(:,:,OP_1),pper79(:,OP_1))
#else
            !Extract appropriate Fourier component of particle contribution
            if (igyroaverage.eq.1) then
               phfac = exp(-rfac*pdata(ipart)%kx(2, ipoint))
            else
               phfac = exp(-rfac*xtemp(2))
            end if
            !setup
            coeffspaf_local(:, itri) = coeffspaf_local(:, itri) + ppar*phfac*wnuhere/4*2
            !coeffspaf_local(:,itri) = coeffspaf_local(:,itri) + ppar*geomterms%g*nrmfac(pdata(ipart)%sps)/4
            coeffspef_local(:, itri) = coeffspef_local(:, itri) + pperp*phfac*wnuhere2/4*2
            !coeffspef_local(:,itri) = coeffspef_local(:,itri) + pperp*geomterms%g*nrmfac(pdata(ipart)%sps)/4
            coeffsdef_local(:, itri) = coeffsdef_local(:, itri) + phfac*wnuhere/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)*2
            coeffsvpf_local(:,itri) = coeffsvpf_local(:,itri) + vpar/(v0_norm/100.)*phfac*wnuhere/4&
               *(b0_norm/1e4)**2/(4*3.14159*1e-7)/(n0_norm*1e6)*2
            !dofspen = dofspen + pperp*phfac*deltaBhere
            !dofspan = intx2(mu79(:,:,OP_1),ppar79(:,OP_1))
            !dofspen = intx2(mu79(:,:,OP_1),pper79(:,OP_1))
#endif
         end if
      end do
      !pdata(ipart)%wt=0.
   end do

   coeffspaf = 0.; coeffspef = 0.
   coeffspai = 0.; coeffspei = 0.
   coeffsdef0 = 0.; coeffsdef0 = 0.
   coeffsdei0 = 0.; coeffsdei0 = 0.
   coeffsdef = 0.; coeffsdef = 0.
   coeffsdei = 0.; coeffsdei = 0.
   coeffsvpi = 0.
#ifdef USECOMPLEX
   call mpi_allreduce(coeffspaf_local, coeffspaf, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffspef_local, coeffspef, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsdef_local, coeffsdef, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsdef0_local, coeffsdef0, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
#else
   call mpi_allreduce(coeffspaf_local, coeffspaf, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffspef_local, coeffspef, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsdef_local, coeffsdef, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsdef0_local, coeffsdef0, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
#ifdef USECOMPLEX
   call mpi_allreduce(coeffspai_local, coeffspai, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffspei_local, coeffspei, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsdei_local, coeffsdei, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsdei0_local, coeffsdei0, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsvpi_local, coeffsvpi, coeffs_per_element*nelms_global, MPI_DOUBLE_COMPLEX,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
#else
   call mpi_allreduce(coeffspai_local, coeffspai, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffspei_local, coeffspei, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsdei_local, coeffsdei, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsdei0_local, coeffsdei0, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
   call mpi_allreduce(coeffsvpi_local, coeffsvpi, coeffs_per_element*nelms_global, MPI_DOUBLE_PRECISION,&
      MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

   p_f_par%vec = 0.; p_f_perp%vec = 0.
   p_i_par%vec = 0.; p_i_perp%vec = 0.
   den_i_1%vec = 0.
   den_f_1%vec = 0.
   den_i_0%vec = 0.
   den_f_0%vec = 0.
   v_i_par%vec = 0.
   do ielm = 1, nelms
      if (iprecompute_metric .eq. 1) then
         cl = ctri(:, :, ielm)
      else
         call local_coeff_vector(ielm, cl)
      end if
#ifdef USE3D
      call m3dc1_ent_getglobalid(3, ielm - 1, ielm_global)
#else
      call m3dc1_ent_getglobalid(2, ielm - 1, ielm_global)
#endif
      ielm_global = ielm_global + 1
      call vector_insert_block(p_f_par%vec, ielm, 1, matmul(cl, coeffspaf(:, ielm_global)), VEC_ADD)
      call vector_insert_block(p_f_perp%vec, ielm, 1, matmul(cl, coeffspef(:, ielm_global)), VEC_ADD)
      call vector_insert_block(den_f_1%vec, ielm, 1, matmul(cl, coeffsdef(:, ielm_global)), VEC_ADD)
      call vector_insert_block(den_f_0%vec, ielm, 1, matmul(cl, coeffsdef0(:, ielm_global)), VEC_ADD)
      call vector_insert_block(den_i_1%vec, ielm, 1, matmul(cl, coeffsdei(:, ielm_global)), VEC_ADD)
      call vector_insert_block(p_i_par%vec, ielm, 1, matmul(cl, coeffspai(:, ielm_global)), VEC_ADD)
      call vector_insert_block(p_i_perp%vec, ielm, 1, matmul(cl, coeffspei(:, ielm_global)), VEC_ADD)
      call vector_insert_block(den_i_0%vec, ielm, 1, matmul(cl, coeffsdei0(:, ielm_global)), VEC_ADD)
      call vector_insert_block(v_i_par%vec, ielm, 1, matmul(cl, coeffsvpi(:,ielm_global)), VEC_ADD)
   end do

end subroutine particle_pressure_rhs
!---------------------------------------------------------------------------
subroutine solve_pi_tensor
   use basic
   use newvar_mod
   use arrays
   use matrix_mod
   implicit none
   integer :: ierr

   !call newvar_solve(p_f_par%vec,  diff_mat)
   !call newvar_solve(p_f_perp%vec,  diff_mat)
   !call newvar_solve(p_i_par%vec,  diff_mat)
   !call newvar_solve(p_i_perp%vec,  diff_mat)
   !call newvar_solve(den_f_1%vec,  diff_mat)
   !call newvar_solve(den_f_0%vec,  diff_mat)
   call sum_shared(p_f_par%vec)
   call newsolve(diff2_mat, p_f_par%vec, ierr)
   call sum_shared(p_f_perp%vec)
   call newsolve(diff2_mat, p_f_perp%vec, ierr)
   call sum_shared(den_f_1%vec)
   call newsolve(diff2_mat, den_f_1%vec, ierr)
   call sum_shared(den_f_0%vec)
   call newsolve(diff2_mat, den_f_0%vec, ierr)

   call sum_shared(den_i_1%vec)
   call newsolve(diff2_mat, den_i_1%vec, ierr)
   call sum_shared(p_i_par%vec)
   call newsolve(diff2_mat, p_i_par%vec, ierr)
   call sum_shared(p_i_perp%vec)
   call newsolve(diff2_mat, p_i_perp%vec, ierr)
   !call sum_shared(v_i_par%vec)
   !call newsolve(diff3_mat, v_i_par%vec, ierr)
   call sum_shared(den_i_0%vec)
   call newsolve(diff2_mat, den_i_0%vec, ierr)
   call sum_shared(v_i_par%vec)
   call newsolve(diff2_mat, v_i_par%vec, ierr)
end subroutine solve_pi_tensor
!---------------------------------------------------------------------------
! Dump particle data for current timeslice using parallel HDF5.
subroutine hdf5_write_particles(ierr)
   use basic
   use hdf5_output
   implicit none

   include 'mpif.h'

   integer, intent(out) :: ierr

   character(LEN=32) :: part_file_name
   real, dimension(:, :), allocatable :: values
   integer, parameter :: pdims = vspdims + 12
   integer(HID_T) :: plist_id, part_file_id, part_root_id, group_id
   integer(HID_T) :: filespace, memspace, dset_id
   integer(HSIZE_T), dimension(2) :: local_dims, global_dims
   integer(HSSIZE_T), dimension(2) :: off_h5
   integer :: info
   integer ielm, poffset, pindex, np, ipart, size
   real, dimension(3) :: xtemp
   real, dimension(vspdims) :: vtemp
   real :: spsq
   type(xgeomterms)   :: geomterms
   integer :: itri

   !Calculate data size
   call MPI_Bcast(nparticles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   !Calculate offset of current process
   if (hostrank == 0) then
      call delete_particle(.true.)
   end if
   call MPI_Bcast(ipart_begin, 1, MPI_INTEGER, 0, hostcomm, ierr)
   call MPI_Bcast(ipart_end, 1, MPI_INTEGER, 0, hostcomm, ierr)
   call MPI_Bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
   call mpi_barrier(mpi_comm_world, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   np = int(nparticles/size) + 1
   !call mpi_scan(locparts, poffset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   !poffset = poffset - locparts
   poffset = np*myrank
   if (poffset + np > nparticles) np = nparticles - poffset

   write (part_file_name, '("ions_",I4.4,".h5")') times_output-1

   !Allocate buffer for element particle data
   allocate (values(pdims, np))

   !Create new file for timeslice
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !Set up the file access property list with parallel I/O
   call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
   info = MPI_INFO_NULL
   call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, ierr)

   !Open the new file
   call h5fcreate_f(part_file_name, H5F_ACC_TRUNC_F, part_file_id, ierr, &
                    access_prp=plist_id)
   if (ierr .lt. 0) then
      if (myrank .eq. 0) &
         print *, "Error: could not open ", part_file_name, &
         " for HDF5 output.  error = ", ierr
      return
   end if

   !Open the root group
   call h5gopen_f(part_file_id, "/", part_root_id, ierr)

   if (myrank .eq. 0 .and. iprint .ge. 1) &
      print *, ' Writing particle time slice file ', part_file_name

   !Write attributes
   call write_real_attr(part_root_id, "time", time, ierr)
   call write_int_attr(part_root_id, "velocity space dims", vspdims, ierr)
   call write_vec_attr(part_root_id, "particle nrmfac", nrmfac/kinetic_nrmfac_scale, 2, ierr)
   !call write_real_attr(part_root_id, "particle delta-t", dt_ion, ierr)

   !Create global dataset
   call h5gcreate_f(part_root_id, "particles", group_id, ierr)
   global_dims(1) = pdims; local_dims(1) = pdims
   global_dims(2) = nparticles
   call h5screate_simple_f(2, global_dims, filespace, ierr)
   if (ierr .ne. 0) then
      print *, myrank, ': error', ierr, ' after h5screate_simple_f'
      call safestop(101)
   end if
   if (idouble_out .eq. 1) then
      call h5dcreate_f(group_id, "data", H5T_NATIVE_DOUBLE, filespace, dset_id, ierr)
   else
      call h5dcreate_f(group_id, "data", H5T_NATIVE_REAL, filespace, dset_id, ierr)
   end if
   if (ierr .ne. 0) then
      print *, myrank, ': error', ierr, ' after h5dcreate_f'
      call safestop(101)
   end if
   call h5sclose_f(filespace, ierr)

   !Add labels, units
   !call write_real_attr(group_id, "atomic number", q_ion/e_mks, ierr)
   !call write_real_attr(group_id, "atomic mass", m_ion/m_proton, ierr)
   !call write_str_attr(group_id, "Col. 1 label", "Global ID", ierr)
   !call write_str_attr(group_id, "Col. 2 label", "R", ierr)
   !call write_str_attr(group_id, "Col. 2 units", "m", ierr)
   !call write_str_attr(group_id, "Col. 3 label", "phi", ierr)
   !call write_str_attr(group_id, "Col. 3 units", "radians", ierr)
   !call write_str_attr(group_id, "Col. 4 label", "z", ierr)
   !call write_str_attr(group_id, "Col. 4 units", "m", ierr)
   !call write_str_attr(group_id, "Col. 5 label", "weight", ierr)
   !if (vspdims.eq.2) then
   !call write_str_attr(group_id, "Col. 6 label", "v_parallel", ierr)
   !call write_str_attr(group_id, "Col. 6 units", "m/s", ierr)
   !call write_str_attr(group_id, "Col. 7 label", "mu/q", ierr)
   !call write_str_attr(group_id, "Col. 7 units", "m**2/s", ierr)
   !else
   !call write_str_attr(group_id, "Col. 6 label", "v_R", ierr)
   !call write_str_attr(group_id, "Col. 6 units", "m/s", ierr)
   !call write_str_attr(group_id, "Col. 7 label", "v_phi", ierr)
   !call write_str_attr(group_id, "Col. 7 units", "m/s", ierr)
   !call write_str_attr(group_id, "Col. 8 label", "v_z", ierr)
   !call write_str_attr(group_id, "Col. 8 units", "m/s", ierr)
   !endif

   !Output the particle data
   off_h5(1) = 0
   !nelms = size(pdata)
   !do ipart=1,nparticles
   !np = pdata(ielm)%np
   if (np .gt. 0) then

      !Select local hyperslab within dataset
      local_dims(2) = np
      call h5screate_simple_f(2, local_dims, memspace, ierr)
      if (ierr .ne. 0) then
         print *, myrank, ': error', ierr, ' after h5screate_simple_f'
         call safestop(102)
      end if
      call h5dget_space_f(dset_id, filespace, ierr)
      if (ierr .ne. 0) then
         print *, myrank, ': error', ierr, ' after h5dget_space_f'
         call safestop(102)
      end if
      off_h5(2) = poffset
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off_h5, &
                                 local_dims, ierr)
      if (ierr .ne. 0) then
         print *, myrank, ': error', ierr, ' after h5sselect_hyperslab_f'
         call safestop(102)
      end if

      !Copy data to buffer
      do ipart = 1, np
         pindex = poffset + ipart
         values(1, ipart) = pdata(pindex)%gid
         values(2, ipart) = pdata(pindex)%x(1)
         values(3, ipart) = pdata(pindex)%x(2)
         values(4, ipart) = pdata(pindex)%x(3)
         values(5, ipart) = pdata(pindex)%wt
         values(6, ipart) = pdata(pindex)%sps
         values(7, ipart) = pdata(pindex)%v(1)
         values(8, ipart) = pdata(pindex)%v(2)
         if (vspdims == 5) then
            values(9, ipart) = pdata(pindex)%v(3)
            values(10, ipart) = pdata(pindex)%v(4)
         end if
         ! values(pdims, ipart) = pdata(pindex)%v(vspdims)
         values(9, ipart) = pdata(pindex)%f0
         values(10, ipart) = pdata(pindex)%x0(1)
         values(11, ipart) = pdata(pindex)%x0(2)
         values(12, ipart) = pdata(pindex)%x0(3)
         values(13, ipart) = pdata(pindex)%v0(1)
         values(14, ipart) = pdata(pindex)%v0(2)
         xtemp = pdata(pindex)%x
         vtemp = pdata(pindex)%v
         itri = pdata(pindex)%jel
         call get_geom_terms(xtemp, itri, geomterms, .false., ierr)
#ifdef USEST
         call update_geom_terms_st(geomterms, elfieldcoefs(itri), .false.)
#endif
         spsq = vtemp(1)**2 + 2.0*qm_ion(pdata(pindex)%sps)*vtemp(2)*pdata(pindex)%B0
         !values(11, ipart) = dot_product(elfieldcoefs(itri)%rho, geomterms%g)
         !values(12, ipart) = vtemp(1)/sqrt(spsq)
         !values(13, ipart) = m_ion(pdata(pindex)%sps)*0.5*spsq
         !values(10, ipart) = vtemp(1)/sqrt(spsq)
         values(10, ipart) = 2.0*qm_ion(pdata(pindex)%sps)*vtemp(2)*1.99/spsq
         values(11, ipart) = m_ion(pdata(pindex)%sps)*0.5*spsq/1.6e-19
         !values(12, ipart) = 1-(dot_product(elfieldcoefs(itri)%psiv0, geomterms%g)-0.19)/0.2569 + 0*(1./qm_ion(pdata(pindex)%sps)) * vtemp(1) *  dot_product(elfieldcoefs(itri)%Bzv0, geomterms%g)/pdata(pindex)%B0
         values(12, ipart) = (dot_product(elfieldcoefs(itri)%psiv0, geomterms%g)+(1./qm_ion(pdata(pindex)%sps)) * vtemp(1) *  dot_product(elfieldcoefs(itri)%Bzv0, geomterms%g)/pdata(pindex)%B0-psibound)/(psimin-psibound)-1
     end do !ipart

      !Write the dataset
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, values, global_dims, ierr, &
                      file_space_id=filespace, mem_space_id=memspace)
      if (ierr .ne. 0) then
         print *, myrank, ': error', ierr, ' after h5dwrite_f'
         call safestop(103)
      end if

      call h5sclose_f(filespace, ierr)
      call h5sclose_f(memspace, ierr)

      !poffset = poffset + np
   end if
   !enddo !ielm

   !Close the particle dataset and group
   call h5dclose_f(dset_id, ierr)
   call h5gclose_f(group_id, ierr)

   !Close the file
   call h5gclose_f(part_root_id, ierr)
   call h5fclose_f(part_file_id, ierr)
   call h5pclose_f(plist_id, ierr)

   !Free the buffer
   deallocate (values)
end subroutine hdf5_write_particles

!---------------------------------------------------------------------------
!Read stored HDF5 particle data in parallel.
subroutine hdf5_read_particles(filename, ierr)
   use basic
   use diagnostics
   use hdf5_output
   implicit none
   include 'mpif.h'

   character(len=*), intent(in) :: filename
   integer, intent(out) :: ierr

   integer, parameter :: chunksize = 1024 ! # of particles to read in one pass
   integer, parameter :: ldim = vspdims + 12 ! Leading dimension of particle data array

   type(particle) :: dpar
   real, dimension(:, :), allocatable :: valbuf
   real :: tstart, tend
   real, dimension(3) :: mmsa
#ifdef USE3D
   real xi, zi
#endif
   integer(HID_T) :: fileid, access_props, part_root_id, group_id
   integer(HID_T) :: dset_id, filespace, memspace
   integer(HSIZE_T), dimension(2) :: pdims, pmaxdims, cdims, off_h5, off_m
   integer vsp_file, datarank, np, ipart
   integer :: ielm, ielmold = 1, size, poffset, isghost
   integer :: info
   logical :: pass = .false.
   integer :: ipart_begin_temp, ipart_end_temp

   locparts = 0

   !Set up the file access property list to use parallel I/O
   call h5pcreate_f(H5P_FILE_ACCESS_F, access_props, ierr)
   info = MPI_INFO_NULL
   call h5pset_fapl_mpio_f(access_props, MPI_COMM_WORLD, info, ierr)

   !Open the file for reading
   call second(tstart)
   call h5fopen_f(filename, H5F_ACC_RDONLY_F, fileid, ierr, access_props)
   if (ierr .lt. 0) then
      if (myrank .eq. 0) &
         print *, "Error: could not open ", filename, &
         " for HDF5 input.  error = ", ierr
      return
   end if

   !Open the root group
   call h5gopen_f(fileid, "/", part_root_id, ierr)

   !Read global attributes
   !call read_int_attr(part_root_id, "velocity space dims", vsp_file, ierr)
   !if (vsp_file.ne.vspdims) then
   !if (myrank.eq.0) then !Abort if full-orbit vs drift kinetic models do not agree
   !print *,'Error: incompatible velocity space dimensions in restart file!'
   !print *,vsp_file,' vs ',vspdims
   !endif
   !call h5gclose_f(part_root_id, ierr)
   !call h5fclose_f(fileid, ierr)
   !call h5pclose_f(access_props, ierr)
   !ierr = -1
   !return
   !endif
   !call read_real_attr(part_root_id, "particle delta-t", dt_ion, ierr)
   !if(myrank.eq.0) print *,'particle dt from file = ',dt_ion
   call read_vec_attr(part_root_id, "particle nrmfac", nrmfac, 2, ierr)

   !Open the particle group for reading
   call h5gopen_f(part_root_id, "particles", group_id, ierr)
   if (ierr .lt. 0) then !group not found
      if (myrank .eq. 0) &
         print *, 'Error ', ierr, ': particle group not found in restart file!'
   else
      !Open the dataset
      call h5dopen_f(group_id, "data", dset_id, ierr)
      if (ierr .lt. 0) then !dataset not found
         if (myrank .eq. 0) &
            print *, 'Error ', ierr, ': dataset not found in restart file!'
      else
         !Get a copy of the dataspace for this dataset
         call h5dget_space_f(dset_id, filespace, ierr)
         if (ierr .lt. 0) then !dataspace not acquired
            if (myrank .eq. 0) &
               print *, 'Error ', ierr, ': dataspace not found in restart file!'
         else
            !Query dataspace dimensions
            call h5sget_simple_extent_ndims_f(filespace, datarank, ierr)
            if (datarank .ne. 2) then
               if (myrank .eq. 0) &
                  print *, 'Error ', ierr, ': dataspace has wrong dimensions in restart file!'
            else
               call h5sget_simple_extent_dims_f(filespace, pdims, pmaxdims, ierr)
               if (pdims(1) .ne. ldim) then
                  if (myrank .eq. 0) &
                     print *, 'Error ', ierr, &
                     ': dataspace has wrong dimensions in restart file!'
               else
                  if (myrank .eq. 0) print *, pdims(2), ' particles in restart file.'
                  nparticles = pdims(2)
                  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
                  np = int(nparticles/size) + 1
                  allocate (valbuf(ldim, np))
                  !call mpi_scan(locparts, poffset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
                  !poffset = poffset - locparts
                  poffset = np*myrank

                  if (poffset + np > nparticles) np = nparticles - poffset
                  cdims(1) = ldim; cdims(2) = np
                  call h5screate_simple_f(2, cdims, memspace, ierr)
                  off_h5(1) = 0; off_h5(2) = poffset
                  off_m = off_h5

                  !Read the particle data, assign particles to processors
                  !do while (np.gt.0)      !For each chunk...
                  !Set location to read from
                  !if (np.lt.chunksize) then
                  !   cdims(2) = np
                  !   call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                  !        off_m, cdims, ierr)
                  !endif
                  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off_h5, &
                                             cdims, ierr)
                  if (ierr .lt. 0) return
                  !off_h5(2) = off_h5(2) + cdims(2)

                  !Read ID, phase space coordinates, weight
                  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, valbuf, cdims, ierr, &
                                 mem_space_id=memspace, file_space_id=filespace)

                  if (ierr .lt. 0) return !Read error

                  !Assign particles to local domain as needed
                  do ipart = 1, cdims(2)
                     !Test for local residence
                     mmsa = valbuf(2:4, ipart)
                     call mesh_search(ielmold, mmsa, ielm)
                     !ielm = ielm + 1
                     !if (ielm.le.0) cycle     !Not in local partition; skip.
                     if (ielm <= 0) then
                        dpar%deleted = .true.
                     else
                        dpar%deleted = .false.
                        ielmold = ielm
                     end if
                     !call m3dc1_ent_isghost(2, ielm-1, isghost)
                     !if (isghost.eq.1) cycle  !In ghost layer; skip.

                     !Particle is local -> add it to pdata.
                     dpar%gid = valbuf(1, ipart)
                     dpar%x(1) = valbuf(2, ipart)
                     dpar%x(2) = valbuf(3, ipart)
                     !dpar%x(2) = 3.1415926/2.
                     dpar%x(3) = valbuf(4, ipart)
                     dpar%wt = valbuf(5, ipart)
                     ! dpar%wt = 0.
                     ! dpar%wt = dpar%wt*1.19
                     ! dpar%wt = dpar%wt*2.666
                     ! dpar%wt = dpar%wt*17.
                     dpar%sps = valbuf(6, ipart)
                     !if (dpar%sps == 1) kinetic_thermal = 1
                     dpar%v(1) = valbuf(7, ipart)
                     !dpar%v(1) = dpar%v(1)*10.
                     !if (mod(dpar%gid,2)==1) dpar%v(1) = -valbuf(7,ipart)
                     dpar%v(2) = valbuf(8, ipart)
                     if (vspdims == 5) then
                        dpar%v(3) = valbuf(9, ipart)
                        dpar%v(4) = valbuf(10, ipart)
                     end if
                     ! dpar%v(vspdims) = valbuf(ldim, ipart)
                     dpar%f0 = valbuf(9, ipart)
                     dpar%x0(1) = valbuf(10, ipart)
                     dpar%x0(2) = valbuf(11, ipart)
                     dpar%x0(3) = valbuf(12, ipart)
                     dpar%v0(1) = valbuf(13, ipart)
                     dpar%v0(2) = valbuf(14, ipart)
                     dpar%B0 = 0.
                     dpar%kx(:, 1) = dpar%x
                     dpar%kx(:, 2) = dpar%x
                     dpar%kx(:, 3) = dpar%x
                     dpar%kx(:, 4) = dpar%x
                     dpar%jel = ielm
                     dpar%kel(:) = ielm
                     !if (dpar%gid==0.3*1e6) then
                     !   dpar%x(1)=0.9
                     !   dpar%x(2)=3.1415926/2.
                     !   dpar%x(3)=0
                     !endif
                     !call add_particle(pdata(ielm), dpar)
                     pdata(ipart + poffset) = dpar
                     locparts = locparts + 1
                  end do
                  ipart_begin = poffset + 1
                  ipart_end = poffset + np
                  !np = np - chunksize
                  !if (np.lt.pdims(2)-10) exit !TMP!!!
                  !enddo !while

                  if (myrank .eq. 0) then
                     call second(tend)
                     write (0, '(I12,A,f9.2,A)') &
                        pdims(2), ' particles read, assigned in ', tend - tstart, ' sec.'
                     !print *,'last =', valbuf(1,cdims(2)),valbuf(7,cdims(2))
                  end if
                  call h5sclose_f(memspace, ierr)

                  print *, myrank, ': ', locparts, ' local particles.'
                  !call mpi_reduce(locparts, nparticles, 1, MPI_INTEGER, MPI_SUM, 0, &
                  !     MPI_COMM_WORLD, ierr)
                  !if (nparticles.ne.pdims(2)) then
                  !   if (myrank.eq.0) print *,&
                  !        "Error: ",nparticles," assigned, ",pdims(2)," expected."
                  !else
                  !   pass = .true.
                  !endif
                  deallocate (valbuf)
               end if !Wrong leading dimension
            end if !Wrong number of dimensions

            !Close the dataspace
            call h5sclose_f(filespace, ierr)
         end if !dataspace not acquired

         !Close the particle dataset
         call h5dclose_f(dset_id, ierr)
      end if !dataset failed to open

      !Close the particle group
      call h5gclose_f(group_id, ierr)
   end if !particle group failed to open

   !Close the file along with other open data structures
   call h5gclose_f(part_root_id, ierr)
   call h5fclose_f(fileid, ierr)
   call h5pclose_f(access_props, ierr)

   if (.not. pass) ierr = -42
   ipart_begin_temp=ipart_begin
   call mpi_allreduce(ipart_begin_temp, ipart_begin, 1, MPI_INTEGER, MPI_MIN, hostcomm, ierr)
   ipart_end_temp=ipart_end
   call mpi_allreduce(ipart_end_temp, ipart_end, 1, MPI_INTEGER, MPI_MAX, hostcomm, ierr)
   if (hostrank == 0) then
      write (0, *) rowrank, ipart_begin, ipart_end
   end if

end subroutine hdf5_read_particles

subroutine set_parallel_velocity

   use mesh_mod
   use basic
   use arrays
   use sparse
   use m3dc1_nint
   use diagnostics
   use boundary_conditions
   use matrix_mod
   use transport_coefficients
   use gyroviscosity
   use runaway_mod
   use auxiliary_fields
   use newvar_mod
   use model

   implicit none

   vectype, dimension(dofs_per_element) :: r4
   integer :: k, itri, izone
   integer :: ieq(3)
   integer, dimension(dofs_per_element) :: imask
   type(vector_type), pointer :: vsource
   type(field_type) ::   p_v
   vectype, dimension(dofs_per_element) :: dofs
   integer :: ierr
   type(field_type) ::   u_v
   type(field_type) ::  vz_v
   type(field_type) :: chi_v
   type(vector_type) :: vel_vec
   type(vector_type) :: b1_vel
   vectype, dimension(dofs_per_element, dofs_per_element) :: tempxx
   logical, save :: first_time = .true.
    

     !call matvecmult(d1_0_mat, vel_vec, b1_vel)
    call create_vector(vel_vec,      numvar)
    call associate_field(u_v,    vel_vec,      1)
    call associate_field(vz_v,  vel_vec,    2)
    call associate_field(chi_v,  vel_vec,    3)
    u_v = u_field(1)
    vz_v = vz_field(1)
    chi_v = chi_field(1)
   call create_vector(b1_vel, numvar)
     
   if (first_time) then
      call set_matrix_index(s1_0_mat, 172)
       call create_mat(s1_0_mat, numvar, numvar, icomplex, 1)
       u_i = 1
       vz_i = 2
       chi_i = 3
       do itri=1,local_elements()
          call define_element_quadrature(itri, int_pts_main, int_pts_tor)
          call define_fields(itri,FIElD_N,1,0)
          call get_vor_mask(itri, imask)
          tempxx = intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DR),r2_79)
          tempxx = tempxx+intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),r2_79)
          call apply_boundary_mask(itri, 0, tempxx, imask)
          call insert_block(s1_0_mat,itri,u_i,u_i,tempxx,MAT_ADD)
          tempxx = intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DZ),ri_79)
          tempxx = tempxx-intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DR),ri_79)
          call apply_boundary_mask(itri, 0, tempxx, imask)
          call insert_block(s1_0_mat,itri,u_i,chi_i,tempxx,MAT_ADD)
          call get_vz_mask(itri, imask)
          tempxx = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_1),r2_79)
          call apply_boundary_mask(itri, 0, tempxx, imask)
          call insert_block(s1_0_mat,itri,vz_i,vz_i,tempxx,MAT_ADD)
          call get_chi_mask(itri, imask)
          tempxx = intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DZ),ri_79)
          tempxx = tempxx-intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DR),ri_79)
          call apply_boundary_mask(itri, 0, tempxx, imask)
          call insert_block(s1_0_mat,itri,chi_i,u_i,tempxx,MAT_ADD)
          tempxx = -intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DR),ri4_79)
          tempxx = tempxx-intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),ri4_79)
          call apply_boundary_mask(itri, 0, tempxx, imask)
          call insert_block(s1_0_mat,itri,chi_i,chi_i,tempxx,MAT_ADD)
          if (isplitstep.eq.0) then
             call get_flux_mask(itri, imask)
             tempxx = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
             !call apply_boundary_mask(itri, 0, tempxx, imask)
             call insert_block(s1_0_mat,itri,psi_i,psi_i,tempxx,MAT_ADD)
             call get_bz_mask(itri, imask)
             tempxx = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
             !call apply_boundary_mask(itri, 0, tempxx, imask)
             call insert_block(s1_0_mat,itri,bz_i,bz_i,tempxx,MAT_ADD)
             call get_pres_mask(itri, imask)
             tempxx = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
             !call apply_boundary_mask(itri, 0, tempxx, imask)
             call insert_block(s1_0_mat,itri,p_i,p_i,tempxx,MAT_ADD)
             tempxx = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
             call insert_block(s1_0_mat,itri,den_i,den_i,tempxx,MAT_ADD)
             !call get_bf_mask(itri, imask)
             !tempxx = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
             !!call apply_boundary_mask(itri, 0, tempxx, imask)
             !call insert_block(s1_0_mat,itri,bf_i,bf_i,tempxx,MAT_ADD)
          endif
       enddo
       call flush(s1_0_mat)

       call boundary_vel(b1_vel, u_v, vz_v, chi_v, s1_0_mat)
       call finalize(s1_0_mat)
       first_time = .false.
    endif

   ieq(1) = u_i
   ieq(2) = vz_i
   ieq(3) = chi_i

   if (isplitstep .ge. 1) then
      vsource => r4_vec
   else
      vsource => q4_vec
   end if
   vsource = 0.
   do itri = 1, local_elements()

      call get_zone(itri, izone)

      call define_element_quadrature(itri, int_pts_main, int_pts_tor)
      call define_fields(itri, FIELD_PSI + FIELD_I + FIELD_PHI + FIELD_V + FIELD_CHI+FIELD_N+FIELD_NI+FIELD_KIN, 1, 0)

      !call define_fields(itri, FIELD_PSI + FIELD_I + FIELD_P + FIELD_PHI + FIELD_V + FIELD_CHI + FIELD_ETA + FIELD_MU + FIELD_N + FIELD_NI+FIELD_B2I+FIELD_KIN, 1, 1)

      if (linear.eq.1) then
         pst79 = ps079
         bzt79 = bz079
         bfpt79 = bfp079
      endif
#if defined(USE3D) || defined(USECOMPLEX)
          temp79a= ((ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(r_79*ph179(:,OP_DR)+ri2_79*ch179(:,OP_DZ)) &
                  +(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-r_79*ph179(:,OP_DZ)+ri2_79*ch179(:,OP_DR)) &
                 + bzt79(:,OP_1)*vz179(:,OP_1)) &
            /(ri2_79* &
            ((pst79(:,OP_DR)-r_79*bfpt79(:,OP_DZ))**2 + (pst79(:,OP_DZ)+r_79*bfpt79(:,OP_DR))**2 + bzt79(:,OP_1)*bzt79(:,OP_1)))
#else
          temp79a= ((ri_79*pst79(:,OP_DR))*(r_79*ph179(:,OP_DR)+ri2_79*ch179(:,OP_DZ)) &
                  +(-ri_79*pst79(:,OP_DZ))*(-r_79*ph179(:,OP_DZ)+ri2_79*ch179(:,OP_DR)) &
                 + bzt79(:,OP_1)*vz179(:,OP_1)) &
            /(ri2_79* &
            ((pst79(:,OP_DR))**2 + (pst79(:,OP_DZ))**2 + bzt79(:,OP_1)*bzt79(:,OP_1)))
#endif
      temp79b=temp79a
      !temp79c=(ps079(:,OP_1)-0.18)/0.02/2+0.1
      temp79c=1.
      !where (real(temp79c(:))>1.0)
      !   temp79c(:)=1.0
      !end where
      !where (real(temp79c(:))<0.1)
      !!!where (real(ps079(:,OP_1))<0.18)
      !!   !temp79a(:)=-temp79b(:)
      !   temp79c(:)=0.1
      !end where
      if (ikinetic_vpar.eq.1) then
         temp79a= (vpar_reduce)*temp79c*(1*vipar79(:,OP_1)+1*(1*vfpar79(:,OP_1)-vfpar079(:,OP_1)*nf79(:,OP_1)))&
            /(nfi079(:,OP_1))/sqrt(ri2_79* &
      ! temp79a= 0.99*vipar79(:,OP_1)/nfi079(:,OP_1)/sqrt(ri2_79* &
      ! temp79a= vipar79(:,OP_1)*ni79(:,OP_1)/sqrt(ri2_79* &
      ((pst79(:,OP_DR)-r_79*bfpt79(:,OP_DZ))**2 + (pst79(:,OP_DZ)+r_79*bfpt79(:,OP_DR))**2 + bzt79(:,OP_1)*bzt79(:,OP_1))) &
         -vpar_reduce*temp79c*temp79b
       else
       temp79a=-vpar_reduce*temp79b
        endif
       !where (real(rhof79(:,OP_1))<0.1)
       !   !temp79a(:)=-temp79b(:)
       !   temp79a(:)=0
       !end where
          !temp79a(:)=-temp79b(:)
          ! temp79a(:)=0.
       temp79c=1.

         do k = 1, numvar
         r4 = 0.
         if (izone .eq. ZONE_PLASMA) then
              select case(k)
              case(1)
                 r4 = intx4(mu79(:,:,OP_DR),ph179(:,OP_DR),r2_79,temp79c)
                 r4 = r4 + intx4(mu79(:,:,OP_DZ),ph179(:,OP_DZ),r2_79,temp79c)
                 r4 = r4 + intx4(mu79(:,:,OP_DR),ch179(:,OP_DZ),ri_79,temp79c)
                 r4 = r4 - intx4(mu79(:,:,OP_DZ),ch179(:,OP_DR),ri_79,temp79c)
#if defined(USE3D) || defined(USECOMPLEX)
                 !r4 = 0.
                 r4 = r4 + intx4(mu79(:,:,OP_DR),pst79(:,OP_DR)-r_79*bfpt79(:,OP_DZ),temp79a,temp79c)
                 r4 = r4 + intx4(mu79(:,:,OP_DZ),pst79(:,OP_DZ)+r_79*bfpt79(:,OP_DR),temp79a,temp79c)
#else
                 r4 = r4 + intx4(mu79(:,:,OP_DR),pst79(:,OP_DR),temp79a,temp79c)
                 r4 = r4 + intx4(mu79(:,:,OP_DZ),pst79(:,OP_DZ),temp79a,temp79c)
#endif
              case(2)
                 !r4 = 0.
                 r4 = intx4(mu79(:,:,OP_1),vz179(:,OP_1),r2_79,temp79c)
                 r4 = r4 + intx4(mu79(:,:,OP_1),bzt79(:,OP_1),temp79a,temp79c)
              case(3)
                 r4 = intx4(mu79(:,:,OP_DR),ph179(:,OP_DZ),ri_79,temp79c)
                 r4 = r4-intx4(mu79(:,:,OP_DZ),ph179(:,OP_DR),ri_79,temp79c)
                 r4 = r4-intx4(mu79(:,:,OP_DR),ch179(:,OP_DR),ri4_79,temp79c)
                 r4 = r4-intx4(mu79(:,:,OP_DZ),ch179(:,OP_DZ),ri4_79,temp79c)
#if defined(USE3D) || defined(USECOMPLEX)
                 !r4=0.
                 r4 = r4+intx5(mu79(:,:,OP_DR),pst79(:,OP_DZ)+r_79*bfpt79(:,OP_DR),temp79a,ri3_79,temp79c)
                 r4 = r4-intx5(mu79(:,:,OP_DZ),pst79(:,OP_DR)-r_79*bfpt79(:,OP_DZ),temp79a,ri3_79,temp79c)
#else
                 r4 = r4+intx5(mu79(:,:,OP_DR),pst79(:,OP_DZ),temp79a,ri3_79,temp79c)
                 r4 = r4-intx5(mu79(:,:,OP_DZ),pst79(:,OP_DR),temp79a,ri3_79,temp79c)
#endif
              end select
            select case (k)
            case (1)
               call get_vor_mask(itri, imask)
            case (2)
               call get_vz_mask(itri, imask)
            case (3)
               call get_chi_mask(itri, imask)
            end select
         end if
         call apply_boundary_mask_vec(itri, 0, r4, imask)
         call vector_insert_block(vsource,itri,ieq(k),r4,VEC_ADD)
      end do

   end do
   call sum_shared(vsource)
     b1_vel=r4_vec
     call boundary_vel(b1_vel, u_v, vz_v, chi_v)
     call newsolve(s1_0_mat, b1_vel, ierr)
     vel_vec = b1_vel

     u_field(1)=u_v
     vz_field(1)=vz_v
     chi_field(1)=chi_v

     call destroy_vector(b1_vel)
     call destroy_vector(vel_vec)
   vsource=0.


     call create_field(p_v)
   p_v=0.
  do itri=1,local_elements()
     call define_element_quadrature(itri, int_pts_main, int_pts_tor)
     call define_fields(itri, FIELD_PSI + FIELD_I + FIELD_PHI + FIELD_V + FIELD_CHI+FIELD_N+FIELD_NI+FIELD_KIN, 1, 0)
     !temp79a=(pipar79(:,OP_1)*0+piper79(:,OP_1)*3)/3.
     !temp79a=nfi79(:,OP_1)*te079(:,OP_1)*2
     !temp79a=nfi79(:,OP_1)*te079(:,OP_1)+p179(:,OP_1)-n179(:,OP_1)*te0
     temp79a= ((ri_79*ps079(:,OP_DR)-bfp079(:,OP_DZ))*(r_79*ph179(:,OP_DR)+ri2_79*ch179(:,OP_DZ)) &
             +(-ri_79*ps079(:,OP_DZ)-bfp079(:,OP_DR))*(-r_79*ph179(:,OP_DZ)+ri2_79*ch179(:,OP_DR)) &
              + bz079(:,OP_1)*vz179(:,OP_1)) &
            /sqrt(ri2_79* &
            ((ps079(:,OP_DR)-r_79*bfp079(:,OP_DZ))**2 + (ps079(:,OP_DZ)+r_79*bfp079(:,OP_DR))**2 + bz079(:,OP_1)*bz079(:,OP_1)))
     ! temp79b=(r_79*ph179(:,OP_DR)+ri2_79*ch179(:,OP_DZ))**2+ &
     !         +(-r_79*ph179(:,OP_DZ)+ri2_79*ch179(:,OP_DR))**2+ (r_79*vz179(:,OP_1))**2
     ! temp79a = sqrt(temp79b-temp79a**2)
             
     ! temp79a= ps179(:,OP_GS)
     ! temp79a= ri2_79*ps179(:,OP_DRP)-ri_79*bz179(:,OP_DZ)-ri_79*bfp179(:,OP_DZP)
      dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(p_v%vec,itri,1,dofs,VEC_ADD)
  end do
  !call newvar_solve(p_v%vec,diff_mat)
  !call sum_shared(p_v%vec)
  !call newsolve(diff3_mat, p_v%vec, ierr)
  call newvar_solve(p_v%vec,mass_mat_lhs)
  !if(calc_matrices.eq.1) then
  !write(0,*) "111111111111111111"
  vparsmooth_field=p_v
  call destroy_field(p_v)


end subroutine set_parallel_velocity

subroutine set_density

   use mesh_mod
   use basic
   use arrays
   use sparse
   use m3dc1_nint
   use diagnostics
   use boundary_conditions
   use matrix_mod
   use transport_coefficients
   use gyroviscosity
   use runaway_mod
   use auxiliary_fields
   use newvar_mod
   use model

   implicit none

   type(field_type) ::   p_v
   vectype, dimension(dofs_per_element) :: dofs
   integer :: k, itri, izone
   integer, dimension(dofs_per_element) :: imask
   integer :: ierr
   type(vector_type) :: temp, temp2
   type(field_type) ::   u_v
   type(field_type) ::  vz_v
   type(field_type) :: chi_v
   type(vector_type) :: den_vec, vel_vec
   logical, save :: first_time = .true.

   call create_field(p_v)
   p_v=0.
   do itri=1,local_elements()
      call get_zone(itri, izone)

      call define_element_quadrature(itri, int_pts_main, int_pts_tor)
      call define_fields(itri, FIELD_PSI+FIELD_P+FIELD_N+FIELD_NI+FIELD_KIN, 1, 0)
      ! temp79a = p179(:,OP_1) + 0.1*dt*nfi79(:,OP_1)*te079(:,OP_1)*2 &
      !   -0.1*dt*n179(:,OP_1)*te079(:,OP_1)*2
      temp79a = n179(:,OP_1) + 0.1*(nfi79(:,OP_1)+nf79(:,OP_1)) &
        -0.1*n179(:,OP_1)
      !where (real(rhof79(:,OP_1))>0.5)
      !    !temp79a(:)=-temp79b(:)
      !    temp79a(:)=n179(:,OP_1)
      ! end where

      !temp79a = n179(:,OP_1) + 0.9*(nfi79(:,OP_1))&
      !          -0.5*n179(:,OP_1)
      !temp79a = n179(:,OP_1)
       where (real(temp79a(:)+n079(:,OP_1))<0.035)
      !   ! temp79a(:)=p179(:,OP_1)
          temp79a(:)=0.035-n079(:,OP_1)
      !   !temp79a(:)=0
      end where
      !temp79a(:)=0.
      dofs = intx2(mu79(:,:,OP_1),temp79a)
      call vector_insert_block(p_v%vec,itri,1,dofs,VEC_ADD)
  enddo
  call newvar_solve(p_v%vec,mass_mat_lhs)
  !call sum_shared(p_v%vec)
  !call newsolve(diff3_mat, p_v%vec, ierr)
  den_field(1) = p_v

  call create_vector(den_vec, 1)
  call associate_field(p_v, den_vec, 1)

  call create_vector(temp, 1)
  call create_vector(temp2, 1)

  call create_vector(vel_vec,      numvar)
  call associate_field(u_v,    vel_vec,      1)
  call associate_field(vz_v,  vel_vec,    2)
  call associate_field(chi_v,  vel_vec,    3)
  u_v = u_field(1)
  vz_v = vz_field(1)
  chi_v = chi_field(1)

  call matvecmult(r8_mat,vel_vec,temp)
  call matvecmult(q8_mat,vel_vec,temp2)
  call add(temp, temp2)

  call matvecmult(d8_mat,den_vec,temp2)
  call add(temp, temp2)
     
  ! Insert boundary conditions
  if(first_time) then
     call boundary_den(temp, p_v, s8_mat)
     call finalize(s8_mat)
     first_time = .false.
  else
     call boundary_den(temp, p_v)
  endif

  call newsolve(s8_mat, temp, ierr)

  den_vec = temp
  call destroy_vector(temp)
  call destroy_vector(temp2)
  den_field(1) = p_v
  call destroy_field(p_v)

  call calculate_ne(1, den_field(1), ne_field(1), eqsubtract)
  if(itemp.eq.0 .and. (numvar.eq.3 .or. ipres.gt.0)) then
     call calculate_temperatures(1, te_field(1), ti_field(1), &
          pe_field(1), p_field(1), ne_field(1), den_field(1), eqsubtract)
  else
     call calculate_pressures(1, pe_field(1), p_field(1), ne_field(1), &
          den_field(1), te_field(1), ti_field(1), eqsubtract)
  endif
end subroutine set_density

subroutine set_den_smooth

   use mesh_mod
   use basic
   use arrays
   use sparse
   use m3dc1_nint
   use diagnostics
   use boundary_conditions
   use matrix_mod
   use transport_coefficients
   use gyroviscosity
   use runaway_mod
   use auxiliary_fields
   use newvar_mod

   implicit none

   vectype, dimension(dofs_per_element) :: dofs
   integer :: k, itri, izone
   integer :: ieq(3)
   integer, dimension(dofs_per_element) :: imask
   type(vector_type), pointer :: vsource
    type(field_type) ::   p_v
   integer :: ierr

    call create_field(p_v)
   p_v=0.
  do itri=1,local_elements()
     call define_element_quadrature(itri,int_pts_main,int_pts_tor)
     call define_fields(itri,FIELD_P+FIELD_TE+FIELD_KIN+FIELD_N+FIELD_NI,1,0)
     temp79a = n179(:,OP_1) + 0.9*(nfi79(:,OP_1)+nf79(:,OP_1)) &
        -0.9*n179(:,OP_1)

     !temp79a=(pipar79(:,OP_1)*0+piper79(:,OP_1)*3)/3.
     !temp79a=nfi79(:,OP_1)*te079(:,OP_1)*2
     !temp79a=nfi79(:,OP_1)*te079(:,OP_1)+p179(:,OP_1)-n179(:,OP_1)*te0
     !temp79a=n179(:,OP_1)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(p_v%vec,itri,1,dofs,VEC_ADD)
  end do
  call sum_shared(p_v%vec)
  call newsolve(diff3_mat, p_v%vec, ierr)
  !call newvar_solve(p_v%vec,mass_mat_lhs)
  !if(calc_matrices.eq.1) then
  !write(0,*) "111111111111111111"
  densmooth_field=p_v
  call destroy_field(p_v)

end subroutine set_den_smooth

#endif

end module particles
