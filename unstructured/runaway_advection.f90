! Kinetic energetic ion module, J. Breslau, 2015
module runaway_advection
   use mesh_mod
   use field
   use matrix_mod
   use newvar_mod
   use gradshafranov
   !use mpi_f08
   use ieee_arithmetic
   use iso_c_binding
   implicit none
   private
   public :: runaway_advection_initialize, runaway_advection_step

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
      vectype, dimension(coeffs_per_element) :: rho
      vectype, dimension(coeffs_per_element) :: jrev0, jrev1
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
      vectype                  :: wt    !Particle weighting in delta-f scheme
      vectype                  :: f0
      integer                  :: gid         !Unique global particle index
      integer                  :: jel         !Predicted element of residence
      logical                  :: deleted, deleted2
   end type particle

   type(particle), dimension(:), pointer :: pdata  !Particle arrays
   integer :: win_pdata
   !type(elfield), dimension(:), allocatable, target :: elfieldcoefs        !Receive buffer for jumping particles
   type(elfield), dimension(:), pointer :: elfieldcoefs       !Receive buffer for jumping particles
!$acc declare link(elfieldcoefs)
   integer :: win_elfieldcoefs
   integer :: linear_particle, psubsteps
!$acc declare create(linear_particle, psubsteps)
   real :: dt_particle, t0_norm_particle, v0_norm_particle, b0_norm_particle
!$acc declare create(dt_particle, t0_norm_particle, v0_norm_particle,b0_norm_particle)
   complex :: rfac_particle
!$acc declare create(rfac_particle)
   real, dimension(:, :, :), pointer :: mesh_coord !Neighbor tracking arrays
!$acc declare link(mesh_coord)
   integer :: win_mesh_coord
   integer, dimension(:, :), pointer :: mesh_nodes !Neighbor tracking arrays
   integer :: win_mesh_nodes
   integer, dimension(:), pointer :: mesh_zone !Neighbor tracking arrays
!$acc declare link(mesh_zone)
   integer :: win_mesh_zone
   integer, dimension(:, :), pointer :: particle_map !Neighbor tracking arrays
!$acc declare link(particle_map)
   integer :: win_particle_map
   integer, dimension(:, :), allocatable :: neighborlist !Neighbor tracking arrays
!$acc declare link(neighborlist)
   integer, dimension(:), allocatable :: localmeshid
   integer :: npar, nparticles, locparts
   integer :: mpi_particle !User-defined MPI datatype for particle communication
   integer :: mpi_elfield
   type(matrix_type) :: diff2_mat
   !type(newvar_matrix) :: diff_mat
   type(field_type) :: nre_vec
   integer :: hostcomm, rowcomm
   integer :: hostrank, rowrank, ncols, nrows
   integer, dimension(:), allocatable :: nelm_row, displs_elm
   integer :: ielm_min, ielm_max, ipart_begin, ipart_end
   integer :: nelms, nelms_global, nnodes_global
   real :: toroidal_period_particle
!$acc declare create(toroidal_period_particle)
contains

!Define MPI datatype for particle communication
! Note: any changes to the "particle" user-defined datatype must be reflected
!       in the definitions of pnvars, pblklen, ptyps, and pdspls below.
subroutine define_mpi_particle(ierr)
   implicit none

   include 'mpif.h'

   integer, intent(out) :: ierr
   integer, parameter :: pnvars = 8
   integer, dimension(pnvars), parameter :: pblklen = (/3, vspdims, 1, 1, 1, 1, 1, 1/)
   integer(kind=MPI_ADDRESS_KIND), dimension(pnvars) :: pdspls
   integer, dimension(pnvars), parameter :: ptyps = (/MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
 MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_INTEGER, MPI_LOGICAL, MPI_LOGICAL/)

   type(particle) :: dum_par

   !Set up component displacements array
   call mpi_get_address(dum_par%x, pdspls(1), ierr)
   call mpi_get_address(dum_par%v, pdspls(2), ierr)
   pdspls(2) = pdspls(2) - pdspls(1)
   call mpi_get_address(dum_par%wt, pdspls(3), ierr)
   pdspls(3) = pdspls(3) - pdspls(1)
   call mpi_get_address(dum_par%f0, pdspls(4), ierr)
   pdspls(4) = pdspls(4) - pdspls(1)
   call mpi_get_address(dum_par%gid, pdspls(5), ierr)
   pdspls(5) = pdspls(5) - pdspls(1)
   call mpi_get_address(dum_par%jel, pdspls(6), ierr)
   pdspls(6) = pdspls(6) - pdspls(1)
   call mpi_get_address(dum_par%deleted, pdspls(7), ierr)
   pdspls(7) = pdspls(7) - pdspls(1)
   call mpi_get_address(dum_par%deleted, pdspls(8), ierr)
   pdspls(8) = pdspls(8) - pdspls(1)
   pdspls(1) = 0

   call mpi_type_create_struct(pnvars, pblklen, pdspls, ptyps, mpi_particle, ierr)
   if (ierr .ne. 0) return
   call mpi_type_commit(mpi_particle, ierr)
end subroutine define_mpi_particle

subroutine define_mpi_elfield(ierr)
   implicit none

   include 'mpif.h'

   integer, intent(out) :: ierr
   integer, parameter :: pnvars = 12
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
   call mpi_get_address(dum_elfield%jrev0, pdspls(11), ierr)
   pdspls(11) = pdspls(11) - pdspls(1)
   call mpi_get_address(dum_elfield%jrev1, pdspls(12), ierr)
   pdspls(12) = pdspls(12) - pdspls(1)
   pdspls(1) = 0

   call mpi_type_create_struct(pnvars, pblklen, pdspls, ptyps, mpi_elfield, ierr)
   if (ierr .ne. 0) return
   call mpi_type_commit(mpi_elfield, ierr)
end subroutine define_mpi_elfield
!---------------------------------------------------------------------------
subroutine runaway_advection_initialize
   use basic
   use diagnostics
   use auxiliary_fields
   implicit none
   include 'mpif.h'

   integer, parameter :: trunit = 120
   real :: tstart, tend
   integer :: ierr

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

   !call update_particle_pressure

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
end subroutine runaway_advection_initialize
!---------------------------------------------------------------------------
real function tri_area(x1, z1, x2, z2, x3, z3)
!$acc routine seq
   implicit none

   real :: x1, z1, x2, z2, x3, z3

   tri_area = x1*(z2 - z3) + x2*(z3 - z1) + x3*(z1 - z2)
end function tri_area
!---------------------------------------------------------------------------
subroutine init_particles(lrestart, ierr)
   use basic
   use arrays
   use m3dc1_nint
   use gradshafranov
   use read_ascii
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
   real    :: maxspeed, lambda_min, lambda_max, B0, Bmax, pphi
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
   real, dimension(3) :: gradpsi
   real :: gradcoef
   real :: rho, f0, fx, T00, f_mesh, df0de
   real, allocatable, dimension(:) :: rho_vals, psi_vals, nf_vals, tf_vals
   integer :: ipoint, dotloc

   real, dimension(3) :: bhat, Bstar, Jcyl, svec, dBdR, dBdphi, dBdz, gradB0, gradB1, BxgrdB, dvdt
   real :: Rinv, Bss, dB1, di, R_axis, di_axis
   integer :: num_f
   real, dimension(:), allocatable :: f_array_temp
   real :: radi, pitch, energy
   integer :: radi_i, pitch_i, energy_i
   real, dimension(:, :, :), allocatable :: f_array2
   real :: f1, f2, f3, f4, f5, f6
   integer :: ran_size
   integer, allocatable :: seed(:)
   integer :: locparts2, izone, ipart

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

   disp_unit = 1
   arraysize = 0
   if (hostrank == 0) arraysize = nelms_global*sizeof(1)
   CALL MPI_Win_allocate_shared(arraysize, disp_unit, MPI_INFO_NULL, &
                                hostcomm, baseptr, win_mesh_zone, ierr)
   if (hostrank /= 0) CALL MPI_Win_shared_query(win_mesh_zone, 0, arraysize, disp_unit, baseptr, ierr)
   !CALL C_F_POINTER(baseptr, mesh_coord_shared, [2,3,nelms])
   CALL C_F_POINTER(baseptr, mesh_zone, [nelms_global])
   !allocate(mesh_coord(2,3,nelms))

   disp_unit = 1
   arraysize = 0
   if (hostrank == 0) arraysize = nelms_global*npoints*sizeof(1)
   CALL MPI_Win_allocate_shared(arraysize, disp_unit, MPI_INFO_NULL, &
                                hostcomm, baseptr, win_particle_map, ierr)
   if (hostrank /= 0) CALL MPI_Win_shared_query(win_particle_map, 0, arraysize, disp_unit, baseptr, ierr)
   !CALL C_F_POINTER(baseptr, mesh_coord_shared, [2,3,nelms])
   CALL C_F_POINTER(baseptr, particle_map, [npoints, nelms_global])
   !allocate(mesh_coord(2,3,nelms))

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
      call get_zone(ielm, izone)
      if (ielm_min > ielm_global) ielm_min = ielm_global
      if (ielm_max < ielm_global) ielm_max = ielm_global
      localmeshid(ielm_global) = ielm
      call get_element_nodes(ielm, nodeids)
      do inode = 1, nodes_per_element
         call m3dc1_ent_getglobalid(0, nodeids(inode) - 1, inode_global)
         mesh_nodes(inode, ielm_global) = inode_global + 1
         call get_node_pos(nodeids(inode), mesh_coord(1,inode,ielm_global), mesh_coord(2,inode,ielm_global), mesh_coord(3,inode,ielm_global))
      end do

      call get_zone(ielm, izone)
      mesh_zone(ielm_global) = izone
   end do
   call mpi_barrier(mpi_comm_world, ierr)
   call mpi_allreduce(ielm_min, ielm_min, 1, MPI_INTEGER, MPI_MIN, hostcomm, ierr)
   call mpi_allreduce(ielm_max, ielm_max, 1, MPI_INTEGER, MPI_MAX, hostcomm, ierr)
   allocate (nelm_row(nrows), displs_elm(nrows))
   nelm_row(rowrank + 1) = ielm_max - ielm_min + 1
   call mpi_allgather(nelm_row(rowrank + 1), 1, MPI_INTEGER, nelm_row, 1, MPI_INTEGER, rowcomm, ierr)
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
      call MPI_ALLGATHERV(mesh_nodes(:,ielm_min:ielm_max),sendcount,MPI_INTEGER, mesh_nodes,recvcounts,displs,MPI_INTEGER, rowcomm,ierr)
      sendcount = nodes_per_element*3*nelm_row(rowrank + 1)
      recvcounts = nodes_per_element*3*nelm_row
      displs = nodes_per_element*3*displs_elm
      call MPI_ALLGATHERV(mesh_coord(:,:,ielm_min:ielm_max),sendcount,MPI_DOUBLE_PRECISION, mesh_coord,recvcounts,displs,MPI_DOUBLE_PRECISION, rowcomm,ierr) 

      sendcount = nelm_row(rowrank + 1)
      recvcounts = nelm_row
      displs = displs_elm
      call MPI_ALLGATHERV(mesh_zone(ielm_min:ielm_max),sendcount,MPI_INTEGER, mesh_zone,recvcounts,displs,MPI_INTEGER, rowcomm,ierr)

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

   npar=nelms_global*npoints

   !Allocate local storage for particle data
   disp_unit = 1
   arraysize = 0
   if (hostrank == 0) then
      arraysize = npar
      arraysize = arraysize*sizeof(dpar)
   end if
   CALL MPI_Win_allocate_shared(arraysize, disp_unit, MPI_INFO_NULL, &
                                hostcomm, baseptr, win_pdata, ierr)
   if (hostrank /= 0) CALL MPI_Win_shared_query(win_pdata, 0, arraysize, disp_unit, baseptr, ierr)
   CALL C_F_POINTER(baseptr, pdata, [npar])

   !Set up 'neighborlist' table of element neighbors for ensemble tracking
   call find_element_neighbors

   dt_particle = dt
   t0_norm_particle = t0_norm
   v0_norm_particle = v0_norm
   b0_norm_particle = b0_norm
   rfac_particle = rfac
   linear_particle = linear
   toroidal_period_particle = toroidal_period
   if (cre==0) cre=3.0e8/(v0_norm/100.)
   !psubsteps = 16*dt*cre
   psubsteps = ra_cyc
   linear_particle=0

   ! if (lrestart) then
   !    write (part_file_name, '("ions_",I4.4,".h5")') times_output
   !    write (0, *) part_file_name
   !    call hdf5_read_particles(part_file_name, ierr)
   !    if (myrank .eq. 0) print *, 'read_parts returned with ierr=', ierr
   ! else
       allocate (pdata_local(npar/maxrank*10))
      !First pass: assign particles to processors, elements
      locparts = 0
      locparts2 = 0
      ipart_begin=0; ipart_end=0

      !allocate(ran(npar*5))
      !call random_number(ran)
      !call random_seed(size=ran_size)
      !allocate (seed(ran_size))
      !call random_seed(get=seed)
      !seed = seed + myrank
      !call random_seed(put=seed)

      do ielm = 1, nelms
         call get_zone(ielm, izone)

         call define_element_quadrature(ielm,int_pts_main,int_pts_tor)
         call define_fields(ielm,FIELD_PSI,1,0)
#ifdef USE3D
         call m3dc1_ent_getglobalid(3, ielm - 1, itri)
#else
         call m3dc1_ent_getglobalid(2, ielm - 1, itri)
#endif
         itri = itri + 1
         !x_min = minval(mesh_coord(1, :, itri))
         !x_max = maxval(mesh_coord(1, :, itri))
         !phi_min = mesh_coord(2, 1, itri)
         !phi_max = mesh_coord(2, 4, itri)
         !if (phi_max == 0) phi_max = toroidal_period
         !z_min = minval(mesh_coord(3, :, itri))
         !z_max = maxval(mesh_coord(3, :, itri))

         !!area=x1*(z2-z3)+x2*(z3-z1)+x3*(z1-z2)
         !area=tri_area(mesh_coord(1,1,itri),mesh_coord(3,1,itri),mesh_coord(1,2,itri),mesh_coord(3,2,itri),mesh_coord(1,3,itri),mesh_coord(3,3,itri))
         !!hou!npar_local=int(npar*area*mesh_coord(1,1,itri)*0.03/nplanes)
         !dpar%x = (mesh_coord(:, 1, itri) + mesh_coord(:, 2, itri) + mesh_coord(:, 3, itri))/3.
         !dpar%v = 0.
         !B0 = 1.
         !call get_geom_terms(dpar%x, itri, geomterms, .false., ierr)
!#ifdef USEST
         !call update_geom_terms_st(geomterms, elfieldcoefs(itri), .false.)
!#endif
         !call evalf0(dpar%x, elfieldcoefs(itri), geomterms, f_mesh, T00, rho, gradcoef)
!#ifdef USEST
         !npar_local=int(npar*(x_max-x_min)*(z_max-z_min)*f_mesh*2/nplanes*38)!fullf
!#else
         !npar_local = int(npar*(x_max - x_min)*(x_max + x_min)/2.*(z_max - z_min)*f_mesh*2/nplanes*6.5)!fullf
!#endif
         do ipar = 1, npoints !Loop over z positions
            !call random_number(ran_temp)
!#ifdef USEST
            !x_temp = ran_temp*(x_max-x_min)+x_min
!#else
            !x_temp = sqrt(ran_temp*(x_max**2 - x_min**2) + x_min**2)
!#endif
            !!dpar%x(1) = x_temp*(x2-x1)+x1
            !call random_number(ran_temp)
            !z_temp = z_min + ran_temp*(z_max - z_min)
            !if ((tri_area(mesh_coord(1, 1, itri), mesh_coord(3, 1, itri), mesh_coord(1, 2, itri), mesh_coord(3, 2, itri), x_temp, z_temp) < 0) &
            !.or.(tri_area(x_temp,z_temp,mesh_coord(1,2,itri),mesh_coord(3,2,itri),mesh_coord(1,3,itri),mesh_coord(3,3,itri))<0) &
            !.or.(tri_area(mesh_coord(1,1,itri),mesh_coord(3,1,itri), x_temp,z_temp,mesh_coord(1,3,itri),mesh_coord(3,3,itri))<0)) cycle
            if (izone.eq.ZONE_PLASMA) then
               dpar%x(1) = x_79(ipar)
               dpar%x(3) = z_79(ipar)
               dpar%x(2) = phi_79(ipar)
               dpar%jel = itri
               dpar%v(1) =  -v0_norm/100*cre
               dpar%v(2) = 0.
               !pdata(ielm)%ion(ip)%v(1) = 100000.                        !v_parallel

               dpar%wt = 0.
               dpar%f0 = 0.
               locparts2 = (itri-1)*npoints+ipar
               dpar%gid = locparts2
               locparts = locparts + 1
               pdata_local(locparts) = dpar
            endif
         end do !iz
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
      if (hostrank == 0) then
         write(0,*) ipart_begin, ipart_end
         do  ipart = ipart_begin, ipart_end
            itri=int((pdata(ipart)%gid-1)/npoints)+1
            ipar=mod((pdata(ipart)%gid-1),npoints)+1
            particle_map(ipar,itri)=ipart
         end do
      end if
      call mpi_barrier(mpi_comm_world, ierr)

   ! if (hostrank == 0) then
   !    call delete_particle(.true.)
   ! end if
#ifdef _OPENACC
   if (hostrank < num_devices) then
      write (0, *) 'num_device', num_devices
!$acc enter data create(mesh_coord,neighborlist,mesh_zone)
!!$acc enter data create(pdata(starty:endy)) async(blocky)
!$acc update device(mesh_coord,neighborlist,mesh_zone)
!!$acc update device(pdata(starty:endy)) async(blocky)
!$acc update device(dt_particle,t0_norm_particle,v0_norm_particle,b0_norm_particle,rfac_particle,linear_particle,toroidal_period_particle,psubsteps)
   end if !hostrank
#endif

   call create_field(nre_vec)
   call set_matrix_index(diff2_mat, 75)
   call create_mat(diff2_mat, 1, 1, icomplex, 1)
   !call create_newvar_matrix(diff_mat, NV_NOBOUND, NV_I_MATRIX, 1)
   do itri=1,nelms
      call define_element_quadrature(itri, int_pts_main, int_pts_tor)
      call define_fields(itri,0,1,0)
      tempxx = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
      tempxx = tempxx + 0.e-9*(intxx2(mu79(:,:,OP_DZZ),nu79(:,:,OP_DZZ))+intxx2(mu79(:,:,OP_DRR),nu79(:,:,OP_DRR)))
#ifdef USE3D
      tempxx = tempxx + 0.e-9*intxx3(mu79(:,:,OP_DPP),nu79(:,:,OP_DPP),ri4_79)
#endif
      call insert_block(diff2_mat,itri,1,1,tempxx,MAT_ADD)
    enddo
   call finalize(diff2_mat)
end subroutine init_particles

subroutine advance_particles(tinc)
   use basic  !For MPI variables
   use omp_lib
   implicit none
   include 'mpif.h'

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
   integer :: i, j, temp
   integer, dimension(1024) :: host_seq
   real :: rand_val
   type(particle), dimension(:), allocatable :: pdata_local  !Particle arrays

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
   do istep = 1, psubsteps
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
!$omp parallel do default(none) shared(pdata,tinc,nparticles,starty,endy,psubsteps,istep) PRIVATE(ierr)
#endif
!$acc parallel loop present(pdata(starty:endy))
      do ipart = starty, endy
         if (pdata(ipart)%deleted) cycle  !Skip if element is empty
         !call rk4(pdata(ielm)%ion(ipart), tinc, itri, ierr)
         call rk4(pdata(ipart), tinc, istep .eq. psubsteps, ierr)
         if (ierr .eq. 1) then ! Particle exited local+ghost domain -> lost
            pdata(ipart)%deleted = .true.
            cycle !Break out of tinc loop, go to next particle.
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

   !ierr = 0
   hh = 0.5*dt
   !B0=part%B0
   itri = part%jel
   ! if (ierr .eq. 1) return
   !1st step
   call fdot(part%x, part%v, part%wt, k1, l1, m1, n1, itri, ierr)
   if (ierr .eq. 1) return
   y1 = part%x + hh*k1; z1 = part%v + hh*l1; w1 = part%wt + hh*m1
   !write(0,*) k1(1),k1(2),k1(3)

   !2nd step
   call fdot(y1, z1, w1, k2, l2, m2, n2, itri, ierr)
   if (ierr .eq. 1) return
   y1 = part%x + hh*k2; z1 = part%v + hh*l2; w1 = part%wt + hh*m2

   !3rd step
   call fdot(y1, z1, w1, k3, l3, m3, n3, itri, ierr)
   if (ierr .eq. 1) return
   y1 = part%x + dt*k3; z1 = part%v + dt*l3; w1 = part%wt + dt*m3

   !4th step
   call fdot(y1, z1, w1, k4, l4, m4, n4, itri, ierr)
   if (ierr .eq. 1) return
   part%x = part%x + onethird*dt*(k2 + k3 + 0.5*(k1 + k4))
   part%wt = part%wt + onethird*dt*(m2 + m3 + 0.5*(m1 + m4))
   !if (part%x(1)<0.9) then
   !   ierr=1
   !   return
   !endif

   if (last_step) then
      !Determine final particle element location
      call get_geom_terms(part%x, itri, geomterms, .false., ierr)
      if (ierr.eq.1) return
      part%f0=dot_product(geomterms%g, elfieldcoefs(itri)%jrev1)
      if (linear_particle.eq.0) part%f0=part%f0+dot_product(geomterms%g, elfieldcoefs(itri)%jrev0)
      part%jel = itri
   endif
end subroutine rk4

subroutine fdot(x, v, w, dxdt, dvdt, dwdt, dEpdt, itri, ierr)
!$acc routine seq
   use basic
   implicit none

   real, dimension(3), intent(inout)                          :: x
   !real, dimension(3)                             :: x
   real, dimension(3), intent(out)                            :: dxdt
   real, dimension(vspdims), intent(in)                       :: v
   real, dimension(vspdims), intent(out)                      :: dvdt
   real, intent(in)                                           :: w
   real, intent(out)                                          :: dwdt
   real, intent(out)                                          :: dEpdt
   !type(elfield), dimension(nneighbors+1), target, intent(in) :: fh
   integer, intent(inout)                                     :: itri
   integer, intent(out)                                       :: ierr
   !integer, intent(in)                                        :: gid
   !real, intent(out)                                        :: deltaB0
   !real, intent(in)                                        :: df0de, df0dpsi
   real, parameter :: g_mks = 9.8067 ! earth avg surf grav accel in m/s/s
   type(elfield), target  :: fh_hop
   real, dimension(3) :: lr, lr2
   integer :: itri2
   type(xgeomterms)   :: geomterms, geomterms2
   real, dimension(3) :: x2, B_cyl, B0_cyl, B_cyl2, Jcyl, BxgrdB, deltaB, deltaB2, E_cyl, E_cyl2, bhat, bhat0, svec, svec0, Bstar, Bstar0
   real, dimension(vspdims)                                  :: v2, vs, vu
   real, dimension(3) :: dBdR, dBdphi, dBdz, dB0dR, dB0dphi, dB0dz, dB1dR, dB1dphi, dB1dz
   real, dimension(3) :: gradB0, gradB, gradB02, gradB1, gradB12, dEdR, dEdphi, dEdz
   real, dimension(3) :: weqv0, weqv1, weqvD, weqvD1, gradpsi, gradf, gradpe
   vectype, dimension(3) :: temp
   real f0, T0, tmp1, tmp2, tmp3, tmp4, tmp5, df0de, df0dxi, spd, gradcoef, dB1, dB12, j0xb, ne0, te0, dBdt, dEdt, dxidt
   real :: Rinv, B0inv, Binv, Bss, Bss0
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
   call get_geom_terms(x, itri, geomterms, .false., ierr)
   if (ierr .eq. 0) then
      if (mesh_zone(itri).ne.ZONE_PLASMA) ierr=1
   endif
   if (ierr .ne. 0) then
      !write(0,*) 'outoutout'
      return
   end if

   !if (itor.eq.1) Rinv = 1.0/x(1)
   Rinv = 1.0/x(1)

   !Calculate time derivatives
   !call getBcylprime(x, elfieldcoefs(itri), geomterms, B0_cyl, deltaB, dB0dR, dB0dphi, dB0dz, dB1dR, dB1dphi, dB1dz)
   call getBcyl(x, elfieldcoefs(itri), geomterms, B0_cyl, deltaB)
   B_cyl = B0_cyl + deltaB
   !call getBcyl(x, elfieldcoefs(itri), geomterms, B_cyl, deltaB, gradB0, gradB1, dB1)
   !write(0,*) B_cyl(1), B_cyl(1), B_cyl(2)
   B0inv = 1.0/sqrt(dot_product(B0_cyl, B0_cyl))  !1/magnitude of B
   bhat0 = B0_cyl*B0inv                         !Unit vector in b direction
   Binv = 1.0/sqrt(dot_product(B_cyl, B_cyl))  !1/magnitude of B
   bhat = B_cyl*Binv                         !Unit vector in b direction

      !call getEcylprime(x, fhptr, geomterms, E_cyl, dEdR, dEdphi, dEdz)
      call getEcyl(x, elfieldcoefs(itri), geomterms, E_cyl)
      !call getBcyl_last(x, fhptr, geomterms, B_cyl2, deltaB_last)
      !temp(1) = dot_product(geomterms%dr, fhptr%pe)
      !temp(3) = dot_product(geomterms%dz, fhptr%pe)
!#ifdef USECOMPLEX
      !temp(2) = dot_product(geomterms%g, fhptr%pe)*rfac_particle/x(1)
      !gradpe=real(temp * exp(rfac_particle*x(2)))
!#else
      !temp(2) = dot_product(geomterms%dphi, fhptr%pe)/x(1)
      !gradpe=temp
!#endif
      !temp(1) = dot_product(geomterms%dr, fhptr%pe0)
      !temp(3) = dot_product(geomterms%dz, fhptr%pe0)
      !temp(2) = 0.
      !j0xb=-dot_product(temp,deltaB)*B0inv
      !if (real(dot_product(geomterms%g, fhptr%psiv0))<0.21) then
      !   gradpe=0.
      !   j0xb=0.

   Bstar0 = B0_cyl
   !Bstar = B_cyl
   Bss0 = dot_product(Bstar0, bhat0)

   svec0 = 0.
   !svec = svec - E_cyl  ! - g_mks/qm_ion
   !svec = v(2)*gradB0  ! - g_mks/qm_ion

   dxdt0(1) = (v(1)*Bstar0(1) + bhat0(2)*svec0(3) - bhat0(3)*svec0(2))/Bss0
   dxdt0(2) = (v(1)*Bstar0(2) + bhat0(3)*svec0(1) - bhat0(1)*svec0(3))/Bss0
   dxdt0(3) = (v(1)*Bstar0(3) + bhat0(1)*svec0(2) - bhat0(2)*svec0(1))/Bss0

   dvdt0(1) = 0.
   !dvdt0(1) = -qm_ion*dot_product(bhat0, svec)
   !dvdt(1) = 0
   dvdt0(2) = 0. !magnetic moment is conserved.

   dxdt0(2) = Rinv*dxdt0(2)  !phi-dot = (v_phi / R) for cylindrical case

   !call getBcylprime(x, elfieldcoefs(itri), geomterms, B_cyl, deltaB, dBdR, dBdphi, dBdz, .false.)
   !call getBcyl(x, fhptr, geomterms, B_cyl, deltaB, gradB0, gradB1,1)

   ! Gradient of B0 = grad(B.B)/(2 B0) = (B . grad B)/B0

   Bstar = B_cyl
   !Bstar = B_cyl
   Bss = dot_product(Bstar, bhat)

   svec = 0.
   svec = svec - E_cyl  ! - g_mks/qm_ion

   if (linear_particle .eq. 1) then
      dxdt(1) = (v(1)*Bstar(1) + bhat0(2)*svec(3) - bhat0(3)*svec(2))/Bss0
      dxdt(2) = (v(1)*Bstar(2) + bhat0(3)*svec(1) - bhat0(1)*svec(3))/Bss0
      dxdt(3) = (v(1)*Bstar(3) + bhat0(1)*svec(2) - bhat0(2)*svec(1))/Bss0

      dvdt(1) = 0.
      dvdt(2) = 0. !magnetic moment is conserved.
   else
      dxdt(1) = (v(1)*Bstar(1) + bhat(2)*svec(3) - bhat(3)*svec(2))/Bss
      dxdt(2) = (v(1)*Bstar(2) + bhat(3)*svec(1) - bhat(1)*svec(3))/Bss
      dxdt(3) = (v(1)*Bstar(3) + bhat(1)*svec(2) - bhat(2)*svec(1))/Bss

      dvdt(1) = 0.
      dvdt(2) = 0. !magnetic moment is conserved.
   end if

   dxdt(2) = Rinv*dxdt(2)  !phi-dot = (v_phi / R) for cylindrical case

   !Weights evolve in delta-f method only.
   ! V1 = (ExB)/(B**2) + U deltaB/B
   !weqv1(1) = ((E_cyl(2)*B_cyl(3) - E_cyl(3)*B_cyl(2))*B0inv + v(1)*deltaB(1))*B0inv
   !weqv1(2) = ((E_cyl(3)*B_cyl(1) - E_cyl(1)*B_cyl(3))*B0inv + v(1)*deltaB(2))*B0inv
   !weqv1(3) = ((E_cyl(1)*B_cyl(2) - E_cyl(2)*B_cyl(1))*B0inv + v(1)*deltaB(3))*B0inv
   !spd = sqrt(v(1)*v(1) + 2.0*qm_ion*v(2)/B0inv)
   weqv1 = dxdt - dxdt0
   !weqa1 = dvdt - dvdt0
   !dBdt = dot_product(weqv1, gradB0)
   !dEdt = m_ion*v(1)*weqa1(1) + q_ion*v(2)*dBdt
   !dxidt = weqa1(1)/spd-v(1)/spd**2*(dEdt/m_ion/spd)

   !! vD = (1/(e B**3))(M_i U**2 + mu B)(B x grad B) + ((M_i U**2)/(eB**2))*J_perp
   !tmp1 = (v(1)*v(1))*(B0inv*B0inv)/qm_ion
   !tmp2 = tmp1*B0inv + v(2)*(B0inv*B0inv)
   !weqvD = tmp2*BxgrdB + tmp1*(Jcyl - dot_product(bhat, Jcyl)*bhat)
   !weqv0 = v(1)*bhat + weqvD
   !tmp2 = (v(4)*v(4))*(B0inv*B0inv)/qm_ion*B0inv
   !tmp2 = tmp1*B0inv
   !weqvD1 = tmp2*BxgrdB + tmp1*(Jcyl - dot_product(bhat, Jcyl)*bhat)

   gradpsi = 0.0
   gradpsi(1) = dot_product(elfieldcoefs(itri)%jrev0, geomterms%dr)
   gradpsi(3) = dot_product(elfieldcoefs(itri)%jrev0, geomterms%dz)
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

   !ne0 = dot_product(geomterms%g, elfieldcoefs(itri)%ne0)
   ! write(0,*) tmp1,(dot_product(-gradpe,bhat)+j0xb)
   !   gradcoef=0.
   !   df0de=0.
   !   df0dxi=0.
   !endif

   ! vD = (1/(e B**3))(M_i U**2 + mu B)(B x grad B) + ((M_i U**2)/(eB**2))*J_perp
   !tmp3 = tmp2*dot_product(gradB0, deltaB)
   !!tmp4 = tmp3 - v(1)*dot_product(Jcyl,E_cyl)/qm_ion*(B0inv*B0inv*B0inv)
   !tmp4 = tmp3
   !tmp5 = -tmp4*dot_product(elfieldcoefs(itri)%Bzv0, geomterms%g)*gradcoef

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
   dwdt = (-1.0)*(1.*dot_product(gradpsi,weqv1))
   !if (deltaB_last(2).ne.0.) write(0,*) deltaB(2),deltaB_last(2)
   !dEpdt = q_ion*(dot_product(weqvD, E_cyl)+v(2)*dB0dt)
   !write(0,*) dt_particle, t0_norm_particle
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
   if (linear_particle .eq. 1) then
      dxdt = dxdt0
      dvdt = dvdt0
   end if

   !dwdt = (1.-w)/(dt_particle*t0_norm_particle)
   !B0=1./B0inv
   !deltaB0=dot_product(deltaB,bhat)
   !deltaB0=f0
end subroutine fdot
!---------------------------------------------------------------------------
subroutine runaway_advection_step(pdt)
   use basic
   use diagnostics
   use auxiliary_fields
   implicit none
   include 'mpif.h'

   real, intent(in) :: pdt

   real    :: tstart, tend
   integer :: istep, ierr, ipart

   !Advance particle positions
   call calculate_electric_fields(linear)
   call get_field_coefs(0)
   call mpi_barrier(mpi_comm_world, ierr)
   !call MPI_Win_fence(0, win_elfieldcoefs)
#ifdef _OPENACC
   if (hostrank < num_devices) then
#else
      if (hostrank < ncols) then
#endif
         call second(tstart)
         call advance_particles(pdt/psubsteps)
         call second(tend)
         write (0, '(A,I7,A,f9.2,A)') 'Particle advance completed', psubsteps, ' steps in', &
            tend - tstart, ' seconds.'
      end if
      call mpi_barrier(mpi_comm_world, ierr)
      !Compute particle pressure tensor components
      call update_particle_pressure
      call mpi_barrier(mpi_comm_world, ierr)
end subroutine runaway_advection_step
!---------------------------------------------------------------------------
subroutine update_particle_pressure
   use basic
   use arrays
   use diagnostics
   implicit none
   include 'mpif.h'

   real    :: tstart, tend
   integer :: ierr

   call MPI_Bcast(ipart_begin, 1, MPI_INTEGER, 0, hostcomm, ierr)
   call MPI_Bcast(ipart_end, 1, MPI_INTEGER, 0, hostcomm, ierr)
   call MPI_Bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
   call mpi_barrier(mpi_comm_world, ierr)
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
         call calcavector(ielm, nre_field(0), elfieldcoefs(ielm_global)%jrev0)
         elfieldcoefs(ielm_global)%Bfpv1 = 0.
         elfieldcoefs(ielm_global)%psiv1 = 0.
         elfieldcoefs(ielm_global)%Bzv1 = 0.
         elfieldcoefs(ielm_global)%jrev1 = 0.
      end if
      call calcavector(ielm, bfp_field(1), elfieldcoefs(ielm_global)%Bfpv1)
      call calcavector(ielm, psi_field(1), elfieldcoefs(ielm_global)%psiv1)
      call calcavector(ielm, bz_field(1), elfieldcoefs(ielm_global)%Bzv1)
      call calcavector(ielm, nre_field(1), elfieldcoefs(ielm_global)%jrev1)
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
call MPI_ALLGATHERV(elfieldcoefs(ielm_min:ielm_max),sendcount,MPI_elfield, elfieldcoefs,recvcounts,displs,MPI_elfield, rowcomm,ierr)
      deallocate (recvcounts)
      deallocate (displs)
   end if    !call m3dc1_ent_isghost(2, ielm-1, isghost)
   call mpi_barrier(mpi_comm_world, ierr)
end subroutine get_field_coefs

!---------------------------------------------------------------------------
subroutine getBcyl(x, fh, gh, Bcyl, deltaB)
!$acc routine seq
   use basic
   implicit none

   real, dimension(3), intent(in) :: x       !Position
   type(elfield), intent(in) :: fh           !Field handle
   type(xgeomterms), intent(in) :: gh        !Geometric terms handle
   real, dimension(3), intent(out) :: Bcyl   !Output total magnetic field
   real, dimension(3), intent(out) :: deltaB !Output perturbed part of Bcyl

   vectype, dimension(3) :: temp
   real :: Rinv
!!$OMP THREADPRIVATE(Rinv)

   !if (itor.eq.1) Rinv = 1.0/x(1)
   Rinv = 1.0/x(1)
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
   if (.not. (linear_particle .eq. 1)) Bcyl = Bcyl + deltaB

   Bcyl = Bcyl*b0_norm_particle/1.e4
   deltaB = deltaB*b0_norm_particle/1.e4
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
   Rinv = 1.0/x(1)
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

   !tempphi(1) = tempphi(1) - dot_product(fh%Bfpv1, gh%drphi)
   !tempphi(3) = tempphi(3) - dot_product(fh%Bfpv1, gh%dzphi)
   !dB1dphi = tempphi

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
! Integrate over elements to compute kinetic ion contributions to RHS vectors
! for parallel and perpendicular components of pressure tensor.
  subroutine particle_pressure_rhs
    use basic
    use arrays
    use math
    use m3dc1_nint
    implicit none
    !include 'mpif.h'
    intrinsic matmul

    real, dimension(dofs_per_element,coeffs_per_element) :: cl
    real, dimension(dofs_per_element) :: wnuhere, deltaBhere
    vectype, dimension(dofs_per_element) :: dofspa, dofspe
!#ifdef USECOMPLEX
    complex, dimension(dofs_per_element) :: dofspan, dofspen
    complex phfac
!#endif
    real, dimension(3) :: B_part, deltaB
    real, dimension(vspdims) :: vperp
    !type(elfield), dimension(nneighbors+1) :: elcoefs
    type(xgeomterms) :: geomterms
    real             :: B0, vpar, ppar, pperp
    integer          :: i,ierr, ielm, ielm_local, ielm_global, ipart, part_index, itri, tridex, isghost
    !integer          :: ibp, iwe, iok
    type(element_data) :: eldat
    integer :: ielm2
    integer, dimension(nodes_per_element) :: nodeids, nodeids2
    real :: x1, phi1, z1, x2, phi2, z2, x3, phi3, z3, area
    real :: x21, phi21, z21, x22, phi22, z22, x23, phi23, z23
    integer :: k
    real :: f0, df0de, newphi
    real, dimension(3) :: gradf0
    real, dimension(201) :: zgrid, phigrid, zloss, philoss, zloss_out, philoss_out
    real :: ratio
    integer :: igrid, izone
    !nelms = size(pdata)
    !nelms = local_elements()
    !elcoefs(:)%itri = 0

    !npart_elm(:)=0
    !do ipart=1,nparticles
    !   ielm=pdata(ipart)%jel
    !   if (localmeshid(ielm)>0) then
    !      ielm_local=localmeshid(ielm)
          !write(0,*) ielm_local
    !      npart_elm(ielm_local)=npart_elm(ielm_local)+1
    !      part_elm(npart_elm(ielm_local),ielm_local)=ipart
    !   endif
    !enddo
    nre_vec=0.

     do ielm=1,nelms
       !call m3dc1_ent_isghost(2, ielm-1, isghost)
       !if(isghost.eq.1) then
          !!if (pdata(ielm)%np.gt.0) print *,myrank,': nonzero ghost particle count!'
          !cycle
       !endif
       
       !if (pdata(ielm)%np.lt.1) cycle !If no particles, then no pressure
      call get_zone(ielm, izone)
      if (izone.ne.ZONE_PLASMA) cycle

      dofspa = 0.

      call define_element_quadrature(ielm,int_pts_main,int_pts_tor)
      call define_fields(ielm,FIELD_PSI+FIELD_RE,1,0)

      !call get_element_nodes(ielm, nodeids)
    
      !call get_node_pos(nodeids(1), x1, phi1, z1)
      !call get_node_pos(nodeids(2), x2, phi2, z2)
      !call get_node_pos(nodeids(3), x3, phi3, z3)
    
  !call define_basis(ielm)

 !Sum over particles within this element
       !ibp = 0;  iwe = 0;  iok = 0
       !write(0,*) ielm
       !do ielm2=1,nelms
       !if (pdata(ielm2)%np.lt.1) cycle !If no particles, then no pressure

       !call get_element_nodes(ielm2, nodeids2)
       !call get_node_pos(nodeids2(1), x21, phi21, z21)
       !call get_node_pos(nodeids2(2), x22, phi22, z22)
       !call get_node_pos(nodeids2(3), x23, phi23, z23)
       !if (.not.((nodeids(1)==nodeids2(1)).or.(nodeids(1)==nodeids2(2)).or.(nodeids(1)==nodeids2(3)) &
          !.or.(nodeids(2)==nodeids2(1)).or.(nodeids(2)==nodeids2(2)).or.(nodeids(2)==nodeids2(3)) &
       !   .or.(nodeids(3)==nodeids2(1)).or.(nodeids(3)==nodeids2(2)).or.(nodeids(3)==nodeids2(3)))) cycle !If no particles, then no pressure
       !Need B at particle locations -> Load scalar fields for this element
              !Get basis function polynomial expansions
       !if (.not.(((x21-x1)**2+(z21-z1)**2<1e-2).or.((x22-x1)**2+(z22-z1)**2<1e-2).or.((x23-x1)**2+(z23-z1)**2<1e-2) &
       !   .or.((x21-x2)**2+(z21-z2)**2<1e-2).or.((x22-x2)**2+(z22-z2)**2<1e-2).or.((x23-x2)**2+(z23-z2)**2<1e-2) &
       !   .or.((x21-x3)**2+(z21-z3)**2<1e-2).or.((x22-x3)**2+(z22-z3)**2<1e-2).or.((x23-x3)**2+(z23-z3)**2<1e-2))) cycle !If no particles, then no pressure
       if (iprecompute_metric.eq.1) then
          cl = ctri(:,:,ielm)
       else
          call local_coeff_vector(ielm, cl)
       endif

       !Need B at particle locations -> Load scalar fields for this element
       !call get_field_coefs()

       do part_index=1,npoints
#ifdef USE3D
          call m3dc1_ent_getglobalid(3,ielm-1,ielm_global)
#else
          call m3dc1_ent_getglobalid(2,ielm-1,ielm_global)
#endif
          ielm_global=ielm_global+1
          !ipart=npoints*(ielm_global-1)+part_index
          ipart=particle_map(part_index,ielm_global)
          if (pdata(ipart)%gid.ne.npoints*(ielm_global-1)+part_index) write(0,*) "Wrong"
          !if (pdata(ipart)%deleted) then
             !igrid=pdata(ipart)%x(2)/(twopi/200)+1
             !ratio=(pdata(ipart)%x(2)-phigrid(igrid))/(phigrid(igrid+1)-phigrid(igrid))
             !philoss(igrid)=philoss(igrid)+nre179(part_index,OP_1)*weight_79(part_index)*(1-ratio)
             !if (igrid<200) then
                 !philoss(igrid+1)=philoss(igrid+1)+nre179(part_index,OP_1)*weight_79(part_index)*ratio
             !else
                 !philoss(1)=philoss(1)+nre179(part_index,OP_1)*weight_79(part_index)*ratio
             !endif
             !igrid=(pdata(ipart)%x(3)+1.5)/(3./200)+1
             !ratio=(pdata(ipart)%x(3)-zgrid(igrid))/(zgrid(igrid+1)-zgrid(igrid))
             !zloss(igrid)=zloss(igrid)+nre179(part_index,OP_1)*weight_79(part_index)*(1-ratio)
             !if (igrid<200) then
                 !zloss(igrid+1)=zloss(igrid+1)+nre179(part_index,OP_1)*weight_79(part_index)*ratio
             !endif
          !endif
          newphi=pdata(ipart)%x(2)
          pdata(ipart)%x(1)=x_79(part_index)
          pdata(ipart)%x(3)=z_79(part_index)
          pdata(ipart)%x(2)=phi_79(part_index)
          pdata(ipart)%jel=ielm_global
          itri = ielm_global
          !call get_geom_terms(pdata(ipart)%x, itri,  &
          !     geomterms, .false., ierr)
          !if (ierr.ne.0) then
             !print *,myrank,': Bad particle in pressure tensor integral; skipping.'
             !write(0,*) itri,pdata(ipart)%x
             !!ibp = ibp + 1
             !cycle !next particle
          !endif
          if (itri.ne.ielm_global) then
             print *,myrank,': Particle in wrong element in pressure tensor integral:',pdata(ipart)%gid
             print *,ipart,': ',itri,'.ne.',ielm
             !print *,myrank,': pcoord ',pdata(ielm)%ion(ipart)%gid,' = ',pdata(ielm)%ion(ipart)%x
             !iwe = iwe + 1
             cycle !next particle
          endif
          !iok = iok + 1
          !call getBcyl(pdata(ipart)%x, elfieldcoefs(ielm_global), geomterms, B_part, &
          !     deltaB)
          !B0 = sqrt(dot_product(B_part(1:3:2), B_part(1:3:2)))
          !pdata(ipart)%dt=pdata(ipart)%dt*0.93247923
          !if (pdata(ipart)%dt<dt*t0_norm/psubsteps*0.5) pdata(ipart)%dt=pdata(ipart)%dt*1.9932235932
          !pdata(ipart)%dt=dt*t0_norm/psubsteps/B0*0.25
          !if (pdata(ipart)%dt>dt*t0_norm/psubsteps*2) pdata(ipart)%dt=dt*t0_norm/psubsteps*2
          if (pdata(ipart)%deleted) then
             !pdata(ipart)%deleted=real((pst79(part_index,OP_1)-psimin)/(psibound-psimin))>1
             pdata(ipart)%deleted=.false.
             pdata(ipart)%wt=0.
             pdata(ipart)%f0=0.
             nre179(part_index,OP_1) = 0
             cycle
          endif
          if (pdata(ipart)%jel/=ielm_global) cycle
          !Calculate B field at particle location         
          !B0 = pdata(ipart)%B0
          !!Use B and v to get parallel and perp components of particle velocity
          !if (vspdims.eq.2) then ! drift-kinetic: v_|| = v(1),  mu = q * v(2)
             !vpar = pdata(ipart)%v(1)
             !pperp = q_ion*pdata(ipart)%v(2)*B0
          !else !full orbit: v_|| = v.B/|B|,  v_perp = v - v_||
             !if (B0.gt.0.0) then !non-degenerate
                !vpar = dot_product(pdata(ipart)%v, B_part(1:vspdims))/B0
                !vperp = pdata(ipart)%v - (vpar/B0)*B_part(1:vspdims)
             !else !degenerate case: no B field, pressure is scalar
                !vpar = 0.0
                !vperp = pdata(ipart)%v
             !endif !degenerate?
             !pperp = 0.5 * m_ion * dot_product(vperp, vperp)
          !endif !full-orbit?
          !ppar = m_ion * vpar**2
          !call evalf0(pdata(ielm)%ion(ipart)%x, pdata(ielm)%ion(ipart)%v, B0, elcoefs(tridex), geomterms, f0, gradf0, df0de)

          !temp79a = exp(-((x_79-pdata(ielm2)%ion(ipart)%x(1))**2+(z_79-pdata(ielm2)%ion(ipart)%x(3))**2)/1.e-3)
          !temp79a = exp(-((x_79-1.1)**2+(z_79-0.)**2)/1.e-4)
          !temp79a = temp79a/int1(temp79a)
          !temp79a = temp79a*ri_79

          !temp79a=3*((x_79-x1)*(z2-z1)-(z_79-z1)*(x2-x1))/((pdata(ielm)%ion(ipart)%x(1)-x1)*(z2-z1)-(pdata(ielm)%ion(ipart)%x(3)-z1)*(x2-x1))
          !temp79b=3*((x_79-x2)*(z3-z2)-(z_79-z2)*(x3-x2))/((pdata(ielm)%ion(ipart)%x(1)-x2)*(z3-z2)-(pdata(ielm)%ion(ipart)%x(3)-z2)*(x3-x2))
          !temp79c=3*((x_79-x3)*(z1-z3)-(z_79-z3)*(x1-x3))/((pdata(ielm)%ion(ipart)%x(1)-x3)*(z1-z3)-(pdata(ielm)%ion(ipart)%x(3)-z3)*(x1-x3))
          !where (real(temp79a)>real(temp79b)) temp79a=temp79b
          !where (real(temp79a)>real(temp79c)) temp79a=temp79c
          !call choosemin(temp79a,temp79b)
          !do k=1,npoints
          !   if (temp79a(k)>temp79b(k)) temp79a(k)=temp79b(k)
          !   temp79a(k)=temp79a(k)*temp79c(k)
          !enddo
          !pper79(:,OP_1) = pper79(:,OP_1)+pperp*temp79a*pdata(ielm2)%ion(ipart)%wt
          !ppar79(:,OP_1) = ppar79(:,OP_1)+ppar*temp79a*pdata(ielm2)%ion(ipart)%wt
          !phfac = exp(-rfac*pdata(ielm2)%ion(ipart)%x(2))
          !pper79(:,OP_1) = pper79(:,OP_1)+pperp*phfac*temp79a*pdata(ielm2)%ion(ipart)%wt
          !ppar79(:,OP_1) = ppar79(:,OP_1)+ppar*phfac*temp79a*pdata(ielm2)%ion(ipart)%wt
          !enddo !ipart
          !enddo
          !do k=1, npoints
              !temp79b=exp(-((x_79-x_79(k))**2+(z_79-z_79(k))**2)/1.e-4)
              !area=real(int1(temp79b))
              !ppar79(k,OP_1)=ppar79(k,OP_1)/area
              !pper79(k,OP_1)=pper79(k,OP_1)/area
          !enddo
          !pper79(:,OP_1) = pper79(:,OP_1)/temp79a
          !ppar79(:,OP_1) = ppar79(:,OP_1)/temp79a
          !Add particle contribution to RHS vector (should vectorize well).
#ifdef DELTA_F
          !wnuhere = pdata(ipart)%wt2 * matmul(cl, geomterms%g)
          !wnuhere = pdata(ipart)%wt * matmul(cl, geomterms%g) * pdata(ipart)%x(1)/10.
          !deltaBhere = pdata(ielm)%ion(ipart)%f0 *dot_product(B_part,deltaB)/B0**2* matmul(cl,geomterms%g)
#else
          !wnuhere = matmul(cl, geomterms%g)
#endif
#ifndef USECOMPLEX
          !if (pdata(ipart)%f0==0) write(0,*) "Wrong",pdata(ipart)%gid
          !write(0,*) x_79(part_index),pdata(ipart)%x(1)
          nre179(part_index,OP_1) = pdata(ipart)%wt+pdata(ipart)%f0
        !call get_geom_terms(pdata(ipart)%x, itri, geomterms, .false., ierr)
        !pdata(ipart)%f0=dot_product(geomterms%g, elfieldcoefs(itri)%jrev1)

!#ifd   dofspa = dofspa + ppar*wnuhere
          pdata(ipart)%wt = 0.
          pdata(ipart)%f0 = 0.
          !dofspa = dofspa + wnuhere !orbit
          !dofspa = dofspa + pdata(ipart)%x(1)*wnuhere
          !dofspe = dofspe + pperp*wnuhere
          !dofspa = intx2(mu79(:,:,OP_1),ppar79(:,OP_1))
          !dofspe = intx2(mu79(:,:,OP_1),pper79(:,OP_1))
#else
          !Extract appropriate Fourier component of particle contribution
          phfac = exp(rfac*newphi)
          nre179(part_index,OP_1) = pdata(ipart)%wt+phfac*pdata(ipart)%f0
          !nre179(part_index,OP_1) = pdata(ipart)%wt
          !nre179(part_index,OP_1) = pdata(ipart)%f0
          pdata(ipart)%wt = 0.
          pdata(ipart)%f0 = 0.
          !dofspan = dofspan + ppar*phfac*wnuhere
          !dofspen = dofspen + pperp*phfac*wnuhere
          !dofspen = dofspen + pperp*phfac*deltaBhere
          !dofspan = intx2(mu79(:,:,OP_1),ppar79(:,OP_1))
          !dofspen = intx2(mu79(:,:,OP_1),pper79(:,OP_1))
#endif
       enddo !ipart
#ifndef USECOMPLEX
       if (linear_particle.eq.0) nre179(:,OP_1)=nre179(:,OP_1)-nre079(:,OP_1)
#endif
          dofspa = intx2(mu79(:,:,OP_1),nre179(:,OP_1))
       !PRINT *,'DB',myrank,ielm,pdata(ielm)%np,ibp,iwe,iok

       !Insert element sums into field data
       ! Note: this is only correct if the local index ielm refers to the
       !  same element in meshes with and without ghost zone layers!
       call vector_insert_block(nre_vec%vec,  ielm, 1, dofspa, VEC_ADD)
       !call vector_insert_block(nre_vec%vec, ielm, 1, dofspe, VEC_ADD)

enddo !ielm
  !call mpi_allreduce(philoss, philoss_out, 201, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  !call mpi_allreduce(zloss, zloss_out, 201, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  !if (myrank.eq.0) write(0,*) "ZLost",zloss_out
  !if (myrank.eq.0) write(0,*) "PhiLost",philoss_out
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
   call sum_shared(nre_vec%vec)
   call newsolve(diff2_mat, nre_vec%vec, ierr)
   nre_field(1)=nre_vec
   ! call sum_shared(p_f_par%vec)
   ! call newsolve(diff2_mat, p_f_par%vec, ierr)
   ! call sum_shared(p_f_perp%vec)
   ! call newsolve(diff2_mat, p_f_perp%vec, ierr)
   ! call sum_shared(den_f_1%vec)
   ! call newsolve(diff2_mat, den_f_1%vec, ierr)
   ! call sum_shared(den_f_0%vec)
   ! call newsolve(diff2_mat, den_f_0%vec, ierr)
end subroutine solve_pi_tensor
end module runaway_advection
