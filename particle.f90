! Kinetic energetic ion module, J. Breslau, 2015
module particles
  use mesh_mod
  implicit none

  real, parameter :: e_mks = 1.6022e-19      !Elementary charge, in Coulombs
  real, parameter :: m_proton = 1.6726e-27   !Proton mass in kg
  real, parameter :: A_deuteron = 1.9990075  !Deuteron/proton mass ratio
  real, parameter :: A_alpha = 3.972599689   !alpha/proton mass ratio
  real, parameter :: qm_proton = e_mks/m_proton  !Proton charge/mass ratio in C/kg
  real, parameter :: vp1eV = 1.3841e+04      !Thermal speed of a 1 eV proton in m/s
  integer, parameter :: vspdims = 2          !Dimensions in velocity space (2=gk; 3=full)
#ifdef USE3D
  integer, parameter :: nneighbors = 5       !Max # of nearest neighbors of a prism
#else
  integer, parameter :: nneighbors = 3       !Max # of nearest neighbors of a triangle
#endif

  type elfield
     vectype, dimension(coeffs_per_element) :: psiv0, psiv1, Bzv0, Bzv1
     vectype, dimension(coeffs_per_element) :: Bfv
     vectype, dimension(coeffs_per_element) :: er, ephi, ez
     integer :: itri
  end type elfield

  type xgeomterms
     real, dimension(coeffs_per_element) :: g, dr, dz
     real, dimension(coeffs_per_element) :: drr, drz, dzz
#ifdef USE3D
     real, dimension(coeffs_per_element) :: dphi, drphi, dzphi
     real, dimension(coeffs_per_element) :: drrphi, drzphi, dzzphi
     real, dimension(coeffs_per_element) :: drphiphi, dzphiphi
#endif
  end type xgeomterms

  type particle
     real, dimension(3)       :: x           !Position in cylindrical coords
     real, dimension(vspdims) :: v           !Velocity
     real                     :: wt = 1.0    !Particle weighting in delta-f scheme
     real                     :: tlast       !Time
     integer                  :: gid         !Unique global particle index
     integer                  :: jel         !Predicted element of residence
  end type particle

  type elplist !Inventory of particles within a finite element
     integer :: np = 0
     type(particle), dimension(:), allocatable :: ion
  end type elplist

  type(elplist), dimension(:), allocatable :: pdata, jmppar  !Particle arrays
  type(particle), dimension(:), allocatable :: jinbuf        !Receive buffer for jumping particles
  real :: m_ion, q_ion, qm_ion, dt_ion
  integer, dimension(:,:), allocatable :: neighborlist, dnbr !Neighbor tracking arrays
  integer, dimension(:), allocatable :: dnlist               !Domain neighbor table
  integer :: nparticles, locparts, ndnbr
  integer :: mpi_particle !User-defined MPI datatype for particle communication

contains

#ifdef USEPARTICLES
!#define JBDEBUG

  !Define MPI datatype for particle communication
  subroutine define_mpi_particle(ierr)
    implicit none

    include 'mpif.h'

    integer, intent(out) :: ierr
    integer, parameter :: pnvars = 6
    integer, dimension(pnvars), parameter :: pblklen=(/3, vspdims, 1, 1, 1, 1/)
    integer(kind=MPI_ADDRESS_KIND), dimension(pnvars) :: pdspls
    integer, dimension(pnvars), parameter :: ptyps = (/MPI_DOUBLE, MPI_DOUBLE, &
         MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER, MPI_INTEGER/)

    type(particle) :: dum_par

    !Set up component displacements array
    call mpi_get_address(dum_par%x, pdspls(1), ierr)
    call mpi_get_address(dum_par%v, pdspls(2), ierr)
    pdspls(2) = pdspls(2) - pdspls(1)
    call mpi_get_address(dum_par%wt, pdspls(3), ierr)
    pdspls(3) = pdspls(3) - pdspls(1)
    call mpi_get_address(dum_par%tlast, pdspls(4), ierr)
    pdspls(4) = pdspls(4) - pdspls(1)
    call mpi_get_address(dum_par%gid, pdspls(5), ierr)
    pdspls(5) = pdspls(5) - pdspls(1)
    call mpi_get_address(dum_par%jel, pdspls(6), ierr)
    pdspls(6) = pdspls(6) - pdspls(1)
    pdspls(1) = 0

    call mpi_type_create_struct(pnvars, pblklen, pdspls, ptyps, mpi_particle, ierr)
    if (ierr.ne.0) return
    call mpi_type_commit(mpi_particle, ierr)
  end subroutine define_mpi_particle

!---------------------------------------------------------------------------
  subroutine particle_test
    use basic
    use diagnostics
    use auxiliary_fields
    implicit none

    include 'mpif.h'

    type(field_type) :: pi_parallel, pi_perp, pi_parallel_n, pi_perp_n
    character(len=32) :: line
    logical, parameter :: loutcart = .false. !Output Cartesian particle coords?
    integer, dimension(3), parameter :: trid = (/ 1979, 3405, 4856 /)
    real, parameter :: JpereV = 1.6022e-19
    real :: pdt, keeV, pphi, tstart, tend
    integer :: ierr, ip, istep=0, itr, trunit=120

    if (myrank.eq.0) then
       print *,'xlim2 = ',xlim2
       print *,'nonrect = ',nonrect
       print *,'iwall_is_limiter = ',iwall_is_limiter
       print *,'ifixedb = ',ifixedb
       print *,'GS magnetic axis at ',xmag,', ',zmag
       print *,'jadv = ',jadv
       print *,'imp_hyper = ',imp_hyper
       !print *,'xzero, zzero = ',xzero,', ',zzero
       print *,'rfac = ',rfac
    endif

    !Precompute electric field components (do it this way for testing only!)
    call calculate_auxiliary_fields(eqsubtract)

    !Particle pressure tensor components
    call create_field(pi_parallel);  call create_field(pi_perp)
    call create_field(pi_parallel_n);  call create_field(pi_perp_n)

    !Initialize particle population
    call second(tstart)
    call init_particles(ierr)
    if (ierr.ne.0) return
    if (myrank.eq.0) then
       call second(tend)
       write(0,'(I12,A,f9.2,A)')nparticles,' particles initialized in',&
            tend-tstart,' seconds.'
    endif

    !Test particle push
    !Create particle trajectory output text file?
    do itr=1,size(trid)
       do ierr=1,size(pdata)
          do ip=1,pdata(ierr)%np
             if (pdata(ierr)%ion(ip)%gid.eq.trid(itr)) then
                print *,myrank,': ielm,ip = ',ierr,ip,': gid = ',pdata(ierr)%ion(ip)%gid
                write(line,'(A,I8.8)') 'ptraj_',trid(itr)
                open(unit=trunit, file=trim(line), status='replace', action='write')
                if (loutcart) then !Cartesian
                   write(trunit,*)'x  y  z  t  KE  Pphi  rank'
                else !Cylindrical
                   write(trunit,*)'R  phi  z  t  KE  Pphi  rank'
                endif
                keeV = getke(pdata(ierr)%ion(ip))/JpereV
                pphi = getPphi(pdata(ierr)%ion(ip))/JpereV
                if (loutcart) then !Cartesian
                   write(trunit,'(3f14.8,3e18.8,I6)') &
                        pdata(ierr)%ion(ip)%x(1)*cos(pdata(ierr)%ion(ip)%x(2)), &
                        pdata(ierr)%ion(ip)%x(1)*sin(pdata(ierr)%ion(ip)%x(2)), &
                        pdata(ierr)%ion(ip)%x(3), pdt*istep, keeV, pphi, myrank
                else !Cylindrical
                   write(trunit,'(3f14.8,3e18.8,I6)') &
                        pdata(ierr)%ion(ip)%x(1), &
                        pdata(ierr)%ion(ip)%x(2), &
                        pdata(ierr)%ion(ip)%x(3), pdt*istep, keeV, pphi, myrank
                endif
                close(trunit)
                print *,myrank,': trajectory file ',trim(line),' created.'
                exit
             endif
          enddo !ip
       enddo !ierr
    enddo !itr

    !Advance particle positions
    call second(tstart)
    pdt = 3.0e-7
    do istep=1,1200
       call advance_particles(pdt)
       if (myrank.eq.0) print *,'particle advance',istep,' complete.'

       do itr=1,size(trid)
          do ierr=1,size(pdata)
             do ip=1,pdata(ierr)%np
                if (pdata(ierr)%ion(ip)%gid.eq.trid(itr)) then
                   write(line,'(A,I8.8)'),'ptraj_',trid(itr)
                   open(unit=trunit, file=trim(line), status='old', action='write', &
                        position='append')
                   keeV = getke(pdata(ierr)%ion(ip))/JpereV
                   pphi = getPphi(pdata(ierr)%ion(ip))/JpereV
                   if (loutcart) then !Cartesian
                      write(trunit,'(3f14.8,3e18.8,I6)') &
                           pdata(ierr)%ion(ip)%x(1)*cos(pdata(ierr)%ion(ip)%x(2)), &
                           pdata(ierr)%ion(ip)%x(1)*sin(pdata(ierr)%ion(ip)%x(2)), &
                           pdata(ierr)%ion(ip)%x(3), pdt*istep, keeV, pphi, myrank
                   else !Cylindrical
                      write(trunit,'(3f14.8,3e18.8,I6)') &
                           pdata(ierr)%ion(ip)%x(1), &
                           pdata(ierr)%ion(ip)%x(2), &
                           pdata(ierr)%ion(ip)%x(3), pdt*istep, keeV, pphi, myrank
                   endif
                   close(trunit)
                   exit
                endif
             enddo !ip
          enddo !ierr
       enddo !itr
    enddo !istep
    if (myrank.eq.0) then
       call second(tend)
       write(0,'(A,I7,A,f9.2,A)')'Particle advance completed',istep-1,' steps in',&
            tend-tstart,' seconds.'
    endif

    !Test particle pressure tensor component calculation
    call second(tstart)
    call particle_pressure(pi_parallel, pi_perp, pi_parallel_n, pi_perp_n)
    if (myrank.eq.0) then
       call second(tend)
       write(0,'(A,f9.2,A)')'Pressure tensor RHS vecs calculated in',tend-tstart,' sec.'
    endif

    !Return to non-ghosted mesh
    call m3dc1_ghost_delete
    call reset_trimats

    !Solve for fields
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call second(tstart)
    call solve_pi_tensor(pi_parallel, pi_perp, pi_parallel_n, pi_perp_n)
    if (myrank.eq.0) then
       call second(tend)
       write(0,'(A,f9.2,A)')'Pressure tensor LHS vecs calculated in',tend-tstart,' sec.'
    endif

    call destroy_field(pi_perp);  call destroy_field(pi_parallel)
    call destroy_field(pi_perp_n);  call destroy_field(pi_parallel_n)

    !Clean up
    call finalize_particles
  end subroutine particle_test

!---------------------------------------------------------------------------
  subroutine init_particles(ierr)
    use basic
    implicit none

    include 'mpif.h'

    integer, intent(out) :: ierr

    real, parameter :: c_mks = 2.9979e+8
    integer, parameter :: nglayers = 2 !# of layers of ghost zones around each mesh domain
    type(particle) :: dpar  !Dummy particle
    type(elfield), dimension(nneighbors+1) :: elcoefs
    type(xgeomterms) :: geomterms
    real, dimension(3) :: Bcyl
#ifndef USE3D
    real, dimension(2) :: mmsa
#else
    real xi, zi
#endif
    real    :: x1, x2, z1, z2, pdx, pdz, pdl, vpar, vperp
    real    :: EmaxeV, A_ion, Z_ion, speed, lambda_min, lambda_max, B0, B1
    real    :: gyroperiod, gfrac=5.0e-3, gkfrac=8.0, dtp, ldtmin
    integer :: npr, npz, npe, npmu, ir, iz, ie, imu, ip
    integer :: nelms, ielm, lc, noc, tridex, itri
    integer :: gfx, isghost, nle=0, nge=0

    !Load a ghost mesh with nglayers layers
    call m3dc1_ghost_load(nglayers)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (myrank.eq.0) print *,'Set up ghost mesh with',nglayers,' layers.'

    !Query ghost mesh dimensions
    call m3dc1_mesh_getnument(2, nelms)
    do gfx=1,nelms
       call m3dc1_ent_isghost(2, gfx-1, isghost)
       if (isghost.eq.1) then
          nge = nge + 1
       else
          nle = nle + 1
       endif
    enddo !gfx
    print *,myrank,':',nle,' local elements;',nge,' ghost elements.'
    nelms = local_elements()
    print *,myrank,': ',nelms,' local+ghost elements.'

    !Recompute triangle coefficients to include ghost elements
    call reset_trimats

    !Set up domain neighbor lookup table for fast particle communication
    call init_ndlookup(ierr)
    if (ierr.ne.0) then
       print *,'Error',ierr,' in init_ndlookup.'
       return
    endif

    !Allocate local storage for particle data
    allocate(pdata(nelms), jmppar(ndnbr), jinbuf(32))

    !Set up 'neighborlist' table of element neighbors for ensemble tracking
    call find_element_neighbors

    npr = 64;  npz = 64;  npe = 1;  npmu = 2

    !Particle spatial ranges
    call get_bounding_box(x1, z1, x2, z2)
    !if (myrank.eq.0) print *,'bb = ',x1,z1,x2,z2
    pdx = (x2 - x1)/real(npr);  pdz = (z2 - z1)/real(npz)

    !Particle velocity space ranges
    EmaxeV = 1000.0                             !Peak ion kinetic energy in eV
    A_ion = 1.0                                 !Ion atomic mass number
    m_ion = A_ion * m_proton                    !Ion mass in kg
    Z_ion = 1.0                                 !Ion atomic number (charge state)
    q_ion = Z_ion * e_mks                       !Ion charge in C
    qm_ion = q_ion / m_ion                      !Ion charge/mass ratio
    speed = sqrt(EmaxeV / A_ion) * vp1eV        !Peak ion speed in m/s
    if (myrank.eq.0) print *,'peak ion speed = ',speed,' m/s = ',speed/c_mks,' c.'

    lambda_min = 1.0e-5                          !Minimum particle pitch angle
    lambda_max = 3.14159                         !Maximum particle pitch angle
    pdl = (lambda_max - lambda_min)/(npmu - 1)


    !First pass: assign particles to processors, elements
    locparts = 0
    dpar%x(2) = 0.0
    do iz=1,npz !Loop over z positions
       dpar%x(3) = z1 + (iz - 0.5)*pdz

       do ir=1,npr !Loop over major radius positions
          dpar%x(1) = x1 + (ir - 0.5)*pdx

          !Check for local residence
#ifdef USE3D
          ielm = 0
          call whattri(dpar%x(1), dpar%x(2), dpar%x(3), ielm, xi, zi)
#else
          mmsa = dpar%x(1:3:2)
          call m3dc1_mesh_search(0, mmsa, ielm)
          ielm = ielm + 1
#endif
          if (ielm.le.0) cycle     !Not in local partition; skip.
          call m3dc1_ent_isghost(2, ielm-1, isghost)
          if (isghost.eq.1) cycle  !In ghost layer; skip.

          do ie=1,npe !Loop over kinetic energies
             dpar%v(1) = sqrt(real(ie)/real(npe)) * speed

             do imu=1,npmu  !Loop over pitch angles
                dpar%v(2) = lambda_min + (imu - 1.0)*pdl
                dpar%gid = npmu*(npe*(npr*(iz-1) + (ir-1)) + (ie-1)) + (imu-1)
                dpar%jel = ielm

                call add_particle(pdata, nelms, ielm, dpar, ierr)
                if (ierr.eq.0) locparts = locparts + 1
             enddo !imu
          enddo !ie
       enddo !ir
#ifdef JBDEBUG
       if (myrank.eq.0) print *,npr*npe*npmu,' particles assessed for row',iz
#endif
    enddo !iz

    write(0,'(I6,A,I9,A,f9.2,A)')myrank,':',locparts,' local particle(s). (avg',&
         locparts/real(nle),' per cell)'
    lc = sum(pdata(:)%np)
    if (lc.ne.locparts) print *,myrank,': mismatch in local particle count.'

    call mpi_allreduce(locparts, nparticles, 1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if (myrank.eq.0) then
       write(0,'(I8,A,I8,A)')nparticles,' particle(s) assigned out of ',&
            npr*npz*npe*npmu,' candidates.'
    endif

    !2nd pass: initialize velocity components
    noc = 0
    ldtmin = +1.0e+4
    elcoefs(:)%itri = 0
    do ielm=1,nelms
       if (pdata(ielm)%np.eq.0) cycle

       !Load scalar fields for this element
       call get_field_coefs(ielm, elcoefs(1), .false.)

       !Loop over particles inside
       do ip=1,pdata(ielm)%np
          vpar = pdata(ielm)%ion(ip)%v(1) * cos(pdata(ielm)%ion(ip)%v(2))
          vperp = pdata(ielm)%ion(ip)%v(1) * sin(pdata(ielm)%ion(ip)%v(2))

          itri = ielm
          call get_geom_terms(pdata(ielm)%ion(ip)%x, itri, elcoefs, tridex, &
               geomterms, .false., ierr)
          if (tridex.ne.1.or.itri.ne.elcoefs(tridex)%itri) ierr = 1
          if (ierr.ne.0) then
             print *,myrank,': get_geom_terms call failed for particle',ip,&
                  ' of element',ielm
             cycle
          endif !ierr

          call getBcyl(pdata(ielm)%ion(ip)%x, elcoefs(1), geomterms, Bcyl)
          B0 = sqrt(dot_product(Bcyl, Bcyl))
          gyroperiod = 6.283185307 / (qm_ion * B0)

          if (vspdims.eq.3) then !full orbit
             B1 = sqrt(Bcyl(1)**2 + Bcyl(2)**2)

             pdata(ielm)%ion(ip)%v(1) = vpar*Bcyl(1)/B0 - vperp*Bcyl(2)/B1  !v_R
             pdata(ielm)%ion(ip)%v(2) = vpar*Bcyl(2)/B0 + vperp*Bcyl(1)/B1  !v_phi
             pdata(ielm)%ion(ip)%v(3) = vpar*Bcyl(3)/B0                     !v_z

             dtp = gfrac * gyroperiod
          else !gyro- or drift kinetic
             pdata(ielm)%ion(ip)%v(1) = vpar                        !v_parallel
             pdata(ielm)%ion(ip)%v(2) = (0.5*vperp**2)/(qm_ion*B0)  !mu/q

             dtp = gkfrac * gyroperiod
          endif !vspdims

          if (ldtmin.gt.dtp) ldtmin = dtp
       enddo !ip

       noc = noc + 1
    enddo !ielm
#ifdef JBDEBUG
    print *,myrank,':',noc,' / ',nle,' elements occupied.'
    print *,myrank,': ldtmin = ',ldtmin
#endif

    call mpi_allreduce(ldtmin, dt_ion, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
    if (myrank.eq.0) print *,'Particle dt = ',dt_ion,' s.'
    if (dt_ion.le.0.0) then
       ierr = 1
       return
    endif

    call define_mpi_particle(ierr)
  end subroutine init_particles

!---------------------------------------------------------------------------
! Initialize lookup table for neighboring domains
  subroutine init_ndlookup(ierr)
    use basic
    implicit none
    include 'mpif.h'

    integer, intent(out) :: ierr
    integer, dimension(:), allocatable :: scratch
    integer :: ngel, nlel, ipe, iel, isghost, gid

    ierr = 0

    !Find number of global elements, set up scratch array for tabulation
    call m3dc1_mesh_getnumglobalent(2, ngel)
    if (ngel.lt.1) then
       ierr = 1;  return
    endif
    allocate(scratch(ngel))

    !Find number of local+ghost elements, set up array for lookup
    call m3dc1_mesh_getnument(2, nlel)
    if (nlel.lt.1) then
       ierr = 2;  return
    endif
    allocate(dnbr(2,nlel), dnlist(maxrank))
    dnbr = -1;  ndnbr = 0

    !Main loop (serial): map ghost zones to neighboring PEs
    do ipe=0,maxrank-1
       !Fill in locally owned segment of table for this PE
       if (myrank.eq.ipe) then
          scratch = 0
          do iel=1,nlel
             call m3dc1_ent_isghost(2, iel-1, isghost)
             if (isghost.ne.1) then !true local element
                call m3dc1_ent_getglobalid(2, iel-1, gid)
                if (gid.lt.0.or.gid.ge.ngel) then
                   ierr = 3;  return
                endif
                scratch(gid+1) = iel
             endif !isghost.ne...
          enddo !iel
       endif !myrank.eq...
       call mpi_barrier(MPI_COMM_WORLD, ierr)

       !Communicate this table to all PEs
       call MPI_Bcast(scratch, ngel, mpi_integer, ipe, MPI_COMM_WORLD, ierr)

       !Update neighbor table
       if (myrank.ne.ipe) then
          do iel=1,nlel
             call m3dc1_ent_isghost(2, iel-1, isghost)
             if (isghost.eq.1) then
                call m3dc1_ent_getglobalid(2, iel-1, gid)
                if (gid.lt.0.or.gid.ge.ngel) then
                   ierr = 4;  return
                endif
                if (scratch(gid+1).gt.0) then
                   dnbr(1,iel) = ipe
                   dnbr(2,iel) = scratch(gid+1)
                endif
             endif !isghost.eq...
          enddo !iel
       endif !myrank.ne...
    enddo !ipe
    deallocate(scratch)

    !Tabulate neighbors
    do iel=1,nlel
       call m3dc1_ent_isghost(2, iel-1, isghost)
       if (isghost.eq.1) then
          if (dnbr(1,iel).lt.0) then
             print *,myrank,': Unmapped ghost zone',iel,' in init_ndlookup.'
             ierr = 5;  return
          else !Zone is mapped; look it up in the neighbor list.
             do ipe=1,ndnbr
                if (dnlist(ipe).eq.dnbr(1,iel)) exit
             enddo !ipe
             if (ipe.gt.ndnbr) then !Not found; add new entry
                ndnbr = ndnbr + 1
                if (ndnbr.gt.maxrank) then
                   ierr = 6;  return
                endif
                dnlist(ndnbr) = dnbr(1,iel)
             endif !ipe.gt...
             dnbr(1,iel) = ipe !Remap from PE index to neighbor list position
          endif !dnbr...
       endif !isghost.eq...
    enddo !iel

    print *,myrank,':',ndnbr,' domain neighbor(s).' !:',dnlist(1:ndnbr)
  end subroutine init_ndlookup

!---------------------------------------------------------------------------
  subroutine add_particle(pbuf, buflen, ient, part, ierr)
    implicit none

    integer, intent(in) :: buflen, ient
    type(elplist), dimension(buflen), intent(inout) :: pbuf
    type(particle), intent(in) :: part
    integer, intent(out) :: ierr

    integer, parameter :: mininc = 128
    type(particle), dimension(:), allocatable :: tmparr
    integer origsize

    ierr = 1

    if (ient.lt.1.or.ient.gt.buflen) return

    if (.not.allocated(pbuf(ient)%ion)) allocate(pbuf(ient)%ion(mininc))

    pbuf(ient)%np = pbuf(ient)%np + 1
    origsize = size(pbuf(ient)%ion)

    if (pbuf(ient)%np.gt.origsize) then !expand particle array
       allocate(tmparr(origsize))
       tmparr = pbuf(ient)%ion
       deallocate(pbuf(ient)%ion)
       allocate(pbuf(ient)%ion(origsize + mininc))
       pbuf(ient)%ion(1:origsize) = tmparr
       deallocate(tmparr)
    endif

    pbuf(ient)%ion(pbuf(ient)%np) = part

    ierr = 0
  end subroutine add_particle

!---------------------------------------------------------------------------
  subroutine delete_particle(pbuf, buflen, ient, ipart, ierr)
    implicit none
    intrinsic cshift

    integer, intent(in) :: buflen, ient, ipart
    type(elplist), dimension(buflen), intent(inout) :: pbuf
    integer, intent(out) :: ierr

    integer np

    !Error checking
    ierr = 1; if (ient.lt.1.or.ient.gt.buflen) return
    np = pbuf(ient)%np
    ierr = 2; if (ipart.lt.1.or.ipart.gt.np) return

    !Use intrinsic circular shift function to get rid of this particle
    pbuf(ient)%ion(ipart:np) = cshift(pbuf(ient)%ion(ipart:np), 1)
    pbuf(ient)%np = np - 1

    ierr = 0
  end subroutine delete_particle

!---------------------------------------------------------------------------
  subroutine finalize_particles
    implicit none
    include 'mpif.h'

    integer :: nelms, ielm

    if (allocated(neighborlist)) deallocate(neighborlist)

    if (allocated(pdata)) then
       nelms = size(pdata)
       do ielm=1,nelms
          if (allocated(pdata(ielm)%ion)) deallocate(pdata(ielm)%ion)
       enddo !ielm

       deallocate(pdata)
    endif

    if (allocated(jmppar)) then
       nelms = size(jmppar)
       do ielm=1,nelms
          if (allocated(jmppar(ielm)%ion)) deallocate(jmppar(ielm)%ion)
       enddo !ielm

       deallocate(jmppar)
    endif

    if (allocated(jinbuf)) deallocate(jinbuf)
    if (allocated(dnbr)) deallocate(dnbr)
    if (allocated(dnlist)) deallocate(dnlist)

    call mpi_type_free(mpi_particle, ielm)
  end subroutine finalize_particles

!---------------------------------------------------------------------------
  subroutine find_element_neighbors
    use basic
    implicit none

#ifdef USE3D
    integer, parameter :: maxconnect = 24 !Max # of faces converging on any mesh node
    integer, parameter :: ifverts = 4     !3 or 4 verts define a face
#else
    integer, parameter :: maxconnect = 8  !Max # of edges converging on any mesh node
    integer, parameter :: ifverts = 2     !Two verts define an edge
#endif

    type d_face
       integer :: el0, side
       integer, dimension(ifverts-1) :: v
    end type d_face

    type face
       type(d_face), dimension(maxconnect) :: o
       integer :: n = 0
    end type face

    type(face), dimension(:), allocatable  :: facelist
    integer, dimension(nodes_per_element)  :: enode
    integer, dimension(ifverts,nneighbors) :: sidevecsub !Vector subscripts for iface IDs
    integer, dimension(ifverts) :: iface
    integer, dimension(1) :: ml
    integer :: nelms, lnodes, ielm, side, v1, ivrt
    logical :: ep

    nelms = local_elements();  lnodes = local_nodes()
    if (nelms.lt.1) return

#ifdef USE3D
    if (nodes_per_element.eq.6) then !prisms
       !Define prism faces in terms of nodes
       sidevecsub(:,1) = (/ 1, 2, 3, 1 /)
       sidevecsub(:,2) = (/ 1, 4, 5, 2 /)
       sidevecsub(:,3) = (/ 2, 5, 6, 3 /)
       sidevecsub(:,4) = (/ 3, 6, 4, 1 /)
       sidevecsub(:,5) = (/ 4, 6, 5, 4 /)
#else
    if (nodes_per_element.eq.3) then !triangles
       !Define triangle edges in terms of nodes
       sidevecsub(:,1) = (/ 1, 2 /)
       sidevecsub(:,2) = (/ 2, 3 /)
       sidevecsub(:,3) = (/ 3, 1 /)
#endif
       !Allocate storage, initialize
       allocate(neighborlist(nneighbors,nelms), facelist(lnodes-1))
       neighborlist = -1

       !Loop over local elements
       do ielm=1,nelms
          call get_element_nodes(ielm, enode)

          !Loop over faces of this element
          do side=1,nneighbors
             iface = enode(sidevecsub(:,side))
             ml = minloc(iface)
             iface = cshift(iface, ml(1)-1)  !Arrange so lowest-indexed node comes 1st
             v1 = iface(1)

             !Search if the face is already present
             ep = .false.
             do ivrt=1,facelist(v1)%n
                if (veceq(facelist(v1)%o(ivrt)%v, iface(2:))) then
                   !Yes, update neighbor table
                   neighborlist(side,ielm) = facelist(v1)%o(ivrt)%el0
                   neighborlist(facelist(v1)%o(ivrt)%side, facelist(v1)%o(ivrt)%el0) = ielm
                   ep = .true.
                   exit
                endif !veceq...
             enddo !ivrt

             if (.not.ep) then !Face was not present; add it.
                facelist(v1)%n = facelist(v1)%n + 1
                if (facelist(v1)%n.gt.maxconnect) then !out of range
                   print *,'Error: too many connections in find_element_neighbors.'
                   deallocate(facelist, neighborlist)
                   return
                endif !n out-of-range
                facelist(v1)%o(facelist(v1)%n)%v = iface(2:)
                facelist(v1)%o(facelist(v1)%n)%el0 = ielm
                facelist(v1)%o(facelist(v1)%n)%side = side
             endif
          enddo !side
       enddo !ielm

       deallocate(facelist)
    else
       if(myrank.eq.0)print *,nodes_per_element,' nodes per element; cannot find neighbors.'
    endif !nodes_per_element...
  end subroutine find_element_neighbors

!---------------------------------------------------------------------------
  logical function veceq(v1, v2)
    use basic
    implicit none

#ifdef USE3D
    integer, dimension(3), intent(in) :: v1, v2
    veceq = (v1(1).eq.v2(3) .and. v1(2).eq.v2(2) .and. v1(3).eq.v2(1))
#else
    integer, dimension(1), intent(in) :: v1, v2
    veceq = (v1(1).eq.v2(1))
#endif
  end function veceq

!---------------------------------------------------------------------------
  subroutine advance_particles(tinc)
    use basic  !For MPI variables
    implicit none
    include 'mpif.h'

    real, intent(in) :: tinc  !Time increment for particle advance

    type(elfield), dimension(nneighbors+1) :: elcoefs
    real    :: dtp, trem
    integer :: nelms, ielm, itri, ipart, ierr
    integer :: nlost, ipe, lunf, gunf, nstep
    integer :: nhop, thop, nreas, treas  !Stats on ptcle movement w/in local domain
    integer :: nrec, isghost, maxin
    integer, dimension(ndnbr) :: npin
    integer, dimension(MPI_STATUS_SIZE) :: status

    !if (myrank.eq.0) print *,'advancing particles by ',tinc

    nelms = size(pdata)

    do ielm=1,nelms
       do ipart=1,pdata(ielm)%np
          pdata(ielm)%ion(ipart)%tlast = 0.0
       enddo
    enddo

    elcoefs(:)%itri = 0

    do !Iterate until all particles are in the correct domain
       thop = 0;  treas = 0
       nlost = 0;  lunf = 0
       jmppar(:)%np = 0  !Clear jumping particle buffer

       !Loop over all local elements (good candidate for OMP parallelization)
!$OMP PARALLEL
       do ielm=1,nelms
          if (pdata(ielm)%np.eq.0) cycle  !Skip if element is empty
          nhop = 0;  nreas = 0

          !Load scalar fields for this element & its nearest neighbors
          call update_coef_ensemble(elcoefs, ielm)

          !Advance particles within this element
          ipart = 1
          do !For each particle ipart
             nstep = 0
             dtp = dt_ion;  itri = ielm

             do !Advance particle by tinc
                trem = tinc - pdata(ielm)%ion(ipart)%tlast  !time remaining to advance
                if (trem.le.0.) exit
                if (dtp.gt.trem) dtp = trem

                call rk4(pdata(ielm)%ion(ipart), dtp, elcoefs, itri, ierr)
                !call rk5ck(pdata(ielm)%ion(ipart), dtp, elcoefs, itri, ierr)

                if (ierr.eq.1) then ! Particle exited local+ghost domain -> lost
                   !print *,myrank,': el',ielm,', p',ipart,pdata(ielm)%ion(ipart)%gid,&
                   !     ' exited domain'
                   call delete_particle(pdata, nelms, ielm, ipart, ierr)
                   if (ierr.ne.0) then
                      print *,myrank,': error',ierr,' deleting lost particle',&
                           pdata(ielm)%ion(ipart)%gid,' from elm',ielm
                   endif
                   nlost = nlost + 1;  ipart = ipart - 1
                   exit !Break out of tinc loop, go to next particle.
                endif

                pdata(ielm)%ion(ipart)%tlast = pdata(ielm)%ion(ipart)%tlast + dtp
                nstep = nstep + 1
 
                if (itri.eq.ielm) cycle !Continue push within element

                !Particle has moved to a new element -> outer loop must repeat.
                lunf = 1

                !Test whether new element is in ghost layer
                call m3dc1_ent_isghost(2, itri-1, isghost)
                if (isghost.eq.1) then !It is -> schedule move to new PE
                   !Add particle to jump list for target PE
                   pdata(ielm)%ion(ipart)%jel = dnbr(2,itri)
                   call add_particle(jmppar, ndnbr, dnbr(1,itri), &
                        pdata(ielm)%ion(ipart), ierr)
                   if (ierr.ne.0) print *,myrank,': error in jmp add_particle!'
                else  !Particle is still on local domain
                   if (ierr.eq.2) then ! Particle exited current element ensemble
                      nhop = nhop + 1
                   else                ! Particle moved within current element ensemble
                      nreas = nreas + 1
                   endif

                   !Add it to the new element
                   pdata(ielm)%ion(ipart)%jel = itri
                   call add_particle(pdata, nelms, itri, pdata(ielm)%ion(ipart), ierr)
                   if (ierr.ne.0) print *,myrank,': error in add_particle!'
                endif

                !Remove particle from the current element
                call delete_particle(pdata, nelms, ielm, ipart, ierr)
                if (ierr.ne.0) print *,myrank,': error',ierr,'in delete_particle!'

                ipart = ipart - 1
                exit !Break out of tinc loop, go to next particle.
             enddo !tinc advance

             ipart = ipart + 1
             if (ipart.gt.pdata(ielm)%np) exit
          enddo !ipart

          thop = thop + nhop
          treas = treas + nreas
       enddo !ielm
!$OMP END PARALLEL

#ifdef JBDEBUG
       print *,myrank,':',nlost,' / ',locparts,' total lost.'
       print *,myrank,':',thop,' / ',locparts,' total hopped to new ensemble.'
       print *,myrank,':',treas,' / ',locparts,' total reassigned within ensemble.'
       print *,myrank,':',sum(jmppar(:)%np),' / ',locparts,' total jumped to new domain.'
       print *,myrank,': jump by dest:',jmppar(:)%np
#endif
       locparts = locparts - nlost

       call mpi_allreduce(lunf, gunf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
       if (gunf.le.0) exit ! All particles have reached target time

       !Reassign transiting particles
       !Tally particles to be received from each neighboring domain
       maxin = 0
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       do ipe=1,ndnbr
          call MPI_SendRecv(jmppar(ipe)%np, 1, MPI_INTEGER, dnlist(ipe), 42, &
                            npin(ipe), 1, MPI_INTEGER, dnlist(ipe), 42, &
                            MPI_COMM_WORLD, status, ierr)
          if (ierr.ne.0) print *,myrank,'MPI_SendRecv error',ierr
          if (npin(ipe).gt.maxin) maxin = npin(ipe)
          locparts = locparts - jmppar(ipe)%np
       enddo !ipe
       call MPI_Barrier(MPI_COMM_WORLD, ierr)

       !Make sure receive buffer for incoming particles is large enough.
       if (size(jinbuf).lt.maxin) then
          deallocate(jinbuf)
          allocate(jinbuf(maxin))
       endif

       !Send/receive particle data
       nrec = 0
       do ipe=1,ndnbr
          call MPI_SendRecv(jmppar(ipe)%ion, jmppar(ipe)%np, mpi_particle, dnlist(ipe), 43, &
               jinbuf, npin(ipe), mpi_particle, dnlist(ipe), 43, &
               MPI_COMM_WORLD, status, ierr)

          do ipart=1,npin(ipe)
             !Add particle to local list
             call add_particle(pdata, nelms, jinbuf(ipart)%jel, jinbuf(ipart), ierr)
             if (ierr.eq.0) then
                nrec = nrec + 1
             else
                print *,myrank,': jump error in add_particle!'
             endif
          enddo !ipart
       enddo !ipe
       call MPI_Barrier(MPI_COMM_WORLD, ierr)

       if (nrec.ne.sum(npin)) then
          print *,myrank,': error:',nrec,'/',sum(npin),' jumps receieved.'
       endif
       locparts = locparts + nrec
    enddo !outer loop

    call mpi_allreduce(locparts, nparticles, 1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if (myrank.eq.0) &
         print *,nparticles,' particle(s) remaining after advance step.'
  end subroutine advance_particles

!---------------------------------------------------------------------------
! 4th-order Runge-Kutta integrator, no adaptive time step control.
!  Four derivative evaluations per step.
! (Adapted from Numerical Recipes in C, 2nd edition, 1994, pp. 712-713).
!
  subroutine rk4(part, dt, fh, itri, ierr)
    implicit none

    type(particle), intent(inout) :: part
    real, intent(in) :: dt
    type(elfield), dimension(nneighbors+1), intent(in) :: fh
    integer, intent(inout) :: itri
    integer, intent(out) :: ierr

    real, parameter :: onethird = 1.0/3.0
    real, dimension(3) :: k1, k2, k3, k4, y1
    real, dimension(vspdims) :: l1, l2, l3, l4, z1
    real :: hh, m1, m2, m3, m4, w1

    ierr = 0
    hh = 0.5*dt

    !1st step
    call fdot(part%x, part%v, part%wt, k1, l1, m1, fh, itri, ierr)
    if (ierr.eq.1) return
    y1 = part%x + hh*k1;  z1 = part%v + hh*l1;  w1 = part%wt + hh*m1

    !2nd step
    call fdot(y1, z1, w1, k2, l2, m2, fh, itri, ierr)
    if (ierr.eq.1) return
    y1 = part%x + hh*k2;  z1 = part%v + hh*l2;  w1 = part%wt + hh*m2

    !3rd step
    call fdot(y1, z1, w1, k3, l3, m3, fh, itri, ierr)
    if (ierr.eq.1) return
    y1 = part%x + dt*k3;  z1 = part%v + dt*l3;  w1 = part%wt + hh*m3

    !4th step
    call fdot(y1, z1, w1, k4, l4, m4, fh, itri, ierr)
    if (ierr.eq.1) return
    part%x  = part%x  + onethird*dt*(k2 + k3 + 0.5*(k1 + k4))
    part%v  = part%v  + onethird*dt*(l2 + l3 + 0.5*(l1 + l4))
    part%wt = part%wt + onethird*dt*(m2 + m3 + 0.5*(m1 + m4))
  end subroutine rk4

!----------------------------------------------------------------------------------
! 5th-order Cash-Karp Runge-Kutta integrator with embedded 4th-order error estimate.
!  Six derivative evaluations per call.
! (Adapted from Numerical Recipes in C, 2nd edition, 1994, pp. 719-720).
!
  subroutine rk5ck(part, dt, fh, itri, ierr)
    implicit none

    type(particle), intent(inout) :: part
    real, intent(in) :: dt
    type(elfield), dimension(nneighbors+1), intent(in) :: fh
    integer, intent(inout) :: itri
    !real, dimension(4+vspdims), intent(out) :: perr
    integer, intent(out) :: ierr

    real, parameter :: b21=0.2, b31=0.075, b32=0.225
    real, parameter :: b41=0.3, b42=-0.9, b43=1.2
    real, parameter :: b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0
    real, parameter :: b61=1631.0/55296.0, b62=175.0/512.0
    real, parameter :: b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0
    real, parameter :: c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0
    real, parameter :: c6=512.0/1771.0
    !real, parameter :: dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0
    !real, parameter :: dc4=c4-13525.0/55296.0, dc5=-277.0/14336.0, dc6=c6-0.25

    real, dimension(3)       :: xdot, y1, ak2, ak3, ak4, ak5, ak6
    real, dimension(vspdims) :: vdot, z1, bk2, bk3, bk4, bk5, bk6
    real                     :: wdot, w1, ck2, ck3, ck4, ck5, ck6

    ierr = 0

    !1st step
    call fdot(part%x, part%v, part%wt, xdot, vdot, wdot, fh, itri, ierr)
    if (ierr.eq.1) return
    y1 = part%x  + b21*dt*xdot
    z1 = part%v  + b21*dt*vdot
    w1 = part%wt + b21*dt*wdot

    !2nd step
    call fdot(y1, z1, w1, ak2, bk2, ck2, fh, itri, ierr)
    if (ierr.eq.1) return
    y1 = part%x  + dt*(b31*xdot + b32*ak2)
    z1 = part%v  + dt*(b31*vdot + b32*bk2)
    w1 = part%wt + dt*(b31*wdot + b32*ck2)

    !3rd step
    call fdot(y1, z1, w1, ak3, bk3, ck3, fh, itri, ierr)
    if (ierr.eq.1) return
    y1 = part%x  + dt*(b41*xdot + b42*ak2 + b43*ak3)
    z1 = part%v  + dt*(b41*vdot + b42*bk2 + b43*bk3)
    w1 = part%wt + dt*(b41*wdot + b42*ck2 + b43*ck3)

    !4th step
    call fdot(y1, z1, w1, ak4, bk4, ck4, fh, itri, ierr)
    if (ierr.eq.1) return
    y1 = part%x  + dt*(b51*xdot + b52*ak2 + b53*ak3 + b54*ak4)
    z1 = part%v  + dt*(b51*vdot + b52*bk2 + b53*bk3 + b54*bk4)
    w1 = part%wt + dt*(b51*wdot + b52*ck2 + b53*ck3 + b54*ck4)

    !5th step
    call fdot(y1, z1, w1, ak5, bk5, ck5, fh, itri, ierr)
    if (ierr.eq.1) return
    y1 = part%x  + dt*(b61*xdot + b62*ak2 + b63*ak3 + b64*ak4 + b65*ak5)
    z1 = part%v  + dt*(b61*vdot + b62*bk2 + b63*bk3 + b64*bk4 + b65*bk5)
    w1 = part%wt + dt*(b61*wdot + b62*ck2 + b63*ck3 + b64*ck4 + b65*ck5)

    !6th step
    call fdot(y1, z1, w1, ak6, bk6, ck6, fh, itri, ierr)
    if (ierr.eq.1) return
    part%x  = part%x  + dt*(c1*xdot + c3*ak3 + c4*ak4 + c6*ak6)
    part%v  = part%v  + dt*(c1*vdot + c3*bk3 + c4*bk4 + c6*bk6)
    part%wt = part%wt + dt*(c1*wdot + c3*ck3 + c4*ck4 + c6*ck6)

    !Error estimate
    !perr(1:3) = dt*(dc1*xdot + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)
    !perr(4:3+vspdims) = dt*(dc1*vdot + dc3*bk3 + dc4*bk4 + dc5*bk5 + dc6*bk6)
    !perr(4+vspdims) = dt*(dc1*wdot + dc3*ck3 + dc4*ck4 + dc5*ck5 + dc6*ck6)
  end subroutine rk5ck

!---------------------------------------------------------------------------
  subroutine fdot(x, v, w, dxdt, dvdt, dwdt, fh, itri, ierr)
    use basic
    implicit none

    real, dimension(3), intent(in)                     :: x
    real, dimension(3), intent(out)                    :: dxdt
    real, dimension(vspdims), intent(in)               :: v
    real, dimension(vspdims), intent(out)              :: dvdt
    real, intent(in)                                   :: w
    real, intent(out)                                  :: dwdt
    type(elfield), dimension(nneighbors+1), intent(in) :: fh
    integer, intent(inout)                             :: itri
    integer, intent(out)                               :: ierr

    real, parameter :: g_mks = 9.8067 ! earth avg surf grav accel in m/s/s
    type(elfield) :: fh_hop
    type(xgeomterms) :: geomterms
    real, dimension(3) :: B_cyl, E_cyl, bhat, svec, Bstar
    real, dimension(3) :: dBdR, dBdphi, dBdz, gradB0
    real :: Rinv = 1.0, B0, Bss
    integer :: tridex

    ierr = 0
    if (itor.eq.1) Rinv = 1.0/x(1)

    !Need terms to compute fields to calculate acceleration
    call get_geom_terms(x, itri, fh, tridex, geomterms, vspdims.eq.2, ierr)
    if (ierr.ne.0) return

    !Get electric field components
    if (tridex.le.0) then !Not part of local ensemble!
       ierr = 2
       call get_field_coefs(itri, fh_hop, .true.) !Load field into temp. buffer
       call getEcyl(fh_hop, geomterms, E_cyl)
    else
       call getEcyl(fh(tridex), geomterms, E_cyl)
    endif

    !Calculate time derivatives
    if (vspdims.eq.3) then !full orbit: ma = q(E + vxB) + mg
       dxdt(1:vspdims) = v

       if (tridex.gt.0) then
          call getBcyl(x, fh(tridex), geomterms, B_cyl)
       else
          call getBcyl(x, fh_hop, geomterms, B_cyl)
       endif

       dvdt(1) = qm_ion*(E_cyl(1) + v(2)*B_cyl(3) - v(3)*B_cyl(2))
       dvdt(2) = qm_ion*(E_cyl(2) + v(3)*B_cyl(1) - v(1)*B_cyl(3))
       dvdt(3) = qm_ion*(E_cyl(3) + v(1)*B_cyl(2) - v(2)*B_cyl(1)) - g_mks

       !if (itor.eq.1) then
       !   dvdt(1) = dvdt(1) + Rinv*v(2)**2  !Centripetal acceleration
       !   dvdt(2) = dvdt(2) - 2.0*Rinv*v(1)*v(2)  !Coriolis effect
       !endif
    else ! Drift-kinetic equation
       if (tridex.gt.0) then
          call getBcylprime(x, fh(tridex), geomterms, B_cyl, dBdR, dBdphi, dBdz)
       else
          call getBcylprime(x, fh_hop, geomterms, B_cyl, dBdR, dBdphi, dBdz)
       endif

       B0 = sqrt(dot_product(B_cyl, B_cyl))  !Magnitude of B
       bhat = B_cyl / B0                     !Unit vector in b direction

       ! Gradient of B0 = grad(B.B)/(2 B0) = (B . grad B)/B0
       gradB0(1) = dot_product(bhat, dBdR)
       gradB0(2) = Rinv*dot_product(bhat, dBdphi)
       gradB0(3) = dot_product(bhat, dBdz)

       ! Curl of bhat = curl(B/B0) = curl(B)/B0 - (grad B0 x B)/(B0**2)
       svec(1) = (Rinv*dBdphi(3) - dBdz(2) + &
            (B_cyl(2)*gradB0(3) - B_cyl(3)*gradB0(2))/B0)/B0
       svec(2) = (dBdz(1) - dBdR(3) + &
            (B_cyl(3)*gradB0(1) - B_cyl(1)*gradB0(3))/B0)/B0
       if (itor.eq.1) then
          svec(3) = (Rinv*B_cyl(2) + dBdR(2) - Rinv*dBdphi(1) + &
               (B_cyl(1)*gradB0(2) - B_cyl(2)*gradB0(1))/B0)/B0
       else
          svec(3) = (dBdR(2) - dBdphi(1) + &
               (B_cyl(1)*gradB0(2) - B_cyl(2)*gradB0(1))/B0)/B0
       endif

       Bstar = B_cyl + (v(1)/qm_ion)*svec
       Bss = dot_product(Bstar, bhat)

       svec = v(2)*gradB0 - E_cyl  ! - g_mks/qm_ion

       dxdt(1) = (v(1)*Bstar(1) + bhat(2)*svec(3) - bhat(3)*svec(2))/Bss
       dxdt(2) = (v(1)*Bstar(2) + bhat(3)*svec(1) - bhat(1)*svec(3))/Bss
       dxdt(3) = (v(1)*Bstar(3) + bhat(1)*svec(2) - bhat(2)*svec(1))/Bss

       dvdt(1) = -qm_ion*dot_product(Bstar, svec)/Bss
       dvdt(2) = 0. !magnetic moment is conserved.
    endif

    dxdt(2) = Rinv*dxdt(2)  !phi-dot = (v_phi / R) for cylindrical case

    dwdt = 0.0*w !Evolution of delta-f weights not yet implemented!
  end subroutine fdot

!---------------------------------------------------------------------------
! Compute terms for a function and its partial derivatives with respect to
!  R and z in the reduced quintic expansion at position x.
  subroutine get_geom_terms(x, ielm, fh, tridex, gh, ic2, ierr)
    implicit none

    real, dimension(3), intent(in) :: x
    integer, intent(inout) :: ielm
    type(elfield), dimension(nneighbors+1), intent(in) :: fh
    type(xgeomterms), intent(out) :: gh  !Geometric terms handle
    logical, intent(in)  :: ic2          !Compute 2nd derivative terms?
    integer, intent(out) :: tridex, ierr

    type(element_data) :: eldat
    real :: xi, zi, eta, dxi, deta
    real :: d2xi, d2eta, dxieta
    integer pp
#ifdef USE3D
    real    :: gtmp, drtmp, dztmp, zpow
    integer :: ii, jj
#else
    real, dimension(2) :: mmsa
    integer ktri
#endif

    ierr = 0;  tridex = -1

#ifdef USE3D
    call whattri(x(1),x(2),x(3),ielm,xi,zi)
#else
    mmsa = x(1:3:2)
    if (ielm.lt.1) ielm = 1
    call m3dc1_mesh_search(ielm-1, mmsa, ktri)
    ielm = ktri + 1
#endif
    if (ielm.le.0) then !The triangle is not in the local partition
       ierr = 1
       return
    endif

    tridex = ensemble_index(fh, ielm)

    call get_element_data(ielm, eldat)
    call global_to_local(eldat, x(1), x(2), x(3), xi, zi, eta)

    !Compute terms for function and 1st derivatives
    gh%dr = 0.;  gh%dz = 0.
    do pp=1,coeffs_per_tri
       gh%g(pp) = xi**mi(pp) * eta**ni(pp)

       if (mi(pp).gt.0) then
          dxi = mi(pp) * xi**(mi(pp)-1) * eta**ni(pp)

          gh%dr(pp) = gh%dr(pp) + eldat%co*dxi
          gh%dz(pp) = gh%dz(pp) + eldat%sn*dxi
       endif

       if (ni(pp).gt.0) then
          deta = xi**mi(pp) * ni(pp)*eta**(ni(pp)-1)

          gh%dr(pp) = gh%dr(pp) - eldat%sn*deta
          gh%dz(pp) = gh%dz(pp) + eldat%co*deta
       endif

#ifdef USE3D
       gtmp = gh%g(pp);  drtmp = gh%dr(pp);  dztmp = gh%dz(pp)
       do ii=1,coeffs_per_dphi
          jj = pp  + (ii - 1)*coeffs_per_tri
          zpow = zi**li(ii)

          gh%g(jj)  = gtmp  * zpow
          gh%dr(jj) = drtmp * zpow
          gh%dz(jj) = dztmp * zpow

          !First toroidal derivative
          if (li(ii).gt.0) then
             zpow = li(ii) * zi**(li(ii) - 1)
             gh%dphi(jj)  = gtmp  * zpow
             gh%drphi(jj) = drtmp * zpow
             gh%dzphi(jj) = dztmp * zpow

             !Second toroidal derivative
             if (li(ii).gt.1) then
                zpow = li(ii) * (li(ii) - 1) * zi**(li(ii) - 2)
                gh%drphiphi(jj) = drtmp * zpow
                gh%dzphiphi(jj) = dztmp * zpow
             else
                gh%drphiphi(jj) = 0.
                gh%dzphiphi(jj) = 0.
             endif
          else
             gh%dphi(jj)     = 0.
             gh%drphi(jj)    = 0.
             gh%dzphi(jj)    = 0.
             gh%drphiphi(jj) = 0.
             gh%dzphiphi(jj) = 0.
          endif
       enddo !ii
#endif
    enddo !pp

    if (ic2) then !2nd derivative terms
       gh%drr = 0.;  gh%drz = 0.;  gh%dzz = 0.
       do pp=1,coeffs_per_tri
          if (mi(pp).gt.0) then
             if (mi(pp).gt.1) then
                d2xi = mi(pp)*(mi(pp)-1)*xi**(mi(pp)-2) * eta**ni(pp)

                gh%drr(pp) = gh%drr(pp) + d2xi*eldat%co**2
                gh%drz(pp) = gh%drz(pp) + d2xi*eldat%co*eldat%sn
                gh%dzz(pp) = gh%dzz(pp) + d2xi*eldat%sn**2
             endif !mi > 1

             if (ni(pp).gt.0) then
                dxieta = mi(pp)*ni(pp) * xi**(mi(pp)-1) * eta**(ni(pp)-1)

                gh%drr(pp) = gh%drr(pp) - 2.0*dxieta*eldat%co*eldat%sn
                gh%drz(pp) = gh%drz(pp) + dxieta*(2.0*eldat%co**2 - 1.0)
                gh%dzz(pp) = gh%dzz(pp) + 2.0*dxieta*eldat%co*eldat%sn
             endif !ni > 0
          endif !mi > 0

          if (ni(pp).gt.1) then
             d2eta = xi**mi(pp) * ni(pp)*(ni(pp)-1)*eta**(ni(pp)-2)

             gh%drr(pp) = gh%drr(pp) + d2eta*eldat%sn**2
             gh%drz(pp) = gh%drz(pp) - d2eta*eldat%co*eldat%sn
             gh%dzz(pp) = gh%dzz(pp) + d2eta*eldat%co**2
          endif !ni > 1
#ifdef USE3D
          gtmp = gh%drr(pp);  drtmp = gh%drz(pp);  dztmp = gh%dzz(pp)
          do ii=1,coeffs_per_dphi
             jj = pp  + (ii - 1)*coeffs_per_tri
             zpow = zi**li(ii)

             gh%drr(jj) = gtmp  * zpow
             gh%drz(jj) = drtmp * zpow
             gh%dzz(jj) = dztmp * zpow

             !First toroidal derivative
             if (li(ii).gt.0) then
                zpow = li(ii) * zi**(li(ii) - 1)
                gh%drrphi(jj) = gh%drr(pp) * zpow
                gh%drzphi(jj) = gh%drz(pp) * zpow
                gh%dzzphi(jj) = gh%dzz(pp) * zpow
             else
                gh%drrphi(jj) = 0.
                gh%drzphi(jj) = 0.
                gh%dzzphi(jj) = 0.
             endif
          enddo !ii
#endif
       enddo !pp
    endif !ic2
  end subroutine get_geom_terms

!---------------------------------------------------------------------------
  subroutine update_coef_ensemble(ensemble, itri)
    implicit none

    type(elfield), dimension(nneighbors+1), intent(inout) :: ensemble
    integer, intent(in) :: itri

    logical, dimension(nneighbors+1) :: ladd, ldel
    integer :: tridex, nbr, jtri

    !Determine which elements need to be added, which can be deleted
    ladd = .false.;  ldel = .true.
    tridex = ensemble_index(ensemble, itri)
    if (tridex.lt.1) then
       ladd(1) = .true.
    else
       ldel(tridex) = .false.
    endif
    do nbr=1,nneighbors
       jtri = neighborlist(nbr,itri) !Look up the jth neighbor of this element
       if (jtri.gt.0) then         !If the neighbor exists,
          tridex = ensemble_index(ensemble, jtri) !See if it is loaded already
          if (tridex.lt.1) then !Not loaded; schedule it for addition
             ladd(nbr+1) = .true.
          else                  !Already loaded; prevent deletion
             ldel(tridex) = .false.
          endif
       endif
    enddo !nbr

    !Load elements as necessary
    tridex = 1
    do !Find first available space
       if (ldel(tridex)) exit
       tridex = tridex + 1
    enddo
    if (ladd(1)) then !Load central element here
       call get_field_coefs(itri, ensemble(tridex), .true.)
       tridex = tridex + 1
    endif !ladd(1)
    do nbr=1,nneighbors !Loop through adjacent elements
       if (ladd(nbr+1)) then
          do !Find next available space
             if (ldel(tridex)) exit
             tridex = tridex + 1
          enddo

          !Load jth neighbor here
          call get_field_coefs(neighborlist(nbr,itri), ensemble(tridex), .true.)
          tridex = tridex + 1
       endif !ladd
    enddo !nbr
  end subroutine update_coef_ensemble

!---------------------------------------------------------------------------
  integer function ensemble_index(ensemble, itri)
    implicit none

    type(elfield), dimension(nneighbors+1), intent(in) :: ensemble
    integer, intent(in) :: itri

    integer :: idx

    ensemble_index = -1

    do idx=1,nneighbors+1
       if (ensemble(idx)%itri.eq.itri) then
          ensemble_index = idx
          return
       endif
    enddo !idx
  end function ensemble_index

!---------------------------------------------------------------------------
  subroutine get_field_coefs(ielm, fh, getE)
    use arrays
    use basic
    use auxiliary_fields
    implicit none

    type(elfield), intent(out) :: fh  !Field handle
    integer, intent(in) :: ielm
    logical, intent(in) :: getE

    logical :: use_f = .false.
#ifdef USECOMPLEX
    use_f = .true.
#endif
#ifdef USE3D
    use_f = .true.
#endif

    !Always get magnetic field components
    call calcavector(ielm, psi_field(0), fh%psiv0)
    call calcavector(ielm, bz_field(0), fh%Bzv0)
    if (use_f) call calcavector(ielm, bf_field(linear), fh%Bfv)
    if (linear.eq.1) then
       call calcavector(ielm, psi_field(1), fh%psiv1)
       call calcavector(ielm, bz_field(1), fh%Bzv1)
    endif !linear

    !Get electric field components if needed
    if (getE) then
       call calcavector(ielm, ef_r, fh%er)
       call calcavector(ielm, ef_phi, fh%ephi)
       call calcavector(ielm, ef_z, fh%ez)
    endif

    fh%itri = ielm
  end subroutine get_field_coefs

!---------------------------------------------------------------------------
  subroutine getBcyl(x, fh, gh, Bcyl)
    use basic
    implicit none

    real, dimension(3), intent(in) :: x      !Position
    type(elfield), intent(in) :: fh          !Field handle
    type(xgeomterms), intent(in) :: gh       !Geometric terms handle
    real, dimension(3), intent(out) :: Bcyl  !Output magnetic field

    vectype, dimension(3) :: temp
    real :: Rinv=1.0

    if (itor.eq.1) Rinv = 1.0/x(1)

    !Total/Equilibrium part
    !B_poloidal_axisymmetric = grad psi x grad phi
    Bcyl(1) = -Rinv*dot_product(fh%psiv0, gh%dz)
    Bcyl(3) =  Rinv*dot_product(fh%psiv0, gh%dr)
#ifdef USE3D
    !Non-axisymmetric B_poloidal term: - grad f'
    Bcyl(1) = Bcyl(1) - dot_product(fh%Bfv, gh%drphi)
    Bcyl(3) = Bcyl(3) - dot_product(fh%Bfv, gh%dzphi)
#endif

    !B_toroidal = B_Z / R
    Bcyl(2) = Rinv*dot_product(fh%Bzv0, gh%g)

    if (linear.eq.1) then
       !Perturbed part
       temp(1) = -Rinv*dot_product(fh%psiv1, gh%dz)
       temp(3) =  Rinv*dot_product(fh%psiv1, gh%dr)
       temp(2) =  Rinv*dot_product(fh%Bzv1,  gh%g)
#ifdef USECOMPLEX
       temp(1) = temp(1) - dot_product(fh%Bfv, gh%dr) * rfac
       temp(3) = temp(3) - dot_product(fh%Bfv, gh%dz) * rfac
       Bcyl = Bcyl + real(temp * exp(rfac*x(2)))
#else
       Bcyl = Bcyl + temp
#endif
    endif !linear
  end subroutine getBcyl

!---------------------------------------------------------------------------
  subroutine getBcylprime(x, fh, gh, Bcyl, dBdR, dBdphi, dBdz)
    use basic
    implicit none

    real, dimension(3), intent(in)  :: x
    type(elfield), intent(in)       :: fh
    type(xgeomterms), intent(in)    :: gh
    real, dimension(3), intent(out) :: Bcyl, dBdR, dBdphi, dBdz

    vectype, dimension(3) :: temp, tempR, tempz
    real :: Rinv=1.0

    if (itor.eq.1) Rinv = 1.0/x(1)

    !Total/Equilibrium part
    !B_poloidal_axisym = grad psi x grad phi
    Bcyl(1) = -Rinv*dot_product(fh%psiv0, gh%dz)
    dBdR(1) = -Rinv*dot_product(fh%psiv0, gh%drz)
    dBdz(1) = -Rinv*dot_product(fh%psiv0, gh%dzz)

    Bcyl(3) =  Rinv*dot_product(fh%psiv0, gh%dr)
    dBdR(3) =  Rinv*dot_product(fh%psiv0, gh%drr)
    dBdz(3) =  Rinv*dot_product(fh%psiv0, gh%drz)

    !B_toroidal = B_Z / R
    Bcyl(2) = Rinv*dot_product(fh%Bzv0, gh%g)
    dBdR(2) = Rinv*dot_product(fh%Bzv0, gh%dr)
    dBdz(2) = Rinv*dot_product(fh%Bzv0, gh%dz)

    if (itor.eq.1) dBdR = dBdR - Rinv*Bcyl

#ifdef USE3D
    !Non-axisymmetric B_poloidal term: - grad f'
    Bcyl(1) = Bcyl(1) - dot_product(fh%Bfv, gh%drphi)
    dBdR(1) = dBdR(1) - dot_product(fh%Bfv, gh%drrphi)
    dBdphi(1) = -Rinv*dot_product(fh%psiv0, gh%dzphi) - dot_product(fh%Bfv, gh%drphiphi)
    dBdz(1) = dBdz(1) - dot_product(fh%Bfv, gh%drzphi)

    Bcyl(3) = Bcyl(3) - dot_product(fh%Bfv, gh%dzphi)
    dBdR(3) = dBdR(3) - dot_product(fh%Bfv, gh%drzphi)
    dBdphi(3) = Rinv*dot_product(fh%psiv0, gh%drphi) - dot_product(fh%Bfv, gh%dzphiphi)
    dBdz(3) = dBdz(3) - dot_product(fh%Bfv, gh%dzzphi)

    dBdphi(2) = Rinv*dot_product(fh%Bzv0, gh%dphi)
#else
    dBdphi = 0.
#endif

    if (linear.eq.1) then
       !Perturbed part
       !B_poloidal = grad psi x grad phi
       temp(1)  = -Rinv*dot_product(fh%psiv1, gh%dz)
       tempR(1) = -Rinv*dot_product(fh%psiv1, gh%drz)
       tempz(1) = -Rinv*dot_product(fh%psiv1, gh%dzz)

       temp(3)  =  Rinv*dot_product(fh%psiv1, gh%dr)
       tempR(3) =  Rinv*dot_product(fh%psiv1, gh%drr)
       tempz(3) =  Rinv*dot_product(fh%psiv1, gh%drz)

       !B_toroidal = B_Z / R
       temp(2)  = Rinv*dot_product(fh%Bzv1, gh%g)
       tempR(2) = Rinv*dot_product(fh%Bzv1, gh%dr)
       tempz(2) = Rinv*dot_product(fh%Bzv1, gh%dz)

       if (itor.eq.1) tempR = tempR - Rinv*temp

#ifdef USECOMPLEX
       temp(1) = temp(1) - dot_product(fh%Bfv, gh%dr) * rfac
       temp(3) = temp(3) - dot_product(fh%Bfv, gh%dz) * rfac
       Bcyl = Bcyl + real(temp * exp(rfac*x(2)))

       tempR(1) = tempR(1) - dot_product(fh%Bfv, gh%drr) * rfac
       tempR(3) = tempR(3) - dot_product(fh%Bfv, gh%drz) * rfac
       dBdR = dBdR + real(tempR * exp(rfac*x(2)))

       tempz(1) = tempz(1) - dot_product(fh%Bfv, gh%drz) * rfac
       tempz(3) = tempz(3) - dot_product(fh%Bfv, gh%dzz) * rfac
       dBdz = dBdz + real(tempz * exp(rfac*x(2)))

       dBdphi = real(temp * rfac * exp(rfac*x(2)))
#else
       Bcyl = Bcyl + temp
       dBdR = dBdR + tempR
       dBdz = dBdz + tempz
#endif
    endif !linear
  end subroutine getBcylprime

!---------------------------------------------------------------------------
  subroutine getEcyl(fh, gh, Ecyl)
    use arrays
    use basic
    use auxiliary_fields
    implicit none

    type(elfield), intent(in) :: fh
    type(xgeomterms), intent(in) :: gh
    real, dimension(3), intent(out) :: Ecyl

    vectype, dimension(3) :: temp

    temp(1) = dot_product(fh%er, gh%g)
    temp(2) = dot_product(fh%ephi, gh%g)
    temp(3) = dot_product(fh%ez, gh%g)

#ifdef USECOMPLEX
    Ecyl = real(temp)
#else
    Ecyl = temp
#endif
  end subroutine getEcyl

!---------------------------------------------------------------------------
! Return particle kinetic energy, in Joules
  real function getke(p)
    use basic
    implicit none

    type(particle), intent(in) :: p

    type(elfield), dimension(nneighbors+1) :: elcoefs
    type(xgeomterms)                       :: geomterms
    real, dimension(3)                     :: B_cyl
    real                                   :: B0
    integer                                :: itri, tridex, ierr

    if (vspdims.eq.3) then
       getke = 0.5*e_mks*dot_product(p%v, p%v)/qm_ion
    else
       itri = p%jel;  elcoefs(:)%itri = 0
       call get_geom_terms(p%x, itri, elcoefs, tridex, geomterms, .false., ierr)
       if (ierr.ne.0) then
          getke = -1.0
          return
       endif

       call get_field_coefs(itri, elcoefs(1), .false.)
       call getBcyl(p%x, elcoefs(1), geomterms, B_cyl)

       B0 = sqrt(dot_product(B_cyl, B_cyl))
       getke = e_mks*(0.5*p%v(1)**2/qm_ion + p%v(2)*B0)
    endif
  end function getke

!---------------------------------------------------------------------------
! Return particle canonical angular momentum in kg-m**2/s
  real function getPphi(p)
    use arrays
    use basic
    implicit none

    type(particle), intent(in) :: p

    vectype                                :: psi
    type(elfield), dimension(nneighbors+1) :: elcoefs
    type(xgeomterms)                       :: geomterms
    real, dimension(3)                     :: B_cyl
    real                                   :: B0
    integer                                :: itri, tridex, ierr

    itri = p%jel;  elcoefs(:)%itri = 0
    call get_geom_terms(p%x, itri, elcoefs, tridex, geomterms, .false., ierr)
    if (ierr.ne.0) then
       getPphi = -1.0
       return
    endif

    call get_field_coefs(itri, elcoefs(1), .false.)

    ! Poloidal magnetic flux
    getPphi = e_mks * dot_product(elcoefs(1)%psiv0, geomterms%g)

    if (linear.eq.1) then
       psi = dot_product(elcoefs(1)%psiv1, geomterms%g)
#ifdef USECOMPLEX
       getPphi = getPphi + e_mks * real(psi * exp(rfac*p%x(2)))
#else
       getPphi = getPphi + e_mks * real(psi)
#endif
    endif !linear

    if (vspdims.eq.3) then
       getPphi = getPphi + (e_mks/qm_ion) * p%v(2) * p%x(1)
    else
       call getBcyl(p%x, elcoefs(1), geomterms, B_cyl)
       B0 = sqrt(dot_product(B_cyl, B_cyl))

       getPphi = getPphi + (e_mks/qm_ion) * p%v(1) * B_cyl(2) * p%x(1) / B0
    endif
  end function getPphi

!---------------------------------------------------------------------------
! Integrate over elements to compute kinetic ion contributions to RHS vectors
! for parallel and perpendicular components of pressure tensor.
  subroutine particle_pressure(p_par_i, p_perp_i, p_par_n, p_perp_n)
    use basic
    use math
    use field
    implicit none
    intrinsic matmul

    type(field_type), intent(inout) :: p_par_i, p_perp_i, p_par_n, p_perp_n

    real, dimension(dofs_per_element,coeffs_per_element) :: cl
    real, dimension(dofs_per_element) :: wnuhere
    vectype, dimension(dofs_per_element) :: dofspa, dofspe
#ifdef USECOMPLEX
    complex, dimension(dofs_per_element) :: dofspan, dofspen
    complex phfac
#endif
    real, dimension(3) :: B_part
    real, dimension(vspdims) :: vperp
    type(elfield), dimension(nneighbors+1) :: elcoefs
    type(xgeomterms) :: geomterms
    real             :: B0, vpar, ppar, pperp
    integer          :: ierr, nelms, ielm, ipart, itri, tridex, isghost

    nelms = size(pdata)
    elcoefs(:)%itri = 0

    p_par_i%vec = 0.;  p_perp_i%vec = 0.

    !Loop over all local+ghost elements to construct RHS; ghosts will be empty.
!$OMP PARALLEL
    do ielm=1,nelms
       call m3dc1_ent_isghost(2, ielm-1, isghost)
       if(isghost.eq.1) then
          if (pdata(ielm)%np.gt.0) print *,myrank,': nonzero ghost particle count!'
          cycle
       endif
       if (pdata(ielm)%np.lt.1) cycle !If no particles, then no pressure

       !Get basis function polynomial expansions
       if (iprecompute_metric.eq.1) then
          cl = ctri(:,:,ielm)
       else
          call local_coeff_vector(ielm, cl)
       endif

       !Need B at particle locations -> Load scalar fields for this element
       call get_field_coefs(ielm, elcoefs(1), .false.)

       dofspa = 0.;  dofspe = 0.
#ifdef USECOMPLEX
       dofspan = 0.;  dofspen = 0.
#endif

       !Sum over particles within this element
       do ipart=1,pdata(ielm)%np

          !Calculate B field at particle location
          itri = ielm
          call get_geom_terms(pdata(ielm)%ion(ipart)%x, itri, elcoefs, tridex, &
               geomterms, .false., ierr)
          if (ierr.ne.0) then
             print *,myrank,': Bad particle in pressure tensor integral; skipping.'
             cycle !next particle
          endif
          if (itri.ne.ielm) then
             print *,myrank,': Particle in wrong element in pressure tensor integral.'
             cycle !next particle
          endif
          call getBcyl(pdata(ielm)%ion(ipart)%x, elcoefs(tridex), geomterms, B_part)
          B0 = sqrt(dot_product(B_part, B_part))

          !Use B and v to get parallel and perp components of particle velocity
          if (vspdims.eq.2) then ! drift-kinetic: v_|| = v(1),  mu = q * v(2)
             vpar = pdata(ielm)%ion(ipart)%v(1)
             pperp = q_ion*pdata(ielm)%ion(ipart)%v(2)*B0
          else !full orbit: v_|| = v.B/|B|,  v_perp = v - v_||
             if (B0.gt.0.0) then !non-degenerate
                vpar = dot_product(pdata(ielm)%ion(ipart)%v, B_part(1:vspdims))/B0
                vperp = pdata(ielm)%ion(ipart)%v - (vpar/B0)*B_part(1:vspdims)
             else !degenerate case: no B field, pressure is scalar
                vpar = 0.0
                vperp = pdata(ielm)%ion(ipart)%v
             endif !degenerate?
             pperp = 0.5 * m_ion * dot_product(vperp, vperp)
          endif !full-orbit?
          ppar = m_ion * vpar**2

          !Add particle contribution to RHS vector (should vectorize well).
          wnuhere = pdata(ielm)%ion(ipart)%wt * matmul(cl, geomterms%g)
          dofspa = dofspa + ppar*wnuhere
          dofspe = dofspe + pperp*wnuhere
#ifdef USECOMPLEX
          !Extract appropriate Fourier component of particle contribution
          phfac = exp(rfac*pdata(ielm)%ion(ipart)%x(2))
          dofspan = dofspan + ppar*phfac*wnuhere
          dofspen = dofspen + pperp*phfac*wnuhere
#endif
       enddo !ipart

       !Insert element sums into field data
       ! Note: this is only correct if the local index ielm refers to the
       !  same element in meshes with and without ghost zone layers!
       call vector_insert_block(p_par_i%vec,  ielm, 1, dofspa, VEC_ADD)
       call vector_insert_block(p_perp_i%vec, ielm, 1, dofspe, VEC_ADD)
#ifdef USECOMPLEX
       call vector_insert_block(p_par_n%vec,  ielm, 1, dofspan, VEC_ADD)
       call vector_insert_block(p_perp_n%vec, ielm, 1, dofspen, VEC_ADD)
#endif
    enddo !ielm
!$OMP END PARALLEL

    !Normalize 2D toroidal integrals
#ifndef USE3D
    call mult(p_par_i, 1.0/twopi)
    call mult(p_perp_i, 1.0/twopi)
#ifdef USECOMPLEX
    call mult(p_par_n, 1.0/pi)
    call mult(p_perp_n, 1.0/pi)
#endif
#endif
  end subroutine particle_pressure

!---------------------------------------------------------------------------
  subroutine solve_pi_tensor(p_par_i, p_perp_i, p_par_n, p_perp_n)
    !use basic
    use field
    use newvar_mod
    implicit none

    type(field_type), intent(inout) :: p_par_i, p_perp_i, p_par_n, p_perp_n

    !print *,myrank,': solving for LHS vector...'
    call newvar_solve(p_par_i%vec,  mass_mat_lhs)
    call newvar_solve(p_perp_i%vec, mass_mat_lhs)
#ifdef USECOMPLEX
    call newvar_solve(p_par_n%vec,  mass_mat_lhs)
    call newvar_solve(p_perp_n%vec, mass_mat_lhs)
#endif
  end subroutine solve_pi_tensor

!---------------------------------------------------------------------------
! Recompute triangle coefficients for new element list
  subroutine reset_trimats
    use basic
    implicit none

    integer nelms

    nelms = local_elements()
#ifdef JBDEBUG
    print *,myrank,': redoing tridef for',nelms,' local elements.'
#endif
    if (nelms.lt.1) return

    deallocate(gtri,htri)
    if(iprecompute_metric.eq.1) deallocate(ctri)

    allocate(gtri(coeffs_per_tri,dofs_per_tri,nelms))
    allocate(htri(coeffs_per_dphi,dofs_per_dphi,nelms))
    if(iprecompute_metric.eq.1) &
         allocate(ctri(dofs_per_element,coeffs_per_element,nelms))

    call tridef
  end subroutine reset_trimats

#endif

end module particles
