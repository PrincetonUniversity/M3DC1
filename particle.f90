! Kinetic energetic ion module, J. Breslau, 2015
module particles
  implicit none

  real, parameter :: csquared = 8.98740441e+16  !Speed of light in m/s, squared
  real, parameter :: m_proton = 1.6726e-27  !Proton mass in kg
  real, parameter :: qm_proton = 9.5791e+07 !Proton charge/mass ratio in C/kg
  real, parameter :: vp1eV = 1.3841e+04     !Thermal speed of a 1 eV proton in m/s
  integer, parameter :: vspdims = 2         !Dimensions in velocity space (2=gk; 3=full)

  type particle
     real, dimension(3)       :: x              !Position in cylindrical coords
     real, dimension(vspdims) :: v              !Velocity
     real                     :: wt             !Particle weighting in delta-f scheme
     real                     :: tlast, dtlast  !Time, time increment
     integer                  :: fh             !file unit handle
  end type particle

  type(particle), dimension(:), allocatable :: ion
  real, dimension(3+vspdims) :: partscale
  real :: qm_ion
  real :: dbtime !For debugging purposes
  integer, dimension(:), allocatable :: llpid, lpcount !Lost particle tracking arrays
  integer :: nparticles, locparts
  integer :: mpi_particle !User-defined MPI datatype for particle communication

contains

  subroutine particle_test
    use basic
    use arrays
    use diagnostics
    use auxiliary_fields
    implicit none

    real, parameter :: JpereV = 1.6022e-19
    real :: fv, fv2, pdt, keeV, pphi
    integer :: itri=0, ierr, istep, ipart
    logical :: lop

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

    itri = 0
    call evaluate(xmag, 0.0, zmag, fv, fv2, bz_field(0), itri, ierr)
    if (ierr.eq.0.and.itri.ge.0) then
       print *,myrank,': bz vals = ',fv, fv2
    endif

    call calculate_auxiliary_fields(eqsubtract)

    !Initialize particle population
    call init_particles

    if (myrank.eq.3 .and. locparts.gt.1) then
       ion(2)%fh = 120
       open(unit=ion(2)%fh, file='ptraj3_2', status='replace', action='write')
       write(ion(2)%fh,*)'x  y  z  t  KE  Pphi'
       !write(ion(2)%fh,*)'x  y  z  t  B0  RB0'
    endif

    !Advance particle positions
    pdt = 9.625e-7
    dbtime = 0.0
    do istep=1,128
       call advance_particles(pdt)
       dbtime = dbtime + pdt
       do ipart=1,locparts
          if (ion(ipart)%fh.ge.0) then
             keeV = getke(ion(ipart))/JpereV
             pphi = getPphi(ion(ipart))/JpereV
             write(ion(ipart)%fh,'(3f14.8,3e18.8)')ion(ipart)%x(1)*cos(ion(ipart)%x(2)), &
                  ion(ipart)%x(1)*sin(ion(ipart)%x(2)), ion(ipart)%x(3), &
                  pdt*istep, keeV, pphi
             if (istep.eq.128) print *,'dtlast =',ion(ipart)%dtlast
          endif
       enddo !ipart
    enddo !istep

    !Close any open particle trajectory output files
    do ipart=1,locparts
       if (ion(ipart)%fh.ge.0) then
          inquire(ion(ipart)%fh, opened=lop)
          if (lop) close(ion(ipart)%fh)
       endif
    enddo

    !Clean up
    call finalize_particles
  end subroutine particle_test

!---------------------------------------------------------------------------
  subroutine init_particles
    use mesh_mod
    use basic
    implicit none

    include 'mpif.h'

    real, parameter :: c_mks = 2.9979e+8
    integer, dimension(6), parameter :: pblklen=(/3, vspdims, 1, 1, 1, 1/)
    integer(kind=MPI_ADDRESS_KIND), dimension(6) :: pdspls
    integer, dimension(6), parameter :: ptyps = (/MPI_DOUBLE, MPI_DOUBLE, &
         MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER/)

    type(particle) :: dumpar
    real, dimension(3) :: x_part, Bcyl
    real, dimension(coeffs_per_element):: gterm, drterm, dzterm
    real, dimension(vspdims) :: v_part
    real :: A_ion, Z_ion, m_ion, speed, vperp, vpar, EeV, B0
    real :: pitch_angle, gfrac=5.0e-3, gkfrac=32.0, Rinv=1.0
    integer :: itri=0, jj, globparts, ierr

    nparticles = 2
    allocate(ion(nparticles))
    allocate(llpid(nparticles))
    allocate(lpcount(maxrank))

    EeV = 100.0                          !Ion kinetic energy in eV
    A_ion = 1.0                         !Ion atomic mass number
    Z_ion = 1.0                         !Ion atomic number (charge state)
    m_ion = A_ion * m_proton            !Ion mass
    qm_ion = Z_ion*qm_proton/A_ion      !Ion charge/mass ratio
    speed = sqrt(EeV / A_ion) * vp1eV   !Ion speed
    if(myrank.eq.0)print *,'ion speed = ',speed,' m/s = ',speed/c_mks,' c.'
    !mu_max = 0.5 * m_ion * speed**2 / B0

    !Particle dimension scale sizes
    if (vspdims.eq.3) then !full orbit
       partscale(1) = speed/(qm_ion * B0)
       partscale(2) = partscale(1)
       partscale(3) = partscale(1)
       partscale(4) = speed
       partscale(5) = speed
       partscale(6) = speed
    else !drift kinetic
       partscale(1) = xlim
       partscale(2) = 6.283
       partscale(3) = partscale(1)
       partscale(4) = speed
       partscale(5) = 1.0
    endif

    locparts = 0
    do jj=1,nparticles
       !Realspace position in cylindrical coordinates
       x_part(1) = xmag + 0.1
       x_part(2) = 0.0
       x_part(3) = zmag + 0.2

       !Find the mesh domain and element containing this point, compute terms
       call get_geom_terms(x_part, itri, gterm, drterm, dzterm)

       if (itri.ge.0) then
          locparts = locparts + 1
          ion(locparts)%x = x_part
          !print *
          !print *,'x = ',x_part,':'
          !print *,'Particle ',jj,' is in element ',itri,' on PE ',myrank
          !print *,x_part(1:3:2),' -> ',si,zi

          !pitch_angle = 1.570796*(jj-1.0)/(nparticles-1.0)
          !pitch_angle = 1.570796e-1*(jj-1.0)/(nparticles-1.0)
          pitch_angle = 0.8 *(jj-1.0)/(nparticles-1.0)
          vpar  = speed * cos(pitch_angle)
          vperp = speed * sin(pitch_angle)

          !Need B to get v components from lambda
          if (itor.eq.1) Rinv = 1.0/x_part(1)
          call getBcyl(itri, gterm, drterm, dzterm, Rinv, x_part(2), Bcyl)
          B0 = sqrt(dot_product(Bcyl, Bcyl))
          !print *,'B0 ~ ',B0,' T.'

          if (vspdims.eq.3) then !full orbit
             v_part(1) = vpar*Bcyl(1)/B0 - vperp*Bcyl(2)/sqrt(Bcyl(1)**2 + Bcyl(2)**2)
             v_part(2) = vpar*Bcyl(2)/B0 + vperp*Bcyl(1)/sqrt(Bcyl(1)**2 + Bcyl(2)**2)
             v_part(3) = vpar*Bcyl(3)/B0
             !PRINT *,'v(R,phi,z) = ',v_part
             !print *,'speedrat = ',sqrt(dot_product(v_part,v_part))/speed
             !print *,'v . bhat = ',dot_product(v_part, Bcyl)/B0
             ion(locparts)%v = v_part

             ion(locparts)%dtlast = gfrac * (6.2831853 / (qm_ion * B0))
             print *,'Estimated gyroperiod = ',1.0e+9*ion(locparts)%dtlast/gfrac,' ns.'
          else !gyro- or drift kinetic
             ion(locparts)%v(1) = vpar                        !v_parallel
             ion(locparts)%v(2) = (0.5*vperp**2)/(qm_ion*B0)  !mu/q

             ion(locparts)%dtlast = gkfrac * (6.2831853 / (qm_ion * B0))
          endif

          ion(locparts)%fh = -1
       endif !itri...
    end do !jj

    print *,myrank,': ',locparts,' local particles.'

    call mpi_reduce(locparts, globparts, 1, MPI_INTEGER, MPI_SUM, &
         0, MPI_COMM_WORLD, ierr)
    if (myrank.eq.0 .and. globparts.ne.nparticles) then
       print *,'Warning: inconsistent particle count in init_particles!'
    endif

    !Define MPI datatype for particle communication
    call mpi_get_address(dumpar%x, pdspls(1), ierr)
    call mpi_get_address(dumpar%v, pdspls(2), ierr)
    pdspls(2) = pdspls(2) - pdspls(1)
    call mpi_get_address(dumpar%wt, pdspls(3), ierr)
    pdspls(3) = pdspls(3) - pdspls(1)
    call mpi_get_address(dumpar%tlast, pdspls(4), ierr)
    pdspls(4) = pdspls(4) - pdspls(1)
    call mpi_get_address(dumpar%dtlast, pdspls(5), ierr)
    pdspls(5) = pdspls(5) - pdspls(1)
    call mpi_get_address(dumpar%fh, pdspls(6), ierr)
    pdspls(6) = pdspls(6) - pdspls(1)
    pdspls(1) = 0
    !if (myrank.eq.0) print *,'pdspls = ',pdspls

    call mpi_type_create_struct(6, pblklen, pdspls, ptyps, mpi_particle, ierr)
    call mpi_type_commit(mpi_particle, ierr)
  end subroutine init_particles

!---------------------------------------------------------------------------
  subroutine finalize_particles
    implicit none

    if (allocated(ion)) deallocate(ion)
    if (allocated(llpid)) deallocate(llpid)
    if (allocated(lpcount)) deallocate(lpcount)
    !Free up particle MPI type?

  end subroutine finalize_particles

!---------------------------------------------------------------------------
  subroutine advance_particles(tinc)
    use basic  !For MPI variables
    implicit none
    include 'mpif.h'

    real, intent(in) :: tinc  !Time increment for particle advance

    real, parameter :: epsilon = 1.0e-6

    type(particle) :: testpart
    real           :: dtp, trem, dtnext, xi, zi
    integer        :: nexit, ipart, nstep, ierr, ipe, itri, globparts, uid
    logical        :: lop

    do ipart=1,locparts
       ion(ipart)%tlast = 0.0
    enddo

    do !Iterate until all particles are in the correct domain

       nexit = 0  !Number of local particles exiting domain during advance

       !Advance particles within local domain
       do ipart=1,locparts !For each local particle...
          nstep = 0
          dtp = ion(ipart)%dtlast
          dtnext = dtp

          do !Advance particle by tinc
             trem = tinc - ion(ipart)%tlast  !time remaining to advance
             if (trem.le.0.) exit
             if (dtp.gt.trem) dtp = trem

             call rk4(ion(ipart), dtp, ierr)
             !call rk5qs(ion(ipart), dtp, epsilon, partscale, dtnext, ierr)

             select case (ierr)
             case (1) ! Particle exited local domain
                !PRINT *,myrank,': Particle ',ipart,' exited at t=',ion(ipart)%tlast

                !Record particle ID and last time in table
                nexit = nexit + 1
                if (nexit.gt.nparticles) stop
                llpid(nexit) = ipart

                exit  !Break out of this loop, go to next particle.

             case (2) ! Adaptive integration step shrank to zero
                PRINT *,myrank,': Integration failed for particle ',ipart
                stop

             case default !ierr==0 -> Success
                ion(ipart)%tlast = ion(ipart)%tlast + ion(ipart)%dtlast
                dtp = dtnext
                nstep = nstep + 1
             end select
          enddo !tinc advance

          ion(ipart)%dtlast = dtp

          !if (myrank.eq.3) print *,myrank,': Advanced particle ',&
          !     ipart,nstep,' steps to t=',ion(ipart)%tlast,'/',tinc,&
          !     ', dtpf = ',dtp
       enddo !ipart

       !Reassign particles that have migrated out of local domain
       if (nexit.gt.0) print *,myrank,': nexit = ',nexit,':',llpid(1:nexit)
       call mpi_allgather(nexit, 1, MPI_INTEGER, lpcount, 1, MPI_INTEGER, &
            MPI_COMM_WORLD, ierr)
       if (SUM(lpcount).eq.0) exit !No more particles to reassign

       !Reassign transiting particles
       if (myrank.eq.0) print *,'reassigning ',SUM(lpcount),' particles.'
       do ipe=0,maxrank-1
          !print *,myrank,': awaiting ',lpcount(ipe+1),' particles from PE ',ipe
          do ipart=lpcount(ipe+1),1,-1
             if (myrank.eq.ipe) then
                !print *,myrank,': sending llpid = ',llpid(ipart)
                testpart = ion(llpid(ipart))
             endif
             call mpi_bcast(testpart, 1, mpi_particle, ipe, MPI_COMM_WORLD, ierr)
             !PRINT *,myrank,': x = ',testpart%x

             if (myrank.eq.ipe) then !Delete particle from local list
                !Close associated open output file, if present
                uid = ion(llpid(ipart))%fh
                if (uid.ge.0) then
                   inquire(uid, opened=lop)
                   if (lop) close(uid)
                endif

                do itri=llpid(ipart),locparts-1
                   ion(itri) = ion(itri+1)
                enddo
                locparts = locparts - 1
                !print *,myrank,': locparts =',locparts
                !llpid(ipart+1:nexit) = llpid(ipart+1:nexit) - 1
             else
                call whattri(testpart%x(1), testpart%x(2), testpart%x(3), &
                     itri, xi, zi)
                if (itri.gt.0) then
                   locparts = locparts + 1
                   if (locparts.gt.nparticles) stop
                   ion(locparts) = testpart
                   print *,'Particle moved to PE ',myrank,':',ion(locparts)%x

                   !Open associated file, if available
                   if (ion(locparts)%fh.ge.0) then
                      print *,myrank,': reopening file for appending...'
                      open(unit=ion(locparts)%fh, file='ptraj3_2', action='write', &
                           position='append')
                   endif
                endif !itri...
             endif !myrank...
          enddo !ipart          
       enddo !ipe

       call mpi_reduce(locparts, globparts, 1, MPI_INTEGER, MPI_SUM, &
            0, MPI_COMM_WORLD, ierr)
       if (myrank.eq.0) &
            print *,globparts,' particle(s) remaining after reassignment.'
    enddo !outer loop

  end subroutine advance_particles

!---------------------------------------------------------------------------
! 4th-order Runge-Kutta integrator, no adaptive time step control.
!  Four derivative evaluations per step.
! (Adapted from Numerical Recipes in C, 2nd edition, 1994, pp. 712-713).
!
  subroutine rk4(part, dt, ierr)
    implicit none

    type(particle), intent(inout) :: part
    real, intent(in) :: dt
    integer, intent(out) :: ierr

    real, parameter :: onethird = 1.0/3.0
    real, dimension(3) :: k1, k2, k3, k4, y1
    real, dimension(vspdims) :: l1, l2, l3, l4, z1
    real :: hh

    ierr = 0
    hh = 0.5*dt

    !1st step
    call fdot(part%x, part%v, k1, l1, ierr)
    if (ierr.ne.0) return
    y1 = part%x + hh*k1;  z1 = part%v + hh*l1

    !2nd step
    call fdot(y1, z1, k2, l2, ierr)
    if (ierr.ne.0) return
    y1 = part%x + hh*k2;  z1 = part%v + hh*l2

    !3rd step
    call fdot(y1, z1, k3, l3, ierr)
    if (ierr.ne.0) return
    y1 = part%x + dt*k3;  z1 = part%v + dt*l3

    !4th step
    call fdot(y1, z1, k4, l4, ierr)
    if (ierr.ne.0) return
    part%x = part%x + onethird*dt*(k2 + k3 + 0.5*(k1 + k4))
    part%v = part%v + onethird*dt*(l2 + l3 + 0.5*(l1 + l4))
    part%dtlast = dt
  end subroutine rk4

!----------------------------------------------------------------------------
! 5th-order Runge-Kutta step with monitoring of truncation error for accuracy,
! step size adjustment.
! (Adapted from Numerical Recipes in C, 2nd edition, 1994, p. 719).
!
  subroutine rk5qs(part0, htry, eps, yscal, hnext, ierr)
    implicit none

    type(particle), intent(inout) :: part0
    real, intent(in) :: htry, eps
    real, dimension(3+vspdims), intent(in) :: yscal
    real, intent(out) :: hnext
    integer, intent(out) :: ierr

    real, parameter :: safety = 0.9, pshrnk = -0.25, pgrow = -0.2
    real, parameter :: errcon = (5.0/safety)**(1.0/pgrow)

    type(particle) :: ptemp
    real, dimension(3+vspdims) :: perr
    real :: h, errmax, htemp

    ierr = 0
    ptemp = part0  !Copy over any particle data aside from x,v

    h = htry  !Set step size to initial trial value
    do
       call rk5ck(part0, h, ptemp, perr, ierr)
       if (ierr.ne.0) return
       errmax = MAXVAL(ABS(perr/yscal))/eps
       if (errmax.le.1.0) exit  !Step succeeded.

       !Error too large -> step failed.
       !Reduce step size by no more than 10x, retry.
       htemp = safety*h*(errmax**pshrnk)
       h = 0.1*h
       if (htemp.gt.h) h = htemp

       if (htry+h.eq.htry) then !Step size underflow condition
          ierr = 2
          return
       endif
    enddo

    !Increase step size by no more than 5x.
    if (errmax.gt.errcon) then
       hnext = safety*h*(errmax**pgrow)
    else
       hnext = 5.0*h
    endif

    part0 = ptemp
    part0%dtlast = h
  end subroutine rk5qs

!----------------------------------------------------------------------------------
! 5th-order Cash-Karp Runge-Kutta integrator with embedded 4th-order error estimate.
!  Six derivative evaluations per call.
! (Adapted from Numerical Recipes in C, 2nd edition, 1994, pp. 719-720).
!
  subroutine rk5ck(part0, dt, part1, perr, ierr)
    implicit none

    type(particle), intent(in) :: part0
    real, intent(in) :: dt
    type(particle), intent(out) :: part1
    real, dimension(3+vspdims), intent(out) :: perr
    integer, intent(out) :: ierr

    real, parameter :: b21=0.2, b31=0.075, b32=0.225
    real, parameter :: b41=0.3, b42=-0.9, b43=1.2
    real, parameter :: b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0
    real, parameter :: b61=1631.0/55296.0, b62=175.0/512.0
    real, parameter :: b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0
    real, parameter :: c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0
    real, parameter :: c6=512.0/1771.0, dc5=-277.0/14336.0
    real, parameter :: dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0
    real, parameter :: dc4=c4-13525.0/55296.0, dc6=c6-0.25

    real, dimension(3)       :: xdot, y1, ak2, ak3, ak4, ak5, ak6
    real, dimension(vspdims) :: vdot, z1, bk2, bk3, bk4, bk5, bk6

    ierr = 0

    !1st step
    call fdot(part0%x, part0%v, xdot, vdot, ierr)
    if (ierr.ne.0) return
    y1 = part0%x + b21*dt*xdot
    z1 = part0%v + b21*dt*vdot

    !2nd step
    call fdot(y1, z1, ak2, bk2, ierr)
    if (ierr.ne.0) return
    y1 = part0%x + dt*(b31*xdot + b32*ak2)
    z1 = part0%v + dt*(b31*vdot + b32*bk2)

    !3rd step
    call fdot(y1, z1, ak3, bk3, ierr)
    if (ierr.ne.0) return
    y1 = part0%x + dt*(b41*xdot + b42*ak2 + b43*ak3)
    z1 = part0%v + dt*(b41*vdot + b42*bk2 + b43*bk3)

    !4th step
    call fdot(y1, z1, ak4, bk4, ierr)
    if (ierr.ne.0) return
    y1 = part0%x + dt*(b51*xdot + b52*ak2 + b53*ak3 + b54*ak4)
    z1 = part0%v + dt*(b51*vdot + b52*bk2 + b53*bk3 + b54*bk4)

    !5th step
    call fdot(y1, z1, ak5, bk5, ierr)
    if (ierr.ne.0) return
    y1 = part0%x + dt*(b61*xdot + b62*ak2 + b63*ak3 + b64*ak4 + b65*ak5)
    z1 = part0%v + dt*(b61*vdot + b62*bk2 + b63*bk3 + b64*bk4 + b65*bk5)

    !6th step
    call fdot(y1, z1, ak6, bk6, ierr)
    if (ierr.ne.0) return
    part1%x = part0%x + dt*(c1*xdot + c3*ak3 + c4*ak4 + c6*ak6)
    part1%v = part0%v + dt*(c1*vdot + c3*bk3 + c4*bk4 + c6*bk6)

    !Error estimate
    perr(1:3) = dt*(dc1*xdot + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)
    perr(4:3+vspdims) = dt*(dc1*vdot + dc3*bk3 + dc4*bk4 + dc5*bk5 + dc6*bk6)
  end subroutine rk5ck

!---------------------------------------------------------------------------
  subroutine fdot(x, v, dxdt, dvdt, ierr)
    use mesh_mod
    use basic
    implicit none

    real, dimension(3), intent(in)        :: x
    real, dimension(3), intent(out)       :: dxdt
    real, dimension(vspdims), intent(in)  :: v
    real, dimension(vspdims), intent(out) :: dvdt
    integer, intent(out)                  :: ierr

    real, parameter :: g_mks = 9.8067 ! earth avg surf grav accel in m/s/s
    real, dimension(coeffs_per_element):: gterm, drterm, dzterm
    real, dimension(coeffs_per_element):: drrterm, drzterm, dzzterm
    real, dimension(3) :: B_cyl, E_cyl, bhat, svec, Bstar
    real, dimension(3) :: dBdR, dBdphi, dBdz, gradB0
    real :: Rinv = 1.0, B0, Bss
    integer :: itri=0

    ierr = 0

    !Need fields to calculate acceleration
    if (vspdims.eq.3) then
       call get_geom_terms(x, itri, gterm, drterm, dzterm)
    else
       call get_geom_terms2(x, itri, gterm, drterm, dzterm,&
            drrterm, drzterm, dzzterm)
    endif
    if (itri.le.0) then
       !print *,myrank,': particle exited local domain!'
       ierr = 1
       return
    endif
    !if (myrank.eq.3) print *,'dbt=',dbtime,': itri = ',itri !DEBUG

    !Get magnetic, electric field components
    if (itor.eq.1) Rinv = 1.0/x(1)
    call getEcyl(itri, gterm, x(2), E_cyl)
    !E_cyl = 0. !DEBUG ONLY!!!

    if (vspdims.eq.3) then !full orbit: ma = q(E + vxB) + mg
       dxdt(1:vspdims) = v

       call getBcyl(itri, gterm, drterm, dzterm, Rinv, x(2), B_cyl)

       dvdt(1) = qm_ion*(E_cyl(1) + v(2)*B_cyl(3) - v(3)*B_cyl(2))
       dvdt(2) = qm_ion*(E_cyl(2) + v(3)*B_cyl(1) - v(1)*B_cyl(3))
       dvdt(3) = qm_ion*(E_cyl(3) + v(1)*B_cyl(2) - v(2)*B_cyl(1)) - g_mks

       !if (itor.eq.1) then
       !   dvdt(1) = dvdt(1) + Rinv*v(2)**2  !Centripetal acceleration
       !   dvdt(2) = dvdt(2) - 2.0*Rinv*v(1)*v(2)  !Coriolis effect
       !endif
    else ! Drift-kinetic equation
       call getBcyl2(itri, gterm, drterm, dzterm, drrterm, drzterm, dzzterm, &
            Rinv, x(2), B_cyl, dBdR, dBdphi, dBdz)

       B0 = sqrt(dot_product(B_cyl, B_cyl))
       bhat = B_cyl / B0

       ! Gradient of mod B
       gradB0(1) = dot_product(bhat, dBdR)
       gradB0(2) = Rinv*dot_product(bhat, dBdphi)
       gradB0(3) = dot_product(bhat, dBdz)

       ! Curl of bhat
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

    dxdt(2) = Rinv*dxdt(2)  !phi-dot = v_phi / R for cylindrical case
  end subroutine fdot

!---------------------------------------------------------------------------
! Compute terms for a function and its partial derivatives with respect to
!  R and z in the reduced quintic expansion at position x.
  subroutine get_geom_terms(x, itri, gterm, drterm, dzterm)
    use mesh_mod
    implicit none

    real, dimension(3), intent(in) :: x
    integer, intent(out) :: itri
    real, dimension(coeffs_per_element), intent(out) :: gterm, drterm, dzterm

    type(element_data) :: eldat
    real :: xi, zi, eta, dxi, deta
    integer :: pp

    itri = 0
    call whattri(x(1),x(2),x(3),itri,xi,zi)
    if (itri.le.0) return

    call get_element_data(itri, eldat)
    call global_to_local(eldat, x(1), x(2), x(3), xi, zi, eta)
    !print *,'Coords in element: xi,eta =',xi,eta

    drterm = 0.;  dzterm = 0.
    do pp=1,coeffs_per_element
       gterm(pp) = xi**mi(pp) * eta**ni(pp)

       if (mi(pp).gt.0) then
          dxi = mi(pp) * xi**(mi(pp)-1) * eta**ni(pp)

          drterm(pp) = drterm(pp) + eldat%co*dxi
          dzterm(pp) = dzterm(pp) + eldat%sn*dxi
       endif

       if (ni(pp).gt.0) then
          deta = xi**mi(pp) * ni(pp)*eta**(ni(pp)-1)

          drterm(pp) = drterm(pp) - eldat%sn*deta
          dzterm(pp) = dzterm(pp) + eldat%co*deta
       endif
    enddo !pp

  end subroutine get_geom_terms

!---------------------------------------------------------------------------
! Compute 2nd partial derivative terms as well.
  subroutine get_geom_terms2(x, itri, gterm, drterm, dzterm, drr, drz, dzz)
    use mesh_mod
    implicit none

    real, dimension(3), intent(in) :: x
    integer, intent(out) :: itri
    real, dimension(coeffs_per_element), intent(out) :: gterm, drterm, dzterm
    real, dimension(coeffs_per_element), intent(out) :: drr, drz, dzz

    type(element_data) :: eldat
    real :: xi, zi, eta, dxi, deta, d2xi, d2eta, dxieta
    integer :: pp

    itri = 0
    call whattri(x(1),x(2),x(3),itri,xi,zi)
    if (itri.le.0) return

    call get_element_data(itri, eldat)
    call global_to_local(eldat, x(1), x(2), x(3), xi, zi, eta)

    drterm = 0.;  dzterm = 0.;  drr = 0.;  drz = 0.;  dzz = 0.

    do pp=1,coeffs_per_element
       gterm(pp) = xi**mi(pp) * eta**ni(pp)

       if (mi(pp).gt.0) then
          dxi = mi(pp) * xi**(mi(pp)-1) * eta**ni(pp)

          drterm(pp) = drterm(pp) + eldat%co*dxi
          dzterm(pp) = dzterm(pp) + eldat%sn*dxi

          if (mi(pp).gt.1) then
             d2xi = mi(pp)*(mi(pp)-1)*xi**(mi(pp)-2) * eta**ni(pp)

             drr(pp) = drr(pp) + d2xi*eldat%co**2
             drz(pp) = drz(pp) + d2xi*eldat%co*eldat%sn
             dzz(pp) = dzz(pp) + d2xi*eldat%sn**2
          endif

          if (ni(pp).gt.0) then
             dxieta = mi(pp)*ni(pp) * xi**(mi(pp)-1) * eta**(ni(pp)-1)

             drr(pp) = drr(pp) - 2.0*dxieta*eldat%co*eldat%sn
             drz(pp) = drz(pp) + dxieta*(2.0*eldat%co**2 - 1.0)
             dzz(pp) = dzz(pp) + 2.0*dxieta*eldat%co*eldat%sn
          endif
       endif

       if (ni(pp).gt.0) then
          deta = xi**mi(pp) * ni(pp)*eta**(ni(pp)-1)

          drterm(pp) = drterm(pp) - eldat%sn*deta
          dzterm(pp) = dzterm(pp) + eldat%co*deta

          if (ni(pp).gt.1) then
             d2eta = xi**mi(pp) * ni(pp)*(ni(pp)-1)*eta**(ni(pp)-2)

             drr(pp) = drr(pp) + d2eta*eldat%sn**2
             drz(pp) = drz(pp) - d2eta*eldat%co*eldat%sn
             dzz(pp) = dzz(pp) + d2eta*eldat%co**2
          endif
       endif
    enddo !pp

  end subroutine get_geom_terms2

!---------------------------------------------------------------------------
  subroutine getBcyl(itri, gterm, drterm, dzterm, Rinv, phi, Bcyl)
    use arrays
    use basic
    implicit none

    integer, intent(in) :: itri
    real, dimension(coeffs_per_element), intent(in) :: gterm, drterm, dzterm
    real, intent(in) :: Rinv, phi
    real, dimension(3), intent(out) :: Bcyl

    vectype, dimension(coeffs_per_element) :: avector
    vectype, dimension(3) :: temp

    if(itri.eq.0) return

    !Equilibrium part
    !B_poloidal = grad psi x grad phi
    call calcavector(itri, psi_field(0), avector)
    Bcyl(1) = -Rinv*dot_product(avector, dzterm)
    Bcyl(3) =  Rinv*dot_product(avector, drterm)

    !B_toroidal = B_Z / R
    call calcavector(itri, bz_field(0), avector)
    Bcyl(2) = Rinv*dot_product(avector, gterm)

    if (linear.eq.1) then
       !Perturbed part

       !B_poloidal = grad psi x grad phi
       call calcavector(itri, psi_field(1), avector)
       temp(1) = -Rinv*dot_product(avector, dzterm)
       temp(3) =  Rinv*dot_product(avector, drterm)

       !B_toroidal = B_Z / R
       call calcavector(itri, bz_field(1), avector)
       temp(2) = Rinv*dot_product(avector, gterm)

#ifdef USECOMPLEX
       Bcyl = Bcyl + real(temp * exp(rfac*phi))
#else
       Bcyl = Bcyl + temp
#endif
    endif !linear
  end subroutine getBcyl

!---------------------------------------------------------------------------
  subroutine getBcyl2(itri, gterm, drterm, dzterm, drr, drz, dzz, Rinv, phi, &
       Bcyl, dBdR, dBdphi, dBdz)
    use arrays
    use basic
    implicit none

    integer, intent(in) :: itri
    real, dimension(coeffs_per_element), intent(in) :: gterm, drterm, dzterm
    real, dimension(coeffs_per_element), intent(in) :: drr, drz, dzz
    real, intent(in) :: Rinv, phi
    real, dimension(3), intent(out) :: Bcyl, dBdR, dBdphi, dBdz

    vectype, dimension(coeffs_per_element) :: avector
    vectype, dimension(3) :: temp, tempR, tempz

    if(itri.eq.0) return

    !Equilibrium part
    !B_poloidal = grad psi x grad phi
    call calcavector(itri, psi_field(0), avector)
    Bcyl(1) = -Rinv*dot_product(avector, dzterm)
    Bcyl(3) =  Rinv*dot_product(avector, drterm)

    dBdR(1) = -Rinv*dot_product(avector, drz)
    dBdR(3) =  Rinv*dot_product(avector, drr)

    dBdz(1) = -Rinv*dot_product(avector, dzz)
    dBdz(3) =  Rinv*dot_product(avector, drz)

    !B_toroidal = B_Z / R
    call calcavector(itri, bz_field(0), avector)
    Bcyl(2) = Rinv*dot_product(avector, gterm)

    dBdR(2) = Rinv*dot_product(avector, drterm)
    dBdz(2) = Rinv*dot_product(avector, dzterm)

    if (itor.eq.1) dBdR = dBdR - Rinv*Bcyl

    dBdphi = 0.

    if (linear.eq.1) then
       !Perturbed part
       !B_poloidal = grad psi x grad phi
       call calcavector(itri, psi_field(1), avector)
       temp(1) = -Rinv*dot_product(avector, dzterm)
       temp(3) =  Rinv*dot_product(avector, drterm)

       tempR(1) = -Rinv*dot_product(avector, drz)
       tempR(3) =  Rinv*dot_product(avector, drr)

       tempz(1) = -Rinv*dot_product(avector, dzz)
       tempz(3) =  Rinv*dot_product(avector, drz)

       !B_toroidal = B_Z / R
       call calcavector(itri, bz_field(1), avector)
       temp(2) = Rinv*dot_product(avector, gterm)

       tempR(2) = Rinv*dot_product(avector, drterm)
       tempz(2) = Rinv*dot_product(avector, dzterm)

       if (itor.eq.1) tempR = tempR - Rinv*temp

#ifdef USECOMPLEX
       Bcyl = Bcyl + real(temp * exp(rfac*phi))
       dBdR = dBdR + real(tempR * exp(rfac*phi))
       dBdz = dBdz + real(tempz * exp(rfac*phi))
       dBdphi = real(temp * rfac * exp(rfac*phi))
#else
       Bcyl = Bcyl + temp
       dBdR = dBdR + tempR
       dBdz = dBdz + tempz
#endif
    endif !linear
  end subroutine getBcyl2

 !---------------------------------------------------------------------------
  subroutine getEcyl(itri, gterm, phi, Ecyl)
    use arrays
    use basic
    use auxiliary_fields
    implicit none

    integer, intent(in) :: itri
    real, dimension(coeffs_per_element), intent(in) :: gterm
    real, intent(in) :: phi
    real, dimension(3), intent(out) :: Ecyl

    vectype, dimension(coeffs_per_element) :: avector
    vectype, dimension(3) :: temp

    if(itri.eq.0) return

    call calcavector(itri, ef_r, avector)
    temp(1) = dot_product(avector, gterm)

    call calcavector(itri, ef_phi, avector)
    temp(2) = dot_product(avector, gterm)

    call calcavector(itri, ef_z, avector)
    temp(3) = dot_product(avector, gterm)

#ifdef USECOMPLEX
    Ecyl = real(temp) ! * exp(rfac*phi))
#else
    Ecyl = temp
#endif
  end subroutine getEcyl

!---------------------------------------------------------------------------
! Return mod B in Tesla
  function getB0(x)
    use basic
    implicit none

    real :: getB0
    real, dimension(3), intent(in) :: x

    real, dimension(coeffs_per_element):: gterm, drterm, dzterm
    real, dimension(3) :: B_cyl
    real :: Rinv = 1.0
    integer :: itri

    call get_geom_terms(x, itri, gterm, drterm, dzterm)
    if (itri.le.0) then
       getB0 = -1.0
       return
    endif

    if (itor.eq.1) Rinv = 1.0/x(1)
    call getBcyl(itri, gterm, drterm, dzterm, Rinv, x(2), B_cyl)
    getB0 = sqrt(dot_product(B_cyl, B_cyl))
  end function getB0

!---------------------------------------------------------------------------
! Return particle kinetic energy, in Joules
  function getke(p)
    use basic
    implicit none

    real :: getke
    type(particle), intent(in) :: p

    real, parameter :: e_p = 1.6022e-19  !Elementary charge, in Coulombs
    real, dimension(coeffs_per_element):: gterm, drterm, dzterm
    real, dimension(3) :: B_cyl
    real :: Rinv = 1.0, B0
    integer :: itri

    if (vspdims.eq.3) then
       getke = 0.5*e_p*dot_product(p%v, p%v)/qm_ion
    else
       call get_geom_terms(p%x, itri, gterm, drterm, dzterm)
       if (itri.le.0) then
          getke = -1.0
          return
       endif

       if (itor.eq.1) Rinv = 1.0/p%x(1)
       call getBcyl(itri, gterm, drterm, dzterm, Rinv, p%x(2), B_cyl)
       B0 = sqrt(dot_product(B_cyl, B_cyl))

       getke = e_p*(0.5*p%v(1)**2/qm_ion + p%v(2)*B0)
    endif
  end function getke

!---------------------------------------------------------------------------
! Return particle canonical angular momentum in kg-m**2/s
  function getPphi(p)
    use arrays
    use basic
    implicit none

    real :: getPphi
    type(particle), intent(in) :: p

    real, parameter :: e_p = 1.6022e-19  !Elementary charge, in Coulombs
    vectype, dimension(coeffs_per_element) :: avector
    vectype :: psi
    real, dimension(coeffs_per_element):: gterm, drterm, dzterm
    real, dimension(3) :: B_cyl
    real :: Rinv=1.0, B0
    integer :: itri

    call get_geom_terms(p%x, itri, gterm, drterm, dzterm)
    if (itri.le.0) then
       getPphi = -1.0
       return
    endif

    ! Poloidal magnetic flux
    call calcavector(itri, psi_field(0), avector)
    getPphi = e_p * dot_product(avector, gterm)

    if (linear.eq.1) then
       call calcavector(itri, psi_field(1), avector)
       psi = dot_product(avector, gterm)
#ifdef USECOMPLEX
       getPphi = getPphi + e_p * real(psi * exp(rfac*p%x(2)))
#else
       getPphi = getPphi + e_p * psi
#endif
    endif !linear
    !if (myrank.eq.3) print *,'psi = ',getPphi/e_p

    if (vspdims.eq.3) then
       getPphi = getPphi + (e_p/qm_ion) * p%v(2) * p%x(1)
    else
       if (itor.eq.1) Rinv = 1.0/p%x(1)
       call getBcyl(itri, gterm, drterm, dzterm, Rinv, p%x(2), B_cyl)
       B0 = sqrt(dot_product(B_cyl, B_cyl))

       getPphi = getPphi + (e_p/qm_ion) * p%v(1) * B_cyl(2) * p%x(1) / B0
    endif
  end function getPphi
end module particles
