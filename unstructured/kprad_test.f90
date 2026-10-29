program main
  use kprad
  
  implicit none
  integer :: zimp, i, j, ntimemax,ntime,nzones, ierr
  real :: dts, facimp, time
  real, allocatable, dimension(:) :: ne, te, den, ti, p, pbrem, dw_brem
  real, allocatable, dimension(:,:) :: imp_rad, dw_rad, dw_ion, dw_reck, &
       dw_recp, nz, source
  
  open(unit=10,file='rad.txt')

  !zimp = 6    ! CARBON
  ZIMP = 10   ! NEON
  !      ZIMP = 18   ! Argon
  
  ikprad = 1
  ikprad_min_option = 1
  ikprad_max_dt = 0
  kprad_max_dt = -1.
  ikprad_evolve_internal = 0

  call kprad_allocate(zimp)
  call kprad_atomic_data_sub(zimp, ierr)
  dts = 1.e-7
  ntimemax = 100
  nzones = 1   ! one zone
  j = 1

  allocate(ne(nzones), te(nzones), den(nzones), ti(nzones), p(nzones))
  ne(j) = 1.e14   ! density in cm-3
  te(j) = 1000.   ! temp in ev
  den(j) = ne(j)
  ti(j) = te(j)
  p(j) = ne(j)*te(j)*1.6022e-12
  allocate(nz(nzones,0:zimp))
  allocate(pbrem(nzones), imp_rad(nzones,zimp+1), &
       dw_brem(nzones), dw_rad(nzones,zimp+1), &
       dw_ion(nzones,0:zimp), dw_reck(nzones,0:zimp), &
       dw_recp(nzones,0:zimp), source(nzones,0:zimp))
  source = 0.

  !     initialize
  facimp = .01
  nz(j,0) = facimp*ne(1)
  do i=1,zimp
     nz(j,i) = 0.0
  enddo
  time = 0.
  
  !     start time loop
  do ntime=1, ntimemax
     time = time + dts

     call kprad_advance_densities(dts, nzones, zimp, p, ne, te, den, ti, &
          nz, dw_rad, dw_brem, dw_ion, dw_reck, &
          dw_recp, source)

     ! calculate power from radiated energy
     pbrem = dw_brem / dts
     imp_rad = dw_rad / dts

     write(10,1000) time,(nz(j,i)*1.e6,i=0,zimp), &
          imp_rad(j,zimp+1)*1.e6,pbrem(j)*1.e6

     dts = dts * 1.1
  enddo
  close(unit=10)

  call kprad_deallocate()

  deallocate(ne, te, den, ti, p, nz, pbrem, imp_rad, dw_brem, dw_rad, &
       dw_ion, dw_reck, dw_recp, source)

  
  stop
1000 format(1p24e12.4)

end program main
