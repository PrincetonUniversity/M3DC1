program main
  use kprad
  
  implicit none
  integer :: zimp, i, j, ntimemax,ntime,nzones
  real :: dts, facimp, time
  real, allocatable, dimension(:) :: ne,te,tekev,pbrem,dw_brem
  real, allocatable, dimension(:,:) :: imp_rad,dw_rad,nz
  
  open(unit=10,file='rad.txt')

  !zimp = 6    ! CARBON
  ZIMP = 10   ! NEON
  !      ZIMP = 18   ! Argon
  
  call kprad_atomic_data_sub(zimp)
  dts = 1.e-7
  ntimemax = 100
  nzones = 1   ! one zone
  j = 1

  allocate(ne(nzones),te(nzones),tekev(nzones))
  ne(j) = 1.e14   ! density in cm-3
  te(j) = 1000.   ! temp in ev
  tekev(j) = 1.e-3*te(j)
  allocate(nz(nzones,0:zimp))
  allocate(pbrem(nzones), imp_rad(nzones,zimp+1), &
       dw_brem(nzones), dw_rad(nzones,zimp+1))

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

     call kprad_advance_densities(dts, nzones, zimp, ne, te, nz, &
          dw_rad, dw_brem)

     ! calculate power from radiated energy
     pbrem = dw_brem / dts
     imp_rad = dw_rad / dts

     write(10,1000) time,(nz(j,i)*1.e6,i=0,zimp), &
          imp_rad(j,zimp+1)*1.e6,pbrem(j)*1.e6

     dts = dts * 1.1
  enddo
  close(unit=10)

  call kprad_deallocate()

  deallocate(nz, pbrem, imp_rad, dw_brem, dw_rad)

  
  stop
1000 format(1p24e12.4)

end program main
