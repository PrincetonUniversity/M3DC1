program main
  use kprad
  
  implicit none
  integer :: ZIMP,i,j,ntimemax,ntime,nzones
  real :: dts,facimp,time
  real, allocatable, dimension(:) :: ne,te,pbrem, &
       aimp,bimp,cimp,dimp,ework,fwork,ne0,tekev
  real, allocatable, dimension(:,:) :: sion,srec, &
       imp_rad,pion,prec,nzeff,nz,nzt
  
  open(unit=10,file='rad.txt')
  smult = 1.
  ZIMP = 6    ! CARBON
  !       ZIMP = 10   ! NEON
  !      ZIMP = 18   ! Argon
  allocate(aimp(0:zimp),bimp(0:zimp),cimp(0:zimp),dimp(0:zimp)  &
       ,ework(0:zimp),fwork(0:zimp))
  
  call kprad_atomic_data_sub(ZIMP)
  dts = 1.e-12
  ntimemax = 1000
  nzones = 1   ! one zone
  j = 1
  allocate(ne(nzones),ne0(nzones),te(nzones), &
       pbrem(nzones),tekev(nzones))
  ne(j) = 1.e14   ! density in cm-3
  te(j) = 1000.   ! temp in ev
  tekev(j) = 1.e-3*te(j)
  ne0(j) = 0.     ! not used ??
  allocate(sion(0:zimp-1,nzones))
  allocate(srec(0:zimp,nzones))
  allocate(imp_rad(nzones,zimp+1),pion(nzones,zimp+1), &
       prec(nzones,zimp+1))
  allocate(nzeff(nzones,2))
  allocate(nz(0:zimp,nzones),nzt(0:zimp,nzones))
  facimp = .01
  !     initialize
  call kprad_ionization_rate(nzones,m2,ne,te,zimp, &
       sion_coeff,sion)
  call kprad_recombination_rate(nzones,ne,zed,te, &
       z_ei,zimp,srec)
  nz(0,j) = facimp*ne(1)
  do i=1,zimp
     nz(i,j) = 0.0
  enddo
  time = 0.
  
  !     start time loop
  do ntime=1,ntimemax
     time = time + dts
     do i=0,zimp
        nzt(i,j) = nz(i,j)
     enddo
     
     aimp(0) = 0.0
     do i = 0,zimp
        if(i.gt.0) aimp(i) = -dts*sion(i-1,j)
        bimp(i) =  1. + dts*srec(i,j)
        if(i.lt.zimp) then
           bimp(i) = bimp(i) + dts*sion(i,j)
           cimp(i) = -dts*srec(i+1,j)
        end if
        dimp(i) = nzt(i,j)
     enddo
     cimp(zimp) = 0.0
     !
     call tridiag(aimp,bimp,cimp,dimp,nz(:,j), &
          ework,fwork,zimp)
     !
     call kprad_energy_loses(nzones,m1,zimp,te, &
          ne,ne0,c,sion,srec,z_ei, &
          nz,nzeff,pion,prec,imp_rad,pbrem)
     !        write output in SI units
     write(10,1000) time,(nz(i,j)*1.e6,i=0,zimp), &
          imp_rad(j,zimp+1)*1.e6,pbrem(j)*1.e6
     dts = dts * 1.02
  enddo
  
  stop
1000 format(1p24e12.4)

end program main
