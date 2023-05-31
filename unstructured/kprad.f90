module kprad

  ! All densities in cm^-3
  ! All temperatures in eV

  implicit none

  real, allocatable, private, dimension(:) :: z_ei, zed
  real, allocatable, private, dimension(:,:) :: c, sion_coeff

  ! mass of chosen impurity species (in amu)
  real :: kprad_mz

  integer :: ikprad_max_dt ! use max dt in KPRAD evolution
  real :: kprad_max_dt
  integer :: ikprad_evolve_internal


  ! polynomial order for evaluating 
  ! radiation and ionization rates, respectively
  integer, private :: m1, m2

  real, private, parameter :: kprad_min_dt_default = 1e10
  real, private :: kprad_min_dt = kprad_min_dt_default     ! minimum final timestep from previous calculation
  real :: kprad_dt = 1e-10      ! kprad integration time step (in seconds)

  ! minimum values for KPRAD evolution
  integer :: ikprad_min_option
  real :: kprad_nemin
  real :: kprad_temin
  real, allocatable, private, dimension(:) :: ne_int, te_int

contains

  subroutine kprad_rebase_dt()
    implicit none

    include 'mpif.h'

    integer :: ier
    real :: temp

    call mpi_allreduce(kprad_min_dt, temp, 1, MPI_DOUBLE_PRECISION, &
         MPI_MIN, MPI_COMM_WORLD, ier)

    kprad_dt = temp
    kprad_min_dt = kprad_min_dt_default
  end subroutine kprad_rebase_dt

  subroutine kprad_deallocate()
    implicit none

    if(allocated(z_ei)) deallocate(z_ei)
    if(allocated(zed)) deallocate(zed)
    if(allocated(c)) deallocate(c)
    if(allocated(sion_coeff)) deallocate(sion_coeff)
    if(allocated(ne_int)) deallocate(ne_int)
    if(allocated(te_int)) deallocate(te_int)

  end subroutine kprad_deallocate

  subroutine kprad_instantaneous_radiation(npts, z, ne, te, nz, prad, pbrem)
    implicit none

    integer, intent(in) :: npts
    integer, intent(in) :: z
    real, intent(in) :: ne(npts)             ! electron density in cm^-3
    real, intent(in) :: te(npts)             ! electron temperature in eV
    real, intent(inout) :: nz(npts,0:z)      ! density
    real, intent(out) :: prad(npts,0:z)    ! energy lost via radiation
    real, intent(out) :: pbrem(npts)       ! energy lost via bremsstrahlung
    
    real, dimension(npts,0:z-1) :: sion
    real, dimension(npts,0:z) :: srec
    real, dimension(npts,z+1) :: pion, preck, precp
    real, dimension(npts, 2) :: nzeff

    ! calculate ionization and recombination rates
    call kprad_ionization_rate(npts,ne,te,z,sion)
    call kprad_recombination_rate(npts,ne,te,z,srec)
    call kprad_energy_losses(npts,z,te, &
         ne,sion,srec,nz,nzeff,pion,preck,precp,prad,pbrem)
           
  end subroutine kprad_instantaneous_radiation


  subroutine kprad_advance_densities(dt, npts, z, p, ne, te, den, ti, nz, dw_rad, dw_brem,&
       dw_ion, dw_reck, dw_recp, source)

    use basic, only : pefac, ipres, itemp

    implicit none

    real, intent(in) :: dt                    ! time step to advance densities
    integer, intent(in) :: npts
    integer, intent(in) :: z
    real, intent(in) :: p(npts)              ! pressure in dyne/cm^2
    real, intent(inout) :: ne(npts)          ! electron density in cm^-3
    real, intent(in) :: te(npts)             ! electron temperature in eV
    real, intent(in) :: den(npts)            ! main ion density in cm^-3
    real, intent(in) :: ti(npts)             ! ion temperature in eV
    real, intent(inout) :: nz(npts,0:z)      ! density
    real, intent(out) :: dw_rad(npts,0:z)    ! energy lost via radiation
    real, intent(out) :: dw_brem(npts)       ! energy lost via bremsstrahlung
    real, intent(out) :: dw_ion(npts,0:z)    ! energy lost via ionization
    real, intent(out) :: dw_reck(npts,0:z)   ! kinetic energy lost via recombination
    real, intent(out) :: dw_recp(npts,0:z)   ! potential energy lost via recombination
    real, intent(in) :: source(npts,0:z)     ! optional density source
    
    real :: t, dts
    integer :: i
    real, dimension(npts) :: ne_old, nzt, delta, ti_over_te
    real, dimension(npts) :: te_int, p_int, dp_int   ! internal Te and p
    real, dimension(npts,0:z) :: nz_old
    real, dimension(npts,0:z-1) :: sion
    real, dimension(npts,0:z) :: srec
    real, dimension(npts) :: pbrem
    real, dimension(npts,z+1) :: imp_rad, pion, preck, precp
    real, dimension(npts, 2) :: nzeff
    real, dimension(npts,0:z) :: aimp, bimp, cimp, dimp, ework, fwork
    real :: max_change
    logical :: last_step
    real :: dts_min, dts_max

    dts_min = dt/1e6
    dts_max = 0.
    if(ikprad_max_dt.eq.1) then
       dts_max = dt/(z+1.0)  ! use one step per charge state
       if(kprad_max_dt.gt.0 .and. kprad_max_dt.lt.dts_max) &
            dts_max = kprad_max_dt
    else
       if(kprad_max_dt.gt.0) dts_max = kprad_max_dt
    end if

    dts = kprad_dt
    t = 0.
    dw_ion = 0.
    dw_reck = 0.
    dw_recp = 0.
    dw_rad = 0.
    dw_brem = 0.

    aimp(:,0) = 0.0
    cimp(:,z) = 0.0

    te_int = te

    if(ikprad_evolve_internal.eq.0) then
       call kprad_ionization_rate(npts, ne, te, z, sion)
       call kprad_recombination_rate(npts, ne, te, z, srec)
    else
       if(ipres.eq.0) then
          ! use total pressure
          p_int = p
          if(itemp.eq.1) then
             ti_over_te = ti/te
             if(any(abs(1. - ti_over_te*size(ti_over_te)/sum(ti_over_te)).gt.1e-3)) then
                print *, "Warning: Ti/Te not constant in kprad_advance_densities"
             end if
          end if
       else if(ipres.eq.1) then
          ! use electron pressure (in erg/cm^3)
          p_int = ne*te*1.6022e-12
       end if
    end if

    ! start time loop
    last_step = .false.
    do while(.not.last_step)
       if(t+dts.ge.dt) then
          if(kprad_min_dt.eq.kprad_min_dt_default .or. dts.lt.kprad_min_dt) &
               kprad_min_dt = dts
          dts = dt - t
          last_step = .true.
       end if

       nz_old = nz
       ne_old = ne

       ! calculate ionization and recombination rates
       if(ikprad_evolve_internal.eq.1) then
          call kprad_ionization_rate(npts, ne, te_int, z, sion)
          call kprad_recombination_rate(npts, ne, te_int, z, srec)
       end if

       do i=0, z
          if(i.gt.0) aimp(:,i) = -dts*sion(:,i-1)
          bimp(:,i) =  1. + dts*srec(:,i)
          if(i.lt.z) then
             bimp(:,i) = bimp(:,i) + dts*sion(:,i)
             cimp(:,i) = -dts*srec(:,i+1)
          end if
          dimp(:,i) = nz(:,i) + source(:,i)*dts
       enddo
     
       call tridiag(aimp,bimp,cimp,dimp,nz, &
            ework,fwork,npts,z)
       
       call kprad_energy_losses(npts,z,te_int, &
            ne,sion,srec,nz,nzeff,pion,preck,precp,imp_rad,pbrem)

       nzt = 0.
       do i=1, z
          ne = ne + i*(nz(:,i) - nz_old(:,i))
          nzt = nzt + nz(:,i)
       end do

       ! change in electron density
       where(ne_old.gt.0.)
          delta = abs(ne - ne_old)/ne_old
       elsewhere
          delta = 0.
       end where

       if(ikprad_evolve_internal.eq.1) then
          ! change in thermal energy
          dp_int = (imp_rad(:,z) + pion(:,z) + pbrem + preck(:,z))*dts * 1.e7
          where(p_int.gt.0) delta = sqrt(delta**2 + (dp_int/p_int)**2)
       end if

       max_change = maxval(delta)

       if((max_change .gt. 0.2) .and. (dts .gt. dts_min)) then

          ! If ne change is > 20%, backtrack changes and reduce time step
          ne = ne_old
          nz = nz_old
          dts = dts / 2.
          last_step = .false.

          if(dts.lt.dts_min) dts = dts_min

       else

          t = t + dts
          dw_ion  = dw_ion + pion*dts
          dw_reck = dw_reck + preck*dts
          dw_recp = dw_recp + precp*dts
          dw_brem = dw_brem + pbrem*dts
          dw_rad  = dw_rad + imp_rad*dts

          if(ikprad_evolve_internal.eq.1) then
             p_int = p_int - dp_int
             if(ipres.eq.0) then
                if(itemp.eq.0) then
                   te_int = pefac*p_int/(ne*1.6022e-12)
                else
                   te_int = p_int/((ne + (den+nzt)*ti_over_te)*1.6022e-12)
                end if
             else
                te_int = p_int/(ne*1.6022e-12)
             end if
          end if

          ! If ne change is < 2%, increase time step
          if(max_change .lt. 0.02) dts = dts * 1.5
          if((dts_max.gt.0.).and.(dts.gt.dts_max)) dts = dts_max

       end if

    enddo

  end subroutine kprad_advance_densities



!-----------------------------------------------------------------------
! kprad_ionization_rates gets the ionization rates for each charge state
!-----------------------------------------------------------------------
  subroutine KPRAD_IONIZATION_RATE(N,NE,TE,Z,sion)              

    !CALCULATE ionization rate for each charge state and both electr
    !       populations in s-1                                              

    implicit none 
                                                                        
    integer, intent(IN) :: N,Z
    real, dimension(N), intent(IN)::ne,te
    real, dimension(N,0:Z-1), intent(OUT) ::sion

    real, dimension(N) :: siont
    integer :: i
    
    do i=0,Z-1 
       call DPOLY_VAL(M2,N,sion_coeff(:,i+1),log10(TE),siont) 
       sion(:,i) = ne*10**siont
       if(ikprad_min_option.eq.2 .or. ikprad_min_option.eq.3) then
          where(ne.lt.kprad_nemin .or. te.lt.kprad_temin) sion(:,i) = 0.
       end if
    enddo
    where(sion.ne.sion) sion = 0.
  
  end subroutine KPRAD_IONIZATION_RATE

  !-----------------------------------------------------------------------
  ! kprad_recombination_rates gets the recomb. rates for each charge state
  !-----------------------------------------------------------------------
  subroutine KPRAD_RECOMBINATION_RATE(N,NE,TE,Z,srec)

    !CALCULATE recombination rate out of each charge state in s-1     

    implicit none 
                                                                            
    integer :: i
    integer, intent(in) :: N,Z
    real, intent(out) :: SREC(N,0:Z)
    real, intent(in) :: NE(N),TE(N)

    if(.not. allocated(ne_int)) allocate(ne_int(N))
    if(.not. allocated(te_int)) allocate(te_int(N))

    if (ikprad_min_option.eq.2) then
       ne_int = merge(ne, kprad_nemin, ne.ge.kprad_nemin)
       te_int = merge(te, kprad_temin, te.ge.kprad_temin)
    else
       ne_int = ne
       te_int = te
    end if

    SREC(:,0) = 0.0
    
    do i=1,Z
       SREC(:,i) = ne_int(:)*5.2E-14*ZED(i+1)*sqrt(Z_EI(i)/     &
            te_int(:))*(0.43+0.5*log(Z_EI(i)/te_int(:)) +           &
            0.469*(Z_EI(i)/te_int(:))**(-0.33))
       if(ikprad_min_option.eq.3) then
          where(ne.lt.kprad_nemin .or. te.lt.kprad_temin) srec(:,i) = 0.
       end if
    end do
    where(srec.ne.srec) srec = 0.

  end subroutine KPRAD_RECOMBINATION_RATE
                                                                        
  !*****************************************************      
  !-----------------------------------------------------------------------
  ! kprad_energy_losses gets the radiated power,ionization power, etc.     
  !-----------------------------------------------------------------------
  subroutine KPRAD_ENERGY_LOSSES(N,Z,TE,NE,SION,             &
       SREC,NZ,nZeff,pion,preck,precp,IMP_RAD,PBREM)

    implicit none 
                                                                        
    integer, intent(in) :: N,Z
    real, intent(out) :: IMP_RAD(N,Z+1),PION(N,Z+1),PRECK(N,Z+1),PRECP(N,Z+1)
    real, intent(out) :: PBREM(N)
    real, intent(in) :: TE(N), NE(N)
    real, intent(in) :: SION(N,0:Z-1),SREC(N,0:Z),NZ(N,0:Z)
    real, intent(out) :: nZeff(N,2)
    
    integer :: L
    real :: SNZ(N), impradt(N)

    SNZ=sum(NZ,DIM=2)
    nZeff(:,1)=0.0 
    nZeff(:,2)=1.0
    IMP_RAD=0.0
    !from all the electrons                          
    do L=1,Z 
       call DPOLY_VAL(M1,N,C(:,L),LOG10(TE*1.0e-3),impradt)
       
       IMP_RAD(:,L) = (10.0**impradt)*(NE/1.0E13)*NZ(:,L-1)
       if(ikprad_min_option.eq.2 .or. ikprad_min_option.eq.3) then
          where(ne.lt.kprad_nemin .or. te.lt.kprad_temin) IMP_RAD(:,L) = 0.
       end if
       
       PION(:,L)= SION(:,L-1)*NZ(:,L-1)*Z_EI(L)*1.6E-19 
       nZeff(:,1)=nZeff(:,1)+real(L)*NZ(:,L-1)/SNZ
       nZeff(:,2)=nZeff(:,2)+real((L**2-L))*NZ(:,L-1)/NE
       
       PRECK(:,L) = SREC(:,L)*NZ(:,L)*TE*1.6E-19 
       PRECP(:,L) = SREC(:,L)*NZ(:,L)*Z_EI(L)*1.6E-19
    end do
    
    PION(:,Z+1)  = sum(PION(:,1:Z),2) 
    !Total recombination losses                                      

    !totals -- recombination                                         
    PRECK(:,Z+1)=sum(PRECK(:,1:Z),2) 
    PRECP(:,Z+1)=sum(PRECP(:,1:Z),2) 

    !sum of all charge states
    IMP_RAD(:,Z+1)=sum(IMP_RAD(:,1:Z),DIM=2)

!       for radiation                                                   
    !CALCULATE instaneous power loss/gain due to ionization/recombinatio
    ! assume radiative recombination is dominant...so we lose electron  
!       energy                                                          
       !to radiation during recombination...this means we are neglecting
       !three-body recombination                                        
       !totals -- ionization                                            
                                                                        
       !CALCULATE average charge state of impurity and Zeff             
       ![ZED,NE]=meshgrid(zed,ne);                                      
                                                                        
       !ZZ(:,1)=sum(ZED*NZ,2)/sum(NZ,DIM=2)                             
    !ZZ(:,2)=1.0+ sum( (ZED**2-ZED)*( NZ/NE ),DIM=2 ) + 30.*0.03*NE0/NE 
                                                                        
       ! note...we have added in a guess at typical carbon density for  
!       initial Zeff ~ 2
       !CALCULATE radiative losses to bremsstrahlung
    ! This appears to be in units of W / cm^3, with ne in cm^-3 (-NF)
    PBREM = 1.69E-32*NE**2.0*SQRT(TE)*nZeff(:,2)
    if(ikprad_min_option.eq.2 .or. ikprad_min_option.eq.3) then
       where(ne.lt.kprad_nemin .or. te.lt.kprad_temin) PBREM = 0.
    end if

  end subroutine KPRAD_ENERGY_LOSSES
                                                                        
  subroutine kprad_atomic_data_sub(Z, ierr)
    implicit none
    
    integer, intent(in) :: Z
    integer, intent(out) :: ierr

    integer :: I
    
    ierr = 0
    select case (Z)
       
    case (1) !DUMMY ARRAYS FOR ZIMP=1
       allocate(C(10,1))
       allocate(SION_COEFF(7,1))
       allocate(Z_EI(1:Z+1))
       allocate(ZED(1:Z+1))
       
       C=transpose(reshape((/          &
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,        &
            0.0, 0.0, 0.0, 0.0/),(/1,10/)))



       SION_COEFF=transpose(reshape((/0.0, 0.0, 0.0, 0.0,         &
            0.0, 0.0, 0.0/),(/1,7/)))
        
       Z_EI = (/0.0, 0.0/)
       ZED=(/(real(I),I=0,Z)/)
       kprad_mz = 2.0
       write(*,*)  'NO IMPURITY SPECIES WITH ZIMP=1.'
       
       ! select IMPURITY SPECIES
    case (2) !SET HELIUM FOR IMPURITY
       allocate(C(10,2))
       allocate(SION_COEFF(7,2))
       allocate(Z_EI(1:Z+1))
       allocate(ZED(1:Z+1))
       m1 = 10
       m2 = 7
        
       
       C=transpose(reshape((/                                   &
            -12.126,   -12.4583,                                       &
            -0.238985,    -0.219838,                                  &
            -0.0808302,   -0.0924304,                            &
            0.0306689,    0.0360498,                            & 
            -0.00490822,  -0.00662168,                           &
            0.0046768,   0.00688757,                            &
            -0.00680186,  -0.00849505,                           &
            -7.11595e-06, -1.82597e-05,                          &
            0.000677444,   0.00078276,                          &
            0.000459187,  0.000537722/),(/2,10/)))



       SION_COEFF=transpose(reshape((/-17.6080,   -27.9344,      &
            17.2317,34.6191,                                           &
            -13.2220, -26.7690,                                        &
            5.74682,11.4064,                                           &
            -1.45270, -2.77375,                                        &
            0.196651, 0.359153,                                     &  
            -0.0109806,   -0.0191983/),(/2,7/)))
        
       Z_EI = (/24.5876,54.416,1.0E6/)
       ZED=(/(real(I),I=0,Z)/)
       kprad_mz = 4.0
       
    case (4) !SET BERYLLIUM FOR IMPURITY
        
       allocate(C(8,4))
       allocate(SION_COEFF(8,4))
       allocate(Z_EI(1:Z+1))
       allocate(ZED(1:Z+1)) 
       m1 = 8
       m2 = 8
        
       C=RESHAPE((/-11.59812535,-0.268122674,-0.056049099,0.016721906,       &
            -0.00383657,-0.003596497,0.00093653,0.001040357,                 &
            -11.75867978,-0.268313536,-0.059673665,0.010949875,0.00534501,   &
            -0.002664926,-0.001076657,0.000641886,                           &
            -12.35611424,-0.081873364,-0.194357749,0.149787442,-0.066529555, &
            -0.002487972,-0.015711093,0.014728804,                           &
            -12.67632568,-0.059377317,-0.215335252,0.145547429,              &
            -0.062877605,   0.012962914,    -0.020546062,                    &
            0.011978393/),(/8,4/))                                                                                  
        
       SION_COEFF=RESHAPE((/-11.76875877,9.848173141,              &
            -9.606463432,5.737738132,                                        &
            -2.171055079,0.497681499,-0.062647402,0.003309997,               &
            -15.61272144,15.64862251,-13.56947136,6.444764614,-1.734611154,  &
            0.246102199,-0.01427075,0.0,-54.81453705,76.31816864,            &
            -54.81626129,   21.62486267,-4.881194115,                        &
            0.591508269,    -0.029873915,   0.0,                             &
            -69.67385101,99.21066284,-70.1591568,27.06711578,-5.950128555,   &
            0.701139331,-0.034435261,0.0/),(/8,4/))
       Z_EI= (/9.3227,   18.21114,   153.89661, 217.71865,1.0E6/)
       ZED=(/(real(I),I=0,Z)/)
       kprad_mz = 9.012
       
    case (5) !SET BORON FOR IMPURITY

       allocate(C(9,5))
       allocate(SION_COEFF(9,5))
       allocate(Z_EI(1:Z+1))
       allocate(ZED(1:Z+1))
       m1 = 9
       m2 = 9
       C=RESHAPE((/ -1.192132e+01, -4.967397e-01,  5.043827e-05,  8.960577e-03, &
                    -8.693015e-03, -7.489224e-04,  1.166458e-03,  5.175516e-04, &
                    -1.962665e-04, -1.208493e+01, -4.815787e-01, -1.395910e-02, &
                    -1.374892e-02, -1.608155e-03,  1.194821e-02, -1.608300e-04, &
                    -1.128394e-03, -3.192678e-04, -1.259333e+01, -4.941146e-01, &
                    -7.179767e-04,  1.163202e-02, -1.347079e-02,  1.174152e-04, &
                     2.328052e-03,  3.779181e-04, -3.167863e-04, -1.252475e+01, &
                    -2.948459e-01, -2.393933e-01,  1.938396e-01, -1.055447e-01, &
                     3.702422e-02, -1.731613e-02,  1.032923e-02, -2.446645e-03, &
                    -1.285255e+01, -2.366438e-01, -3.069907e-01,  2.439879e-01, &
                    -1.318474e-01,  4.906848e-02, -2.254865e-02,  1.183835e-02, &
                    -2.623528e-03/), (/9,5/))
       SION_COEFF=RESHAPE((/ -1.152575e+01,  8.755813e+00, -8.383173e+00,  5.520414e+00, &
                             -2.551776e+00,  7.772134e-01, -1.459613e-01,  1.520845e-02, &
                             -6.704090e-04, -1.900611e+01,  2.421069e+01, -2.471980e+01, &
                              1.557728e+01, -6.381716e+00,  1.690138e+00, -2.782783e-01, &
                              2.584187e-02, -1.032499e-03, -2.444181e+01,  3.411982e+01, &
                             -3.320268e+01,  1.931898e+01, -7.189498e+00,  1.725666e+00, &
                             -2.589333e-01,  2.211362e-02, -8.208679e-04, -1.013409e+02, &
                              1.776986e+02, -1.564890e+02,  8.132509e+01, -2.698106e+01, &
                              5.805324e+00, -7.868877e-01,  6.119722e-02, -2.084705e-03, &
                             -1.266743e+02,  2.228490e+02, -1.934931e+02,  9.902706e+01, &
                             -3.235052e+01,  6.857843e+00, -9.166120e-01,  7.036009e-02, &
                              -2.367897e-03/), (/9,5/))

       Z_EI= (/8.2980, 25.1548, 37.93064, 259.37521, 340.22580,1.0E6/)

       ZED=(/(real(I),I=0,Z)/)
       kprad_mz = 10.81

    case (6) !CARBON
       
       allocate(C(8,6))
       allocate(SION_COEFF(10,6))
       allocate(Z_EI(1:Z+1))
       allocate(ZED(1:Z+1))
       m1 = 8
       m2 = 10
       
       C=RESHAPE((/-11.658639,-0.24560987,-0.044048871,0.01451389,           &
            -0.022141926,-0.005714527,0.005500657,0.003130559,               &
            -11.649103,     -0.24190558,    -0.037042534,   0.010883051,     &
            -0.031552067,   -0.003814105,   0.007096615,    0.003875159,     &
            -11.70192,      -0.2377956,     -0.022745361,   0.010836725,     &
            -0.048891777,   -0.006297751,   0.012263192,       0.005220342,  &
            -11.913389,     -0.2128092,     -0.045676775,   0.001569224,     &
            -0.020538656,   -0.004588293,   0.006418909,       0.003218322,  &
            -12.550067,     0.1304338,      -0.40753917,        0.30143108,  &
            -0.14696592,    0.054084502,-0.046054565,0.020007227,            &
            -12.868148,0.15009505,  -0.43177885,        0.32590553,          &
            -0.15323027, 0.049376739,-0.047432245,  0.022250319/),(/8,6/))

       SION_COEFF=transpose(RESHAPE((/-12.7340,-17.4107,-25.1086,-31.1445,             &
            -73.8292,-86.1572,                                               &
            10.4813,   18.3657 ,29.8344,39.4998, 83.4077 ,                   &
            99.2425,                                                         &
            -8.38791,-14.6705 ,-22.2500,-29.8922,  -43.7563                  &
            ,-52.2210,                                                       &
            3.75569,6.38180,8.98325  ,12.2766 ,11.5784                       &
            ,13.8259,                                                        &
            -0.960082, -1.56452 ,-2.04632,-2.83479, -1.53664,                &
            -1.83214,                                                        &
            0.128895,0.201356 ,0.246368,0.344735,0.0813417                   &
            ,0.0967250,                                                      &
            -0.00702319, -0.0105572, -0.0121805, -0.0171659,                 &
            0.00000,0.00000,                                                 &
            0.00000,  0.00000,0.00000, 0.00000 ,0.00000                      &
            ,0.00000,                                                        &
            0.00000,  0.00000 ,0.00000,  0.00000,0.00000                     &
            ,0.00000,                                                        &
            0.00000, 0.00000 ,0.00000, 0.00000,0.00000                       &
            ,0.00000/),                                                      &
            (/6,10/)))

       Z_EI= (/11.26, 24.384, 47.888, 64.5, 392.1, 490.0, 1.0E6/)
       ZED=(/(real(I),I=0,Z)/)
       kprad_mz = 12.0
       
    case (10)
       
       allocate(C(10,10))
       allocate(SION_COEFF(10,10))
       allocate(Z_EI(1:Z+1))
       allocate(ZED(1:Z+1))
       m1 = 10
       m2 = 10
       
       C=reshape((/ -11.707716, -0.18676494,  -0.13528950,                   &
            0.025605723,                                                     &
            0.088009309,    0.0084664872,                                    &
            -0.069115677,  -0.0054429311,  0.013196953, 0.0031395287,        &
            -11.707716, -0.18676494 ,-0.13528950,0.025605723,                &
            0.088009309,  0.0084664872,                                      &
            -0.069115677,  -0.0054429311,0.013196953,0.0031395287,           &
            -11.723273, -0.15958272  ,-0.050532946, -0.033958697,            &
            -0.088516297, 0.073820408,                                       &
            0.036099811,  -0.013526015, -0.0085533573, -0.00070188989,       &
            -11.748191,  -0.15814151    ,-0.054800120, 0.010859518,          &
            -0.10618337,  0.037611955,                                       &
            0.049460301, -0.0046416701  ,-0.010120509, -0.0014656924,        &
            -11.767116 ,-0.15555648 ,-0.074620798,0.028127561,               &
            -0.074889102,   0.018080982,                                     &
            0.033130626,-0.00047289192,-0.0070322880, -0.0010321265,         &
            -11.830380,-0.14634359,-0.092327940,0.037655326,                 &
            -0.042212739,  0.0035730051,                                     &
            0.015450174,0.0018603897,-0.0034545906,-0.00041907844,           &
            -11.918433, -0.14148970 ,-0.094057644  ,0.024527836 ,            &
            -0.0092384784,-0.0021089001,                                     &
            0.00029707864, 0.0021866082, -0.00051323880,                     &
            8.7974239e-05,                                                   &
            -12.174002, -0.088504778, -0.082410590, -0.028255051,            &
            0.019315606, 0.0073984201,                                       &
            -0.0084316763,-0.00030755515, 0.00097443423,                     &
            0.00036464413,                                                   &
            -12.961943  ,0.76030676,-1.0994758,0.85890118,                   &
            -0.52044842, 0.22250819,                                         &
            -0.063636605,0.023730279,-0.015142696, 0.0042680263,             &
            -13.281673 ,0.79175595, -1.1343042, 0.88148321 ,                 &
             -0.52384009,0.21989544,                                         &
            -0.067363452,0.027969167,-0.015372851,                           &
            0.0038092913/),(/10,10/))

       sion_coeff=transpose(reshape((/-15.5946,   -23.3026,     -28.0792,    &
            -35.5425,                                                        &
            -42.8937,     -52.0208,     -49.8740,                            &
            -58.4576,     -81.0648 ,-90.7343,                                &
            14.2879,    26.3362  ,34.5052,                                   &
            44.8789,      55.5430,      68.4173,                             &
            56.6446,                                                         &
            66.4594,    68.1825  ,77.7918,                                   &
            -11.2628,   -19.1807 ,-25.5282,                                  &
            -31.4065,     -38.2860,     -46.1427,                            &
            -32.1844,                                                        &
            -36.1988, -24.7538 ,-28.4056,                                    &
            5.06897,  7.70118  ,10.3490,                                     &
            11.9438,  14.3466,      16.9024,    9.46790,                     &
            9.88111,  4.01827  ,4.62336,                                     &
            -1.30957,   -1.76208 ,-2.37449,                                  &
            -2.56998,   -3.04571,     -3.51118,     -1.48702,                &
            -1.34501,  -0.245672  ,-0.282649,                                &
            0.177775,   0.213985 ,0.288012,                                  &
            0.293157,     0.343692,   0.388672,     0.113323,                &
            0.0726658,    0.00000,  0.00000,                                 &
            -0.00977938,   -0.0106844,   -0.0143345,                         &
           -0.0137726,   -0.0160197, -0.0178211,  -0.00296459,               &
            0.00000  ,0.00000,    0.00000,                                   &
            0.00000  ,0.00000,    0.00000,     0.00000   ,                   &
            0.00000,      0.00000,      0.00000,                             &
            0.00000  ,0.00000,    0.00000,                                   &
            0.00000  ,0.00000,    0.00000,      0.00000 ,                    &
            0.00000,      0.00000,      0.00000,                             &
            0.00000  ,0.00000,    0.00000,                                   &
            0.00000  ,0.00000,    0.00000,      0.00000,                     &
            0.00000,    0.00000,    0.00000,                                 &
            0.00000  ,0.00000,    0.00000/),(/10,10/)))
       
       z_ei =                          &
            (/21.6,41.0,63.5,97.0,126.3,157.9,207.2 ,239.0,1195.0,1362.3,1.0e6/)
       ZED=(/(real(I),I=0,Z)/)
       kprad_mz = 20.0
               
       
    case (18)
       
       allocate(C(8,18))
       allocate(SION_COEFF(10,18))
       allocate(Z_EI(1:Z+1))
       allocate(ZED(1:Z+1))
       m1 = 8
       m2 = 10
       
       C=reshape((/                                                          &
            -11.115350,-0.22648251, -0.072592330, 0.0052344029,              &
            0.0096323252,  0.0031468124, -0.0034768159,                      &
            -0.00026370319,                                                  &
            -11.202989   ,-0.19262990,-0.064467383, 0.015013939,             &
            -0.022544940 ,-0.0040357671,0.0046553203,0.0037237526,           &
            -11.234236   ,-0.19759543,  -0.033967713   ,0.013103205,         &
            -0.052866921 ,-0.0077906587,   0.013168522,0.0062897673,         &
            -11.301820   ,-0.17872160,  -0.032038500   ,0.020457332,         &
            -0.068808417,-0.012267377,   0.016797858,0.0086463004,           &
            -11.364419   ,-0.16439383,  -0.047202305   ,0.023570925,         &
            -0.057902687,-0.013666551,   0.014636065,0.0080280030,           &
            -11.452980   ,-0.13082109,  -0.065423357   ,0.027492870,         &
            -0.055904037,-0.015748508,   0.014373749,0.0087232750,           &
            -11.478044   ,-0.12660441,  -0.075988354   ,0.016616691,         &
            -0.035326777,-0.012410099,  0.0097844802,0.0070637667,           &
            -11.644427  ,-0.042213033,   -0.15096961    ,0.0057700113,       &
            0.015007880 ,-0.0082830257,-0.00046463904,0.0039770109,          &
            -11.912902,0.11835730,   -0.38404234  ,0.32084511,               &
            -0.14930369, 0.020998977,-0.035623931, 0.023762242,              &
            -11.924930,0.17001001,   -0.56077370  ,0.24104250,               &
            0.16047043,-0.069081786,-0.044303359,0.0086861172,               &
            -11.971962,0.19447692,   -0.50929221  ,0.10429281,               &
            0.19770146,-0.017291468,-0.054363083,-0.00010139356,             &
            -12.033731,0.19331504,   -0.42677735   ,0.032466075,             &
            0.17391070,0.0022071038,-0.047505560, -0.0013848471,             &
            -12.065764,0.12104935,   -0.26799862  ,0.0021716118,             &
            0.078794484,0.0070322263,-0.026860932,  0.0032693364,            &
            -12.175932,0.11579738,   -0.21762619 ,-0.0049662849,             &
            0.047193618,0.0071477932,-0.020270099,  0.0049998252,            &
            -12.271151 ,0.046395265,   -0.10401582,-0.016001716,             &
            -0.0050251942,0.0050219330 ,-0.0074081161,  0.0064065231,        &
            -12.578830 ,0.060816027,   0.049031648  ,-0.059414309,           &
            -0.077380778 ,0.00034639424, 0.014759651,  0.0087798639,         &
            -14.085028, 2.9982267,    -3.6844645, 2.7995513,                 &
            -1.6019418,  0.77264625, -0.29292775,   0.056285854,             &
            -14.411208,3.0439809,-3.7234131,2.8177722,                       &
            -1.6175255, 0.79514635, -0.30843034 ,0.060052147/),              &
            (/8,18/))                                              

       sion_coeff=transpose(reshape((/                                       &
            -14.7011,   -17.7744,   -22.1440,   -26.737,                     &
            -30.2985,     -28.4935,     -32.9706,                            &
            -40.1073, -70.4135,   -73.5380, -82.5107,                        &
            -91.5811,   -95.7802,   -109.215,                                &
            -126.953, -92.6558,   -166.579, -165.765,                        &
            14.7364,  19.7206,    26.3383,  33.1043,                         &
            36.8346,27.9872,  30.7708,                                       &
            43.8658,  78.1980,    80.6788,  90.9403,                         &
            103.262,  107.459,    124.130,                                   &
            148.970,  92.4737,    141.628,  139.648,                         &
            -12.0714, -15.7357,   -20.1592, -24.7612,                        &
            -25.9344 ,-15.1888, -14.9428,                                    &
            -24.5349, -40.1045,   -41.0685, -45.7779,                        &
            -52.6188 ,-54.3072, -62.8024,                                    &
            -76.7611, -41.3994,   -48.8054, -47.7125,                        &
            5.45205,  6.81562,    8.32276,  10.0316,                         &
            9.99728,4.09221,  3.54245,                                       &
            6.85033,  10.3529,    10.5907,  11.6410,                         &
            13.5382,13.8601,  16.0478,                                       &
            19.9138,  9.21174,    7.53545,  7.28987,                         &
            -1.39128, -1.66178,   -1.93193, -2.29483,                        &
            -2.21018,  -0.550309,-0.413532,                                  &
            -0.950696, -1.33956,   -1.37514, -1.48795,                       &
            -1.74855,   -1.77653, -2.06017,                                  &
            -2.58857, -1.01471,  -0.438891,-0.419560,                        &
            0.185806, 0.212769,   0.236282, 0.277732,                        &
            0.262330,  0.0294354,0.0190083,                                  &
            0.0521924,0.0691496,  0.0714507,0.0760625,                       &
            0.0901641,  0.0909846, 0.105755,                                 &
            0.134293,0.0439989,    0.00000,  0.00000,                        &
            -0.0100628 ,-0.0111056, -0.0118376 ,-0.0138125,                  &
            -0.0129407,    0.00000,    0.00000,                              &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,  0.00000,  0.00000,                                     &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,    0.00000,    0.00000,                                 &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,  0.00000,  0.00000,                                     &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,0.00000,  0.00000,                                       &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,0.00000,  0.00000,                                       &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,0.00000,    0.00000,                                     &
            0.00000,  0.00000,    0.00000,  0.00000,                         &
            0.00000,0.00000,    0.00000,                                     &
            0.00000,  0.00000,    0.00000,  0.00000/),(/18,10/)))
        
       z_ei = (/15.76,27.63,40.74,59.81,75.02,91.0,124.324,143.5,  &
            422.5,478.7,618.3,538.96,686.11,755.75,854.78,918.05,   &  
            4120.87,4426.24,1.0e6/)
       ZED=(/(real(I),I=0,Z)/)
       kprad_mz = 40.0
         
    case DEFAULT
       write(*,*) 'NO DATA FOR THIS ELEMENT EXISTS!'
       ierr = 1
       
    end select

  end subroutine kprad_atomic_data_sub


  !-----------------------------------------------------------------------
  ! dpoly_val evaluates polynomials                                       
  !-----------------------------------------------------------------------
  subroutine DPOLY_VAL( N, M,COEFFS, X, Y )
  
    !***********************************************************************
    !                                                                       
    !! DPOLY_VAL evaluates a polynomial using Horner's method.              
    !                                                                       
    !  Parameters:                                                          
    !                                                                       
    !    Input, integer N, the degree of the polynomial.                    
    !                                                                       
    !    Input, integer M, the size of X                                    
    !                                                                       
    !    Input, real C(0:N), the polynomial coefficients.                   
    !    C(I) is the coefficient of X**I.                                   
    !                                                                       
    !    Input, real X(0:M), the points at which the polynomial             
    !    is to be evaluated.                                                
    !                                                                       
    !    Output, real CX(0:M), the value of the polynomial at X(0:M).       

                                                                        
    implicit none 

    integer :: I,N,M
    
    real :: COEFFS(N)
    real :: Y(1:M)
    real :: X(1:M)
    
    Y = COEFFS(1) 
    do I = 2,N 
       Y=Y+COEFFS(I)*(X**real((I-1)))
    enddo
    
  end subroutine dpoly_val


  subroutine tridiag(a,b,c,d,x,e,f,n,m)
    ! solves n simultaneous tri-diagonal matrices of rank m

    implicit none

    integer, intent(in) :: m,n
    real, intent(in), dimension(n,0:m) :: a,b,c,d
    real, intent(out), dimension(n,0:m) :: x,e,f

    integer :: i, j
    
    e(:,0) = b(:,0)
    f(:,0) = d(:,0)/e(:,0)
    do i = 1, m
       e(:,i) = b(:,i)-a(:,i)*c(:,i-1)/e(:,i-1)
       f(:,i) = (d(:,i)-a(:,i)*f(:,i-1))/e(:,i)
    end do
    
    x(:,m) = f(:,m)
    do i = 2, m+1
       j = m-i+1
       x(:,j) = f(:,j)-c(:,j)*x(:,j+1)/e(:,j)
    end do
  end subroutine tridiag

end module kprad
