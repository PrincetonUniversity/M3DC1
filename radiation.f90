module radiation
  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: dp = selected_real_kind(15,307)
contains

  ! get_Prad_simple calls get_Prad
  ! assuming one impurity species
  ! having same temperature as electrons
  subroutine get_Prad_simple(prad,Te,nz,Z,ne,ierr)
    implicit none
    vectype, intent(out) :: prad
    vectype, intent(in) :: Te, nz, ne
    integer, intent(in) :: Z
    integer, intent(out) :: ierr

    real(dp), allocatable :: pradv(:), zeffv(:)
    real(dp), dimension(1) :: tev, nzv
    integer, dimension(1) :: Zv
    real(dp) :: nev

    tev(1) = te
    nzv(1) = nz
    Zv(1) = Z
    nev = ne

    if(tev(1).le.0.03) tev(1) = 0.03

!    write(*,'(A,2E12.4,I5,E12.4)') 'calling get_prad', tev, nzv, zv, nev
    call get_Prad(pradv, zeffv, tev, nzv, zv, nev,ierr)

    prad = pradv(1)
    deallocate(pradv, zeffv)
  end subroutine get_Prad_simple

    
  subroutine get_Prad(Prad,Zeff,Te,nZ,Z,ne,ierr)
    ! Inputs: Te: Temperature in keV, vector
    !         nZ,ne:  Density of Impurities/Electron in m^-3, vector/scalar
    !         Z: Atomic number of Impurities, vector
    ! Outputs: Prad: Radiation power in W*m^-3, vector
    !          Zeff: Effective chage state, vector
    !Example: call get_Prad(Prad,Zeff,(/1._dp,2._dp,3._dp/),(/1.E20_dp,1.E20_dp/),(/6,26/),1.E20_dp)
      implicit none
      real(dp), intent(out),allocatable :: Prad(:),Zeff(:)
      real(dp),intent(in) :: Te(:),nZ(:),ne
      integer :: N,i
      integer,intent(in) :: Z(:)
      integer, intent(out) :: ierr
      real(dp),allocatable :: Lz(:), Zave(:), Zsqave(:), LD(:)
      ierr = 0
      N=size(Te)
      allocate(Prad(N),Zeff(N),LD(N))
      Zeff=(/(1,i=1,N)/)
      Prad=(/(1,i=1,N)/)     
      LD=5.35E-37*sqrt(Te)
      do i=1,size(Z)
         call get_Lz(Lz,Zave,Zsqave,Te,Z(i),ierr)
         if(ierr.ne.0) goto 100
         Zeff=Zeff+nZ(i)/ne*(Zave**2-Zave)
         Prad=Prad-nZ(i)/ne*Zave+nZ(i)/ne*Lz/LD
      end do
      Prad=Prad*ne**2*LD
100   deallocate(Lz,Zave,Zsqave,LD)
    end subroutine get_Prad
  subroutine get_Lz(Lz,Zave,Zsqave,Te,Z, ierr)
    ! Inputs: Te: Temperature in keV, vector
    !         Z:  The atomic number of the element, scalar
    ! Outputs: Lz: the cooling rate in 1e13erg*cm^3/sec(W*m^3), vector
    !          Zave: average ionization <Z>, vector
    !          Zsqave: <Z^2>, vector
    !Example: call Lz(Lz,Zave,Zsqave,(/1,2,3/),6)
    real(dp),intent(out),allocatable :: Lz(:), Zave(:), Zsqave(:)
    real(dp),intent(in) :: Te(:)
    integer, intent(in) :: Z
    integer, intent(out) :: ierr
    integer :: N,i,j
    real(dp),allocatable :: T(:,:), A(:,:), B(:,:), C(:,:)

    ierr = 0
    N=size(Te)
    allocate(Lz(N),Zave(N),Zsqave(N))
    select case (Z)
    case (6)
       !C
       allocate(T(5,2),A(5,6),B(5,6),C(5,6))
       T=reshape((/3.e-3,2.e-2,&
            2e-2,2e-1,&
            2e-1,2e0,&
            2e0,2e1,&
            2e1,1e2/),(/5,2/),(/0./),(/2,1/))
       A=reshape((/1.965300E03, 4.572039E03, 4.159590E03, 1.871560E03, 4.173889E02, 3.699382E01,&
            7.467599E01, 4.549038E02, 8.372937E02, 7.402515E02, 3.147607E02, 5.164578E01,&
            -2.120151E01, -3.668933E-01, 7.295099E-01, -1.944827E-01, -1.263576E-01, -1.491027E-01,&
            -2.121979E01, -2.346986E-01, 4.093794E-01, 7.874548E-02, -1.841379E-01, 5.590744E-02,&
            -2.476796E01, 9.408181E00, -9.657446E00, 4.999161E00, -1.237382E00, 1.160610E-01/),&
            (/5,6/),(/0./),(/2,1/))
       B=reshape((/2.093736E03, 5.153766E03, 5.042105E03, 2.445345E03, 5.879207E02, 5.609128E01,&
            7.286051E00, 4.506650E01, 1.595483E02, 2.112702E02, 1.186473E02, 2.400040E01,&
            5.998782E00, 6.808206E-03, -3.096242E-02, -4.794194E-02, 1.374732E-01, 5.157571E-01,&
            5.988153E00, 7.335329E-02, -1.754858E-01, 2.034126E-01, -1.145930E-01, 2.517150E-02,&
            -3.412166E01, 1.250454E02, -1.550822E02, 9.568297E01, -2.937297E01, 3.589667E00/),&
            (/5,6/),(/0./),(/2,1/))
       C=reshape((/1.569604E04, 3.872013E04, 3.793754E04, 1.842954E04, 4.438563E03, 4.241526E02,&
            -8.567197E01, -1.867449E02, 4.249500E02, 1.072139E03, 7.377224E02, 1.646775E02,&
            3.598505E01, 7.432950E-02, -2.552972E-01, -4.325508E-01, 8.446683E-01, 4.719095E00,&
            3.588152E01, 7.335339E-01, -1.754860E00, 2.034127E00, -1.145931E00, 2.517152E-01,&
            -2.047300E02, 7.502727E02, -9.304934E02, 5.740977E02, -1.762379E02, 2.153801E01/),&
            (/5,6/),(/0./),(/2,1/))
    case (26)
       !Fe
       allocate(T(4,2),A(4,6),B(4,6),C(4,6))
       T=reshape((/2e-2,2e-1,&
            2e-1,2e0,&
            2e0,2e1,&
            2e1,1e2/),(/4,2/),(/0./),(/2,1/))
       A=reshape((/-2.752599e1,-3.908228e1,-6.469423e1,-5.555048e1,-2.405568e1,-4.093160e0,&
            -1.834973e1,-1.252028e0,-7.533115e0,-3.289693e0,2.866739e1,2.830249e1,&
            -1.671042e1,-1.646143e1,3.766238e1,-3.944080e1,1.918529e1,-3.509238e0,&
            -2.453957e1,1.795222e1,-2.356360e1,1.484503e1,-4.542323e0,5.477462e-1/),&
            (/4,6/),(/0./),(/2,1/))
       B=reshape((/8.778318e0,-5.581412e1,-1.225124e2,-1.013985e2,-3.914244e1,-5.851445e0,&
            1.959496e1,1.997852e1,9.593373e0,-7.609962e1,-1.007190e2,-1.144046e1,&
            1.478150e1,6.945311e1,-1.991202e2,2.653804e2,-1.631423e2,3.770026e1,&
            2.122081E1,6.891607E0,-4.076853E0,1.577171E0,-5.139468E-1,8.934176E-2/),&
            (/4,6/),(/0./),(/2,1/))
       C=reshape((/-1.629057E00, -2.038262E03, -4.720863E03, -4.299215E03, -1.797781E03, -2.879963E02,&
            3.885212E02, 7.909700E02, 6.077126E02, -2.903589E03, -5.163387E03, -1.865747E03,&
            1.181937E02, 3.464816E03, -9.953434E03, 1.328083E04, -8.173660E03, 1.891447E03,&
            -3.243616E02, 2.476134E03, -2.551817E03, 1.348646E03, -3.619529E02, 3.921327E01/),&
            (/4,6/),(/0./),(/2,1/))
    case (18)
       !Argon
       allocate(T(4,2),A(4,6),B(4,6),C(4,6))
       T=reshape((/3E-2, 2E-1,&
            2E-1, 2E0,&
            2E0, 2E1,&
            2E1, 1E2/),(/4,2/),(/0./),(/2,1/))
       A=reshape((/-2.053043E01, -2.834287E00, 1.506902E01, 3.517177E01, 2.400122E01, 5.072723E00,&
            -1.965204E01, -1.172763E-01, 7.833220E00, -6.351577E00, -3.058849E01, -1.528534E01,&
            -1.974883E01, 2.964839E00, -8.829391E00, 9.791004E00, -4.960018E00, 9.820032E-01,&
            -2.117935E01, 5.191481E00, -7.439717E00, 4.969023E00, -1.553180E00, 1.877047E-01/),&
            (/4,6/),(/0./),(/2,1/))
       B=reshape((/-6.351877E01, -4.145188E02, -8.502200E02, -8.074783E02, -3.621262E02, -6.187830E01,&
            1.591067E01, -7.886774E-01, 2.874539E00, 3.361189E01, -3.306891E01, -7.162601E01,&
            1.296383E01, 1.833252E01, -2.834795E01, 2.267165E01, -9.219744E00, 1.507133E00,&
            -7.890005E01, 2.939922E02, -3.551084E02, 2.133685E02, -6.376566E01, 7.582572E00/),&
            (/4,6/),(/0./),(/2,1/))
       C=reshape((/-9.617104E02, -6.478980E03, -1.394668E04, -1.372414E04, -6.350906E03, -1.120351E03,&
            2.534142E02, -2.346306E01, 7.075633E01, 1.014115E03, -7.342520E02, -1.963298E03,&
            1.582535E02, 5.763311E02, -8.295298E02, 5.954534E02, -2.057947E02, 2.602888E01,&
            2.959939E02, 6.590942E01, -7.358733E01, 4.647609E01, -1.567553E01, 2.167101E00/),&
            (/4,6/),(/0./),(/2,1/))

    case default
       !  print *, 'Z = ', Z, ' not implemented in Prad module.'
       !  print *, "Z's implemented: 6, 18, 26"
       ierr = 1
       return

    end select
    do i=1,N
       if (Te(i)<T(1,1) .or. Te(i)>T(size(T,1),2)) then
          Lz(i)=-1.
          Zave(i)=-1.
          Zsqave(i)=-1.
          return
       end if
       j=1
       do while (.not.(Te(i)>=T(j,1).and.Te(i)<=T(j,2)))
          j=j+1
       end do
       Lz(i)=1.E-13*10**(A(j,1)+A(j,2)*log10(Te(i))+A(j,3)*log10(Te(i))**2+A(j,4)*log10(Te(i))**3+&
            A(j,5)*log10(Te(i))**4+A(j,6)*log10(Te(i))**5)
       Zave(i)=(B(j,1)+B(j,2)*log10(Te(i))+B(j,3)*log10(Te(i))**2+B(j,4)*log10(Te(i))**3+&
            B(j,5)*log10(Te(i))**4+B(j,6)*log10(Te(i))**5)
       Zsqave(i)=(C(j,1)+C(j,2)*log10(Te(i))+C(j,3)*log10(Te(i))**2+C(j,4)*log10(Te(i))**3+&
            C(j,5)*log10(Te(i))**4+C(j,6)*log10(Te(i))**5)
    end do
  end subroutine get_Lz

end module radiation
