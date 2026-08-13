module physconst
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: e = 1.602176634d-19
  real(dp), parameter :: epsilon0 = 8.8541878128d-12
  real(dp), parameter :: me = 9.1093837015d-31
  real(dp), parameter :: c = 2.99792458d8
  real(dp), parameter :: pi = 3.14159265358979323846d0
end module physconst

module HT_funcs
    use physconst
    implicit none
contains

    function getf(v, n0, vT0, tau) result(f)
        real(dp), intent(in) :: v, n0, vT0, tau
        real(dp) :: f
        real(dp) :: prefactor, expo

        prefactor = n0 / (pi**1.5d0 * vT0**3)
        expo = -((v**3 / vT0**3) + 3.d0*tau)**(2.d0/3.d0)
        f = prefactor * exp(expo)
    end function getf 

    function getfm(v, n0, vt) result(fm)
        real(dp), intent(in) :: v, n0, vt
        real(dp) :: fm
        real(dp) :: prefactor, expo

        prefactor = n0 / (pi**1.5d0 * vt**3)
        expo = -(v**2 / vt**2)
        fm = prefactor * exp(expo)
    end function getfm

    ! Temp evo used to recreate figure 
    function temperature(t, T0, Tf, ti) result(Tout)
        real(dp), intent(in) :: t, T0, Tf, ti
        real(dp) :: Tout

        Tout = Tf + (T0 - Tf) * exp(-t / ti)
    end function temperature

    function getvT0(T0) result(vT0)
        use physconst
        implicit none
        real(dp), intent(in) :: T0      
        real(dp) :: vT0                 

        vT0 = sqrt(2.0_dp * T0 * e / me)
    end function getvT0


    function getnu0(n, T0) result(nu0)
        real(dp), intent(in) :: n, T0
        real(dp) :: nu0
        real(dp) :: vT0, lnL

        vT0 = sqrt(2.d0 * T0 * e / me)
        lnL = 15.d0

        nu0 = n * e**4 * lnL / &
            (4.d0 * pi * epsilon0**2 * me**2 * vT0**3)
    end function getnu0

    function getTau(nu0, t, t0) result(tau)
        real(dp), intent(in) :: nu0, t, t0
        real(dp) :: tau

        if (t > 0.d0) then
            tau = nu0 * t
        else
            tau = 0.d0
        end if
    end function getTau

    ! Needed the Dreicer FIeld function so added here

    function getDreicerField(ne, lnA, T_eV) result(Ed)
        real(dp), intent(in) :: ne, lnA, T_eV       
        real(dp) :: Ed                  ! Dreicer field [V/m]
        real(dp) :: T_J

        T_J = T_eV * e


        Ed = e**3 * ne * lnA / (T_J * 4.0_dp * pi * epsilon0**2)

    end function getDreicerField

    function get_vc(T_eV, Efield, ne) result(vc)
        real(dp), intent(in) :: T_eV, Efield, ne    
        real(dp) :: vc
        real(dp) :: vc2, Ed
        !real(dp) :: e_charge

        Ed = getDreicerField(ne, 15.0_dp, T_eV) ! columb log is assumed 15
        !print *, "D field", Ed

        !e_charge = 1.602176634d-19

        vc2 = ((T_eV * e) / me) * (Ed / Efield)

        !print *, "first half:", (T_eV * e / me)
        !print *, "2nd half:", (Ed / Efield)
        !print *, "temp times e", (T_eV * e)
        !print *, "e= ", e

        if (vc2 > 0.d0) then
            vc = sqrt(vc2)
        else
            vc = 0.d0
        end if
    end function get_vc

    ! https://fortranwiki.org/fortran/show/integration
    pure function integrate(x, y) result(r)
        !! Calculates the integral of an array y with respect to x using the trapezoid
        !! approximation. Note that the mesh spacing of x does not have to be uniform.
        real(dp), intent(in)  :: x(:)         !! Variable x
        real(dp), intent(in)  :: y(size(x))   !! Function y(x)
        real(dp)              :: r            !! Integral ∫y(x)·dx

        ! Integrate using the trapezoidal rule
        associate(n => size(x))
        r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
        end associate
    end function



    ! Tried using aforementioned trapazoid integration method
    function getHT(ne, T0, Efield, temp, t) result(nHT)

        real(dp), intent(in) :: ne     ! electron density [m^-3]
        real(dp), intent(in) :: T0    ! initial temperature [eV]
        real(dp), intent(in) :: Efield  ! instantaneous electric field at time t [V/m]
        real(dp), intent(in) :: temp  ! temperature at time t [eV]
        real(dp), intent(in) :: t   ! time [s]
        !real(dp), intent(in) :: td ! usually 0 but added here
        real(dp) :: nHT

        integer, parameter :: N = 50000
        real(dp), parameter :: xmax = 50.0_dp
        real(dp), parameter :: lnA = 15.0_dp

        ! put other stuff used here
        real(dp) :: vT0, vT
        real(dp) :: nu0, tau
        real(dp) :: vc
        real(dp) :: dx
        real(dp), dimension(N) :: x, v, ftot, integrand
        integer :: i, i0

        ! thermal velocities
        vT0 = sqrt(2.0_dp * T0 * e / me)
        vT  = sqrt(2.0_dp * temp  * e / me)

        ! col freq
        nu0 = getnu0(ne, T0)
        tau = getTau(nu0, t, 0.0_dp)

    
        dx = xmax / real(N-1, dp) ! for the linspace

        do i = 1, N
            x(i) = dx * real(i-1, dp)
            v(i) = sqrt(x(i)) * vT0
        end do

        ! get the total distribution 
        do i = 1, N
!~             ftot(i) = getfm(v(i), ne, vT) + &
!~                     getf(v(i), ne, vT0, tau) ! cold + hot pop. 
            ftot(i) = getf(v(i), ne, vT0, tau) ! cold + hot pop. 
        end do

        ! critical v at the instantaneous electric field
        vc = get_vc(temp, Efield, ne)
        !print *, "vc =", vc


        ! integral doesn't start from 0 so find the first index of v over v_critical

        i0 = N + 1

        do i = 1, N
            if (v(i) >= vc) then
                i0 = i
                exit
            end if
        end do

        ! array for integrand  
       ! do i = i0, N
            !integrand(i) = 4.0_dp*pi * (v(i)**2 - vc**2) * ftot(i)
        !    integrand(i) = 4.0_dp*pi * (v(i)-vc)*(v(i)+vc) * ftot(i)
       ! end do

       
       if (i0 <= N) then
            integrand = 0.0_dp

            do i = i0, N
                !integrand(i) = 4.0_dp*pi * (v(i)**2 - vc**2) * ftot(i)
                integrand(i) = 4.0_dp*pi * max(0.0_dp, v(i)**2 - vc**2) * ftot(i)
            end do

            nHT = integrate(v(i0:N), integrand(i0:N))

        else
            nHT = 0.0_dp

        end if
        !print *, nHT

    end function getHT

    function get_HT_rate(ne, temp, t, T0, Efield, nht, dt, nht_new) result(rate)

        real(dp), intent(in)  :: ne   ! electron density [m^-3]
        real(dp), intent(in)  :: temp ! temperature at time t [eV]
        real(dp), intent(in)  :: t ! time [s]
        real(dp), intent(in)  :: T0  ! initial temperature [eV]
        real(dp), intent(in)  :: Efield ! instantaneous electric field at time t [V/m]
        real(dp), intent(in)  :: nht ! previous HT density [m^-3]
        real(dp), intent(in)  :: dt ! time step [s]
        !real(dp), intent(in)  :: td

        real(dp), intent(out) :: nht_new
        real(dp) :: rate
        real(dp) :: nht_current

        ! current hot-tail density from distribution
        nht_current = getHT(ne, T0, Efield, temp, t)

        
        if (nht_current > nht) then
            nht_new = nht_current
        else
            nht_new = nht
        end if

        
        rate = (nht_new - nht) / dt
        !print *, rate

    end function get_HT_rate


end module HT_funcs
