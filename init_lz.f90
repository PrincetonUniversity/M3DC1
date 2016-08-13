module LZeqbm
  IMPLICIT NONE

  !Parameters
  REAL, PARAMETER :: aspect_ratio = 18.0  ! R_major / r_wall_b
  REAL, PARAMETER :: B0 = 1.0             ! Toroidal field at r=0
  REAL, PARAMETER :: dfrac = 0.1          ! Bndry layer width / r_plas_a
  REAL, PARAMETER :: Trat = 100.0         ! Ratio of plas to vac temperature
  REAL, PARAMETER :: nrat = 1.0           ! Ratio of plas to vac density
  REAL, PARAMETER :: q_interior = 1.0     ! Safety factor across bulk plasma
  REAL, PARAMETER :: q_a = 0.75           ! Safety factor on plasma surface
  REAL, PARAMETER :: r_plas_a = 0.6       ! Plasma minor radius
  REAL, PARAMETER :: r_wall_b = 1.0       ! Minor radius of ideal wall
  REAL, PARAMETER :: r_tile_c = 0.7       ! Minor radius of tile surface

  !Derived quantities
  REAL R_major                            ! Major radius of magnetic axis
  REAL delta                              ! Width of boundary layer
  REAL r_interior                         ! Radius to interior of bndry layer
  REAL prat                               ! Ratio pf plas to vac pressure
  REAL p_interior                         ! Pressure in bulk plasma

  REAL p_vacuum                           ! Pressure in vacuum region
  REAL q0R0                               ! Central safety factor * major radius
  REAL Jdenom                             ! Intermediate derived quantity
  REAL Bbz                                ! Toroidal field in bndry layer
  REAL gamma                              ! Pressure gradient in bndry layer
 
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Profiles for a pseudo-VMEC equilibrium structure resembling L. Zakharov's    !
! circular cross-section, straight-cylinder, flat q, ideal equilibrium with a  !
! surface current.                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setupLZeqbm

    IMPLICIT NONE

    delta = dfrac * r_plas_a
    r_interior = r_plas_a - delta
    R_major = aspect_ratio * r_wall_b
    q0R0 = q_interior * R_major
    Jdenom = q0R0**2 + r_interior**2
    Bbz = q0R0**2 * B0 / Jdenom
    gamma = (1.5 * (q0R0**2 * B0**2 * r_interior**4 / Jdenom**2 - &
         r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / &
         (r_interior**3 - r_plas_a**3)) ! /aspect_ratio
    prat = nrat * Trat
    p_vacuum = gamma*delta/(prat - 1.0)
    p_interior = p_vacuum + gamma*delta
  end subroutine setupLZeqbm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function p_LZ(psi)

    implicit none
    REAL p_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       p_LZ = p_interior
    ELSE IF (rmin < r_plas_a) THEN
       p_LZ = p_interior - gamma*(rmin - r_interior)
    ELSE
       p_LZ = p_vacuum
    endIF
  end function p_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pprime_LZ(psi)

    implicit none
    REAL pprime_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, dpdr

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       pprime_LZ = 0.0
       RETURN
    endIF

    IF (rmin < r_plas_a) THEN
       dpdr = -gamma
    ELSE
       dpdr = 0.0
    endIF

    pprime_LZ = 0.5*dpdr/rmin
  end function pprime_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function q_LZ(psi)

    implicit none
    intrinsic sqrt
    REAL q_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, Btheta


    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       q_LZ = q_interior
    ELSE IF (rmin < r_plas_a) THEN
       Btheta = SQRT(2.0*gamma*(rmin**3 - r_plas_a**3)/3.0 + &
            r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / rmin
       q_LZ = rmin * Bbz / (R_major * Btheta)
    ELSE
       q_LZ = q_a * (rmin/r_plas_a)**2
    endIF
  end function q_LZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function qprime_LZ(psi)

    implicit none
    REAL qprime_LZ
    REAL, INTENT(IN) :: psi  ! square of minor radius
    REAL rmin, Btheta, Bthetaprime, dqdr

    rmin = SQRT(psi)

    IF (rmin <= r_interior) THEN
       qprime_LZ = 0.0
       RETURN
    endIF

    IF (rmin < r_plas_a) THEN
       Btheta = SQRT(2.0*gamma*(rmin**3 - r_plas_a**3)/3.0 + &
            r_plas_a**4 * Bbz**2 / (q_a**2 * R_major**2)) / rmin
       dqdr = Bbz * (2.0 - gamma*rmin/Btheta**2) / (R_major * Btheta)
    ELSE
       dqdr = 2.0*rmin*q_a / r_plas_a**2
    endIF

    qprime_LZ = 0.5*dqdr/rmin
  end function qprime_LZ
end module LZeqbm
