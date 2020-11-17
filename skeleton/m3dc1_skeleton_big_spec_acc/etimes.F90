module etimes

  implicit none

  real :: etime_define_fields      = 0.
  real :: etime_eval_ops           = 0.
  real :: etime_matdefphi          = 0.
  real :: etime_precalculate_terms = 0.

contains

  subroutine get_etime(t)

    implicit none

    real, intent(out) :: t
    integer*8 kount, kountrate, kountmax

    call system_clock(kount, kountrate, kountmax)

    t = real(kount) / real(kountrate)

  end subroutine get_etime

end module etimes
