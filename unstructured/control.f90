subroutine control_pid

  use basic
  use diagnostics

  implicit none

  real :: error_p, error_d
  real, save :: error_old = 0.
  real, save :: error_i = 0.

  error_p = totcur - tcuro
  error_i = error_i + error_p*dt
  error_d = (error_p - error_old)/dt

  vloop = vloop*(1.- control_p*(error_p + control_i*error_i + control_d*error_d))

  error_old = error_p

end subroutine control_pid
