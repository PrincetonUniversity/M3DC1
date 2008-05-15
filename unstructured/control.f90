module pid_controller

  type :: pid_control
     real :: p, i, d
     real :: err_p_old = 0.
     real :: err_i = 0.
     real :: target_val
  end type pid_control
  
contains

  subroutine control(val, control_param, pid, dt)
    implicit none

    real, intent(in) :: val
    real, intent(inout) :: control_param
    type(pid_control), intent(inout) :: pid
    real, intent(in) :: dt

    real :: err_p, err_d

    if(dt.eq.0.) return

    err_p = val - pid%target_val
    pid%err_i = pid%err_i + err_p*dt
    err_d = (err_p - pid%err_p_old)/dt

    control_param = control_param - control_param*dt* &
         (err_p*pid%p + pid%err_i*pid%i + err_d*pid%d)
    
    pid%err_p_old = err_p

  end subroutine control

end module pid_controller
