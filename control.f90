module pid_controller

  type :: pid_control
     real :: p, i, d
     real :: err_p_old
     real :: err_i
     real :: target_val
     integer :: icontrol_type
  end type pid_control
  
contains

  subroutine control(val, control_param, pid, dt)
    implicit none

    real, intent(in) :: val
    real, intent(inout) :: control_param
    type(pid_control), intent(inout) :: pid
    real, intent(in) :: dt

    real :: err_p, err_d

      err_p = val - pid%target_val
      pid%err_i = pid%err_i + err_p*dt
      if(dt.gt.0) then
         err_d = (err_p - pid%err_p_old)/dt
      endif

    select case (pid%icontrol_type)
    case(0)   ! this was the original coding and is the default
      if(pid%target_val.gt.0) then
        control_param = control_param - control_param*dt* &
            (err_p*pid%p + pid%err_i*pid%i + err_d*pid%d)
      else
        control_param = control_param + control_param*dt* &
            (err_p*pid%p + pid%err_i*pid%i + err_d*pid%d)
      endif
!
!
    case(1)   ! this is the "standard" PID controller

      control_param = - (err_p*pid%p + pid%err_i*pid%i + err_d*pid%d)
    end select


    if(dt.eq.0) return
    pid%err_p_old = err_p

  end subroutine control

end module pid_controller
