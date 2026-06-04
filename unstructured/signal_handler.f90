module signal_handler
  use iso_c_binding
  implicit none

  ! Use integer, not logical, for signal safety
  integer(c_int), volatile :: timeout_flag = 0

  private
  public :: install_signal_handler, timeout_flag

  ! C struct sigset_t (size is platform dependent)
  ! On Linux x86_64 it is 128 bytes
  integer, parameter :: SIGSET_SIZE = 128

  type, bind(C) :: sigset_t
     integer(c_int8_t) :: val(SIGSET_SIZE)
  end type sigset_t

  type, bind(C) :: sigaction_t
     type(c_funptr)    :: sa_handler
     type(sigset_t)    :: sa_mask
     integer(c_int)    :: sa_flags
     type(c_funptr)    :: sa_restorer  ! not used but required on glibc
  end type sigaction_t

  interface
     function c_sigaction(signum, act, oldact) bind(C, name="sigaction")
       import :: c_int, sigaction_t
       integer(c_int), value :: signum
       type(sigaction_t), intent(in)  :: act
       type(sigaction_t), intent(out) :: oldact
       integer(c_int) :: c_sigaction
     end function c_sigaction

     function c_sigemptyset(set) bind(C, name="sigemptyset")
       import :: sigset_t, c_int
       type(sigset_t), intent(out) :: set
       integer(c_int) :: c_sigemptyset
     end function c_sigemptyset
  end interface

contains

  subroutine handle_signal(sig) bind(C)
    use iso_c_binding
    integer(c_int), value :: sig

    ! ONLY set a flag — nothing else!
    timeout_flag = 1
  end subroutine handle_signal


  subroutine install_signal_handler()
    type(sigaction_t) :: act, oldact
    integer(c_int) :: ierr

    integer(c_int), parameter :: SIGUSR1 = 10_c_int

    ! Set handler
    act%sa_handler = c_funloc(handle_signal)
    act%sa_flags   = 0_c_int
    act%sa_restorer = c_null_funptr

    ! Clear mask
    ierr = c_sigemptyset(act%sa_mask)

    ! Install handler
    ierr = c_sigaction(SIGUSR1, act, oldact)

  end subroutine install_signal_handler

end module signal_handler
