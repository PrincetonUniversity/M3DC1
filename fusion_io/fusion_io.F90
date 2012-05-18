module fusion_io

#include "fusion_io_defs.h"


contains

  subroutine fio_add_field(icfield, ifield, iop, fac, ierr)
    integer, intent(in) :: icfield
    integer, intent(in) :: ifield
    integer, intent(in) :: iop
    real, intent(in) :: fac
    integer, intent(out) :: ierr
    call add_field(icfield, ifield, iop, fac, ierr)
  end subroutine fio_add_field

  subroutine fio_close_field(ifield, ierr)
    integer, intent(in) :: ifield
    integer, intent(out) :: ierr
    call close_field(ifield, ierr)
  end subroutine fio_close_field

  subroutine fio_close_source(isrc, ierr)
    integer, intent(in) :: isrc
    integer, intent(out) :: ierr
    call close_source(isrc, ierr)
  end subroutine fio_close_source

  subroutine fio_create_compound_field(ifield, ierr)
    integer, intent(out) :: ifield
    integer, intent(out) :: ierr
    call create_compound_field(ifield, ierr)
  end subroutine fio_create_compound_field

  subroutine fio_eval_field(ifield, x, v, ierr)
    integer, intent(in) :: ifield
    real, intent(in), dimension(*) :: x
    real, intent(out), dimension(*) :: v
    integer, intent(out) :: ierr
    call eval_field(ifield, x, v, ierr)
  end subroutine fio_eval_field

  subroutine fio_get_field(isrc, itype, ifield, ierr)
    integer, intent(in) :: isrc, itype
    integer, intent(out) :: ifield, ierr
    call get_field(isrc, itype, ifield, ierr)
  end subroutine fio_get_field

  subroutine fio_get_options(isrc, ierr)
    integer, intent(in) :: isrc
    integer, intent(out) :: ierr
    call get_options(isrc, ierr)
  end subroutine fio_get_options
  
  subroutine fio_open_source(type, filename, isrc, ierr)
    integer, intent(in) :: type
    character(len=*), intent(in) :: filename
    integer, intent(out) :: isrc, ierr   
    call open_source(type, filename//char(0), isrc, ierr)
  end subroutine fio_open_source

  subroutine fio_set_int_option(iopt, val, ierr)
    integer, intent(in) :: iopt
    integer, intent(in) :: val
    integer, intent(out) :: ierr
    call set_int_option(iopt, val, ierr)
  end subroutine fio_set_int_option

  subroutine fio_set_real_option(iopt, val, ierr)
    integer, intent(in) :: iopt
    real, intent(in) :: val
    integer, intent(out) :: ierr
    call set_real_option(iopt, val, ierr)
  end subroutine fio_set_real_option

  subroutine fio_set_str_option(iopt, val, ierr)
    integer, intent(in) :: iopt
    character(len=*), intent(in) :: val
    integer, intent(out) :: ierr
    call set_str_option(iopt, val//char(0), ierr)
  end subroutine fio_set_str_option
end module fusion_io
