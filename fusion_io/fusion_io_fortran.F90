module fusion_io

#include "fusion_io_defs.h"


contains

  subroutine fio_add_field_f(icfield, ifield, iop, fac, ierr)
    integer, intent(in) :: icfield
    integer, intent(in) :: ifield
    integer, intent(in) :: iop
    real, intent(in) :: fac
    integer, intent(out) :: ierr
    call fio_add_field(icfield, ifield, iop, fac, ierr)
  end subroutine fio_add_field_f

  subroutine fio_close_field_f(ifield, ierr)
    integer, intent(in) :: ifield
    integer, intent(out) :: ierr
    call fio_close_field(ifield, ierr)
  end subroutine fio_close_field_f

  subroutine fio_close_source_f(isrc, ierr)
    integer, intent(in) :: isrc
    integer, intent(out) :: ierr
    call fio_close_source(isrc, ierr)
  end subroutine fio_close_source_f

  subroutine fio_create_compound_field_f(ifield, ierr)
    integer, intent(out) :: ifield
    integer, intent(out) :: ierr
    call fio_create_compound_field(ifield, ierr)
  end subroutine fio_create_compound_field_f

  subroutine fio_eval_field_f(ifield, x, v, ierr)
    integer, intent(in) :: ifield
    real, intent(in), dimension(*) :: x
    real, intent(out), dimension(*) :: v
    integer, intent(out) :: ierr
    call fio_eval_field(ifield, x, v, ierr)
  end subroutine fio_eval_field_f

  subroutine fio_get_field_f(isrc, itype, ispec, ifield, ierr)
    integer, intent(in) :: isrc, itype, ispec
    integer, intent(out) :: ifield, ierr
    call fio_get_field(isrc, itype, ispec, ifield, ierr)
  end subroutine fio_get_field_f

  subroutine fio_get_options_f(isrc, ierr)
    integer, intent(in) :: isrc
    integer, intent(out) :: ierr
    call fio_get_options(isrc, ierr)
  end subroutine fio_get_options_f
  
  subroutine fio_open_source_f(type, filename, isrc, ierr)
    integer, intent(in) :: type
    character(len=*), intent(in) :: filename
    integer, intent(out) :: isrc, ierr   
    call fio_open_source(type, filename//char(0), isrc, ierr)
  end subroutine fio_open_source_f

  subroutine fio_set_int_option_f(iopt, val, ierr)
    integer, intent(in) :: iopt
    integer, intent(in) :: val
    integer, intent(out) :: ierr
    call fio_set_int_option(iopt, val, ierr)
  end subroutine fio_set_int_option_f

  subroutine fio_set_real_option_f(iopt, val, ierr)
    integer, intent(in) :: iopt
    real, intent(in) :: val
    integer, intent(out) :: ierr
    call fio_set_real_option(iopt, val, ierr)
  end subroutine fio_set_real_option_f

  subroutine fio_set_str_option_f(iopt, val, ierr)
    integer, intent(in) :: iopt
    character(len=*), intent(in) :: val
    integer, intent(out) :: ierr
    call fio_set_str_option(iopt, val//char(0), ierr)
  end subroutine fio_set_str_option_f
end module fusion_io
