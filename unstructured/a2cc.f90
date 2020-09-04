program a2cc
  use eqdsk_a

  implicit none

  integer :: args
  character(len=256) :: filename

  args = command_argument_count()
  if(args.lt.1) then
     write(0,*) 'Usage: a2cc <aeqdsk>'
     stop
  end if
     
  call getarg(1, filename)

  call load_eqdsk_a(filename)

end program a2cc
