module ouput

contains

  subroutine write_boundary_bn
    use basic
    use arrays
    use vacuum_interface
    implicit none

#include "mpif.h"

    integer :: i, ierr
    complex :: bn(nodes), buff(nodes)
    real :: xbuff(nodes), zbuff(nodes)
    vectype, dimension(dofs_per_node) :: psi_data, f_data
    logical :: is_boundary
    integer :: izone, izonedim
    real ::  normal(2), curv, x(nodes), z(nodes)
    integer, parameter :: bnout=77
    character(len=*), parameter :: bnfilename = "bn.out"

    bn = 0
    do i=1, nodes
       if(local_id(i).le.0) cycle

       call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x(i),z(i))

       call get_node_data(psi_field(1), i, psi_data)
       
       bn(i) = (normal(1)*psi_data(3) - normal(2)*psi_data(2))/x(i)
#ifdef USECOMPLEX
       call get_node_data(bf_field(1), i, f_data)
       bn(i) = bn(i) - rfac*(normal(1)*f_data(2) + normal(2)*f_data(3))
#endif
    end do

    call mpi_reduce(bn, buff, 2*nodes, MPI_DOUBLE_PRECISION, MPI_SUM, &
         0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(x, xbuff, nodes, MPI_DOUBLE_PRECISION, MPI_SUM, &
         0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(z, zbuff, nodes, MPI_DOUBLE_PRECISION, MPI_SUM, &
         0, MPI_COMM_WORLD, ierr)

    if(myrank.eq.0) then
       open(unit=bnout, file=bnfilename, status='unknown') 
       write(bnout, '(I5)') nodes
       do i=1, nodes
          write(bnout, '(4f12.6)') xbuff(i), zbuff(i), buff
       end do
       close(bnout)
    endif

  end subroutine write_boundary_bn

end module ouput
