module vacuum_interface

  implicit none

  integer :: nodes
  integer, allocatable :: node_id(:)
  real, allocatable :: normal(:,:)
  real, allocatable :: zgrbth(:,:), zgrbph(:,:)

  vectype, allocatable :: psi_mat(:,:), bf_mat(:,:)

contains 

  subroutine populate_boundary_matrices(ierr)
    implicit none

    integer, intent(out) :: ierr

    integer :: i, j, i1, j1
    real :: dtheta
    real :: ident

    ierr = 0

    allocate(psi_mat(6*nodes,6*nodes), bf_mat(6*nodes,6*nodes))

    psi_mat = 0.
    bf_mat = 0.

    dtheta = 2.*3.14159265/nodes
    do i=1,nodes
       i1 = 6*(i-1) + 1

       do j=1,i
          j1 = 6*(j-1) + 1

          if(i .eq. j) then
             ident = 1.
          else 
             ident = 0.
          end if
          psi_mat(i1  ,j1+1) = ident - dtheta*zgrbth(i,j)*normal(1,j)
          psi_mat(i1  ,j1+2) = ident - dtheta*zgrbth(i,j)*normal(2,j)

          psi_mat(i1+1,j1+3) = ident - dtheta*zgrbth(i,j)*normal(1,j)
          psi_mat(i1+1,j1+4) = ident - dtheta*zgrbth(i,j)*normal(2,j)

          psi_mat(i1+2,j1+4) = ident - dtheta*zgrbth(i,j)*normal(1,j)
          psi_mat(i1+2,j1+5) = ident - dtheta*zgrbth(i,j)*normal(2,j)
       end do
    end do
    
  end subroutine populate_boundary_matrices

  subroutine load_boundary_nodes(ierr)
    implicit none
    
    integer, intent(out) :: ierr

    character(len=*), parameter :: filename = 'ordered.points'
    integer, parameter :: ifile = 1
    integer :: i
    real :: dummy

    ierr = 0

    open(unit=ifile,file=filename,action='read',status='unknown')

    ! read total number of nodes
    read(ifile, '(I8)') nodes

    allocate(node_id(nodes))
    allocate(normal(2,nodes))

    ! write node entry for first node
    do i=1, nodes
       read(ifile, '(I8,4f12.8)') node_id(i), dummy, dummy, &
            normal(1,i), normal(2,i)
    end do

    close(ifile)
    
  end subroutine load_boundary_nodes


  subroutine load_response_matrix(ierr)
    implicit none

    integer, intent(out) :: ierr
    
    character(len=*), parameter :: filename = 'RESPONSE-M3DC1'
    integer, parameter :: ifile = 2
    integer :: idum, i, j

    ierr = 0

    open(unit=ifile,file=filename,action='read',status='unknown')

    read(ifile, '(5x, "Response Matrix in theta space for M3D_C1" )')

    read(ifile, '(/,1x, "No of internal points used = ", i5 )') idum 
    write(*, '(/,1x, "No of internal points used = ", i5 )') idum 

    read(ifile, '(/,1x, "Min Fourier mode used = ", i5 )') idum
    write(*, '(/,1x, "Min Fourier mode used = ", i5 )') idum

    read(ifile, '(  1x, "Max Fourier mode used = ", i5 )' ) idum
    write(*, '(  1x, "Max Fourier mode used = ", i5 )' ) idum

    read(ifile, '(/,1x, "Toroidal mode no. = ", i5 )' ) idum
    write(*, '(/,1x, "Toroidal mode no. = ", i5 )' ) idum

    read(ifile, &
         '(/,1x, "Number of surface points on the closed domain = ", i5 )' ) &
         idum

    if(idum .ne. nodes+1) then 
       print *, 'Error: nodes in response file different from nodes in ordered.points', idum-1, nodes
       ierr = 1
       goto 100
    endif

    allocate(zgrbph(nodes+1,nodes+1),zgrbth(nodes+1,nodes+1))

    read(ifile, & 
         '(/,1x, "B_theta response Matrix, Rth(obs,srce):" )' )

    do i=1, nodes+1
       read(ifile, '(/, 1x, "i_obs = ", i5 )' ) idum
       read(ifile, '( (1x, 8es14.6) )' ) (zgrbth(idum,j), j=1, nodes+1)
    END DO
    
    read(ifile, '(/,1x, "IMAG. B_phi response Matrix, Rph(obs,src):" )')
    
    do i=1, nodes+1
       read(ifile, '(/, 1x, "i_obs = ", i5 )' ) idum
       read(ifile, '( (1x, 8es14.6) )' ) (zgrbph(idum,j), j=1, nodes+1)
    end do
    
100 continue

    close(ifile)
    
  end subroutine load_response_matrix

  subroutine load_vacuum_data(ierr)
    implicit none

    integer, intent(out) :: ierr

    call load_boundary_nodes(ierr)
    if(ierr .ne. 0) return
    
    call load_response_matrix(ierr)
    if(ierr .ne. 0) return

    call populate_boundary_matrices(ierr)
    if(ierr .ne. 0) return
       
  end subroutine load_vacuum_data

  subroutine unload_vacuum_data
    implicit none

    print *, 'unloading vacuum data'
    if(allocated(node_id)) deallocate(node_id)
    if(allocated(normal)) deallocate(normal)
    if(allocated(zgrbth)) deallocate(zgrbth)
    if(allocated(zgrbph)) deallocate(zgrbph)
    if(allocated(psi_mat)) deallocate(psi_mat)
    if(allocated(bf_mat)) deallocate(bf_mat)
  end subroutine unload_vacuum_data


end module vacuum_interface
