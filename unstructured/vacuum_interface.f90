module vacuum_interface

  implicit none

  integer :: nodes, boundary_edges
  integer, allocatable :: global_id(:)    ! global dof of boundary node i
  integer, allocatable :: local_id(:)     ! local dof of boundary node i
  integer, allocatable :: edge_elm(:)     ! local elm of boundary node i
  integer, allocatable :: edge_nodes(:,:) ! boundary nodes of boundary edge j
  integer, allocatable :: edge_edge(:)    ! elm node of boundary edge i
  real, allocatable :: xnode(:),znode(:)   ! position of boundary node i
  real, allocatable :: nxnode(:),nznode(:) ! normal vector at boundary node i
  complex, allocatable :: zgrbth(:,:), zgrbph(:,:)   ! response matrices
  complex, allocatable :: zgrbthp(:,:), zgrbphp(:,:) ! tangential derivs

contains 

  subroutine load_boundary_nodes(ierr)
    use mesh_mod
    implicit none
        
    integer, intent(out) :: ierr

    character(len=*), parameter :: filename = 'ordered.points'
    integer, parameter :: ifile = 1
    integer :: i, j, k, numnodes, num_local, n, numelms, itri
    integer :: inode(nodes_per_element)
    real :: x, z
    real, parameter :: tol = 1e-5
    logical :: is_edge(3)  ! is inode on boundary
    real :: norm(2,3)
    integer :: iedge, idim(3)
    
    ierr = 0

    open(unit=ifile,file=filename,action='read',status='unknown')

    ! read toroidal mode number
    read(ifile, '(I8)') n

    ! read total number of nodes
    read(ifile, '(I8)') nodes

    allocate(global_id(nodes))
    allocate(local_id(nodes))
    allocate(xnode(nodes+1),nxnode(nodes+1))
    allocate(znode(nodes+1),nznode(nodes+1))

    ! read list of nodes
    do i=1, nodes
       read(ifile, '(I8,4f12.6)') global_id(i), xnode(i), znode(i), &
            nxnode(i), nznode(i)
    end do

    close(ifile)

    xnode(nodes+1) = xnode(1)
    znode(nodes+1) = znode(1)
    
    ! allocate response matrices
    allocate(zgrbph (nodes+1,nodes+1),zgrbth (nodes+1,nodes+1))
    allocate(zgrbphp(nodes+1,nodes+1),zgrbthp(nodes+1,nodes+1))


    ! find local and global ids of each node
    numnodes = local_nodes()
    num_local = 0
    do i=1, nodes
       local_id(i) = -1
       do j=1, numnodes
          call get_node_pos(j,x,z)
          if((abs(x - xnode(i)).lt.tol) .and. (abs(z-znode(i)).lt.tol)) then
             local_id(i) = j
             itri = global_node_id(local_id(i))
             if(itri.ne.global_id(i)) then
                print *, "Error: global ids don't match: ", local_id(i), itri, global_id(i), x, z
                call safestop(0)
             endif
             num_local = num_local + 1
          endif
       end do
    end do
    print *, 'Found ', num_local, ' local boundary nodes.'

    ! compile list of boundary edges
    allocate(edge_elm(nodes),edge_nodes(2,nodes),edge_edge(nodes))
    boundary_edges = 0
    numelms = local_elements()
    do itri=1, numelms
       call boundary_edge(itri, is_edge, norm, idim)
     
       do iedge=1,3
          if(.not.is_edge(iedge)) cycle

          call get_element_nodes(itri, inode)

          boundary_edges = boundary_edges + 1
          edge_elm(boundary_edges) = itri
          edge_nodes(1,boundary_edges) = inode(iedge)
          edge_nodes(2,boundary_edges) = inode(mod(iedge,nodes_per_element)+1)
          edge_edge(boundary_edges) = iedge
       end do
    end do
    ! convert from local node numbering to boundary numbering
    do i=1, boundary_edges
       do j=1,nodes_per_edge
          do k=1, nodes
             if(edge_nodes(j,i).eq.local_id(k)) then
                edge_nodes(j,i) = k
                exit
             endif
          end do
       end do
    end do
    print *, 'Found ', boundary_edges, ' local boundary edges.'
    
  end subroutine load_boundary_nodes


  subroutine load_response_matrix(n, ierr)
    use math

    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: ierr
    
    character(len=*), parameter :: filename = 'RESPONSE-M3DC1'
    integer, parameter :: ifile = 2
    integer :: idum, i, j
    real :: dtheta

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
       close(ifile)
       call safestop(6)
    endif

    read(ifile, '(/,1x, "B_theta response Matrix, Rth(obs,srce):" )' )

    do i=1, nodes+1
       read(ifile, '(/, 1x, "i_obs = ", i5 )' ) idum
       read(ifile, '( (1x, 8es14.6) )' ) (zgrbth(idum,j), j=1, nodes+1)
    end do
    
    read(ifile, '(/,1x, "IMAG. B_phi response Matrix, Rph(obs,src):" )')
    
    do i=1, nodes+1
       read(ifile, '(/, 1x, "i_obs = ", i5 )' ) idum
       read(ifile, '( (1x, 8es14.6) )' ) (zgrbph(idum,j), j=1, nodes+1)
    end do

    dtheta = -twopi/nodes
!    dtheta = twopi/nodes

    ! multiply by dtheta/twopi
    zgrbth = zgrbth/nodes
    zgrbph = zgrbph/nodes

    zgrbth = -zgrbth

!    zgrbph = -zgrbph

    ! calculate derivatives (wrt i)
    zgrbthp(1,:) = 0.5*(zgrbth(2,:) - zgrbth(nodes+1,:))/dtheta
    zgrbphp(1,:) = 0.5*(zgrbph(2,:) - zgrbph(nodes+1,:))/dtheta
    do i=2,nodes
       zgrbthp(i,:) = 0.5*(zgrbth(i+1,:) - zgrbth(i-1,:))/dtheta
       zgrbphp(i,:) = 0.5*(zgrbph(i+1,:) - zgrbph(i-1,:))/dtheta
    end do
    
    close(ifile)
  end subroutine load_response_matrix

  subroutine load_vacuum_data(n, ierr)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: ierr

    call load_boundary_nodes(ierr)
    if(ierr .ne. 0) return
    
    call load_response_matrix(n, ierr)
    if(ierr .ne. 0) return
      
  end subroutine load_vacuum_data

  subroutine unload_vacuum_data
    implicit none

    print *, 'unloading vacuum data'
    if(allocated(global_id)) deallocate(global_id)
    if(allocated(local_id)) deallocate(local_id)
    if(allocated(zgrbth)) deallocate(zgrbth)
    if(allocated(zgrbph)) deallocate(zgrbph)
    if(allocated(zgrbthp)) deallocate(zgrbthp)
    if(allocated(zgrbphp)) deallocate(zgrbphp)
    if(allocated(xnode)) deallocate(xnode)
    if(allocated(znode)) deallocate(znode)
    if(allocated(nxnode)) deallocate(nxnode)
    if(allocated(nznode)) deallocate(nznode)
    if(allocated(edge_elm)) deallocate(edge_elm)
    if(allocated(edge_nodes)) deallocate(edge_nodes)
    if(allocated(edge_edge)) deallocate(edge_edge)

  end subroutine unload_vacuum_data


end module vacuum_interface
