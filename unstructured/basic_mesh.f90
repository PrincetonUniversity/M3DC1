module basic_mesh_mod
  use element

  implicit none

  logical :: is_rectilinear

  real :: xzero, zzero
  real :: tiltangled
  integer :: iper, jper
  integer :: icurv
  integer :: nplanes

  integer, parameter :: maxi = 20

  logical, private :: has_bounding_box = .false.
  real, private :: bb(4)

  logical, parameter :: IGNORE_PHI = .true.

  type, private :: node_type
     integer :: izone, idim
     real :: R, Phi, Z
     real :: normal(2), curv(3)
  end type node_type

  type, private :: element_type
     integer :: inode(nodes_per_element)
  end type element_type

  integer, private :: num_global_nodes
  integer, private :: num_local_nodes     ! local nodes, excluding ghosts
  integer, private :: total_local_nodes   ! local nodes, including ghosts
  integer, private :: num_local_elements
  integer, private :: nodes_per_plane

  type(node_type), allocatable, private :: local_node(:)
  type(element_type), allocatable, private :: local_elm(:)
  integer, private :: num_ghost_nodes
  integer, allocatable, target, private :: ghost_node(:)
  integer, allocatable, private :: global_id(:)
  integer, allocatable :: num_adjacent(:)

  integer, private :: myplane ! toroidal plane of current process

contains

  !======================================================================
  ! load_sms_mesh
  ! ~~~~~~~~~~~~~
  ! load mesh from symmetrix mesh file
  ! if nplanes > 1, create 3D mesh using nplanes copies of 2D mesh
  !======================================================================
  subroutine load_sms_mesh(nplanes, inodes, ielms, nodes, elms)
    use math

    implicit none

    integer, intent(in) :: nplanes
    integer, intent(out) :: inodes, ielms
    type(node_type), allocatable, intent(inout) :: nodes(:)
    type(element_type), allocatable, intent(inout) :: elms(:)

    integer :: i, j
    integer, parameter :: ifile = 87
    character(len=128) :: buffer
    integer :: inodespp, ielmspp, iedges, dum
    integer, allocatable :: edges(:,:)
    integer, dimension(3) :: ed

    real :: x0, x1, x2, z0, z1, co, sn, cs, ss, l0
    
    if(is_rectilinear) then
       open(unit=ifile, file='struct-dmg.sms', status='old')
    else
       open(unit=ifile, file='struct-curveDomain.sms', status='old')
    end if

    read(ifile,*) buffer
    read(ifile,*) dum, ielmspp, iedges, inodespp, dum
    write(*,*) 'Reading file with '
    write(*,*) ' elements = ', ielmspp
    write(*,*) ' edges = ', iedges
    write(*,*) ' nodes = ', inodespp
    write(*,*) ' nplanes = ', nplanes
    if(nplanes.gt.1) then
       ielms = ielmspp*nplanes
       inodes = inodespp*nplanes
    else
       ielms = ielmspp
       inodes = inodespp
    end if

    allocate(nodes(inodes), elms(ielms), edges(iedges,2))

    do i=1, inodespp
       read(ifile,*) dum, nodes(i)%idim, dum
       read(ifile,*) nodes(i)%R, nodes(i)%Z, nodes(i)%Phi
    end do
    do i=1, iedges
       read(ifile,*) dum, dum, edges(i,1), edges(i,2), dum, dum
    end do
    do i=1, ielmspp
       read(ifile,*) dum, dum, dum, ed(1), ed(2), ed(3), dum
       do j=1, 3
          if(ed(j).gt.0) then 
             elms(i)%inode(j) = edges(ed(j),1)
          else
             elms(i)%inode(j) = edges(-ed(j),2)
          end if
       end do
    end do

    close(ifile)

    deallocate(edges)

    if(is_rectilinear) then
       do i=1, inodes
          nodes(i)%R = nodes(i)%R + xzero
          nodes(i)%Z = nodes(i)%Z + zzero
          if(nodes(i)%R.lt.bb(1) .or. i.eq.1) bb(1) = nodes(i)%R
          if(nodes(i)%R.gt.bb(3) .or. i.eq.1) bb(3) = nodes(i)%R
          if(nodes(i)%Z.lt.bb(2) .or. i.eq.1) bb(2) = nodes(i)%Z
          if(nodes(i)%Z.gt.bb(4) .or. i.eq.1) bb(4) = nodes(i)%Z
       end do
       do i=1, inodes
          if(nodes(i)%idim .ge.2) cycle
          
          nodes(i)%curv = 0
          if(nodes(i)%R.eq.bb(1)) then
             nodes(i)%normal(1) = -1.
             nodes(i)%normal(2) = 0.
          else if(nodes(i)%R.eq.bb(3)) then
             nodes(i)%normal(1) = 1.
             nodes(i)%normal(2) = 0.
          else if(nodes(i)%Z.eq.bb(2)) then
             nodes(i)%normal(1) = 0.
             nodes(i)%normal(2) = -1.
          else if(nodes(i)%Z.eq.bb(4)) then
             nodes(i)%normal(1) = 0.
             nodes(i)%normal(2) = 1.
          else
             print *, 'Error! boundary undetermined.'
          end if
       end do
    else
       ! Calculate curv and normal from modelfile
       open(unit=ifile, file='AnalyticModel', status='old')
       read(ifile, *) x0, x1, x2, z0, z1
       close(ifile)

       ! R(theta) = x0 + x1 cos(theta + x2 sin(theta))
       ! Z(theta) = z0 + z1 sin(theta)
       do i=1, inodespp
          if(nodes(i)%idim .ge.2) cycle

          sn = (nodes(i)%Z - z0)/z1    ! sin(theta)
          co = sqrt(1. - sn*sn)        ! cos(theta)
          if(nodes(i)%R .lt. x0) co = -co
          cs = (nodes(i)%R - x0)/x1    ! cos(theta + x2*sin(theta))
          ss = sqrt(1. - cs*cs)        ! sin(theta + x2*sin(theta))
          if(nodes(i)%Z .lt. z0) ss = -ss

          l0 = sqrt((z1*co)**2 + (x1*(1.+x2*co)*ss)**2)
          nodes(i)%normal(1) = z1*co/l0
          nodes(i)%normal(2) = x1*(1.+x2*co)*ss/l0
          nodes(i)%curv(1) = x1*z1*(co*(1.+x2*co)**2*cs + sn*ss)/l0**3
!          write(*,'(6f12.4)') (nodes(i)%R-x0)/x1, (nodes(i)%Z-z0)/z1, &
!               nodes(i)%normal(1), nodes(i)%normal(2), nodes(i)%curv
       end do
    end if

#ifdef USE3D
    ! Create nodes in other planes
    do i=2, nplanes
       nodes((i-1)*inodespp+1:i*inodespp) = nodes(1:inodespp)
       nodes((i-1)*inodespp+1:i*inodespp)%Phi = (i-1.)*twopi/nplanes
    end do

    ! Set "ghost" plane of elements in first plane
    do i=1, ielmspp
       elms(i)%inode(4:6) = elms(i)%inode(1:3) + inodespp
    end do

    ! Set nodes in other planes
    do i=2, nplanes
       do j=1, ielmspp
          elms((i-1)*ielmspp+j)%inode(:) = elms(j)%inode(:) + (i-1)*inodespp
       end do
    end do

    ! Enforce periodic boundary conditions
    do i=1, ielmspp
       elms((nplanes-1)*ielmspp+i)%inode(4:6) = elms(i)%inode(1:3)
    end do
#endif

  end subroutine load_sms_mesh

  !======================================================================
  ! create_rectangular_mesh
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  ! creates a rectangular mesh with m x n nodes.
  ! physical dimensions are set by bb array.
  !======================================================================
  subroutine create_rectangular_mesh(m,n, nplanes,inodes, ielms, nodes, elms)
    use math
    implicit none

    integer, intent(in) :: m, n, nplanes
    integer, intent(out) :: inodes, ielms
    type(node_type), allocatable, intent(inout) :: nodes(:)
    type(element_type), allocatable, intent(inout) :: elms(:)

    real :: dx, dz, dphi
    integer :: i, j, k, l
    
    dx = (bb(3) - bb(1))/(m-1)
    dz = (bb(4) - bb(2))/(n-1)
    dphi = twopi/nplanes

    inodes = m*n*nplanes
    ielms  = 2*(m-1)*(n-1)*nplanes

    print *, 'allocating', inodes, ielms
    allocate(nodes(inodes), elms(ielms))

    print *, 'initializing'
    ! initialize nodes
    ! ~~~~~~~~~~~~~~~~
    k = 1
    do l=1, nplanes
       do i=1,m
          do j=1,n
             nodes(k)%R =   dx  *(j-1) + bb(1)
             nodes(k)%Z =   dz  *(i-1) + bb(2)
             nodes(k)%Phi = dphi*(l-1)
             nodes(k)%idim = 2
             if(j.eq.1 .or. j.eq.n) then
                nodes(k)%idim = nodes(k)%idim - 1
                nodes(k)%normal(2) = 0
                nodes(k)%curv = 0
                if(j.eq.1) nodes(k)%normal(1) = -1
                if(j.eq.n) nodes(k)%normal(1) = 1
             endif
             if(i.eq.1 .or. i.eq.m) then
                nodes(k)%idim = nodes(k)%idim - 1
                nodes(k)%normal(1) = 0
                nodes(k)%curv = 0
                if(i.eq.1) nodes(k)%normal(2) = -1
                if(i.eq.m) nodes(k)%normal(2) = 1
             endif
             k = k + 1
          end do
       end do
    end do

    ! initialize elements
    ! ~~~~~~~~~~~~~~~~~~~
    k = 1
    do l=1, nplanes
       do i=1,m-1
          do j=1,n-1
             elms(k)%inode(1) = (i-1)*n + (j  ) + (l-1)*n*m
             elms(k)%inode(2) = (i  )*n + (j+1) + (l-1)*n*m
             elms(k)%inode(3) = (i  )*n + (j  ) + (l-1)*n*m
#ifdef USE3D
             elms(k)%inode(4) = (i-1)*n + (j  ) + mod(l,nplanes)*n*m
             elms(k)%inode(5) = (i  )*n + (j+1) + mod(l,nplanes)*n*m
             elms(k)%inode(6) = (i  )*n + (j  ) + mod(l,nplanes)*n*m
#endif
             k = k + 1
             elms(k)%inode(1) = (i-1)*n + (j  ) + (l-1)*n*m
             elms(k)%inode(2) = (i-1)*n + (j+1) + (l-1)*n*m
             elms(k)%inode(3) = (i  )*n + (j+1) + (l-1)*n*m
#ifdef USE3D
             elms(k)%inode(4) = (i-1)*n + (j  ) + mod(l,nplanes)*n*m
             elms(k)%inode(5) = (i-1)*n + (j+1) + mod(l,nplanes)*n*m
             elms(k)%inode(6) = (i  )*n + (j+1) + mod(l,nplanes)*n*m
#endif
             k = k + 1
          end do
       end do
    end do
  end subroutine create_rectangular_mesh

  subroutine load_mesh

    use math

    implicit none

#include "finclude/petsc.h"
#include "finclude/petscis.h90"

    integer :: num_global_elements
    type(node_type), allocatable :: global_node(:)
    type(element_type), allocatable :: global_elm(:)
    Mat :: connectivity, adjacency
    MatPartitioning :: partitioning
    integer :: i, j, k, ierr, rank, size, itri
    integer :: min_node, max_node
    integer, allocatable :: local_id(:)
    logical, allocatable :: is_ghost(:)
    PetscInt, pointer :: idx(:)
    integer, allocatable :: nodes_per_proc(:), nodes_per_proc_local(:)
    integer :: new_id(1)
    IS :: global_numbering, node_distribution
    AO :: ordering
    PetscViewer :: pv

    bb(1) = 0.185
    bb(2) = -1.45
    bb(3) = 1.785
    bb(4) = 1.45
!    call create_rectangular_mesh(21,21,nplanes, &
!         num_global_nodes,num_global_elements, &
!         global_node, global_elm)
    call load_sms_mesh(nplanes, num_global_nodes, num_global_elements, &
         global_node, global_elm)

    call MPI_comm_size(PETSC_COMM_WORLD,size,ierr)
    call MPI_comm_rank(PETSC_COMM_WORLD,rank,ierr)

    ! Partition nodes
    ! ~~~~~~~~~~~~~~~
    if(nplanes.gt.1) then
       nodes_per_plane = num_global_nodes/nplanes
       myplane = (nplanes/size)*rank + 1
       min_node = (nplanes/size)*rank * nodes_per_plane + 1
       if(rank.eq.size-1) then
          max_node = num_global_nodes
       else
          max_node = min_node + (nplanes/size) * nodes_per_plane - 1
       endif
       num_local_nodes = max_node - min_node + 1
       print *, 'nodes per plane = ', nodes_per_plane
       print *, 'rank, myplane, min_node, max_node ', &
            rank, myplane, min_node, max_node
    else
       min_node = (num_global_nodes/size) * rank + 1
       if(rank.eq.size-1) then
          max_node = num_global_nodes
       else
          max_node = min_node + num_global_nodes/size - 1
       endif
       num_local_nodes = max_node - min_node + 1

       ! Create Connectivity Matrix
       call MatCreate(PETSC_COMM_WORLD, connectivity, ierr)
       call MatSetSizes(connectivity, num_local_nodes, num_local_nodes, &
            num_global_nodes, num_global_nodes, ierr)
       call MatSetFromOptions(connectivity, ierr)

       print *, 'populating connectivity matrix'
       do i=1, num_global_elements
          do j=1, nodes_per_element
             if(global_elm(i)%inode(j).lt.min_node) cycle
             if(global_elm(i)%inode(j).gt.max_node) cycle
             do k=1, nodes_per_element
                if(j.eq.k) cycle
                call MatSetValue(connectivity, global_elm(i)%inode(j)-1, &
                     global_elm(i)%inode(k)-1, 1., INSERT_VALUES, ierr)
             end do
          end do
       end do
       
       call MatAssemblyBegin(connectivity, MAT_FINAL_ASSEMBLY, ierr)
       call MatAssemblyEnd(connectivity, MAT_FINAL_ASSEMBLY, ierr)
       
       call MatConvert(connectivity, MATMPIADJ, MAT_INITIAL_MATRIX, &
            adjacency, ierr)
       call MatDestroy(connectivity, ierr)
       
       call MatPartitioningCreate(PETSC_COMM_WORLD, partitioning, ierr)
       call MatPartitioningSetAdjacency(partitioning, adjacency, ierr)
       call MatPartitioningSetFromOptions(partitioning, ierr)
       call MatPartitioningApply(partitioning, node_distribution, ierr)
       call MatPartitioningDestroy(partitioning, ierr)
       call MatDestroy(adjacency, ierr)
       call ISPartitioningToNumbering(node_distribution,global_numbering,ierr)
       
       call AOCreateBasicIS(global_numbering, PETSC_NULL, ordering, ierr)
       
       call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'node_distribution',pv,ierr)
       call ISView(node_distribution, pv, ierr)
       call PetscViewerDestroy(pv, ierr)
       
       call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'global_numbering',pv,ierr)
       call ISView(global_numbering, pv, ierr)
       call PetscViewerDestroy(pv, ierr)
       
       call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'ordering', pv, ierr)
       call AOView(ordering, pv, ierr)
       call PetscViewerDestroy(pv, ierr)

       ! Determine the number of nodes per process
       ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       allocate(nodes_per_proc_local(size), nodes_per_proc(size))
       nodes_per_proc_local = 0

       call ISGetIndicesF90(node_distribution, idx, ierr)
       do i=1, num_local_nodes
          nodes_per_proc_local(idx(i)+1) = nodes_per_proc_local(idx(i)+1) + 1
       end do
       call ISRestoreIndicesF90(node_distribution, idx, ierr)
       
       call mpi_allreduce(nodes_per_proc_local, nodes_per_proc, size, &
            MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

       min_node = 1
       do i=1, rank
          min_node = min_node + nodes_per_proc(i)
       end do
       max_node = min_node + nodes_per_proc(i) - 1
       num_local_nodes = nodes_per_proc(rank+1)
       
       deallocate(nodes_per_proc_local, nodes_per_proc)

       ! Re-number element node data
       ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
       do itri=1, num_global_elements
          global_elm(itri)%inode = global_elm(itri)%inode - 1
          call AOPetscToApplication(ordering, nodes_per_element, &
               global_elm(itri)%inode, ierr)
          global_elm(itri)%inode = global_elm(itri)%inode + 1
       end do
    endif

    print *, 'rank, min_node, max_node, local_nodes'
    print *, rank, min_node, max_node, num_local_nodes


    ! Determine local elements
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! take ownership of an element if I own inode(1) of the element
    
    ! calculate number of owned elements
    num_local_elements = 0
    do itri=1, num_global_elements
       if((global_elm(itri)%inode(1) .ge. min_node) .and. &
          (global_elm(itri)%inode(1) .le. max_node)) then
          num_local_elements = num_local_elements + 1
       end if
    end do

    ! create local element array
    allocate(local_elm(num_local_elements))
    k = 0
    do itri=1, num_global_elements
       if((global_elm(itri)%inode(1) .ge. min_node) .and. &
          (global_elm(itri)%inode(1) .le. max_node)) then
          k = k + 1
          local_elm(k) = global_elm(itri)
       end if
    end do


    ! Determine ghost nodes
    ! ~~~~~~~~~~~~~~~~~~~~~
    ! for each owned element, nodes that I don't own are ghosts

    ! calculate number of ghosts
    allocate(is_ghost(num_global_nodes))
    is_ghost = .false.
    num_ghost_nodes = 0
    do itri=1, num_local_elements
       do i=2, nodes_per_element
          if((local_elm(itri)%inode(i) .lt. min_node) .or. &
             (local_elm(itri)%inode(i) .gt. max_node)) then
             if(.not.is_ghost(local_elm(itri)%inode(i))) then
                is_ghost(local_elm(itri)%inode(i)) = .true.
                num_ghost_nodes = num_ghost_nodes + 1
             endif
          endif
       end do
    end do

    ! populate ghost list
    allocate(ghost_node(num_ghost_nodes))
    k = 1
    do i=1, num_global_nodes 
       if(is_ghost(i)) then
          ghost_node(k) = i
          k = k + 1
       end if
    end do
    deallocate(is_ghost)


    ! Create local node array
    ! ~~~~~~~~~~~~~~~~~~~~~~~
    total_local_nodes = num_local_nodes + num_ghost_nodes
    allocate(local_node(total_local_nodes),local_id(num_global_nodes))
    allocate(global_id(total_local_nodes))
    k = 1
    do i=1, num_local_nodes
       global_id(k) = i + min_node - 1
       local_id(global_id(k)) = k
       new_id(1) = global_id(k)-1
       if(nplanes.eq.1) then
          ! the 'global_node' array uses the original ordering
          call AOApplicationToPetsc(ordering, 1, new_id, ierr)
       end if
       local_node(k) = global_node(new_id(1)+1)
       k = k + 1
    end do
    do i=1, num_ghost_nodes
       global_id(k) = ghost_node(i)
       local_id(ghost_node(i)) = k
       new_id(1) = ghost_node(i)-1
       if(nplanes.eq.1) then
          ! the 'global_node' array uses the original ordering
          call AOApplicationToPetsc(ordering, 1, new_id, ierr)
       end if
       local_node(k) = global_node(new_id(1)+1)
       k = k + 1
    end do


    ! Re-number element nodes
    ! ~~~~~~~~~~~~~~~~~~~~~~~
    do itri=1, num_local_elements
       do i=1, nodes_per_element
          local_elm(itri)%inode(i) = local_id(local_elm(itri)%inode(i))
       end do
    end do
    
    deallocate(local_id)

    deallocate(global_node)
    deallocate(global_elm)

    if(nplanes.eq.1) then
       call AODestroy(ordering, ierr)
    end if

    print *, 'global nodes, elements', num_global_nodes, num_global_elements
    print *, 'local nodes, elements', num_local_nodes, num_local_elements    
    print *, 'ghost nodes', num_ghost_nodes

  end subroutine load_mesh

  subroutine unload_mesh
    implicit none

    if(allocated(ghost_node)) deallocate(ghost_node)
    if(allocated(local_node)) deallocate(local_node)
    if(allocated(local_elm)) deallocate(local_elm)
    if(allocated(global_id)) deallocate(global_id)
  end subroutine unload_mesh

  !======================================================================
  ! local_plane
  ! ~~~~~~~~~~~
  ! returns toroidal plane associated with current process
  !======================================================================
  integer function local_plane()
    implicit none
    local_plane = myplane
  end function local_plane

  !==============================================================
  ! global_node_id
  ! ~~~~~~~~~~~~~~
  ! returns the global node id of local node inode
  !==============================================================
  integer function global_node_id(inode)
    implicit none
    integer, intent(in) :: inode
    global_node_id = global_id(inode)
  end function global_node_id


  !==============================================================
  ! local_nodes
  ! ~~~~~~~~~~~
  ! returns the number of elements local to this process
  !==============================================================
  integer function local_elements()
    implicit none

    local_elements = num_local_elements
  end function local_elements


  !==============================================================
  ! local_nodes
  ! ~~~~~~~~~~~
  ! returns the number of nodes local to this process
  !==============================================================
  integer function local_nodes()
    implicit none

    local_nodes = total_local_nodes
  end function local_nodes

  !==============================================================
  ! local_dofs
  ! ~~~~~~~~~~
  ! returns the number of dofs local to this process
  !==============================================================
  integer function local_dofs()
    implicit none

    local_dofs = total_local_nodes*dofs_per_node
  end function local_dofs

  !==============================================================
  ! owned_nodes
  ! ~~~~~~~~~~~
  ! returns the number of dofs local to this process
  !==============================================================
  integer function owned_nodes()
    implicit none

    owned_nodes = num_local_nodes
  end function owned_nodes

  !==============================================================
  ! owned_dofs
  ! ~~~~~~~~~~
  ! returns the number of dofs local to this process
  !==============================================================
  integer function owned_dofs()
    implicit none

    owned_dofs = num_local_nodes*dofs_per_node
  end function owned_dofs
  

  !==============================================================
  ! global_dofs
  ! ~~~~~~~~~~~
  ! returns the global number of dofs
  !==============================================================
  integer function global_dofs()
    implicit none

    global_dofs = num_global_nodes*dofs_per_node
  end function global_dofs


  !==============================================================
  ! get_element_nodes
  ! ~~~~~~~~~~~~~~~~~
  ! returns the indices of each element node
  !==============================================================
  subroutine get_element_nodes(itri,n)
    implicit none
    integer, intent(in) :: itri
    integer, intent(out), dimension(nodes_per_element) :: n

    n = local_elm(itri)%inode
  end subroutine get_element_nodes

  subroutine get_ghost_nodes(ids, num)
    implicit none
    integer, intent(out), pointer :: ids(:)
    integer, intent(out) :: num

    num = num_ghost_nodes
    ids => ghost_node
  end subroutine get_ghost_nodes


  !==============================================================
  ! get_bounding_box
  ! ~~~~~~~~~~~~~~~~
  ! returns the lower-left and upper-right coordinates of the 
  ! mesh bounding box
  !==============================================================
  subroutine get_bounding_box(x1, z1, x2, z2)
    implicit none
    real, intent(out) :: x1, z1, x2, z2
    
    x1 = bb(1)
    z1 = bb(2)
    x2 = bb(3)
    z2 = bb(4) 
  end subroutine get_bounding_box


  !==============================================================
  ! get_bounding_box
  ! ~~~~~~~~~~~~~~~~
  ! returns the width and heightof the
  ! mesh bounding box
  !==============================================================
  subroutine get_bounding_box_size(alx, alz)
    implicit none

    real, intent(out) :: alx, alz
    real :: x1, z1, x2, z2
    call get_bounding_box(x1,z1,x2,z2)
    alx = x2 - x1
    alz = z2 - z1
  end subroutine get_bounding_box_size


  !============================================================
  ! get_node_pos
  ! ~~~~~~~~~~~~
  ! get the global coordinates (R,Phi,Z) of node inode
  !============================================================
  subroutine get_node_pos(inode,R,Phi,Z)
    implicit none
    
    integer, intent(in) :: inode
    real, intent(out) :: R, Phi, Z
    
    R = local_node(inode)%R
    Phi = local_node(inode)%Phi
    Z = local_node(inode)%Z
    
  end subroutine get_node_pos


  !============================================================
  ! get_element_data
  ! ~~~~~~~~~~~~~~~~
  ! calculates element parameters a, b, c, and theta from
  ! node coordinates
  !============================================================
  subroutine get_element_data(itri, d)
    use math
    implicit none
    integer, intent(in) :: itri
    type(element_data), intent(out) :: d
    real :: x2, x3, z2, z3, x2p, x3p, z2p, z3p, hi, phi2, phi3
    integer :: nodeids(nodes_per_element)
    
    d%itri = itri

    call get_element_nodes(itri, nodeids)
    
    call get_node_pos(nodeids(1), d%R, d%Phi, d%Z)
    call get_node_pos(nodeids(2), x2, phi2, z2)
    call get_node_pos(nodeids(3), x3, phi3, z3)

    x2 = x2 - d%R
    z2 = z2 - d%Z
    x3 = x3 - d%R
    z3 = z3 - d%Z
    
    hi = 1./sqrt(x2**2 + z2**2)
    d%co = x2*hi
    d%sn = z2*hi

    x2p =  d%co * x2 + d%sn * z2
    z2p = -d%sn * x2 + d%co * z2
    x3p =  d%co * x3 + d%sn * z3
    z3p = -d%sn * x3 + d%co * z3

    d%a = x2p-x3p
    d%b = x3p
    d%c = z3p

#ifdef USE3D
    call get_node_pos(nodeids(4), x2, phi2, z2)
    d%d = phi2 - d%Phi
    if(d%d .le. 0.) d%d = d%d + twopi
#else
    d%d = 0.
#endif
    
  end subroutine get_element_data


  !============================================================
  ! whattri
  ! ~~~~~~~
  ! Gets the element itri in which global coordinates (x,z)
  ! fall.  If no element on this process contains (x,z) then
  ! itri = -1.  Otherwise, xref and zref are global coordinates
  ! of the first node of element itri.
  !============================================================
  subroutine whattri(x,phi,z,ielm,xref,zref,nophi)
    implicit none
    
    real, intent(in) :: x, phi, z
    integer, intent(inout) :: ielm
    real, intent(out) :: xref, zref
    logical, intent(in), optional :: nophi
    
    integer :: nelms, i
    type(element_data) :: d

    ! first, try ielm
    if(ielm.gt.0) then
       call get_element_data(ielm,d)
       if(is_in_element(d,x,phi,z,nophi)) then
          xref = d%R
          zref = d%Z
          return
       endif
    endif

    ! try all elements
    ielm = -1
    nelms = local_elements()
    do i=1, nelms
       call get_element_data(i,d)
       if(is_in_element(d,x,phi,z,nophi)) then
          xref = d%R
          zref = d%Z
          ielm = i
          return
       end if
    end do
  end subroutine whattri

  
  !=========================================
  ! is_boundary_node
  ! ~~~~~~~~~~~~~~~~
  ! returns true if node lies on boundary
  !=========================================
  logical function is_boundary_node(inode)
    implicit none
    integer, intent(in) :: inode
    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv(3), x, phi,  z

    call boundary_node(inode, is_boundary, izone, izonedim, normal, curv, &
         x, phi, z)
    is_boundary_node = is_boundary
  end function is_boundary_node

  !======================================================================
  ! boundary_node
  ! ~~~~~~~~~~~~~
  ! determines if node is on boundary, and returns relevant info
  ! about boundary surface
  !======================================================================
  subroutine boundary_node(inode,is_boundary,izone,izonedim,normal,curv, &
       x,phi,z)

    implicit none
    
    integer, intent(in) :: inode              ! node index
    integer, intent(out) :: izone,izonedim    ! zone type/dimension
    real, intent(out) :: normal(2), curv(3)
    real, intent(out) :: x,z,phi              ! coordinates of inode
    logical, intent(out) :: is_boundary       ! is inode on boundary

    is_boundary = local_node(inode)%idim.lt.2

    if(is_boundary) then
       izone = local_node(inode)%izone
       izonedim = local_node(inode)%idim
       normal = local_node(inode)%normal
       curv = local_node(inode)%curv
       x = local_node(inode)%R
       phi = local_node(inode)%Phi
       z = local_node(inode)%Z
    endif
  end subroutine boundary_node


  subroutine boundary_edge(itri, is_edge, normal, idim)

    implicit none
    integer, intent(in) :: itri
    integer, intent(out) :: is_edge(3)
    real, intent(out) :: normal(2,3)
    integer, intent(out) :: idim(3)
    
    integer :: inode(nodes_per_element), izone, i, j
    real :: x, phi, z, c(3)
    logical :: is_bound(3)
    
    call get_element_nodes(itri,inode)
    
    do i=1,3
       call boundary_node(inode(i),is_bound(i),izone,idim(i), &
            normal(:,i),c(i),x,phi,z)
    end do
    
    do i=1,3
       j = mod(i,3) + 1
       is_edge(i) = 0
       
       ! skip edges not having both points on a boundary
       if((.not.is_bound(i)).or.(.not.is_bound(j))) cycle
       
       ! skip edges cutting across corners
       if(is_bound(1) .and. is_bound(2) .and. is_bound(3)) then
          if(idim(i).ne.0 .and. idim(j).ne.0) cycle
       endif
       
       ! skip suspicious edges (edges w/o corner point where normal changes
       ! dramatically)
       if(idim(i).eq.1 .and. idim(j).eq.1 .and. idim(mod(i+1,3)+1).eq.2) then
          if(normal(1,i)*normal(1,j) + normal(2,i)*normal(2,j) .lt. .5) cycle
       end if
       
       is_edge(i) = 1
    end do
  end subroutine boundary_edge
  
end module basic_mesh_mod
