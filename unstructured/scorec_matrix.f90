module scorec_matrix_mod

  use element

  implicit none

  type scorec_matrix
     integer :: imatrix
     integer :: icomplex
     integer :: isize, m, n
     integer :: lhs
  end type scorec_matrix
  integer, parameter :: MAT_SET = 0
  integer, parameter :: MAT_ADD = 1

  integer, parameter :: maxnumofsolves = 4
  real, allocatable:: kspits(:)

#ifdef REORDERED
  integer :: is1_agg_blk_cnt = 1
  integer :: is1_agg_scp = 0
#endif

  interface clear_mat
     module procedure scorec_matrix_clear
  end interface

  interface create_mat
     module procedure scorec_matrix_create
  end interface

  interface destroy_mat
     module procedure scorec_matrix_destroy
  end interface

  interface finalize
     module procedure scorec_matrix_finalize
  end interface

  interface zero_mat
     module procedure scorec_matrix_zero
  end interface

  interface update_mat
     module procedure scorec_matrix_update
  end interface

  interface flush
     module procedure scorec_matrix_flush
  end interface

  interface get_global_node_indices
     module procedure scorec_matrix_get_global_node_indices
  end interface

  interface get_element_indices
     module procedure scorec_matrix_get_element_indices
  end interface

  interface get_node_indices
     module procedure scorec_matrix_get_node_indices
  end interface

  interface insert
     module procedure scorec_matrix_insert_real
#ifdef USECOMPLEX
     module procedure scorec_matrix_insert_complex
#endif
  end interface

  interface insert_block
     module procedure scorec_matrix_insert_block
  end interface

  interface insert_global
     module procedure scorec_matrix_insert_global_real
#ifdef USECOMPLEX
     module procedure scorec_matrix_insert_global_complex
#endif
  end interface

  interface matvecmult
     module procedure scorec_matrix_matvecmult
  end interface

  interface newsolve
     module procedure scorec_matrix_solve
  end interface

  interface newsolve_with_guess
     module procedure scorec_matrix_solve_with_guess
  end interface

  interface set_matrix_index
     module procedure scorec_matrix_set_index
  end interface

  interface write_matrix
     module procedure scorec_matrix_write
  end interface

  interface allocate_kspits
     module procedure scorec_allocate_kspits
  end interface

contains


  !====================================================================
  ! -memory_view - Print memory usage at end of run
  ! -log_view_memory - Display memory information for each logged event
  ! -malloc_view - Print usage of PetscMalloc() in PetscFinalize()
  ! PetscMemorySetGetMaximumUsage

  ! PETSc provides a number of tools to aid in detection of problems with memory allocation, 
  ! including leaks and use of uninitialized space. We briefly describe these below.

  !• The PETSc memory allocation (which collects statistics and performs error checking), 
  !  is employed by default for codes compiled in a debug-mode (configured with --with-debugging=1). 
  !  PETSc memory allocation can be activated for optimized-mode (configured with --with-debugging=0)
  !  using the option -malloc_debug. The option -malloc_debug 0 forces the use of conventional memory
  !  allocation when debugging is enabled. When running timing tests, one should build libraries in optimized mode.

  !• When the PETSc memory allocation routines are used, the option -malloc_dump will print a list of
  !  unfreed memory at the conclusion of a program. If all memory has been freed, only a message stating the
  !  maximum allocated space will be printed. However, if some memory remains unfreed, this information
  !  will be printed. Note that the option -malloc_dump activates a call to PetscMallocDump() during
  !  PetscFinalize() the user can also call PetscMallocDump() elsewhere in a program.

  !• Another useful option for use with PETSc memory allocation routines is -malloc_view, 
  !  which activates logging of all calls to malloc and reports memory usage, including all Fortran arrays. 
  !  This option provides a more complete picture than -malloc_dump for codes that employ Fortran with hardwired arrays. 
  !  The option -malloc_view activates logging by calling PetscMallocViewSet() in PetscInitialize() and then 
  !  prints the log by calling PetscMallocView() in PetscFinalize(). 
  !  The user can also call these routines elsewhere in a program. When finer granularity is desired,
  !  the user should call PetscMallocGetCurrentUsage() and PetscMallocGetMaximumUsage() for memory 
  !  allocated by PETSc, or PetscMemoryGetCurrentUsage() and PetscMemoryGetMaximumUsage() for the total memory 
  !  used by the program. Note that PetscMemorySetGetMaximumUsage() must be called before 
  !  PetscMemoryGetMaximumUsage() (typically at the beginning of the program).

  !• When running with -log_view the additional option -log_view_memory causes the display of additional
  !  columns of information about how much memory was allocated and freed during each logged event.
  !  This is useful to understand what phases of a computation require the most memory.
  !====================================================================
  subroutine printMemStatTotalStart
    integer :: ierr
    call PetscMemorySetGetMaximumUsage(ierr);
  end subroutine printMemStatTotalStart
  subroutine printMemStatTotal
    !call cc_printMemStatTotal
    use petsc
    implicit none
    real :: mem, mem_max
    integer :: myrank, ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call PetscMemoryGetCurrentUsage(mem,ierr)
    call PetscMemoryGetMaximumUsage(mem_max,ierr);
    if(mod(myrank,10)==0) &
    write(*,'(A,2x,i4,2x,F14.2,2x, F14.2)') "Memory usage (MB) reported by Petsc:",myrank,mem/1e6,mem_max/1e6
  end subroutine printMemStatTotal

  !====================================================================
  ! set_matrix_index
  ! ~~~~~~~~~~~~~~~~
  ! sets the index of a scorec matrix
  ! this must be called before the matrix is created
  !====================================================================
  subroutine scorec_matrix_set_index(mat, imat)
    implicit none
    type(scorec_matrix) :: mat
    integer, intent(in) :: imat

    mat%imatrix = imat
  end subroutine scorec_matrix_set_index

  !====================================================================
  ! create
  ! ~~~~~~
  ! creates a scorec solve matrix with index imat and size isize
  !====================================================================
#ifdef REORDERED
  subroutine scorec_matrix_create(mat, m, n, icomplex, lhs, agg_blk_cnt_opt, agg_scp_opt)
#else
  subroutine scorec_matrix_create(mat, m, n, icomplex, lhs)
#endif

#ifndef M3DC1_TRILINOS

#if PETSC_VERSION >= 38
    use petsc
    implicit none
#elif PETSC_VERSION >= 36
    implicit none
#include "petsc/finclude/petsc.h"
#else
    implicit none
#include "finclude/petsc.h"
#endif

#endif

    type(scorec_matrix) :: mat
    integer, intent(in) :: m, n, icomplex
    integer, intent(in) :: lhs

#ifdef REORDERED
    integer, intent(in), optional :: agg_blk_cnt_opt
    integer, intent(in), optional :: agg_scp_opt

    integer :: agg_blk_cnt
    integer :: agg_scp

    if(.not. present(agg_blk_cnt_opt)) then
       agg_blk_cnt = 1
    else
       agg_blk_cnt = agg_blk_cnt_opt
    end if

    if(.not. present(agg_scp_opt)) then
       agg_scp = 0
    else
       agg_scp = agg_scp_opt
    end if
#endif

    if(mat%imatrix .le. 0) then
       print*, 'Error: scorec matrix index not set!'
       return
    endif

    mat%m = m
    mat%n = n
    mat%isize = max(m, n)
    mat%icomplex = icomplex
    mat%lhs = lhs
!#ifndef M3DC1_TRILINOS
!    if(lhs .eq. 1) then
!       call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc',flg_petsc,ierr)
!       call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ierr)
!       if(flg_solve2.eq.PETSC_TRUE) flg_petsc=PETSC_TRUE
!    endif
!#endif
#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_create(mat%imatrix, lhs, mat%icomplex, mat%isize)
#else
#ifdef REORDERED
    call m3dc1_matrix_create(mat%imatrix, lhs, mat%icomplex, mat%isize, agg_blk_cnt, agg_scp)
#else
    call m3dc1_matrix_create(mat%imatrix, lhs, mat%icomplex, mat%isize)
#endif !reordered
#endif
  end subroutine scorec_matrix_create


  !====================================================================
  ! clear
  ! ~~~~~
  ! zeroes out a scorec solve matrix
  !====================================================================
  subroutine scorec_matrix_clear(mat)
!    use basic, ONLY: myrank
    implicit none

    type(scorec_matrix), intent(inout) :: mat
#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_delete (mat%imatrix)
#else
    call m3dc1_matrix_delete (mat%imatrix)
#endif
    call create_mat(mat,mat%m,mat%n,mat%icomplex,mat%lhs)
!    if(myrank==0) &
!    write(*,'(A,2x,i4)') "matrix deleted and re-created:",mat%imatrix
  end subroutine scorec_matrix_clear


  !====================================================================
  ! destroy
  ! ~~~~~~~
  ! destroys and scorec solve matrix
  !====================================================================
  subroutine scorec_matrix_destroy(mat)
!    use basic, ONLY: myrank
    implicit none
    
    type(scorec_matrix) :: mat
#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_delete (mat%imatrix)
#else
    call m3dc1_matrix_delete (mat%imatrix)
#endif
!    if(myrank==0) &
!    write(*,'(A,2x,i4)') "matrix deleted:",mat%imatrix
  end subroutine scorec_matrix_destroy


  !====================================================================
  ! matvecmult
  ! ~~~~~~~~~~
  ! matrix vector multiply with scorec data structures
  !====================================================================
  subroutine scorec_matrix_matvecmult(mat,vin,vout)

    use vector_mod

    implicit none

    type(scorec_matrix), intent(in) :: mat
    type(vector_type), intent(in), target :: vin
    type(vector_type), intent(inout), target :: vout

    type(vector_type), pointer :: vin_ptr, vout_ptr
    type(vector_type), target :: temp_in, temp_out

    if(mat%m .ne. vout%isize) then
       print *, 'Error: col sizes for matvecmult do not conform', &
            mat%m, vout%isize
       return
    endif
    if(mat%n .ne. vin%isize) then
       print *, 'Error: row sizes for matvecmult do not conform', &
            mat%n, vin%isize
       return
    endif

    ! scorec matvecmult only works for square matrices
    if(vin%isize .eq. mat%isize) then
       vin_ptr => vin
    else
       call create_vector(temp_in, mat%isize)
       temp_in = vin
       vin_ptr => temp_in
    endif
    if(vout%isize .eq. mat%isize) then
       vout_ptr => vout
    else
       call create_vector(temp_out, mat%isize)
       vout_ptr => temp_out
    endif

#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_multiply(mat%imatrix,vin_ptr%id,vout_ptr%id)
#else
    call m3dc1_matrix_multiply(mat%imatrix,vin_ptr%id,vout_ptr%id)
#endif
    nullify(vin_ptr, vout_ptr)

    if(vin%isize .ne. mat%isize) call destroy_vector(temp_in)
    if(vout%isize .ne. mat%isize) then
       vout = temp_out
       call destroy_vector(temp_out)
    endif
     
  end subroutine scorec_matrix_matvecmult


  !====================================================================
  ! insert_real
  ! ~~~~~~~~~~~
  ! inserts (or increments) a matrix element with real value
  !====================================================================
  subroutine scorec_matrix_insert_real(mat,val,i,j,iop)
    implicit none

    type(scorec_matrix), intent(in) :: mat
    real, intent(in) :: val
    integer, intent(in) :: i, j, iop

    print *, "not implement scorec_matrix_insert_real"
    call abort() 
  end subroutine scorec_matrix_insert_real


#ifdef USECOMPLEX
  !====================================================================
  ! insert_complex
  ! ~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with complex value
  !====================================================================
  subroutine scorec_matrix_insert_complex(mat,val,i,j,iop)
    implicit none

    type(scorec_matrix), intent(in) :: mat
    complex, intent(in) :: val
    integer, intent(in) :: i, j, iop

    if(val.eq.0.) then
       if(i.ne.j) return
    endif
    print *, "not implement scorec_matrix_insert_complex"
    call abort()
  end subroutine scorec_matrix_insert_complex
#endif

  !====================================================================
  ! insert_global_real
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with real value
  ! using global indices
  !====================================================================
  subroutine scorec_matrix_insert_global_real(mat,val,i,j,iop)

    implicit none

    type(scorec_matrix), intent(in) :: mat
    real, intent(in) :: val
    integer, intent(in) :: i, j, iop
#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_insert(mat%imatrix, val, i-1, j-1, 0,val, iop)
#else
    call m3dc1_matrix_insert(mat%imatrix, val, i-1, j-1, 0,val, iop)
#endif
  end subroutine scorec_matrix_insert_global_real


#ifdef USECOMPLEX
  !====================================================================
  ! insert_global_complex
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with complex value
  ! using local indices
  !====================================================================
  subroutine scorec_matrix_insert_global_complex(mat,val,i,j,iop)

    implicit none

    type(scorec_matrix), intent(in) :: mat
    complex, intent(in) :: val
    integer, intent(in) :: i, j, iop
#ifdef M3DC1_TRILINOS
    call m3dc1_epertra_insert(mat%imatrix, val, i-1, j-1, 1,val, iop)
#else
    call m3dc1_matrix_insert(mat%imatrix, val, i-1, j-1, 1,val, iop)
#endif
  end subroutine scorec_matrix_insert_global_complex
#endif


  !====================================================================
  ! solve
  ! ~~~~~
  ! linear matrix solve with scorec data structures
  !====================================================================
  subroutine scorec_matrix_solve(mat, v, ierr)
    use vector_mod

#if PETSC_VERSION >= 38
    use petsc
    implicit none
#elif PETSC_VERSION >= 36
    implicit none
#include "petsc/finclude/petsc.h"
#else
    implicit none
#include "finclude/petsc.h"
#endif
    
    type(scorec_matrix), intent(in) :: mat
    type(vector_type), intent(inout) :: v
    integer, intent(out) :: ierr
    real :: t1,t2
    integer :: myrank

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    t1=MPI_WTIME() 
#ifdef M3DC1_TRILINOS
    call m3dc1_solver_aztec(mat%imatrix,v%id,v%id,num_iter,&
         solver_tol,krylov_solver,preconditioner,sub_dom_solver,&
         subdomain_overlap,graph_fill,ilu_drop_tol,ilu_fill,&
         ilu_omega,poly_ord)    
    !! call m3dc1_solver_aztec(mat%imatrix,v%id,v%id,num_iter,&
    !!     solver_tol,aztec_options)    

#else
    call m3dc1_matrix_solve(mat%imatrix,v%id,solver_type,solver_tol)
#endif
    ierr = 0
    t2=MPI_WTIME()

#ifdef M3DC1_TRILINOS
    call m3dc1_solver_getnumiter(mat%imatrix,num_iter)
#else
       !2013-jan-17 only for hopper
    call m3dc1_matrix_getiternum(mat%imatrix,num_iter)
!    if(mat%imatrix==5) call m3dc1_matrix_getinfo(mat%imatrix,v%id)
#endif
   if(myrank==0) write(*,'(A,I4,A,I5,A,ES11.3,A)') &
           "For mat%imatrix=",mat%imatrix,", newsolve takes num_iter=",num_iter,", time=", t2-t1, " secs"
   if (num_iter.ge.10000) then
      print*, 'Error: solver', mat%imatrix, 'diverged. JOB stopped'
      call safestop(10000)
   endif
   if(mat%imatrix== 5) kspits(1)=num_iter
   if(mat%imatrix== 1) kspits(2)=num_iter
   if(mat%imatrix==17) kspits(3)=num_iter
   if(mat%imatrix== 6) kspits(4)=num_iter
  end subroutine scorec_matrix_solve

  subroutine scorec_matrix_solve_with_guess(mat, v, x_guess, ierr)
   use vector_mod
#if PETSC_VERSION >= 38
   use petsc
   implicit none
#elif PETSC_VERSION >= 36
   implicit none
#include "petsc/finclude/petsc.h"
#else
   implicit none
#include "finclude/petsc.h"
#endif
   
   type(scorec_matrix), intent(in) :: mat
   type(vector_type), intent(inout) :: v, x_guess
   integer, intent(out) :: ierr
   real :: t1,t2
   integer :: myrank

   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
   t1=MPI_WTIME()
#ifdef M3DC1_TRILINOS
   call m3dc1_solver_aztec(mat%imatrix,v%id,v%id,num_iter,&
        solver_tol,krylov_solver,preconditioner,sub_dom_solver,&
        subdomain_overlap,graph_fill,ilu_drop_tol,ilu_fill,&
        ilu_omega,poly_ord)    
  
#else
   call m3dc1_matrix_solve_with_guess(mat%imatrix,v%id, x_guess%id, solver_type,solver_tol)
#endif
   ierr = 0
   t2=MPI_WTIME()

#ifdef M3DC1_TRILINOS
   call m3dc1_solver_getnumiter(mat%imatrix,num_iter)
#else
   call m3dc1_matrix_getiternum(mat%imatrix,num_iter)
#endif

  if(myrank==0) write(*,'(A,I4,A,I5,A,ES11.3,A)') &
          "For mat%imatrix=",mat%imatrix,", newsolve_with_guess takes num_iter=",num_iter,", time=", t2-t1, " secs"
  if (num_iter.ge.10000) then
     print*, 'Error: solver', mat%imatrix, 'diverged. JOB stopped'
     call safestop(10000)
  endif
  if(mat%imatrix== 5) kspits(1)=num_iter
  if(mat%imatrix== 1) kspits(2)=num_iter
  if(mat%imatrix==17) kspits(3)=num_iter
  if(mat%imatrix== 6) kspits(4)=num_iter
 end subroutine scorec_matrix_solve_with_guess


  !====================================================================
  ! finalize
  ! ~~~~~~~~
  ! finalizes matrix 
  !====================================================================
  subroutine scorec_matrix_finalize(mat)
    implicit none
    type(scorec_matrix) :: mat   
#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_freeze(mat%imatrix)
#else
    call m3dc1_matrix_freeze(mat%imatrix)
#endif
  end subroutine scorec_matrix_finalize


  !====================================================================
  ! zero: mat=0, to be update at ntime+1
  ! ~~~~~~~~
  ! zero matrix 
  !====================================================================
  subroutine scorec_matrix_zero(mat)
    implicit none
    type(scorec_matrix) :: mat   
     call m3dc1_matrix_reset(mat%imatrix)
  end subroutine scorec_matrix_zero


  !====================================================================
  ! update: mat=0, to be update at ntime+1
  ! ~~~~~~~~
  ! update solver matrix and preconditioner matrix
  !====================================================================
  subroutine scorec_matrix_update(mat)
    implicit none
    type(scorec_matrix) :: mat
     call m3dc1_matrix_update(mat%imatrix)
  end subroutine scorec_matrix_update

  !====================================================================
  ! flush
  ! ~~~~~
  ! flushes matrix 
  !====================================================================
  subroutine scorec_matrix_flush(mat)
    implicit none
    type(scorec_matrix) :: mat
  end subroutine scorec_matrix_flush

  !======================================================================
  ! get_node_indices
  ! ~~~~~~~~~~~~~~~~
  ! given a matrix mat and node inode, returns:
  ! irow(i,j): the local row index associated with dof j associated 
  !           with field i
  ! icol(i,j): the local column index associated with dof j associated 
  !           with field i
  !======================================================================
  subroutine scorec_matrix_get_node_indices(mat, inode, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: inode
    integer, intent(out), dimension(mat%m,dofs_per_node) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_node) :: icol

    integer, dimension(mat%isize,dofs_per_node) :: ind
    integer :: i

    call get_node_indices(mat%isize,inode,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do    
  end subroutine scorec_matrix_get_node_indices

  !======================================================================
  ! get_global_node_indices
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  ! given a matrix mat and node inode, returns:
  ! irow(i,j): the global row index associated with dof j associated 
  !           with field i
  ! icol(i,j): the global column index associated with dof j associated 
  !           with field i
  !======================================================================
  subroutine scorec_matrix_get_global_node_indices(mat, inode, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: inode
    integer, intent(out), dimension(mat%m,dofs_per_node) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_node) :: icol

    integer, dimension(mat%isize,dofs_per_node) :: ind
    integer :: i

    call get_global_node_indices(mat%isize,inode,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do
  end subroutine scorec_matrix_get_global_node_indices

  !======================================================================
  ! get_element_indices
  ! ~~~~~~~~~~~~~~~~~~~
  ! given a matrix mat and element itri, returns:
  ! irow(i,j): the local row index associated with dof j associated 
  !           with field i
  ! icol(i,j): the local column index associated with dof j associated 
  !           with field i
  !======================================================================
  subroutine scorec_matrix_get_element_indices(mat, itri, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: itri
    integer, intent(out), dimension(mat%m,dofs_per_element) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_element) :: icol

    integer, dimension(mat%isize,dofs_per_element) :: ind
    integer :: i
    call get_element_indices(mat%isize,itri,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do
  end subroutine scorec_matrix_get_element_indices

  subroutine scorec_matrix_insert_block(mat, itri, m, n, val, iop)
    implicit none
    type(matrix_type) :: mat
    integer, intent(in) :: itri, m, n, iop
    vectype, intent(in), dimension(dofs_per_element,dofs_per_element) :: val

    integer, dimension(mat%m,dofs_per_element) :: irow
    integer, dimension(mat%n,dofs_per_element) :: icol

    if (iop .ne. MAT_ADD) then
       print *, " mat block insert only supports MAT_ADD"
       call abort()
    end if
    call get_element_indices(mat, itri, irow, icol)
#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_addblock(mat%imatrix, itri-1, m-1, n-1, TRANSPOSE(val));
#else
    call m3dc1_matrix_insertblock(mat%imatrix, itri-1, m-1, n-1, TRANSPOSE(val));
#endif
  end subroutine scorec_matrix_insert_block

  subroutine identity_row(mat, irow)
    implicit none
    type(scorec_matrix) :: mat
    integer, intent(in) :: irow

#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_setbc (mat%imatrix, irow-1)
#else
    call m3dc1_matrix_setbc (mat%imatrix, irow-1)
#endif    
  end subroutine identity_row

  subroutine set_row_vals(mat, irow, ncols, icols, vals)
    use vector_mod

    type(scorec_matrix) :: mat
    integer, intent(in) :: irow, ncols
    integer, intent(in), dimension(ncols) :: icols
    vectype, intent(in), dimension(ncols) :: vals
#ifdef M3DC1_TRILINOS
    call m3dc1_epetra_setlaplacebc(mat%imatrix, irow-1, ncols, icols-1, vals)
#else
    call m3dc1_matrix_setlaplacebc(mat%imatrix, irow-1, ncols, icols-1, vals)
#endif
  end subroutine set_row_vals

  subroutine scorec_matrix_write(mat, file)
    type(scorec_matrix), intent(in) :: mat
    character(len=*) :: file
    call m3dc1_matrix_write(mat%imatrix, file, 0);
  end subroutine scorec_matrix_write

  subroutine scorec_allocate_kspits
    implicit none
  if(.not.allocated(kspits)) allocate(kspits(1:maxnumofsolves))
  end subroutine scorec_allocate_kspits


#ifdef TODO
  !====================================================================
  ! m3dc1_matrix_getinfo
  ! ~~~~~~~~
  ! get matrix info
  !====================================================================
  subroutine m3dc1_matrix_getinfo(imatrix,ivec)
    use vector_mod
    use petsc
    implicit none

    integer :: ierr, myrank
    integer :: imatrix,ivec
    integer :: num_local_dof, num_global_dof
    integer :: num_local_dof_f, num_own_dof_f, num_global_dof_f, num_ghost_dof_f
    integer :: num_local_ent, num_own_ent, num_global_ent, num_ghost_ent
    integer :: num_local_elements, num_own_elements, num_global_elements, num_ghost_elements, tot_elms

    call m3dc1_matrix_getlocalnumdof(imatrix,num_local_dof)
    call m3dc1_matrix_getglobalnumdof(imatrix,num_global_dof)

    call m3dc1_field_getnumlocaldof(ivec, num_local_dof_f)
    call m3dc1_field_getnumowndof(ivec, num_own_dof_f)
    call m3dc1_field_getnumglobaldof(ivec, num_global_dof_f)
    call m3dc1_field_getnumghostdof(ivec, num_ghost_dof_f)

    call m3dc1_mesh_getnument(0, num_local_ent)
    call m3dc1_mesh_getnumownent(0, num_own_ent)
    call m3dc1_mesh_getnumglobalent(0, num_global_ent)
    call m3dc1_mesh_getnumghostent(0, num_ghost_ent)

    call m3dc1_mesh_getnument(3, num_local_elements)
    call m3dc1_mesh_getnumglobalent(3, num_global_elements)

    call mpi_allreduce(num_local_elements, tot_elms, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)  

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    if(myrank==0) then
       write(*,'(a,I8)') "CJ num_local_dof_mat5=", num_local_dof
       write(*,'(a,I8)') "CJ num_global_dof_mat5=", num_global_dof

       write(*,'(a,I8)') "CJ num_local_dof_f=", num_local_dof_f
       write(*,'(a,I8)') "CJ num_own_dof_f=", num_own_dof_f       !=num_local_dof
       write(*,'(a,I8)') "CJ num_global_dof_f=", num_global_dof_f !=num_global_dof

       write(*,'(a,I8)') "CJ num_local_ent=", num_local_ent
       write(*,'(a,I8)') "CJ num_global_ent=", num_global_ent

       write(*,'(a,I8)') "CJ num_local_elements=", num_local_elements
       write(*,'(a,I8,2x,I8)') "CJ num_global_elements=", num_global_elements, tot_elms
    endif
  end subroutine m3dc1_matrix_getinfo
#endif

end module scorec_matrix_mod
