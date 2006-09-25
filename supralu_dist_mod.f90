! Attempt at a distributed-memory parallel version of A. Pletzer's supralu
! module to interface with SuperLU_DIST v2.0. J. Breslau, 7/30/04

module supralu_dist_mod
  implicit none

  integer, parameter :: r8 = selected_real_kind(12,100)
  
!integer, parameter :: i8 = selected_int_kind(12)

  integer SLUprocsize  !Number of active MPI processes
  integer*8 SLUprocgrid  !Opaque handle to SuperLU process grid

  type sparseR8d_obj
     !Size specifiers
     integer base        !0=C-style; 1=Fortran-style
     integer n           !Linear dimension of matrix
     integer nnz         !Number of nominal non-zero entries
     integer nnzfilled   !Number of entries set so far
     integer colfilled   !Number of columns set so far

     !Shape & data holders
     integer, pointer :: rowind(:), colptr(:)
     real(r8), pointer :: values(:)

     !Opaque SuperLU handles
     integer*8 options          !SuperLU options settings for this matrix
     integer*8 ScalePermstruct  !Permutations & pivot scalings done to matrix
     integer*8 LUstruct         !Factorization of matrix (L & U)
     integer*8 SOLVEstruct      !Data related to parallel solution
     integer*8 SM               !The SuperLU distributed supermatrix itself
     integer*8 stat             !Solver statistics

     !Initialization flag
     logical initialized      !Has the rmatrix been created?
  end type sparseR8d_obj

contains
!============================================================
subroutine SLUD_init  !Initialize SuperLU process grid on global communicator
  implicit none
  include 'mpif.h'

  integer nprow, npcol, ierr

  !Divide matrix among processors by rows
  call MPI_Comm_size(MPI_COMM_WORLD, SLUprocsize, ierr)
      nprow = 1
!     if(SLUprocsize .gt. 1) nprow = 2
      npcol = SLUprocsize/nprow

  call f_create_gridinfo_handle(SLUprocgrid)
  call f_superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, SLUprocgrid)
end subroutine SLUD_init

!============================================================
subroutine SLUD_exit  !Finalize SuperLU process grid
  implicit none

  call f_superlu_gridexit(SLUprocgrid)
  call f_destroy_gridinfo_handle(SLUprocgrid)
end subroutine SLUD_exit

!============================================================
subroutine sparseR8d_init(this, n, nnz, base, ier)

  ! Constructor: sets up data structure/container for distributed
  !  column-compressed sparse matrix, allocates storage, but does
  !  not call any SuperLU routines
  !
  ! this: handle of the container
  ! n: linear dimension of the (square) matrix
  ! nnz: number of non-zero elements in the matrix
  ! base: indexing base (0 like C or 1 like Fortran)

  implicit none
  type(sparseR8d_obj) :: this
  integer, intent(in) :: n, nnz, base
  integer, intent(out) :: ier

  ier = 0

  !Define matrix size
  this % base = base
  this % n = n
  this % nnz = nnz
  this % nnzfilled = 0
  this % colfilled = 0

  !Allocate memory for matrix entries, shape descriptors
  allocate(this%values(this%nnz), stat=ier)
  if (ier /=0) return
  allocate(this%rowind(this%nnz), stat=ier)
  if (ier /=0) return
  allocate(this%colptr(this%n+1), stat=ier)
  if (ier /=0) return

  !Initialize values to zero
  this%values(:) = 0._r8
  this%rowind(:) = 0
  this%colptr(:) = 0

  !Allocate SuperLU handles
  call f_create_options_handle(this % options)
  call f_create_ScalePermstruct_handle(this % ScalePermstruct)
  call f_create_LUstruct_handle(this % LUstruct)
  call f_create_SOLVEstruct_handle(this % SOLVEstruct)
  call f_create_SuperMatrix_handle(this % SM)
  call f_create_SuperLUStat_handle(this % stat)

  this % initialized = .FALSE.
end subroutine sparseR8d_init

!============================================================
subroutine sparseR8d_free(this, ier)

  implicit none
  type(sparseR8d_obj) :: this
  integer, intent(out) :: ier

  ier = 0

  !Deallocate data storage
  deallocate(this%values, stat=ier)
  deallocate(this%rowind, stat=ier)
  deallocate(this%colptr, stat=ier)

  !Deallocate SuperLU handles
  call f_destroy_options_handle(this % options)
  call f_destroy_ScalePermstruct_handle(this % ScalePermstruct)
  call f_destroy_LUstruct_handle(this % LUstruct)
  call f_destroy_SOLVEstruct_handle(this % SOLVEstruct)
  call f_destroy_SuperMatrix_handle(this % SM)
  call f_destroy_SuperLUStat_handle(this % stat)
end subroutine sparseR8d_free

!============================================================
subroutine sparseR8d_set_next(this, row, newcol_flag, value, ier)

  implicit none
  type(sparseR8d_obj) :: this
  integer, intent(in) :: row 
  logical, intent(in) :: newcol_flag ! .TRUE. if new column, .FALSE. otherwise
  real(r8), intent(in) :: value
  integer, intent(out) :: ier

  ier = 0

  !Increment the counter of filled-in nonzeroes
  this%nnzfilled = this%nnzfilled + 1
  if (this%nnzfilled > this%nnz) then !Error: too many entries!
     ier = 3
     return
  endif

  !Set the row number of the new element
  this%rowind(this%nnzfilled) = row

  !Set the value of the new element
  this%values(this%nnzfilled) = value

  !Increment the filled-column counter if the column is new
  if (newcol_flag) then
     this%colfilled = this%colfilled + 1
     if (this%colfilled > this%n) then !Error: too many columns!
        ier = 4
        return
     endif

     !Set correct value for ptr to 1st row in column to the left
     this%colptr(this%colfilled) = this%nnzfilled - 1 + this%base
  endif

  !Set final column pointer to upper limit for entry list
  this%colptr(this%colfilled + 1) = this%nnzfilled + this%base
end subroutine sparseR8d_set_next

!============================================================
!Create a distributed SuperLU "SuperMatrix" attached to the
! object provided.
subroutine sparseR8d_new(this, ier)
  use superlu_mod
  implicit none
  type(sparseR8d_obj)     :: this
  integer, intent(out)    :: ier
  real(r8) :: rtmp
  real :: stmp

  ier = 0

  !Check that the matrix has been properly initialized
  if (this%colfilled < this%n) then
     write(0,*)'Error: nonfilled matrix supplied to sparseR8d_new.'
     write(0,*)this%colfilled,' of ',this%n,' columns filled in.'
     ier = 1
     return
  endif
!  if (this%nnzfilled < this%nnz) then
!     write(0,*)'Warning: nonfilled matrix supplied to sparseR8d_new.'
!     write(0,*)'All ',this%n,' columns filled in.'
!     write(0,*)this%nnzfilled,' of ',this%nnz,' entries filled in.'
!  endif

!     write (*,*) 'debug sgi mpi - 7-8-1'
      rtmp = 1._r8
      stmp = 1.
!     write(*,*) 'debug sgi mpi - 7-8-2' , EPSILON(rtmp), EPSILON(stmp)
  !Distribute the matrix to the process grid
  call f_dcreate_dist_matrix(this%SM, this%n, this%n, this%nnzfilled, &
       this%values(1), this%rowind(1), this%colptr(1), SLUprocgrid)

!     write (*,*) 'debug sgi mpi - 7-8-3'
  !Initialize options
  !Set up the default input options:
  ! Fact = DOFACT : Factorize the matrix from scratch, as if brand new.
  ! Equil = YES   : Equilibrate the system to get a large diagonal
  ! ColPerm = MMD_AT_PLUS_A : column ordering is minimum degree on A^T+A
  ! RowPerm = LargeDiag : Use weighted bipartite matching algorithm for perm
  ! ReplaceTinyPivot = YES : tiny diags -> sqrt(eps)||A|| during factorization
  ! Trans = NOTRANS
  ! IterRefine = DOUBLE : accumulate double precision residual in refinement
  ! SolveInitialized = NO : Has initialization been done to triang. solve?
  ! RefineInitialized = NO : Has initialized been done to matvecmult routine?
  ! PrintStat = NO : Print solver statistics? (may have to modify SRC/util.c.)
  call f_set_default_options(this%options)

!     write (*,*) 'debug sgi mpi - 7-8-3'
  !Customize column permutation option: 0=Natural; 1=MMD_ATA;
  ! 2=MMD_AT_PLUS_A (default); 3=COLAMD)
! call set_superlu_options(this%options, ColPerm=3)

!     write (*,*) 'debug sgi mpi - 7-8-4'
  !Allocate storage for ScalePermStruct and LUstruct
  call f_ScalePermstructInit(this%n, this%n, this%ScalePermstruct)
  call f_LUstructInit(this%n, this%n, this%LUstruct)

!     write (*,*) 'debug sgi mpi - 7-8-5'
  !Initialize statistics variables
  call f_PStatInit(this%stat)

  !Set flags
  this % initialized = .TRUE.
end subroutine sparseR8d_new

!============================================================
!Deallocate superLU structures associated with a sparse matrix
! object.
subroutine sparseR8d_del(this, ier)
  implicit none
  type(sparseR8d_obj)     :: this
  integer, intent(out)    :: ier

  ier = 0

  if (this % initialized) then
     !Statistics structure
     call f_PStatFree(this%stat)

     !SuperMatrix
     call f_Destroy_CompRowLoc_Matrix_dist(this%SM)

     !Column permutation, pivot scaling structure
     call f_ScalePermstructFree(this%ScalePermstruct)

     !LU factorization
     call f_Destroy_LU(this%n, SLUprocgrid, this%LUstruct)

     this % initialized = .FALSE.
  else
     ier = 1
  endif
end subroutine sparseR8d_del

!============================================================
subroutine sparseR8d_solve(this, rhs, ier)
  use superlu_mod
  implicit none
  include 'mpif.h'

  type(sparseR8d_obj)     :: this
  real(r8), intent(inout) :: rhs(:)
  integer, intent(out)    :: ier

  integer fstrow, info, ldb, nrow, ncol, nnzloc
  real(r8)              :: berr(1)  !Backward error holder
  real(r8), allocatable :: b(:)     !RHS/solution vector holder

  ier = 0

  !Error checking
  if (.NOT. this%initialized) then
     write(0,*)'Error: sparseR8d_solve called for unallocated matrix.'
     ier = 1
     return
  endif

  !Set up the right hand side by finding local piece of global vector
  ! Important args: ldb    = number of local rows.
  !                 fstrow = global index of 1st local row.
  call f_get_CompRowLoc_Matrix(this%SM, nrow, ncol, nnzloc, ldb, fstrow)
! write(*,*) fstrow, info, ldb, nrow, ncol, nnzloc
  allocate(b(ldb))
  b(1:ldb) = rhs(fstrow+1:fstrow+ldb)

  !!Double-check factorization
   call get_superlu_options(this%options, Fact=info)
!  write(*,*)'Factorization state =',info
!  call SLUD_printopts(this)

  !Call the linear equation "expert" solver routine
  call f_pdgssvx(this%options, this%SM, this%ScalePermstruct, b, ldb, &
       1, SLUprocgrid, this%LUstruct, this%SOLVEstruct, berr, this%stat, &
       info)

  !Check info
  if (info == 0) then !Solve was successful
!     write (*,*) 'Backward error: ',berr(1)

     !Copy (distribute) the local solution back into the global vector
     !do info=1,ldb
     !   print *,'x[',info+fstrow,'] = ',b(info) 
     !enddo
     if (SLUprocsize == 1) then !copy answer
        rhs(1:ldb) = b(1:ldb)
     else                       !redistribute answer
        !ALLTOALL should work but doesn't, so do it this way:
        call MPI_Gather(b, ldb, MPI_DOUBLE_PRECISION, &
             rhs, ldb, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
        call MPI_BCast(rhs, this%n, MPI_DOUBLE_PRECISION, 0, &
             MPI_COMM_WORLD, ier)
     endif

     !Tell options the matrix is factored, so only back-substitution needs
     ! to be done next time. (3=Factored)
     call set_superlu_options(this%options, Fact=3)

  else !Solve failed
!    write (*,*) 'INFO from f_pdgssvz = ', info
  endif

  deallocate(b)
end subroutine sparseR8d_solve

!============================================================
!Solve with a matrix and rhs that will only be used once...
subroutine sparseR8d_driver(this, rhs, ier)
  implicit none
  type(sparseR8d_obj)     :: this
  real(r8), intent(inout) :: rhs(:)
  integer, intent(out)    :: ier

  !Create the supermatrix
  call sparseR8d_new(this, ier)

  !Solve
  call sparseR8d_solve(this, rhs, ier)

  !Free the supermatrix
  call sparseR8d_del(this, ier)
end subroutine sparseR8d_driver
 
!============================================================ 
subroutine sparseR8d_A_dot_x(this, x, res, ier) 
  implicit none
  type(sparseR8d_obj) :: this
  real(r8), intent(in) :: x(:) ! vector
  real(r8), intent(out) :: res(:) ! result
  integer, intent(out) :: ier

  integer row, nzrow, col

  ier = 0
  res(:) = 0

  do col=1,this%colfilled
     do nzrow=this%colptr(col)+1-this%base, this%colptr(col+1)-this%base
        row = this%rowind(nzrow) + 1 - this%base
        res(row) = res(row) + this%values(nzrow)*x(col)
     enddo !nzrow
  enddo !col
end subroutine sparseR8d_A_dot_x

!============================================================
subroutine sparseR8d_dump(this)
  implicit none
  type(sparseR8d_obj) :: this

  integer col, elm

  PRINT *,'Target is ',this%n,' by ',this%n,' with ',this%nnzfilled, &
       ' nonzero entries.'

  do col=1,this%colfilled
     PRINT *,'Column ',col,' begins with element ',this%colptr(col)
     do elm=this%colptr(col)+1-this%base,this%colptr(col+1)-this%base
        PRINT *,' row ',this%rowind(elm)+this%base,': ',this%values(elm)
     enddo
  enddo
end subroutine sparseR8d_dump

!============================================================
subroutine SLUD_printopts(this)
  implicit none
  include 'mpif.h'
  type(sparseR8d_obj) :: this

  integer Fact, Trans, Equil, RowPerm, ColPerm, myrank
  integer ReplaceTiny, IterRefine, SolveInit, RefineInit

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, Fact)

  call f_get_superlu_options(this%options, Fact, Trans, Equil, RowPerm, &
       ColPerm, ReplaceTiny, IterRefine, SolveInit, RefineInit)

  if (myrank == 0) then
     write(0,*)myrank,': Fact =',Fact
     write(0,*)myrank,': ColPerm =',ColPerm
     write(0,*)myrank,': SolveInit =',SolveInit
     write(0,*)myrank,': RefineInit =',RefineInit
  endif
end subroutine SLUD_printopts

subroutine printvalues(this)
  implicit none
  include 'mpif.h'
  type(sparseR8d_obj) :: this
  integer i, j, counter
  
  counter = 1
  
  do i=1,this%colfilled
!     write(*,*) 'columns ',this%colptr(i+1), this%colptr(i)
     do j=1,this%colptr(i+1)-this%colptr(i)
        if(this%values(counter) .ne. 0.) then
           write(*,*) i-1,this%rowind(counter), this%values(counter)
        endif
        counter = counter+1
     enddo
  enddo
  
  call printarray(this%values,this%nnz,0)
end subroutine printvalues

end module supralu_dist_mod
