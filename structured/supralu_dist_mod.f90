! Attempt at a distributed-memory parallel version of A. Pletzer's supralu
! module to interface with SuperLU_DIST v2.0. J. Breslau, 7/30/04
! Modified 3/9/06 to use distributed compressed row storage to parallelize
!  matrix construction. J. Breslau.

module supralu_dist_mod
  implicit none

  integer, parameter :: r8 = selected_real_kind(12,100)
  
  integer SLUprocsize  !Number of active MPI processes
  integer*8 SLUprocgrid  !Opaque handle to SuperLU process grid

  type sparseR8d_obj
     !Shape & data holders
     real(r8), pointer :: values(:) !Storage for local nonzero matrix elements
     integer, pointer  :: colind(:) !Column indices of the nonzeroes
     integer, pointer  :: rowptr(:) !Beginnings of rows in values, colind
     integer, pointer  :: recvcounts(:) !number of local rows, by processor
     integer, pointer  :: displs(:)  !offset of rows, by processor

     !Opaque SuperLU handles
     integer*8 SM               !The SuperLU distributed supermatrix itself
     integer*8 options          !SuperLU options settings for this matrix
     integer*8 ScalePermstruct  !Permutations & pivot scalings done to matrix
     integer*8 LUstruct         !Factorization of matrix (L & U)
     integer*8 SOLVEstruct      !Data related to parallel solution
     integer*8 stat             !Solver statistics

     !Size specifiers
     integer n           !Linear dimension of matrix
     integer nnz_loc     !Number of nominal local non-zero entries
     integer nnzfilled   !Number of entries set so far
     integer m_loc       !Number of local rows
     integer fst_row     !Row number of 1st row in local submatrix
     integer rowsfilled  !Number of rows set so far

     !Initialization flag
     logical alloc            !Has storage been set aside?
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
  npcol = 1
  if(SLUprocsize.eq.2) npcol = 2
  if(SLUprocsize.eq.4) npcol = 2
  if(SLUprocsize.eq.6) npcol = 3
  if(SLUprocsize.eq.8) npcol = 4
  if(SLUprocsize.eq.12) npcol = 4
  if(SLUprocsize.eq.16) npcol = 4
  if(SLUprocsize.eq.24) npcol = 6
  if(SLUprocsize.eq.32) npcol = 8
  nprow = SLUprocsize/npcol

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
subroutine sparseR8d_init(this, n, nnz_loc, m_loc, fst_row, ier)

  ! Constructor: sets up data structure/container for distributed
  !  column-compressed sparse matrix, allocates storage, but does
  !  not call any SuperLU routines
  !
  ! this: handle of the container
  ! n: linear dimension of the (square) matrix
  ! nnz: number of non-zero elements in the matrix
  ! base: indexing base (0 like C or 1 like Fortran)

  implicit none
  include 'mpif.h'

  type(sparseR8d_obj)  :: this
  integer, intent(in)  :: n, nnz_loc, m_loc, fst_row
  integer, intent(out) :: ier
  integer irow

  ier = 0

  !An error check
  call MPI_Reduce(m_loc, irow, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, &
       ier)
  if ((fst_row.eq.0).and.(irow.ne.n)) then
     print *,'WARNING: local rows add up to ',irow,', not ',n
     ier = 1
  endif

  !Set up communication info for matrix vector products
  allocate(this%recvcounts(SLUprocsize), this%displs(SLUprocsize), stat=ier)
  if (ier /=0) return
  call MPI_AllGather(m_loc, 1, MPI_INTEGER, this%recvcounts, 1, MPI_INTEGER, &
       MPI_COMM_WORLD, ier)
  this%displs(1) = 0
  do irow=2,SLUprocsize
     this%displs(irow) = this%displs(irow-1) + this%recvcounts(irow-1)
  end do

  !Define matrix size, reset counters
  this % n = n
  this % nnz_loc = nnz_loc
  this % nnzfilled = 0
  this % m_loc = m_loc
  this % fst_row = fst_row
  this % rowsfilled = 0

  !Allocate memory for matrix entries, shape descriptors
  allocate(this%values(this%nnz_loc))!, stat=ier)
  if (ier /=0) return
  allocate(this%colind(this%nnz_loc), stat=ier)
  if (ier /=0) return
  allocate(this%rowptr(this%m_loc + 1), stat=ier)
  if (ier /=0) return

  !Initialize values to zero
  this%values(:) = 0._r8
  this%colind(:) = 0
  this%rowptr(:) = 0

  this % initialized = .FALSE.
end subroutine sparseR8d_init

!=======================================================================
subroutine sparseR8d_reset(this, ier)
  implicit none

  type(sparseR8d_obj) :: this
  integer, intent(out) :: ier

  ier = 0

  this % nnzfilled = 0
  this % rowsfilled = 0
  this%values(1:this%nnz_loc) = 0._r8
!  this%colind(1:this%nnz_loc) = 0
!  this%rowptr(1:this%m_loc + 1) = 0

!  this % initialized = .FALSE.
end subroutine sparseR8d_reset

!============================================================
subroutine sparseR8d_free(this, ier)

  implicit none
  type(sparseR8d_obj) :: this
  integer, intent(out) :: ier

  ier = 0

  !Deallocate data storage
  if (associated(this%values)) then
     deallocate(this%recvcounts, stat=ier)
     deallocate(this%displs, stat=ier)
     deallocate(this%values, stat=ier)
     deallocate(this%colind, stat=ier)
     deallocate(this%rowptr, stat=ier)
  end if

end subroutine sparseR8d_free

!============================================================
subroutine sparseR8d_set_next(this, col, newrow_flag, value, ier)

  implicit none
  type(sparseR8d_obj)  :: this
  integer, intent(in)  :: col
  logical, intent(in)  :: newrow_flag ! .TRUE. if new column, .FALSE. otherwise
  real(r8), intent(in) :: value
  integer, intent(out) :: ier

  ier = 0

  !Increment the counter of filled-in nonzeroes
  this%nnzfilled = this%nnzfilled + 1

  !Check for overflow
  if (this%nnzfilled > this%nnz_loc) then !Error: too many entries!
     ier = 3
     return
  endif

  !Set the value of the new element
  this%values(this%nnzfilled) = value

  !Set the column number of the new element
  this%colind(this%nnzfilled) = col

  !Check whether we're starting a new row
  if (newrow_flag) then
     !Increment the filled-row counter
     this%rowsfilled = this%rowsfilled + 1

     !Error check
     if (this%rowsfilled > this%m_loc) then !Error: too many columns!
        ier = 4
        return
     endif

     !Set correct value for ptr to 1st column in this row
     this%rowptr(this%rowsfilled) = this%nnzfilled - 1
  endif

  !Set final column pointer to upper limit for entry list
  this%rowptr(this%rowsfilled + 1) = this%nnzfilled
end subroutine sparseR8d_set_next

!=======================================================================
subroutine sparseR8d_set_next_block(this, blocksize, col0, newrow_flag, &
     values, ier)
  implicit none
  type(sparseR8d_obj)  :: this
  integer, intent(in)  :: blocksize
  integer, intent(in)  :: col0
  logical, intent(in)  :: newrow_flag
  real(r8), intent(in) :: values(blocksize)
  integer, intent(out) :: ier

  integer col

  ier = 0

  if (this%initialized) then !Error: already set up
     print *,'Error: tried to add entry to initialized matrix'
     ier = 4
     return
  endif

  !Increment the counter of filled-in nonzeroes
  this%nnzfilled = this%nnzfilled + blocksize

  !Check lfor overflow
  if (this%nnzfilled > this%nnz_loc) then !Error: too many entries!
     ier = 3
     return
  endif

  !Set the values of the new elements
  this%values(this%nnzfilled-blocksize+1:this%nnzfilled) = values

  !Set the column numbers of the new elements
  do col=1,blocksize
     this%colind(this%nnzfilled-blocksize+col) = col0 + col - 1
  end do

  !Check whether we're starting a new row
  if (newrow_flag) then
     !Increment the filled-row counter
     this%rowsfilled = this%rowsfilled + 1

     !Error check
     if (this%rowsfilled > this%m_loc) then !Error: too many columns!
        ier = 4
        return
     endif

     !Set correct value for ptr to 1st column in this row
     this%rowptr(this%rowsfilled) = this%nnzfilled - blocksize
  endif

  !Set final column pointer to upper limit for entry list
  this%rowptr(this%rowsfilled + 1) = this%nnzfilled
end subroutine sparseR8d_set_next_block

!============================================================
!Create a distributed SuperLU "SuperMatrix" attached to the
! object provided.
subroutine sparseR8d_new(this, ier)
  use superlu_mod
  implicit none
  include 'mpif.h'

  type(sparseR8d_obj)     :: this
  integer, intent(out)    :: ier

  ier = 0

  !Allocate SuperLU handles
  call f_create_options_handle(this % options)
  call f_create_ScalePermstruct_handle(this % ScalePermstruct)
  call f_create_LUstruct_handle(this % LUstruct)
  call f_create_SOLVEstruct_handle(this % SOLVEstruct)
  call f_create_SuperMatrix_handle(this % SM)
  call f_create_SuperLUStat_handle(this % stat)

  !Check that the matrix has been properly initialized
  if (this%initialized) then
     print *,'Error: attempted to initialize existing matrix'
     ier = 2
     return
  end if

  !Create distributed SuperLU matrix with data in "this"
  call f_dCreate_CompRowLoc_Matrix_dist(this%SM, this%n, this%n, &
       this%nnz_loc, this%m_loc, this%fst_row, &
       this%values, this%colind, this%rowptr, &
       SLU_NR_loc, SLU_D, SLU_GE)

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
  call printoffc(this%options) !JB 3/6/06 turn off stat printing
  

  !Customize column permutation option: 0=Natural; 1=MMD_ATA;
  ! 2=MMD_AT_PLUS_A (default); 3=COLAMD)
!  call set_superlu_options(this%options, ColPerm=3)
  call set_superlu_options(this%options, RowPerm=0, ColPerm=2)

  !Allocate storage for ScalePermStruct and LUstruct
  call f_ScalePermstructInit(this%n, this%n, this%ScalePermstruct)
  call f_LUstructInit(this%n, this%n, this%LUstruct)

  !Initialize statistics variables
  call f_PStatInit(this%stat)

  !Set flags
  this % initialized = .TRUE.

  call MPI_Barrier(MPI_COMM_WORLD, ier)
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
     call f_Destroy_SuperMatrix_Store_dist(this%SM)

     !Column permutation, pivot scaling structure
     call f_ScalePermstructFree(this%ScalePermstruct)

     !LU factorization
     call f_Destroy_LU(this%n, SLUprocgrid, this%LUstruct)
     call f_LUstructFree(this%LUstruct)

     !Deallocate SuperLU handles
     call f_destroy_options_handle(this % options)
     call f_destroy_ScalePermstruct_handle(this % ScalePermstruct)
     call f_destroy_LUstruct_handle(this % LUstruct)
     call f_destroy_SOLVEstruct_handle(this % SOLVEstruct)
     call f_destroy_SuperMatrix_handle(this % SM)
     call f_destroy_SuperLUStat_handle(this % stat)

     this % initialized = .FALSE.
  else
     print *,'deleting noninitialized matrix!'
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

  real(r8)              :: berr(1)  !Backward error holder
  real(r8), allocatable :: b(:)     !RHS/solution vector holder
  integer fstrow, info, ldb, nrow, ncol, nnzloc

  ier = 0

  !Error checking
!  if (.NOT. this%initialized) then
!     write(0,*)'Error: sparseR8d_solve called for unallocated matrix.'
!     ier = 1
!     return
!  endif

  !Set up the right hand side by finding local piece of global vector
  ! Important args: ldb    = number of local rows.
  !                 fstrow = global index of 1st local row.
  call f_get_CompRowLoc_Matrix(this%SM, nrow, ncol, nnzloc, ldb, fstrow)
  allocate(b(ldb))
  b(1:ldb) = rhs(fstrow+1:fstrow+ldb)

  if(ldb.ne.this%m_loc) write(*,*) "Nate: ldb, m_loc", ldb, this%m_loc
  if(fstrow.ne.this%fst_row) write(*,*) "Nate: fstrow, fst_row", &
       fstrow , this%fst_row

  !Double-check factorization
  !  call get_superlu_options(this%options, Fact=info)
  !  print *,'Factorization state =',info
  !  call SLUD_printopts(this)

  !Call the linear equation "expert" solver routine
  call f_pdgssvx(this%options, this%SM, this%ScalePermstruct, b, ldb, &
       1, SLUprocgrid, this%LUstruct, this%SOLVEstruct, berr, this%stat, &
       info)

  !Check info
  if (info == 0) then !Solve was successful

     !Tell options the matrix is factored, so only back-substitution needs
     ! to be done next time. (3=Factored)
     call set_superlu_options(this%options, Fact=3)

     !Copy (distribute) the local solution back into the global vector
     if (SLUprocsize == 1) then !copy answer
        rhs(1:ldb) = b(1:ldb)
     else                       !redistribute answer
        call MPI_AllGatherv(b, ldb, MPI_DOUBLE_PRECISION, &
             rhs, this%recvcounts, this%displs, MPI_DOUBLE_PRECISION, &
             MPI_COMM_WORLD, ier)
     endif

  else !Solve failed
    write (*,*) 'INFO from f_pdgssvz = ', info
  endif

  deallocate(b)
end subroutine sparseR8d_solve

!============================================================
subroutine sparseR8d_solve_part(this, b, sol, ier)
  use superlu_mod

  implicit none
  include 'mpif.h'

  type(sparseR8d_obj)     :: this
  real(r8), intent(inout) :: b(:)
  real(r8), intent(out) :: sol(:)   !RHS/solution vector holder
  integer, intent(out)    :: ier

  real(r8) :: berr(1)  !Backward error holder
  integer :: info

  ier = 0

  !Call the linear equation "expert" solver routine
  call f_pdgssvx(this%options, this%SM, this%ScalePermstruct, b, this%m_loc, &
       1, SLUprocgrid, this%LUstruct, this%SOLVEstruct, berr, this%stat, &
       info)

  !Check info
  if (info == 0) then !Solve was successful

     !Tell options the matrix is factored, so only back-substitution needs
     ! to be done next time. (3=Factored)
     call set_superlu_options(this%options, Fact=3)

     !Copy (distribute) the local solution back into the global vector
     if (SLUprocsize == 1) then !copy answer
        sol = b
     else                       !redistribute answer
        call MPI_AllGatherv(b, this%m_loc, MPI_DOUBLE_PRECISION, &
             sol, this%recvcounts, this%displs, MPI_DOUBLE_PRECISION, &
             MPI_COMM_WORLD, ier)
     endif

  else !Solve failed
     write (*,*) 'INFO from f_pdgssvz = ', info
  endif

end subroutine sparseR8d_solve_part

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
  include 'mpif.h'

  type(sparseR8d_obj) :: this
  real(r8), intent(in) :: x(this%n) ! global input vector
  real(r8), intent(out) :: res(this%n) ! global result
  integer, intent(out) :: ier

  real(r8), dimension(:), allocatable, save :: locres
  integer, save :: lastsize=-1
  integer row, inz

  ier = 0

  !Reallocate local storage, as seldom as possible
  if (this%m_loc > lastsize) then
     if (lastsize > 0) deallocate(locres)
     allocate(locres(this%m_loc))
     lastsize = this%m_loc
  end if

  !Compute local piece of matrix-vector product
  do row = 1, this%m_loc
     locres(row) = 0._r8
     do inz = this%rowptr(row)+1, this%rowptr(row+1)
        locres(row) = locres(row) + &
             this%values(inz) * x(this%colind(inz) + 1)
     end do !inz
  end do !row

  !Gather local pieces together on all processors
  if (SLUprocsize > 1) then
     call MPI_AllGatherv(locres, this%m_loc, MPI_DOUBLE_PRECISION, &
          res, this%recvcounts, this%displs, MPI_DOUBLE_PRECISION, &
          MPI_COMM_WORLD, ier)
  else
     res(1:this%m_loc) = locres(1:this%m_loc)
  end if
end subroutine sparseR8d_A_dot_x

!============================================================
subroutine sparseR8d_dump(this)
  implicit none
  type(sparseR8d_obj) :: this

  integer row, elm

  PRINT *,'Target is ',this%n,' by ',this%n,' with ',this%nnzfilled, &
       ' nonzero entries.'

  do row=1,this%rowsfilled
     PRINT *,'Row ',row,' begins with element ',this%rowptr(row)
     do elm=this%rowptr(row)+1,this%rowptr(row+1)
        PRINT *,' col ',this%colind(elm),': ',this%values(elm)
     enddo !elm
  enddo !row
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
     print *,myrank,': Fact =',Fact
     print *,myrank,': ColPerm =',ColPerm
     print *,myrank,': SolveInit =',SolveInit
     print *,myrank,': RefineInit =',RefineInit
  endif
end subroutine SLUD_printopts
!============================================================
subroutine sparseR8d_setup(this, numvard, jbl)
  use basic
  use t_data

  implicit none


  type(sparseR8d_obj), intent(inout) :: this
  integer, intent(in) :: numvard
  integer, intent(out), dimension(m*n, 9) :: jbl
  
  integer :: ibig, jbig, blocksize, il, jl, globsize, index, itype, icol
  integer :: degen

  blocksize = 6*numvard

  this % rowptr(1) = 0
  index = 1
  globsize = m*n
  do il=1,this%m_loc
     ibig = (this%fst_row + il - 1)/blocksize + 1

     degen = 0

     ! Inner loop over local columns
     do jbig=1,globsize
        do itype=1,itypemax(ibig)
           if(incb(ibig,itype) + ibig .eq. jbig) go to 102
        enddo
        cycle
102     icol = (jbig-1)*blocksize
        degen = degen + 1
        jbl(ibig,degen) = jbig 

        do jl=1,blocksize         
           ! define colind
           this%colind(index) = icol

           icol = icol + 1  
           index = index + 1
        enddo            !jl
     enddo               !jbig

     ! define rowptr
     this%rowptr(il+1) = this%rowptr(il) + itypemax(ibig)*blocksize
     
     ! clear the rest of jbl
     if(degen.lt.9) then
        do itype=degen+1,9
           jbl(ibig,itype) = 0
        enddo
     endif

  enddo ! il

end subroutine sparseR8d_setup
!===========================================================
subroutine sparseR8d_increment(this, ig, jg, val, numvard, ovwrt)
  use t_data
  use basic

  implicit none

  type(sparseR8d_obj), intent(inout) :: this
  integer, intent(in) :: ig, jg, numvard
  real, intent(in) :: val
  integer, intent(in) :: ovwrt
  
  integer :: index, blocksize, ibig, jbig, il, jl, itype, ji


  il = ig - this%fst_row

  if(il.lt.1 .or. il.gt.this%m_loc) then
     write(*,*) "Error: i1 not in local submatrix. ig, fst_row, mloc=", &
          ig, this%fst_row, this%m_loc
     return
  end if

  blocksize = 6*numvard      
  ibig = (ig - 1) / blocksize + 1
  
  index = this%rowptr(ig-this%fst_row)+1

  do ji = 1, 9
     jbig = jbiglist(ibig,ji)
     if(jbig.eq.0) return
     do itype=1,itypemax(ibig)
        if(incb(ibig,itype) + ibig .eq. jbig) goto 102
     enddo
     cycle
102  jl = jg - (jbig-1)*blocksize - 1

     if(jl.ge.0 .and. jl.lt.blocksize) then
        if(ovwrt.eq.1) then
           this%values(index+jl) = val
        else 
           this%values(index+jl) = this%values(index+jl) + val
        endif
     endif

     index = index + blocksize
  enddo
end subroutine sparseR8d_increment
!===========================================================
subroutine sparseR8d_increment_array(this, ig, jg, val, numvard)
  use t_data
  use basic

  implicit none

  type(sparseR8d_obj), intent(inout) :: this
  integer, intent(in) :: ig, jg, numvard
  real, intent(in), dimension(3,3) :: val
  
  integer :: rowindex, index, blocksize, ibig, jbig, il, jl, itype
  integer :: ji, offset

  il = ig - this%fst_row

  if(il.lt.1 .or. il.gt.this%m_loc) then
     write(*,*) "Error: i1 not in local submatrix. ig, fst_row, mloc=", &
          ig, this%fst_row, this%m_loc
     return
  end if

  blocksize = 6*numvard      
  ibig = (ig - 1) / blocksize + 1
  
  offset = 0
  rowindex = this%rowptr(ig-this%fst_row)+1
      
  do ji = 1, 9
     jbig = jbiglist(ibig,ji)
     if(jbig.eq.0) return
     do itype=1,itypemax(ibig)
        if(incb(ibig,itype) + ibig .eq. jbig) goto 102
     enddo
     cycle
102  jl = jg - (jbig-1)*blocksize - 1
     
     if(jl.ge.0 .and. jl.lt.blocksize) then
        index = rowindex + offset + jl
        this%values(index) = this%values(index) + val(1,1)
        if(numvard.ge.2) then 
           this%values(index+6) = this%values(index+6) + val(1,2)

           index = this%rowptr(ig-this%fst_row+ 6) + 1 + offset + jl
           this%values(index  ) = this%values(index  ) + val(2,1)
           this%values(index+6) = this%values(index+6) + val(2,2)
        endif
        if(numvard.ge.3) then
           index = this%rowptr(ig-this%fst_row   ) + 1 + offset + jl
           this%values(index+12) = this%values(index+12) + val(1,3)
           
           index = this%rowptr(ig-this%fst_row+ 6) + 1 + offset + jl
           this%values(index+12) = this%values(index+12) + val(2,3)
           
           index = this%rowptr(ig-this%fst_row+12) + 1 + offset + jl
           this%values(index   ) = this%values(index   ) + val(3,1)
           this%values(index+ 6) = this%values(index+ 6) + val(3,2)
           this%values(index+12) = this%values(index+12) + val(3,3)
        endif
     endif
     
     offset = offset + blocksize
  enddo

end subroutine sparseR8d_increment_array
!===========================================================
subroutine sparseR8d_zero_row(this, ig)
  use t_data
  use basic

  implicit none

  type(sparseR8d_obj), intent(inout) :: this
  integer, intent(in) :: ig
  
  integer :: ifirst, ilast
 
  ifirst = this%rowptr(ig-this%fst_row)+1
  ilast = this%rowptr(ig-this%fst_row+1)

  this%values(ifirst:ilast) = 0.

end subroutine sparseR8d_zero_row
!===========================================================
integer function sparseR8d_is_local_row(this, i1)

  implicit none

  type(sparseR8d_obj), intent(in) :: this
  integer, intent(in) :: i1

!     check if this node is assigned to this processor
  if(((i1-1).lt.this%fst_row) .or. &
       ((i1-1).ge.(this%fst_row + this%m_loc))) then
     sparseR8d_is_local_row = 0
  else 
     sparseR8d_is_local_row = 1
  endif
end function sparseR8d_is_local_row
 

end module supralu_dist_mod



module superlu
  use p_data
  use supralu_dist_mod
  type(sparseR8d_obj) :: s1matrix_lu, d1matrix_lu
  type(sparseR8d_obj) :: s2matrix_lu, d2matrix_lu
  type(sparseR8d_obj) :: s3matrix_lu, d3matrix_lu
  type(sparseR8d_obj) :: s4matrix_lu, d4matrix_lu
  type(sparseR8d_obj) :: s5matrix_lu, d5matrix_lu
  type(sparseR8d_obj) :: s6matrix_lu, samatrix_lu
  type(sparseR8d_obj) :: s7matrix_lu, d7matrix_lu
  type(sparseR8d_obj) :: s8matrix_lu, d8matrix_lu
  type(sparseR8d_obj) :: s9matrix_lu, d9matrix_lu, r9matrix_lu, q9matrix_lu
  type(sparseR8d_obj) :: r1matrix_lu, r2matrix_lu, q2matrix_lu, smmatrix_lu
  type(sparseR8d_obj) :: r8matrix_lu, q8matrix_lu
  integer :: base
      
end module superlu
