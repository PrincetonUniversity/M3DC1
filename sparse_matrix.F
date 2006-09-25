c     this module has routines to save a matrix in the harwell-boeing
c     sparse matrix stiffness storage format.  note that the values
c     are stored by column NOT row.  this module also assumes that the
c     sparseness pattern of the matrix is symmetric


      module sparse_matrix
#ifdef mpi
      use supralu_dist_mod
#else
      use superlu_mod
#endif
      use sparse_params
      implicit none
      
      type sparse_matrix_obj
! status specifies what form matrix is in: 
!     0 = uninitialized
!     1 = initialized and some/all values may be set but no dirichlet BC's set
!     2 = values set and some/all dirichlet BC's applied
!     3 = matrix finalized and stored on all processes 
!     4 = matrix stored in superlu matrix
      integer status, which     
      type (sparse_params_obj) , pointer ::  spo_ptr ! pointer to its sparse_params_obj
      real , allocatable:: nzvals(:)                ! array of non-zero values of stiffness matrix
#ifdef mpi
c     superlu_mat is the superlu matrix
      type (sparseR8d_obj) :: superlu_mat
#else
      use supralu_mod
      type (sparseR8_obj) :: superlu_mat
#endif
      end type sparse_matrix_obj
      contains

      subroutine init_smo(this,spo,which_in) 
      implicit none
      type(sparse_params_obj), target :: spo
      type(sparse_matrix_obj) this
      integer which_in
      
      this%which = which_in
      this%spo_ptr => spo
      allocate(this%nzvals(spo%nnz))
      this%nzvals = 0
      this%status = 1
      end subroutine init_smo

      subroutine free_smo(this)
      implicit none
      type(sparse_matrix_obj) this
      integer ier

      if(this%status .eq. 1 .or. this%status .eq. 2 .or. 
     &     this%status .eq. 3) then
         deallocate(this%nzvals)
      else if(this%status.eq. 4) then
#ifdef mpi
         call sparseR8d_del(this%superlu_mat,ier)
         call sparseR8d_free(this%superlu_mat,ier)
#else
         call sparseR8_free(this%superlu_mat,ier)
#endif
      endif
      this%status=0
      end subroutine free_smo
      
      subroutine insert_val(this, val, i, j, operation)
!     do matrix(i,j) = operation*matrix(i,j) + val
!     note that if operation is 1 that the value is added
!     and if operation is 0 that the value is set
      implicit none
      type(sparse_matrix_obj) this
      real val
      integer i, j, rowbegin, rowend, k, done, operation

      if(this%status .ne. 1) call safestop(673)
      rowbegin = this%spo_ptr%colptr(j)
      rowend = this%spo_ptr%colptr(j+1)
      k=rowbegin
      done = 0
      do while(k .lt. rowend .and. done .eq. 0)
         if(i .eq. this%spo_ptr%rowind(k)) then
            if(k .eq. 0 .or. k .gt. this%spo_ptr%nnz) then
               write(*,*) 'wrong k', k, this%spo_ptr%nnz, this%which,i,j
               call safestop(88)
            endif
            this%nzvals(k) = this%nzvals(k)*operation + val
!            write(*,*) 'k is ',k,i,j, val, this%nzvals(k)
            done = 1
         endif
         k = k + 1
      enddo
      if( done.eq. 0) then
         write(*,*) 'k is ',k,rowbegin,rowend,i,j,
     &        this%spo_ptr%numcols,this%which
         call safestop(6689)
      endif

      end subroutine insert_val
      
      subroutine insert_mult_vals(this, vals, iarray, niarray, jarray, 
     &     njarray, operation)
!     do matrix(iarray(i),jarray(j)) = matrix(iarray(i),jarray(j)) + vals(i,j)
!     note that if operation is 1 that the value is added
!     and if operation is 0 that the value is set
      implicit none
      type(sparse_matrix_obj) this
      real vals(niarray,njarray)
      integer i, j, iarray(*), jarray(*), niarray, njarray, ia, ja,
     &     rowbegin, rowend, k, done, operation
      
      if(this%status .ne. 1) call safestop(674)
      do ja=1,njarray
         rowbegin = this%spo_ptr%colptr(jarray(ja))
         rowend = this%spo_ptr%colptr(jarray(ja))
         do ia=1,niarray
            k=rowbegin
            done = 0
            do while(k .le. rowend .and. done .eq. 0)
               if(iarray(ia) .eq. 
     &              this%spo_ptr%rowind(k)) then
                  this%nzvals(k) = this%nzvals(k)*operation+vals(ia,ja)
                  done = 1
               endif
               k = k + 1
            enddo
         enddo
      enddo
         
      end subroutine insert_mult_vals
      
      subroutine set_diri_bc(this,irow)
c     this routine zeroes out all non-diagonal terms
c     and sets the diagonal to 1 for row "irow"
      implicit none
      type(sparse_matrix_obj) this
      integer irow, rowbegin, rowend, k

      if(this%status .ne. 1 .and. this%status .ne. 2) call safestop(675)
      if(irow .gt. this%spo_ptr%numcols .or. irow .lt. 0) then
         call safestop(676)
      endif
      rowbegin = this%spo_ptr%colptr(irow)
      rowend = this%spo_ptr%colptr(irow+1)
      do k=rowbegin,rowend-1
         if(irow.ne.this%spo_ptr%rowind(k)) then
c     since this is the column, now find the row and zero out the value
c     note that this assumes that the matrix structure is symmetric!!!
            this%status = 1 ! temporarily set status to 1
            call insert_val(this, 0., irow, this%spo_ptr%rowind(k),0)
            this%status = 2 ! put status back to 2 for bc
         else
            this%nzvals(k) = 1.
         endif
      enddo
      this%status = 2
      end subroutine set_diri_bc
      
      subroutine zero_array(this, spo_ptr, which_in)
      use sparse_params
      implicit none
      type(sparse_matrix_obj) this
      type(sparse_params_obj), pointer:: spo_ptr
      integer ier, which_in
 
      this%spo_ptr => spo_ptr
      if(this%status .eq. 4) then
         allocate(this%nzvals(this%spo_ptr%nnz))
#ifdef mpi
         call sparseR8d_free(this%superlu_mat,ier)
#else
         call sparseR8_free(this%superlu_mat,ier)
#endif
         if(ier .ne. 0) then
            write(*,*) 'problem freeing superlu matrix, error=',ier
            call safestop(4477)
         endif
      else if(this%status .eq. 0) then
         allocate(this%nzvals(this%spo_ptr%nnz))
      endif
      this%nzvals = 0.
      this%status = 1
      this%which = which_in
      end subroutine zero_array
      
      subroutine finalize_array(this)
!     after all calls to insert_val and set_diri_bc subroutines have been done 
!     we need to do mpi_allreduce on nzval if we are in parallel. then the
!     matrix is converted to a superlu matrix and "my" matrix is deallocated
      implicit none
      type(sparse_matrix_obj) this
      integer ier
#ifdef mpi
      real , allocatable:: nzvals_copy(:)
      include 'mpif.h'
      allocate(nzvals_copy(this%spo_ptr%nnz))
      call mpi_allreduce(this%nzvals, nzvals_copy, this%spo_ptr%nnz,
     &     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
      this%nzvals = nzvals_copy
      deallocate(nzvals_copy)
#endif
      this%status = 3
      end subroutine finalize_array
      
      subroutine solve(this, rhs_sol, ier)
! solve for rhs and replace with solution
      implicit none
      type(sparse_matrix_obj) this
      integer ier, i, j, count
      logical newcol_flag
      real rhs_sol(this%spo_ptr%numcols)
! probably need to check to see if matrix has been given to superlu and factored
      if(this%status .ne. 3 .and. this%status.ne. 4) then
         write(*,*) 'problem with matrix ', this%which, this%status
         call safestop(9933)
      endif

!     now create the superlu matrix if it hasn't been created yet
      if(this%status .eq. 3) then
#ifdef mpi
         call sparseR8d_init(this%superlu_mat, this%spo_ptr%numcols, 
     &        this%spo_ptr%nnz, 0, ier)
         count = 1
         do i=1,this%spo_ptr%numcols
            newcol_flag = .TRUE.
            do j=this%spo_ptr%colptr(i),this%spo_ptr%colptr(i+1)-1
               call sparseR8d_set_next(this%superlu_mat,
     &              this%spo_ptr%rowind(count)-1,newcol_flag,
     &              this%nzvals(count),ier)
               newcol_flag = .FALSE.
               count = count + 1
            enddo
         enddo
c     call printvalues(this%superlu_mat)
         call sparseR8d_new(this%superlu_mat, ier)
#else 
         call safestop(111)     ! put something in here later for serial
#endif
c     deallocate "my" matrix
         deallocate(this%nzvals)
      endif ! this%status = 3
      this%status=4

#ifdef mpi
      call sparseR8d_solve(this%superlu_mat,rhs_sol,ier)
#else
      call dsupralu_solve(this%superlu_mat,rhs_sol,ier)
#endif
      end subroutine solve
      
      subroutine matrix_vector_mult(this, x, b, nsize)   
! matrix * x = b where matrix = this%nzvals
      implicit none
      type(sparse_matrix_obj) this
      real x(*), b(nsize)
      integer nsize, i, j, counter
      
      counter = 1
      if(this%status .ne. 3) call safestop(4896)
      if(nsize .ne. this%spo_ptr%numcols) call safestop(4897)
      do i=1,nsize
         b(i) = 0.
      enddo
      do i=1,this%spo_ptr%numcols
         do j=1,this%spo_ptr%colptr(i+1)-this%spo_ptr%colptr(i)
            b(this%spo_ptr%rowind(counter)) =
     &           b(this%spo_ptr%rowind(counter)) + 
     &           this%nzvals(counter)*x(i)
            counter = counter+1
         enddo
      enddo
      
      
      end subroutine matrix_vector_mult
      
      subroutine getSparseMatrix(this, superlu_mat)
      implicit none
      type (sparse_matrix_obj) this
#ifdef mpi
      type (sparseR8d_obj) :: superlu_mat
#else
      type (sparseR8_obj) :: superlu_mat
#endif
      superlu_mat = this%superlu_mat
      end subroutine getSparseMatrix

      subroutine printsparsearray(this)
      implicit none
      type (sparse_matrix_obj) this
      integer i, j, counter
      real sum, sum2
      sum = 0
      sum2 = 0
      counter = 1
      
      do i=1,this%spo_ptr%numcols
         do j=1,this%spo_ptr%colptr(i+1)-this%spo_ptr%colptr(i)
            if(this%nzvals(counter) .ne. 0.) then
               write(*,100) this%spo_ptr%rowind(counter)-1,i-1,
     &              this%nzvals(counter)
               if(this%nzvals(counter) > 0) then
                  sum = sum + this%nzvals(counter)
               else
                  sum = sum - this%nzvals(counter)
               endif
               sum2 = sum2 + this%nzvals(counter)*this%nzvals(counter)
            endif
            counter = counter+1
         enddo
      enddo
      write(*,*) 'the sum is',sum, sum2
      call safestop(1112)
 100  format(i4, 1x, i4, e12.4)

      end subroutine printsparsearray

      subroutine checksparsearray(this,array, nsize)
      implicit none
      type (sparse_matrix_obj) this
      integer i, j, counter, nsize
      real array(nsize,nsize), diff
      
      counter = 1
      if(nsize .ne. this%spo_ptr%numcols) call safestop(3388)
      
      do i=1,this%spo_ptr%numcols
         do j=1,this%spo_ptr%colptr(i+1)-this%spo_ptr%colptr(i)
            diff = this%nzvals(counter) - 
     &           array(this%spo_ptr%rowind(counter),i)
            if(diff .lt. 0.) diff = -diff
            if(diff .gt. 0.) then
               write(*,*) 'diff is ',this%nzvals(counter), 
     &              array(this%spo_ptr%rowind(counter),i),
     &              this%spo_ptr%rowind(counter), i,counter
               call safestop(9090)
            endif
            counter = counter+1
         enddo
      enddo
      end subroutine checksparsearray
      
      end module sparse_matrix
