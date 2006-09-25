      module linklist
      implicit none
      type linklistobj

      integer icol
      type (linklistobj), POINTER :: next
      end type linklistobj

      
      end module linklist

      module llptr
      use linklist
      implicit none
      type llptrobj

      type (linklistobj), pointer :: ptr
      end type llptrobj
      contains
      subroutine null(this, icol)
      implicit none
      integer icol
      type (llptrobj) this
      type (linklistobj), pointer :: new
c      allocate(new)
c      new%icol = icol
c      nullify(new%next)
c      this%ptr => new
      nullify (this%ptr)
      end subroutine null
      end module llptr

      subroutine listinsert(in, inum)
      use linklist
      use llptr
      implicit none
      type (linklistobj), pointer :: new, prev, next
      type (llptrobj) in
      integer inum, i
      
      if(ASSOCIATED(in%ptr)) then
         if(in%ptr%icol .gt. inum) then
            allocate(new)
            new%next => in%ptr
            new%icol = inum
            in%ptr => new
         else
            prev => in%ptr
            next => prev%next
            do while(associated(next) .and. next%icol .lt. inum )
               prev => next
               next => next%next
            enddo
            if(associated(next)) then
               if(next%icol .gt. inum .and. prev%icol .ne. inum) then
                  allocate(new)
c     if(inum .eq. 1) call safestop(1)
                  new%icol = inum
                  prev%next => new
                  new%next => next
               endif
            else if(prev%icol .ne. inum) then
               allocate(new)
               new%icol = inum
               prev%next => new
               nullify(new%next)
            endif
         endif
      else
         allocate(new)
         new%icol = inum
         in%ptr => new
         nullify(new%next)
      endif

      return
      end
      
      subroutine listdelete(in)
      use linklist
      use llptr
      implicit none
      type (llptrobj) in
      type (linklistobj), pointer :: curr, next
      next => in%ptr
      do while(associated(next))
         curr => next
         next => next%next
         deallocate(curr)
         nullify(curr)
      enddo
      nullify(in%ptr)
      return
      end

      subroutine listlength(in, ilength)
      use linklist
      use llptr
      implicit none
      type (llptrobj) in
      type (linklistobj), pointer :: next
      integer ilength

      ilength = 0
      next => in%ptr
      do while(associated(next))
c         write(*,*) 'icol ', next%icol
         ilength = ilength + 1
         next => next%next
      enddo

      return
      end


      module sparse_params
      implicit none

      type sparse_params_obj
c     this assumes a fortran ordering
      integer numcols       ! matrix rank
      integer nnz          ! number of non-zeroes in sparse stiffness matrix
      integer initialized  ! specifies whether values have been set
      integer , allocatable:: rowind(:)    ! row index of Harwell-Boeing sparse matrix - dimension is nnz
      integer , allocatable:: colptr(:)    ! column pointer of Harwell-Boeing sparse matrix - dimension is numcols+1

      end type sparse_params_obj

      contains

      subroutine init_spo(this,numvarl)
! allocates rowind and colptr and sets numcols, nnz, initialized=1, colptr and rowind
      use linklist
      use llptr
      implicit none
      integer numvarl,i,j,k,numedges,count,maxnumedges,     
     &     idofs, jdofs, ndof, nodeid, ier, ibegin(3), iendplusone(3),
     &     numfaces, nodes(4), iibegin, iiendplusone, irow, icol
      integer , allocatable :: edgeids(:), itemp(:)
      type (sparse_params_obj) this
      type (linklistobj), pointer :: curr
      type (llptrobj), allocatable :: llist(:)
      include 'mpif.h'

      call numfac(numfaces)
      call numdofs(numvarl, ndof)
      allocate(llist(ndof))
      do i=1,ndof
         call null(llist(i), i)
      enddo
      write(*,*) 'there are this many dofs',ndof

      do i=1,numfaces
         call nodfac(i, nodes)
         do j=1,3
            call entdofs(numvarl, nodes(j), 0, ibegin(J), 
     &           iendplusone(J))
         enddo
         do j=1,3
            do irow=ibegin(j),iendplusone(j)-1
               do k=1,3
                  do icol=ibegin(k),iendplusone(k)-1
                     call listinsert(llist(irow), icol)
                  enddo
               enddo
            enddo
         enddo
      enddo
                        ! first calculate nnz and find the maximum number of edges
      j=0
      do i=1,ndof
c         write(*,*) 'lenghting ',i
         call listlength(llist(i), k)
         j = j + k
      enddo
      this%numcols = ndof
      this%nnz = j
      
c     now allocate rowind and colptr and fill in their values
      allocate(this%rowind(this%nnz),this%colptr(this%numcols+1))
      this%rowind = 0
      this%colptr = 0
      count = 0
      allocate(edgeids(maxnumedges),itemp((maxnumedges+1)))
      this%colptr(1) = 1
      count = 0                 ! used to keep track of what was put into rowind
      do i=1,ndof
         call listlength(llist(i), k)
         this%colptr(i+1) = k + this%colptr(i)
         curr => llist(i)%ptr
         do while(associated(curr))
            count = count + 1
            this%rowind(count) = curr%icol
            curr => curr%next
         enddo
      enddo

      if(count .ne. this%nnz) then
         write(*,*) 'wrong nnz ', this%nnz, count, numvarl
         call safestop (3324)
      endif
      
      if(this%colptr(ndof+1) .ne. this%nnz+1) then
         write(*,*) 'colptr wrong ',
     %        this%colptr(ndof+1), this%nnz
         call safestop(3325)
      endif
      
      this%initialized = 1
      do i=1,ndof
         call listdelete(llist(i))
      enddo
      deallocate(llist)
      end subroutine init_spo
      
      subroutine free_spo(this)
                                ! deallocates all memory
      implicit none
      type (sparse_params_obj) this
      deallocate(this%rowind,this%colptr)
      this%initialized = 0
      end subroutine free_spo
      
      end module sparse_params

