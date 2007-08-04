module comparedofarrays
!  compares two arrays of unknown/rhs vectors of the same numvar size
!  the first vector is stored in tempcompare and after calling
!  the comparevecs subroutine, the other passed in vector minus the first
!  vector is also stored in tempcompare
  use arrays
  integer :: mynumberingid, init

  data init /0/
  
contains
  subroutine savefirstdofarray(vec, numberingid)
    implicit none
    integer, intent(in) :: numberingid
    integer :: i, ndof
    real, intent(in), allocatable:: vec(:)

    if(init .eq. 1) then
       call deletevec(tempcompare)
    endif

    call createvec(tempcompare, numberingid)
    mynumberingid = numberingid
    init = 1
    tempcompare = vec

  end subroutine savefirstdofarray


  subroutine comparevecs(vec, numberingid)
    implicit none
    integer, intent(in) :: numberingid
    real, intent(in), allocatable:: vec(:)
    integer :: i, ndofs

    if(init .eq. 0) then
       print *, 'WARNING: need to call savefirstdofarray before calling comparevecs'
       return;
    endif

    if(numberingid .ne. mynumberingid) then
       print *, 'WARNING: need to compare vecs of the same numvar type'
       return
    endif

    call numdofs(numberingid, ndofs)
    do i=1,ndofs
       tempvar(i) = vec(i) - tempvar(i)
    enddo

  end subroutine comparevecs

  subroutine deletecomparevec()
    implicit none
    
    if(init .eq. 1) then
       call deletevec(tempcompare)
       init = 0
    endif
  end subroutine deletecomparevec
    
end module comparedofarrays




!==================================
!#define DRIVER
#ifdef DRIVER
program driver
  
  implicit none
  integer numvar
  
  numvar = 1
  call compload()
  
  call comptecplot(numvar)
  call compfree()
  
  stop
end program driver
#endif

!===================

