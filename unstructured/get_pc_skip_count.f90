  subroutine get_pc_skip_count(MAX_SAME_PC_COUNT)
    use basic, ONLY: iskippc, myrank
    implicit none
   
    integer :: MAX_SAME_PC_COUNT

    MAX_SAME_PC_COUNT = iskippc
    if(myrank.eq.0) &
    write(*,*) "MAX_SAME_PC_COUNT=", MAX_SAME_PC_COUNT
    return
  end subroutine get_pc_skip_count
