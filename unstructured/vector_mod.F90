! vector_mod.F90 — thin wrapper that selects the vector backend at compile time.
! CMake defines exactly one of the following:
!   M3DC1_VECTOR_IS_SCOREC  → use scorec_vector_mod
!   M3DC1_VECTOR_IS_PETSC   → use petsc_vector_mod

module vector_mod
#ifdef M3DC1_VECTOR_IS_SCOREC
  use scorec_vector_mod
#elif defined(M3DC1_VECTOR_IS_PETSC)
  use petsc_vector_mod
#else
  use scorec_vector_mod
#endif
end module vector_mod
