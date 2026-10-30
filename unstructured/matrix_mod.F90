! matrix_mod.F90 — thin wrapper that selects the matrix backend at compile time.
! CMake defines exactly one of the following:
!   M3DC1_MATRIX_IS_SCOREC  → use scorec_matrix_mod
!   M3DC1_MATRIX_IS_PETSC   → use petsc_matrix_mod

module matrix_mod
#ifdef M3DC1_MATRIX_IS_SCOREC
  use scorec_matrix_mod
#elif defined(M3DC1_MATRIX_IS_PETSC)
  use petsc_matrix_mod
#else
  use scorec_matrix_mod
#endif
end module matrix_mod
