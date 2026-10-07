! mesh_mod.F90 — thin wrapper that selects the mesh backend at compile time.
! CMake defines exactly one of the following:
!   M3DC1_MESH_IS_SCOREC   → use scorec_mesh_mod
!   M3DC1_MESH_IS_BASIC    → use basic_mesh_mod

module mesh_mod
#ifdef M3DC1_MESH_IS_SCOREC
  use scorec_mesh_mod
#elif defined(M3DC1_MESH_IS_BASIC)
  use basic_mesh_mod
#else
  use scorec_mesh_mod
#endif
end module mesh_mod
