set(SKELETON_OMP_SOURCES
  etimes.F90
  math.F90
  element.F90
  mesh.F90
  nintegrate.F90
  field.F90
  m3dc1_data.F90
  m3dc1_nint.F90
  matdef.F90
  main.F90
)

find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(m3dc1_skeleton_omp ${SKELETON_OMP_SOURCES})
target_compile_definitions(m3dc1_skeleton_omp PRIVATE USE3D MAX_PTS=60)
target_link_libraries(m3dc1_skeleton_omp PRIVATE
  OpenMP::OpenMP_Fortran
  LAPACK::LAPACK
)

# .F90 files are preprocessed automatically by CMake (capital F)

install(TARGETS m3dc1_skeleton_omp
  RUNTIME DESTINATION bin
)