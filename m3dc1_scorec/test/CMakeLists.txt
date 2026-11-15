set(_data "${PROJECT_SOURCE_DIR}/unstructured/regtest/RMP_nonlin/base")
set(_model "${_data}/circle-0.10-0.0-0.0.txt")
set(_mesh "${_data}/circle-0.10-0.0-0.0-1K.smb")   # circle-0.10-0.0-0.0-1K0.smb
set(_mesh4 "${CMAKE_CURRENT_BINARY_DIR}/circle-4part/part.smb")

add_executable(m3dc1_test_solver solver/main.cc)
target_link_libraries(m3dc1_test_solver PRIVATE m3dc1_scorec)
target_compile_definitions(m3dc1_test_solver PRIVATE M3DC1_PETSC)

add_executable(m3dc1_test_sum_vec_planes sum_vec_planes/main.cc)
target_link_libraries(m3dc1_test_sum_vec_planes PRIVATE m3dc1_scorec)

# Four-part copy of the mesh for the parallel tests.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/circle-4part")
mpi_test(unit_split_mesh 4 split_smb ${_mesh} ${_mesh4} 4)
set_tests_properties(unit_split_mesh PROPERTIES FIXTURES_SETUP circle_4part)

mpi_test(unit_solver_serial 1 m3dc1_test_solver ${_model} ${_mesh} 1)
mpi_test(unit_solver_parallel 4 m3dc1_test_solver ${_model} ${_mesh4} 1)
set_tests_properties(unit_solver_parallel PROPERTIES FIXTURES_REQUIRED circle_4part)
mpi_test(unit_solver_planes 2 m3dc1_test_solver ${_model} ${_mesh} 2
  -pc_type bjacobi -pc_bjacobi_blocks 2 -sub_pc_type lu
  -sub_pc_factor_mat_solver_type superlu_dist -sub_ksp_type preonly
  -ksp_type fgmres -ksp_rtol 1e-9)

mpi_test(unit_sum_vec_planes 2 m3dc1_test_sum_vec_planes ${_model} ${_mesh} 2)
mpi_test(unit_sum_vec_planes_parallel 8 m3dc1_test_sum_vec_planes ${_model} ${_mesh4} 2)
set_tests_properties(unit_sum_vec_planes_parallel PROPERTIES FIXTURES_REQUIRED circle_4part)

set_tests_properties(
  unit_split_mesh unit_solver_serial unit_solver_parallel unit_solver_planes
  unit_sum_vec_planes unit_sum_vec_planes_parallel
  PROPERTIES LABELS unit TIMEOUT 300)
