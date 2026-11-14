find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
  message(WARNING "python3 not found; regression tests are not registered")
  return()
endif()

find_package(SCOREC CONFIG REQUIRED)

set(M3DC1_REGTEST_NPROCS "" CACHE STRING
  "Rank count for the 2D regression tests (empty: the batch script values)")

# m3dc1_add_regtest(<name> EXE <target> MESH <mesh.smb> PARTS <n> NPROCS <n>
#                   [ARGS ...] [RESTART] [COLLAPSE_MODEL <model.dmg>]
#                   [LABELS ...] [TIMEOUT <seconds>])
function(m3dc1_add_regtest name)
  cmake_parse_arguments(REGTEST "RESTART" "EXE;MESH;PARTS;NPROCS;COLLAPSE_MODEL;TIMEOUT"
    "ARGS;LABELS" ${ARGN})
  if(NOT TARGET ${REGTEST_EXE})
    return()
  endif()

  set(parts ${REGTEST_PARTS})
  set(nprocs ${REGTEST_NPROCS})
  if(M3DC1_REGTEST_NPROCS AND parts EQUAL nprocs)
    set(parts ${M3DC1_REGTEST_NPROCS})
    set(nprocs ${M3DC1_REGTEST_NPROCS})
  endif()
  if(NOT REGTEST_TIMEOUT)
    set(REGTEST_TIMEOUT 1800)
  endif()

  set(test regtest_${name})
  set(base_dir ${CMAKE_CURRENT_SOURCE_DIR}/${name}/base)
  set(work_dir ${CMAKE_CURRENT_BINARY_DIR}/mpi_test_runs/${test})

  # Copy the inputs into the work directory.
  file(GLOB inputs CONFIGURE_DEPENDS ${base_dir}/*)
  foreach(input_file IN LISTS inputs)
    get_filename_component(file_name ${input_file} NAME)
    if(NOT file_name STREQUAL "C1ke" AND NOT file_name MATCHES "^batchjob\\.")
      configure_file(${input_file} ${work_dir}/${file_name} COPYONLY)
    endif()
  endforeach()

  mpi_test(${test}_split ${parts} split_smb ${REGTEST_MESH} part.smb ${parts})
  set(steps ${test}_split)
  if(REGTEST_RESTART)
    # Run to the restart point, then restart from it.
    foreach(run 1 2)
      add_test(NAME ${test}_input${run}
        COMMAND ${CMAKE_COMMAND} -E copy C1input.${run} C1input)
      mpi_test(${test}_run${run} ${nprocs} ${REGTEST_EXE} ${REGTEST_ARGS})
      list(APPEND steps ${test}_input${run} ${test}_run${run})
    endforeach()
  else()
    mpi_test(${test}_run ${nprocs} ${REGTEST_EXE} ${REGTEST_ARGS})
    list(APPEND steps ${test}_run)
  endif()
  if(REGTEST_COLLAPSE_MODEL)
    mpi_test(${test}_collapse ${nprocs} SCOREC::collapse
      ${REGTEST_COLLAPSE_MODEL} ts0-adapted.smb mesh.smb ${nprocs})
    list(APPEND steps ${test}_collapse)
  endif()
  add_test(NAME ${test}_compare
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare.py
            ${base_dir}/C1ke C1ke)
  list(APPEND steps ${test}_compare)

  set(previous_step "")
  foreach(step IN LISTS steps)
    set_tests_properties(${step} PROPERTIES
      WORKING_DIRECTORY ${work_dir}
      TIMEOUT ${REGTEST_TIMEOUT}
      LABELS "regression;${REGTEST_LABELS}"
      FIXTURES_SETUP ${step})
    if(previous_step)
      set_tests_properties(${step} PROPERTIES FIXTURES_REQUIRED ${previous_step})
    endif()
    set(previous_step ${step})
  endforeach()
endfunction()

# 2D
m3dc1_add_regtest(KPRAD_2D EXE m3dc1_2d MESH analytic-2K.smb PARTS 48 NPROCS 48)
m3dc1_add_regtest(KPRAD_restart EXE m3dc1_2d MESH analytic-2K.smb PARTS 48 NPROCS 48
  RESTART)
m3dc1_add_regtest(adapt EXE m3dc1_2d MESH circle-0.02-0.0-0.0-2K.smb PARTS 16 NPROCS 16
  COLLAPSE_MODEL circle-0.02-0.0-0.0.dmg)
m3dc1_add_regtest(RMP EXE m3dc1_2d_complex MESH diiid-0.02-2.5-4.0-4K.smb
  PARTS 16 NPROCS 16 ARGS -pc_factor_mat_solver_type mumps)

# 3D: ranks = planes (set in C1input) x mesh parts.
m3dc1_add_regtest(RMP_nonlin EXE m3dc1_3d MESH circle-0.10-0.0-0.0-1K.smb
  PARTS 16 NPROCS 64 ARGS -ipetsc -options_file options_bjacobi.type_superludist
  LABELS long TIMEOUT 3600)
m3dc1_add_regtest(pellet EXE m3dc1_3d MESH analytic-2K.smb
  PARTS 16 NPROCS 64 ARGS -ipetsc -options_file options_bjacobi.type_superludist
  LABELS long TIMEOUT 3600)
m3dc1_add_regtest(NCSX EXE m3dc1_3d_st MESH circle-858.smb
  PARTS 8 NPROCS 64 ARGS -ipetsc -options_file options_bjacobi1
  LABELS long TIMEOUT 3600)
