set(M3DC1_TEST_OMP_NUM_THREADS
    "1"
    CACHE STRING "OMP_NUM_THREADS for each rank of the MPI tests")

function(set_mpi_test_environment TESTNAME)
  set_property(
    TEST ${TESTNAME}
    APPEND
    PROPERTY ENVIRONMENT "OMP_NUM_THREADS=${M3DC1_TEST_OMP_NUM_THREADS}")
endfunction()

# run tests in independent directories to avoid input/output file collisions
function(isolate_mpi_test TESTNAME NPROCS)
  set(dir ${CMAKE_CURRENT_BINARY_DIR}/mpi_test_runs/${TESTNAME})
  file(MAKE_DIRECTORY ${dir})
  set_tests_properties(${TESTNAME} PROPERTIES WORKING_DIRECTORY ${dir}
                                              PROCESSORS ${NPROCS})
endfunction()

# expand the full path of the EXE
function(resolve_test_exe EXE_VAR)
  if(TARGET ${${EXE_VAR}})
    set(${EXE_VAR}
        $<TARGET_FILE:${${EXE_VAR}}>
        PARENT_SCOPE)
  endif()
endfunction()

function(mpi_test TESTNAME PROCS EXE)
  resolve_test_exe(EXE)
  add_test(NAME ${TESTNAME}
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS}
                   ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN})
  set_mpi_test_environment(${TESTNAME})
  isolate_mpi_test(${TESTNAME} ${PROCS})
endfunction(mpi_test)
