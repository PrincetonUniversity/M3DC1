# Builds libm3dc1_cudss: the reusable cuDSS block-Jacobi solver
# (petsc_cudss_solve() + cudss_block) consumed by m3dc1_scorec via
# matrix_solve::solve_cudss.
#
# Normally entered via add_subdirectory() from the repo root (ENABLE_CUDSS=ON)
# or from m3dc1_scorec; can also be configured standalone.  Only the library
# is built here -- the file-based benchmark executables (petsc_cudss_solve,
# petsc_cudss_solve_mg) are still built by build.sbatch, which remains the
# self-contained path on Perlmutter.
cmake_minimum_required(VERSION 3.18)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(m3dc1_cudss LANGUAGES C CXX CUDA)
else()
  enable_language(CUDA)
endif()

if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80) # Perlmutter A100, matching build.sbatch
endif()

find_package(CUDAToolkit REQUIRED)
find_package(MPI REQUIRED)

# Compile-only PETSc dependency (headers; libs come from the final link).
# Same PETSC_INCLUDE_DIR cache variable the m3dc1_scorec build uses.
find_path(PETSC_INCLUDE_DIR NAMES petsc.h
  HINTS ${PETSC_INCLUDE_DIR}
        $ENV{PETSC_DIR}/include
        $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include)
if (NOT EXISTS "${PETSC_INCLUDE_DIR}")
  message(FATAL_ERROR "petsc.h not found; set -DPETSC_INCLUDE_DIR=<petsc>/include")
endif()

# cuDSS: defaults to the copy vendored in deps/cudss
set(CUDSS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/deps/cudss CACHE PATH
    "cuDSS install directory (contains include/ and lib/)")
find_library(CUDSS_LIBRARY NAMES cudss libcudss.so.0 PATHS ${CUDSS_DIR}/lib)
if (NOT EXISTS "${CUDSS_LIBRARY}")
  message(FATAL_ERROR "cuDSS library not found under ${CUDSS_DIR}/lib (set CUDSS_DIR)")
endif()

# cudss_block.cu includes nccl.h (MGMN comm-layer types)
find_path(NCCL_INCLUDE_DIR nccl.h HINTS $ENV{NCCL_DIR}/include $ENV{NCCL_HOME}/include)
find_library(NCCL_LIBRARY nccl HINTS $ENV{NCCL_DIR}/lib $ENV{NCCL_HOME}/lib)
if (NOT NCCL_INCLUDE_DIR OR NOT NCCL_LIBRARY)
  message(FATAL_ERROR "NCCL not found (cudss_block.cu includes nccl.h); set NCCL_DIR")
endif()

add_library(m3dc1_cudss STATIC
  src/petsc_cudss_solve.c
  src/cudss_block.cu)

# library build: compile out the file-based benchmark main()
set_source_files_properties(src/petsc_cudss_solve.c
  PROPERTIES COMPILE_DEFINITIONS "CUDSS_NO_MAIN")

target_include_directories(m3dc1_cudss
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/src
  PRIVATE ${CUDSS_DIR}/include ${NCCL_INCLUDE_DIR} ${PETSC_INCLUDE_DIR})

# consumers (m3dc1_scorec) compile their cuDSS call sites behind this define
target_compile_definitions(m3dc1_cudss INTERFACE M3DC1_CUDSS)

target_link_libraries(m3dc1_cudss
  PUBLIC  ${CUDSS_LIBRARY} ${NCCL_LIBRARY}
          CUDA::cudart CUDA::cusparse CUDA::cublas
  PRIVATE MPI::MPI_C)

install(TARGETS m3dc1_cudss ARCHIVE DESTINATION lib)
install(FILES src/petsc_cudss_solve.h src/cudss_block.h DESTINATION include)
