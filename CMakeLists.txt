cmake_minimum_required(VERSION 3.15)
project(m3dc1 LANGUAGES Fortran C CXX)

option(M3DC1_ENABLE_OPENMP "Enable OpenMP support" OFF)
option(M3DC1_ENABLE_COMPLEX "Build complex version" OFF)
option(M3DC1_ENABLE_3D "Build 3D version" OFF)
option(M3DC1_ENABLE_ST "Build stellarator version " OFF)
option(M3DC1_ENABLE_PARTICLE "Enable kinetic particle module" OFF)
option(M3DC1_ENABLE_GPU "Use GPU/OpenACC Fortran kernels" OFF)
option(M3DC1_ENABLE_ADAS "Enable ADAS radiation data support" OFF)
option(M3DC1_ENABLE_TRILINOS "Enable Trilinos support" OFF)
option(M3DC1_ENABLE_PSPLINE "Enable PSPLINE support for readgato/read_jsolver" OFF)
option(M3DC1_DOWNLOAD_ADAS_DATA "Download OpenADAS ADF11 data as part of the default build" OFF)
set(M3DC1_ADAS_DATA_DIR "${CMAKE_BINARY_DIR}/adas/adf11"
  CACHE PATH "Path to OpenADAS ADF11 data directory")

cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake"
                              "${CMAKE_CURRENT_LIST_DIR}/m3dc1_scorec/cmake")
find_package(PETSc REQUIRED)

# Automatically set exe linker flags from PETSc pkg-config data.
# Construct -L<dir> flags from PETSC_LIBRARY_DIRS.
string(REPLACE ";" " -L" _m3dc1_petsc_libdirs_flags "${PETSC_LIBRARY_DIRS}")
if(_m3dc1_petsc_libdirs_flags)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${_m3dc1_petsc_libdirs_flags}")
endif()
# Append PETSc's extra LDFLAGS (which include rpath from pkg-config).
if(PETSC_LDFLAGS_OTHER)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PETSC_LDFLAGS_OTHER}")
endif()

add_subdirectory(m3dc1_scorec)
add_subdirectory(unstructured)
# add_subdirectory(skeleton/m3dc1_skeleton_big_spec_omp)
# add_subdirectory(skeleton/m3dc1_skeleton_big_spec_acc)
