cmake_minimum_required(VERSION 4.1)
project(m3dc1 LANGUAGES Fortran C CXX)

option(M3DC1_ENABLE_OPENMP "Enable OpenMP support" OFF)
option(M3DC1_ENABLE_COMPLEX "Build complex version" OFF)
option(M3DC1_ENABLE_3D "Build 3D version" OFF)
option(M3DC1_ENABLE_ST "Build stellarator version " OFF)
option(M3DC1_ENABLE_PARTICLE "Enable kinetic particle module" OFF)
option(M3DC1_ENABLE_GPU "Use GPU/OpenACC Fortran kernels" OFF)
option(M3DC1_ENABLE_ADAS "Enable ADAS radiation data support" OFF)
option(M3DC1_DOWNLOAD_ADAS_DATA "Download OpenADAS ADF11 data as part of the default build" OFF)
set(M3DC1_ADAS_DATA_DIR "${CMAKE_BINARY_DIR}/adas/adf11"
  CACHE PATH "Path to OpenADAS ADF11 data directory")

cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

add_subdirectory(m3dc1_scorec)
add_subdirectory(unstructured)
