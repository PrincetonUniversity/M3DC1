cmake_minimum_required(VERSION 3.8)
project(m3dc1 LANGUAGES Fortran C CXX)

# macro(checkSetParam varName isRequired)
#   if("${${varName}}" STREQUAL "") # if varName is empty
#     set(${varName} $ENV{${varName}})
#     if(("${${varName}}" STREQUAL "") AND ${isRequired}) # if ENV{varName} is empty also
#       message(FATAL_ERROR
#         "${varName} is not defined or in the environment, please specify ${varName}")
#     endif()
#   endif()
# endmacro(checkSetParam varName)

option(ENABLE_OPENMP "Enable OpenMP support" OFF)
option(ENABLE_COMPLEX "Build complex version" OFF)
option(ENABLE_3D "Build 3D version" OFF)
option(ENABLE_ST "Build stellarator version " OFF)
option(ENABLE_PARTICLE "Enable kinetic particle module" OFF)

cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")

# # rhel/centos place x86-64 libs in lib64 instead of lib
# # this *shouldn't* effect searching the lib directory
# # first on other distros where the default install location
# # is lib instead of lib64 (like debian derivatives)
# set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS TRUE)

# prefer static libraries unless we explicitly specify no to
# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
# if(BUILD_SHARED_LIBS)
#   set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
# endif()

# set the build type logically
 # default is Debug if .git directory exists
 # otherwise default is Release
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
# include(BuildType)

# # add C++ compiler flags
# if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
#   # todo :  set cray-specific cxx flags
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
#   # todo :  check for regular xl vs bg/q xl compiler
#   #         currently assuming bg/q
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qflag=w -glanglvl=extended0x")
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#   # todo :  set intel-specific cxx flags
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#   # todo :  version checks to see which c++ standards are available
#   #         and warn/error if the version is incompatible
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wextra")
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#   # todo :  version checks to see which c++ standards are available
#   #         and warn/error if the version is incompatible
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -Wall -Wextra")
# endif()

# # add C compiler flags
# if("${CMAKE_C_COMPILER_ID}" STREQUAL "Cray")
#   # todo :  set cray-specific c flags
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
# elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "XL")
#   # todo :  set xl-specific c flags
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
# elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
#   # todo :  set intel-specific c flags
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
# elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
#   # todo :  version checks to see which c standards are available
#   #         and warn/error if the version is incompatible
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#   # todo :  version checks to see which c++ standards are available
#   #         and warn/error if the version is incompatible
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
# endif()

# # add Fortran compiler flags
# if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
#   # todo : set cray-specific fortran flags
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
# elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
#   # todo : set xl-specific fortran flags
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
# elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
#   # todo : set intel-specific fortran flags
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g noarg_temp_created")
# elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -Wall -cpp")
# elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang")
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -Wall -cpp")
# endif()

message(STATUS "Configuring m3dc1_scorec")
add_subdirectory(m3dc1_scorec EXCLUDE_FROM_ALL)
message(STATUS "Configuring unstructured")
add_subdirectory(unstructured)

