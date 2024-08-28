# - Try to find PETSc and SuperLU libraries
# Once done this will define
#  PETSC_FOUND - System has PETSc
#  PETSC_INCLUDE_DIRS - The PETSC include directories
#  PETSC_LIBRARIES - The libraries needed to use PETSC
#  PETSC_DEFINITIONS - Compiler switches required for using PETSC
#
# This implementation assumes a PETSC install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(petscLibCheck libs isRequired)
  foreach(lib ${libs}) 
    unset(petsclib CACHE)
    find_library(petsclib "${lib}" PATHS ${PETSC_LIB_DIR})
    if(petsclib MATCHES "^petsclib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "PETSC library ${lib} not found in ${PETSC_LIB_DIR}")
      else()
        message("PETSC library ${lib} not found in ${PETSC_LIB_DIR}")
      endif()
    else()
      set("PETSC_${lib}_FOUND" TRUE CACHE INTERNAL "PETSC library present")
      set(PETSC_LIBS ${PETSC_LIBS} ${petsclib})
    endif()
  endforeach()
endmacro(petscLibCheck)

set(PETSC_LIBS "")
set(PETSC_LIB_NAMES
  petsc
  parmetis
  metis
)

petscLibCheck("${PETSC_LIB_NAMES}" TRUE)

find_path(PETSC_INCLUDE_DIR 
  NAMES petsc.h 
  PATHS ${PETSC_INCLUDE_DIR})
if(NOT EXISTS "${PETSC_INCLUDE_DIR}")
  message(FATAL_ERROR "PETSC include dir not found")
endif()

set(PETSC_LIBRARIES ${PETSC_LIBS} ${GFORTRAN_LIBRARY})
set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR} )

# if (NOT EXISTS "${MPI_DIR}")
#   set(MPI_DIR "/usr/local/openmpi/latest")
# endif()
find_package(MPI REQUIRED)
set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIRS};${MPI_C_INCLUDE_DIRS})
# set(PETSC_LIBRARIES ${PETSC_LIBRARIES} ${MPI_DIR}/lib/libmpi.a ${MPI_DIR}/lib/libmpi_f90.a ${MPI_DIR}/lib/libmpi_f77.a ${MPI_DIR}/lib/libmpi_cxx.a)

# if (NOT EXISTS "${GCC_DIR}")
#   set (GCC_DIR "/usr/lib/gcc/x86_64-linux-gnu/4.4.5")
# endif()
# find_library(${GCC_DIR} gfortran)

# set(PETSC_LIBRARIES ${PETSC_LIBRARIES} ${GCC_DIR}/libgcc_s.a ${GCC_DIR}/libgfortran.a ${GCC_DIR}/libstdc++.a)

string(REGEX REPLACE 
  "/include$" "" 
  PETSC_INSTALL_DIR
  "${PETSC_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PETSC  DEFAULT_MSG
                                  PETSC_LIBS PETSC_INCLUDE_DIR)

mark_as_advanced(PETSC_INCLUDE_DIR PETSC_LIBS)

set(PETSC_LINK_LIBS "")
foreach(lib ${PETSC_LIB_NAMES})
  set(PETSC_LINK_LIBS "${PETSC_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${PETSC_INSTALL_DIR}")
set(includedir "${PETSC_INCLUDE_DIR}")
configure_file(
  "cmake/libPetsc.pc.in"
  "${CMAKE_BINARY_DIR}/libPetsc.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libPetsc.pc" DESTINATION lib/pkgconfig)

