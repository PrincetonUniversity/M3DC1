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
  pthread 
  ssl 
  crypto 
  m 
  numa 
  rt 
  nsl  
  util 
  pthread 
  dl)

petscLibCheck("${PETSC_LIB_NAMES}" TRUE)

find_path(PETSC_INCLUDE_DIR 
  NAMES petsc.h 
  PATHS ${PETSC_INCLUDE_DIR})
if(NOT EXISTS "${PETSC_INCLUDE_DIR}")
  message(FATAL_ERROR "PETSC include dir not found")
endif()

set(PETSC_LIBRARIES ${PETSC_LIBS} )
set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR} )

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
  "${CMAKE_HOME_DIRECTORY}/cmake/libPetsc.pc.in"
  "${CMAKE_BINARY_DIR}/libPetsc.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libPetsc.pc" DESTINATION lib/pkgconfig)

