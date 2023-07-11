# - Try to find SIMMETRIX libraries
# Once done this will define
#  SIMMETRIX_FOUND - System has SIMMETRIX
#  SIMMETRIX_INCLUDE_DIRS - The SIMMETRIX include directories
#  SIMMETRIX_LIBRARIES - The libraries needed to use SIMMETRIX
#  SIMMETRIX_DEFINITIONS - Compiler switches required for using SIMMETRIX
#
# This implementation assumes a SIMMETRIX install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(simmetrixLibCheck libs isRequired)
  foreach(lib ${libs}) 
    unset(simmetrixlib CACHE)
    find_library(simmetrixlib "${lib}" PATHS ${SIMMETRIX_LIB_DIR})
    if(simmetrixlib MATCHES "^simmetrixlib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "SIMMETRIX library ${lib} not found in ${SIMMETRIX_LIB_DIR}")
      else()
        message("SIMMETRIX library ${lib} not found in ${SIMMETRIX_LIB_DIR}")
      endif()
    else()
      set("SIMMETRIX_${lib}_FOUND" TRUE CACHE INTERNAL "SIMMETRIX library present")
      set(SIMMETRIX_LIBS ${SIMMETRIX_LIBS} ${simmetrixlib})
    endif()
  endforeach()
endmacro(simmetrixLibCheck)

set(SIMMETRIX_LIBS "")
set(SIMMETRIX_LIB_NAMES
  SimLicense  #-- valid for PPPL
  SimPartitionedMesh #-mpi
  SimMeshing
  SimMeshTools
  SimModel
  SimPartitionWrapper #-${SIM_MPI}
  SimAdvMeshing
  #tirpc #for Stellar or MIT Rhel7
  #SimField -- not valid for PPPL
)

simmetrixLibCheck("${SIMMETRIX_LIB_NAMES}" TRUE)

find_path(SIMMETRIX_INCLUDE_DIR 
  NAMES MeshSim.h 
  PATHS ${SIMMETRIX_INCLUDE_DIR})
if(NOT EXISTS "${SIMMETRIX_INCLUDE_DIR}")
  message(FATAL_ERROR "SIMMETRIX include dir not found")
endif()

string(REGEX REPLACE 
  "/include$" "" 
  SIMMETRIX_INSTALL_DIR
  "${SIMMETRIX_INCLUDE_DIR}")

set(SIMMETRIX_LIBRARIES ${SIMMETRIX_LIBS} )
set(SIMMETRIX_INCLUDE_DIRS ${SIMMETRIX_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SIMMETRIX  DEFAULT_MSG
                                  SIMMETRIX_LIBS SIMMETRIX_INCLUDE_DIR)

mark_as_advanced(SIMMETRIX_INCLUDE_DIR SIMMETRIX_LIBS)

set(SIMMETRIX_LINK_LIBS "")
foreach(lib ${SIMMETRIX_LIB_NAMES})
  set(SIMMETRIX_LINK_LIBS "${SIMMETRIX_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${SIMMETRIX_INSTALL_DIR}")
set(includedir "${SIMMETRIX_INCLUDE_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libSimmetrix.pc.in"
  "${CMAKE_BINARY_DIR}/libSimmetrix.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libSimmetrix.pc" DESTINATION lib/pkgconfig)

