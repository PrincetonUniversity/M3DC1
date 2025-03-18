# - Try to find Simmetrix SimModSuite
# Once done this will define
#  SIMMODSUITE_FOUND - System has SimModSuite
#  SIMMODSUITE_INCLUDE_DIR - The SimModSuite include directories
#  SIMMODSUITE_LIBS - The libraries needed to use SimModSuite
#  SIMMODSUITE_<library>_FOUND - System has <library>
#  SIMMODSUITE_MAJOR_VERSION - the leading integer of the version string
#  SIMMODSUITE_MINOR_VERSION - the date code from the version string
#
# Based on input variables:
#  SIM_MPI
#  SIMMETRIX_LIB_DIR
#  SIMMETRIX_INCLUDE_DIR
# And environment variable:
#  CMAKE_PREFIX_PATH
#
# This implementation assumes a simmetrix install has the following structure
# VERSION/
#         include/*.h
#         lib/ARCHOS/*.a

macro(simmetrixLibCheck libs isRequired)
  foreach(lib ${libs}) 
    unset(simmetrixlib CACHE)
    find_library(simmetrixlib "${lib}" PATHS ${SIMMETRIX_LIB_DIR})
    if(simmetrixlib MATCHES "^simmetrixlib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "SIMMODSUITE library ${lib} not found in ${SIMMETRIX_LIB_DIR}")
      else()
        message("SIMMODSUITE library ${lib} not found in ${SIMMETRIX_LIB_DIR}")
      endif()
    else()
      set("SIMMODSUITE_${lib}_FOUND" TRUE CACHE INTERNAL "SIMMODSUITE library present")
      set(SIMMODSUITE_LIBS ${SIMMODSUITE_LIBS} ${simmetrixlib})
    endif()
  endforeach()
endmacro(simmetrixLibCheck)

set(SIMMODSUITE_LIB_NAMES
  SimLicense  #-- valid for PPPL
  SimPartitionedMesh #-mpi
  SimMeshing
  SimMeshTools
  SimModel
  SimPartitionWrapper #-${SIM_MPI}
  SimAdvMeshing
  #tirpc #for Stellar, MIT Rhel7
  #SimField -- not valid for PPPL
)

simmetrixLibCheck("${SIMMODSUITE_LIB_NAMES}" TRUE)

find_path(SIMMODSUITE_INCLUDE_DIR 
  NAMES MeshSim.h
  PATHS ${SIMMETRIX_INCLUDE_DIR})
if(NOT EXISTS "${SIMMODSUITE_INCLUDE_DIR}")
  message(FATAL_ERROR "SIMMODSUITE include dir not found")
endif()

string(REGEX REPLACE 
  "/include$" "" 
  SIMMODSUITE_INSTALL_DIR
  "${SIMMODSUITE_INCLUDE_DIR}")

set(SIMMODSUITE_LIBRARIES ${SIMMODSUITE_LIBS} ${TIRPC_LIBRARY})
set(SIMMODSUITE_INCLUDE_DIRS ${SIMMODSUITE_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SIMMODSUITE DEFAULT_MSG
                                  SIMMODSUITE_LIBS SIMMODSUITE_INCLUDE_DIR)

mark_as_advanced(SIMMODSUITE_INCLUDE_DIR SIMMODSUITE_LIBS ${TIRPC_LIBRARY})

set(SIMMODSUITE_LINK_LIBS "")
foreach(lib ${SIMMODSUITE_LIB_NAMES})
  set(SIMMODSUITE_LINK_LIBS "${SIMMODSUITE_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${SIMMODSUITE_INSTALL_DIR}")
set(includedir "${SIMMODSUITE_INCLUDE_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libSimmetrix.pc.in"
  "${CMAKE_BINARY_DIR}/libSimmetrix.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libSimmetrix.pc" DESTINATION lib/pkgconfig)
