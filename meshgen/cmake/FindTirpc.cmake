# - Try to find TIRPC libraries
# Once done this will define
#  TIRPC_FOUND - System has TIRPC 
#  TIRPC_LIBRARIES - The libraries needed to use TIRPC 
#  TIRPC_DEFINITIONS - Compiler switches required for using TIRPC 
#
# This implementation assumes a TIRPC install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(tirpcLibCheck libs isRequired)
  foreach(lib ${libs}) 
    unset(tirpclib CACHE)
    find_library(tirpclib "${lib}" HINTS ${TIRPC_LIB_DIR})
    if(tirpclib MATCHES "^tirpclib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "TIRPC library ${lib} not found in ${TIRPC_LIB_DIR}")
      else()
        message("TIRPC library ${lib} not found in ${TIRPC_LIB_DIR}")
      endif()
    else()
      set("TIRPC_${lib}_FOUND" TRUE CACHE INTERNAL "TIRPC library present")
      set(TIRPC_LIBS ${TIRPC_LIBS} ${tirpclib})
    endif()
  endforeach()
endmacro(tirpcLibCheck)

set(TIRPC_LIBS "")
set(TIRPC_LIB_NAMES
  tirpc
)

tirpcLibCheck("${TIRPC_LIB_NAMES}" TRUE)

set(TIRPC_LIBRARIES ${TIRPC_LIBS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TIRPC  DEFAULT_MSG
                                  TIRPC_LIBS)

mark_as_advanced(TIRPC_LIBS)

set(TIRPC_LINK_LIBS "")
foreach(lib ${TIRPC_LIB_NAMES})
  set(TIRPC_LINK_LIBS "${TIRPC_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${TIRPC_INSTALL_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libTirpc.pc.in"
  "${CMAKE_BINARY_DIR}/libTirpc.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libTirpc.pc" DESTINATION  lib/pkgconfig)

