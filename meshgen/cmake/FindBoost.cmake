# - Try to find BOOST libraries
# Once done this will define
#  BOOST_FOUND - System has BOOST 
#  BOOST_LIBRARIES - The libraries needed to use BOOST 
#  BOOST_DEFINITIONS - Compiler switches required for using BOOST 
#
# This implementation assumes an BOOST install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(boostLibCheck libs isRequired)
  foreach(lib ${libs}) 
    unset(boostlib CACHE)
    find_library(boostlib "${lib}" HINTS ${BOOST_LIB_DIR})
    if(boostlib MATCHES "^boostlib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "BOOST library ${lib} not found in ${BOOST_LIB_DIR}")
      else()
        message("BOOST library ${lib} not found in ${BOOST_LIB_DIR}")
      endif()
    else()
      set("BOOST_${lib}_FOUND" TRUE CACHE INTERNAL "BOOST library present")
      set(BOOST_LIBS ${BOOST_LIBS} ${boostlib})
    endif()
  endforeach()
endmacro(boostLibCheck)

set(BOOST_LIBS "")
set(BOOST_LIB_NAMES
  boost_mpi #the dependent libraries should follow
)

boostLibCheck("${BOOST_LIB_NAMES}" TRUE)

set(BOOST_LIBRARIES ${BOOST_LIBS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(BOOST  DEFAULT_MSG
                                  BOOST_LIBS)

mark_as_advanced(BOOST_LIBS)

set(BOOST_LINK_LIBS "")
foreach(lib ${BOOST_LIB_NAMES})
  set(BOOST_LINK_LIBS "${BOOST_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${BOOST_INSTALL_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libBoost.pc.in"
  "${CMAKE_BINARY_DIR}/libBoost.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libBoost.pc" DESTINATION  lib/pkgconfig)

