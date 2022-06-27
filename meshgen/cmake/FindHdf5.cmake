# - Try to find HDF5 libraries
# Once done this will define
#  HDF5_FOUND - System has HDF5
#  HDF5_INCLUDE_DIRS - The HDF5 include directories
#  HDF5_LIBRARIES - The libraries needed to use HDF5
#  HDF5_DEFINITIONS - Compiler switches required for using HDF5
#
# This implementation assumes a HDF5 install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(hdf5LibCheck libs isRequired)
  foreach(lib ${libs}) 
    unset(hdf5lib CACHE)
    find_library(hdf5lib "${lib}" HINTS ${HDF5_LIB_DIR})
    if(hdf5lib MATCHES "^hdf5lib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "HDF5 library ${lib} not found in ${HDF5_LIB_DIR}")
      else()
        message("HDF5 library ${lib} not found in ${HDF5_LIB_DIR}")
      endif()
    else()
      set("HDF5_${lib}_FOUND" TRUE CACHE INTERNAL "HDF5 library present")
      set(HDF5_LIBS ${HDF5_LIBS} ${hdf5lib})
    endif()
  endforeach()
endmacro(hdf5LibCheck)

set(HDF5_LIBS "")
set(HDF5_LIB_NAMES
 hdf5hl_fortran
 hdf5_fortran
 hdf5_hl
 hdf5
 z)

hdf5LibCheck("${HDF5_LIB_NAMES}" TRUE)

find_path(HDF5_INCLUDE_DIR 
  NAMES hdf5.h
  PATHS ${HDF5_INCLUDE_DIR})
if(NOT EXISTS "${HDF5_INCLUDE_DIR}")
  message(FATAL_ERROR "HDF5 include dir not found")
endif()

set(HDF5_LIBRARIES ${HDF5_LIBS} )
set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR} )

string(REGEX REPLACE 
  "/include$" "" 
  HDF5_INSTALL_DIR
  "${HDF5_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(HDF5  DEFAULT_MSG
                                  HDF5_LIBS HDF5_INCLUDE_DIR)

mark_as_advanced(HDF5_INCLUDE_DIR HDF5_LIBS)

set(HDF5_LINK_LIBS "")
foreach(lib ${HDF5_LIB_NAMES})
  set(HDF5_LINK_LIBS "${HDF5_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${HDF5_INSTALL_DIR}")
set(includedir "${HDF5_INCLUDE_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libHdf5.pc.in"
  "${CMAKE_BINARY_DIR}/libHdf5.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libHdf5.pc" DESTINATION  lib/pkgconfig)

