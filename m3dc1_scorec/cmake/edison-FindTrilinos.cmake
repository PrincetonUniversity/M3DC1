# - Try to find Trilinos libraries
# Once done this will define
#  TRILINOS_FOUND - System has TRILINOS
#  TRILINOS_INCLUDE_DIRS - The TRILINOS include directories
#  TRILINOS_LIBRARIES - The libraries needed to use TRILINOS
#  TRILINOS_DEFINITIONS - Compiler switches required for using TRILINOS
#
# This implementation assumes a TRILINOS install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(trilinosLibCheck libs isRequired)
  foreach(lib ${libs}) 
    unset(trilinoslib CACHE)
    find_library(trilinoslib "${lib}" PATHS ${TRILINOS_LIB_DIR})
    if(trilinoslib MATCHES "^trilinoslib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "TRILINOS library ${lib} not found in ${TRILINOS_LIB_DIR}")
      else()
        message("TRILINOS library ${lib} not found in ${TRILINOS_LIB_DIR}")
      endif()
    else()
      set("TRILINOS_${lib}_FOUND" TRUE CACHE INTERNAL "TRILINOS library present")
      set(TRILINOS_LIBS ${TRILINOS_LIBS} ${trilinoslib})
    endif()
  endforeach()
endmacro(trilinosLibCheck)

find_library(NETCDF_LIBRARY netcdf)
if (NOT EXISTS "${NETCDF_LIBRARY}")
  message(FATAL ERROR "NETCDF library not found")
endif()

find_library(STDCPP_LIBRARY stdc++)
if (NOT EXISTS "${STDCPP_LIBRARY}")
  message(FATAL ERROR "stdc++ library not found")
endif()

set(TRILINOS_LIBS "")
set(TRILINOS_LIB_NAMES
  amesos
  tpetra
  kokkosnodeapi
  tpi
  aztecoo
  epetra
  sacado
  teuchosparameterlist
  teuchoscomm
  teuchoscore
  teuchosnumerics
  teuchosremainder
)

trilinosLibCheck("${TRILINOS_LIB_NAMES}" TRUE)

find_path(TRILINOS_INCLUDE_DIR 
  NAMES Epetra_config.h
  PATHS ${TRILINOS_INCLUDE_DIR})
if(NOT EXISTS "${TRILINOS_INCLUDE_DIR}")
  message(FATAL_ERROR "TRILINOS include dir not found")
endif()

set(TRILINOS_LIBRARIES ${TRILINOS_LIBS} ${NETCDF_LIBRARY} ${STDCPP_LIBRARY})
set(TRILINOS_INCLUDE_DIRS ${TRILINOS_INCLUDE_DIR})

string(REGEX REPLACE 
  "/include$" "" 
  TRILINOS_INSTALL_DIR
  "${TRILINOS_INCLUDE_DIR}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TRILINOS  DEFAULT_MSG
                                  TRILINOS_LIBS TRILINOS_INCLUDE_DIR)

mark_as_advanced(TRILINOS_INCLUDE_DIR TRILINOS_LIBS ${NETCDF_LIBRARY})

set(TRILINOS_LINK_LIBS "")
foreach(lib ${TRILINOS_LIB_NAMES})
  set(TRILINOS_LINK_LIBS "${TRILINOS_LINK_LIBS} -l${lib}")
endforeach()
set(TRILINOS_LINK_LIBS "-lstdc++ ${TRILINOS_LINK_LIBS}")
message("TRILINOS_LINK_LIBS:" ${TRILINOS_LINK_LIBS})

#pkgconfig  
set(prefix "${TRILINOS_INSTALL_DIR}")
set(includedir "${TRILINOS_INCLUDE_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libTrilinos.pc.in"
  "${CMAKE_BINARY_DIR}/libTrilinos.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libTrilinos.pc" DESTINATION lib/pkgconfig)

