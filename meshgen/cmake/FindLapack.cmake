# - Try to find LAPACK libraries
# Once done this will define
#  LAPACK_FOUND - System has LAPACK 
#  LAPACK_LIBRARIES - The libraries needed to use LAPACK 
#  LAPACK_DEFINITIONS - Compiler switches required for using LAPACK 
#
# This implementation assumes an LAPACK install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(lapackLibCheck libs isRequired)
  foreach(lib ${libs}) 
    unset(lapacklib CACHE)
    find_library(lapacklib "${lib}" HINTS ${LAPACK_LIB_DIR})
    if(lapacklib MATCHES "^lapacklib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "LAPACK library ${lib} not found in ${LAPACK_LIB_DIR}")
      else()
        message("LAPACK library ${lib} not found in ${LAPACK_LIB_DIR}")
      endif()
    else()
      set("LAPACK_${lib}_FOUND" TRUE CACHE INTERNAL "LAPACK library present")
      set(LAPACK_LIBS ${LAPACK_LIBS} ${lapacklib})
    endif()
  endforeach()
endmacro(lapackLibCheck)

set(LAPACK_LIBS "")
if (ENABLE_MIT)
set(LAPACK_LIB_NAMES
  mkl_lapack95_lp64 #lapack should come first since it references routines in the BLAS library
  mkl_blas95_lp64
)	
elseif (ENABLE_SDUMONT)
set(LAPACK_LIB_NAMES
  lapack #lapack should come first since it references routines in the BLAS library
  blas
)
else()
set(LAPACK_LIB_NAMES
  flapack #lapack should come first since it references routines in the BLAS library
  fblas
)
endif()

lapackLibCheck("${LAPACK_LIB_NAMES}" TRUE)

set(LAPACK_LIBRARIES ${LAPACK_LIBS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LAPACK  DEFAULT_MSG
                                  LAPACK_LIBS)

mark_as_advanced(LAPACK_LIBS)

set(LAPACK_LINK_LIBS "")
foreach(lib ${LAPACK_LIB_NAMES})
  set(LAPACK_LINK_LIBS "${LAPACK_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${LAPACK_INSTALL_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libLapack.pc.in"
  "${CMAKE_BINARY_DIR}/libLapack.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libLapack.pc" DESTINATION  lib/pkgconfig)

