# - Try to find SCOREC PUMI libraries
# Once done this will define
#  SCOREC_FOUND - System has SCOREC
#  SCOREC_INCLUDE_DIRS - The SCOREC include directories
#  SCOREC_LIBRARIES - The libraries needed to use SCOREC
#  SCOREC_DEFINITIONS - Compiler switches required for using SCOREC
#
# This implementation assumes a SCOREC install has the following structure
# VERSION/
#         include/*.h
#         lib/*.a

macro(scorecLibCheck libs isRequired)
  foreach(lib ${libs})
    unset(scoreclib CACHE)
    find_library(scoreclib "${lib}" PATHS ${SCOREC_LIB_DIR})
    if(scoreclib MATCHES "^scoreclib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "SCOREC library ${lib} not found in ${SCOREC_LIB_DIR}")
      else()
        message("SCOREC library ${lib} not found in ${SCOREC_LIB_DIR}")
      endif()
    else()
      set("SCOREC_${lib}_FOUND" TRUE CACHE INTERNAL "SCOREC library present")
      set(SCOREC_LIBS ${SCOREC_LIBS} ${scoreclib})
    endif()
  endforeach()
endmacro(scorecLibCheck)

find_package(ZOLTAN REQUIRED)
find_package(PARMETIS REQUIRED)

set(SCOREC_LIBS "")
if(ENABLE_SIMMETRIX)
set(SCOREC_LIB_NAMES
  pumi
  ma
  mds
  apf
  apf_sim
  apf_zoltan
  parma
  dsp
  gmi
  gmi_sim
  mth
  pcu
  sam
  spr
  crv
  lion
  ph 
  size
  )
else()
set(SCOREC_LIB_NAMES
pumi
ma
crv
ph
sam
spr
apf_zoltan
parma
mds
apf
lion
gmi
mth
pcu
  )
endif()
scorecLibCheck("${SCOREC_LIB_NAMES}" TRUE)

find_path(SCOREC_INCLUDE_DIR 
  NAMES apf.h PCU.h ma.h 
  PATHS ${SCOREC_INCLUDE_DIR})
if(NOT EXISTS "${SCOREC_INCLUDE_DIR}")
  message(FATAL_ERROR "SCOREC include dir not found")
endif()

string(REGEX REPLACE 
  "/include$" "" 
  SCOREC_INSTALL_DIR
  "${SCOREC_INCLUDE_DIR}")

set(SCOREC_LIBRARIES ${SCOREC_LIBS} ${ZOLTAN_LIBRARIES} ${PARMETIS_LIBRARIES})
set(SCOREC_INCLUDE_DIRS ${SCOREC_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments
find_package_handle_standard_args(SCOREC  DEFAULT_MSG
                                  SCOREC_LIBRARIES SCOREC_INCLUDE_DIR)

mark_as_advanced(SCOREC_INCLUDE_DIR SCOREC_LIBRARIES)

set(SCOREC_LINK_LIBS "")
foreach(lib ${SCOREC_LIB_NAMES})
  set(SCOREC_LINK_LIBS "${SCOREC_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${SCOREC_INSTALL_DIR}")
set(includedir "${SCOREC_INCLUDE_DIR}")
configure_file(
  "${CMAKE_SOURCE_DIR}/cmake/libScorec.pc.in"
  "${CMAKE_BINARY_DIR}/libScorec.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libScorec.pc" DESTINATION lib/pkgconfig)

