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

find_library(ZOLTAN_LIBRARY zoltan)
if (NOT EXISTS "${ZOLTAN_LIBRARY}")
  message(FATAL ERROR "ZOLTAN library not found")
endif()

find_library(PARMETIS_LIBRARY parmetis)
if (NOT EXISTS "${PARMETIS_LIBRARY}")
  message(FATAL ERROR "PARMETIS library not found")
endif()

find_library(METIS_LIBRARY metis)
if (NOT EXISTS "${METIS_LIBRARY}")
  message(FATAL ERROR "METIS library not found")
endif()

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

set(SCOREC_LIBRARIES ${SCOREC_LIBS} ${ZOLTAN_LIBRARY} ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
set(SCOREC_INCLUDE_DIRS ${SCOREC_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SCOREC  DEFAULT_MSG
                                  SCOREC_LIBS SCOREC_INCLUDE_DIR)

mark_as_advanced(SCOREC_INCLUDE_DIR SCOREC_LIBS ${ZOLTAN_LIBRARY} ${PARMETIS_LIBRARY} ${METIS_LIBRARY})

set(SCOREC_LINK_LIBS "")
foreach(lib ${SCOREC_LIB_NAMES})
  set(SCOREC_LINK_LIBS "${SCOREC_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${SCOREC_INSTALL_DIR}")
set(includedir "${SCOREC_INCLUDE_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libScorec.pc.in"
  "${CMAKE_BINARY_DIR}/libScorec.pc"
  @ONLY)

INSTALL(FILES "${CMAKE_BINARY_DIR}/libScorec.pc" DESTINATION lib/pkgconfig)

