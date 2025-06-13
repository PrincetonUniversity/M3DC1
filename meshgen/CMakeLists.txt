#This is the top M3DC1_SCOREC CMakeList File for the Build

#Setting Version Number, Project Name
cmake_minimum_required (VERSION 3.0)
project (m3dc1_scorec)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
#get_filename_component (default_prefix ".." ABSOLUTE)
#set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
#      "Choose the installation directory; by default it installs in the NORMA directory."
#      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message (Fortran_COMPILER_NAME = ${Fortran_COMPILER_NAME})

#unless building shared libs, then select static libs 
# if both static and shared libs are available 
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so") 
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

#Settings options for testing
enable_testing()
include(CTest)
#This will be set to ON by the CTest driver script (and only by that)
option(ENABLE_TESTING "Build for CTest" OFF)
set(MPIRUN "mpirun"
    CACHE string 
    "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
    CACHE string 
    "the command line flag to give process count to MPIRUN")

#Doxygen generation system
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
               ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_scorecConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/m3dc1_scorecConfig.cmake")

#Source and header files
set(SOURCES
  api/m3dc1_scorec.cc
  src/BSpline.cc
  src/Expression.cc
  src/m3dc1_ls.cc
  src/m3dc1_mesh.cc
  src/m3dc1_ghost.cc
  src/m3dc1_sizeField.cc
  src/PolyNomial.cc
  src/ReducedQuinticImplicit.cc
  src/CMODGeoExpression.cc
  src/m3dc1_matrix.cc
  src/m3dc1_model.cc
  src/m3dc1_slnTransfer.cc
  src/ReducedQuinticExplicit.cc
  src/slntransferUtil.cc 
  fortran/c1_12_element.F90
  fortran/explicit_element.F90
  fortran/explicit_element_interface.F90)

set(HEADERS
  api/name_convert.h
  api/m3dc1_scorec.h 
  include/BSpline.h
  include/m3dc1_field.h
  include/m3dc1_ls.h
  include/m3dc1_mesh.h
  include/m3dc1_ghost.h
  include/m3dc1_model.h
  include/m3dc1_slnTransfer.h
  include/PolyNomial.h
  include/ReducedQuinticImplicit.h
  include/CMODGeoExpression.h  
  include/m3dc1_matrix.h
  include/m3dc1_sizeField.h
  include/ReducedQuinticExplicit.h
  include/slntransferUtil.h
  include/Expression.h

  include/ReducedQuintic.h)

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/api)
set(TARGET_LIB_NAME m3dc1_scorec)

find_package(Scorec QUIET REQUIRED)
include_directories(${SCOREC_INCLUDE_DIRS})
set(DEP_LIBS ${DEP_LIBS} ${SCOREC_LIBRARIES} )

if (ENABLE_PPPL)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPPPL")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPPPL")
endif()

if (ENABLE_MIT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMIT")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMIT")
endif()

if (ENABLE_PU)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSTELLAR")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTELLAR")
endif()

if (ENABLE_SDUMONT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSDUMONT")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSDUMONT")
endif()

set(ENABLE_SIMMETRIX "ON")
set(ENABLE_PETSC "OFF")
set(ENABLE_TRILINOS "OFF")

if (ENABLE_SCOREC)
  find_package(Simmetrix QUIET REQUIRED)
else()
  if (ENABLE_SIMLICENSE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLICENSE")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLICENSE")
  endif()
  find_package(SimModSuite QUIET REQUIRED)
endif()

include_directories(${SIMMODSUITE_INCLUDE_DIR})
set(DEP_LIBS ${DEP_LIBS} ${SIMMODSUITE_LIBS} )
set(TARGET_LIB_NAME ${TARGET_LIB_NAME}_simmetrix)
find_package(Lapack QUIET REQUIRED)
set(DEP_LIBS ${DEP_LIBS} ${LAPACK_LIBRARIES} )

  add_executable(m3dc1_meshgen test/m3dc1_meshgen.cc)
  target_include_directories(m3dc1_meshgen PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(m3dc1_meshgen ${DEP_LIBS})
  target_link_libraries(m3dc1_meshgen ${TARGET_LIB_NAME})
  INSTALL(TARGETS m3dc1_meshgen RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
  
  add_executable(m3dc1_mfmgen test/multi_meshgen.cc)
  target_include_directories(m3dc1_mfmgen PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(m3dc1_mfmgen ${DEP_LIBS})
  target_link_libraries(m3dc1_mfmgen ${TARGET_LIB_NAME})
  INSTALL(TARGETS m3dc1_mfmgen RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

  add_executable(polar_meshgen test/polar_meshgen.cc)
  target_include_directories(polar_meshgen PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(polar_meshgen ${DEP_LIBS})
  target_link_libraries(polar_meshgen ${TARGET_LIB_NAME})
  INSTALL(TARGETS polar_meshgen RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

  add_executable(convert_sms test/convert_sms.cc)
  target_include_directories(convert_sms PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(convert_sms ${DEP_LIBS})
  target_link_libraries(convert_sms ${TARGET_LIB_NAME})
  INSTALL(TARGETS convert_sms RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

  add_executable(simToM3dc1 test/simToM3dc1.cc)
  target_include_directories(simToM3dc1 PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(simToM3dc1 ${DEP_LIBS})
  target_link_libraries(simToM3dc1 ${TARGET_LIB_NAME})
  INSTALL(TARGETS simToM3dc1 RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_library(${TARGET_LIB_NAME} ${SOURCES})
target_link_libraries(${TARGET_LIB_NAME} ${DEP_LIBS})

INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(TARGETS ${TARGET_LIB_NAME}
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)

#binary distribution package
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION "1.0.1")
include(CPack)

message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("CMAKE_Fortran_COMPILER_INIT = ${CMAKE_Fortran_COMPILER_INIT}")
message("CMAKE_Fortran_COMPILER_FULLPATH =
${CMAKE_Fortran_COMPILER_FULLPATH}")
message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
message("CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
message("CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
message("DEP_LIBS = ${DEP_LIBS}")
