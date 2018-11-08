cmake_minimum_required(VERSION 3.8)
project(m3dc1_scorec)
enable_language(Fortran)
enable_language(C)
enable_language(CXX)

cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# prefer static libraries unless we explicitly specify no to
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

# set the build type logically
 # default is Debug if .git directory exists
 # otherwise default is Release
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#include(BuildType)

# add C++ compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
  # todo :  set cray-specific cxx flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  # todo :  check for regular xl vs bg/q xl compiler
  #         currently assuming bg/q
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qflag=w -glanglvl=extended0x")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # todo :  set intel-specific cxx flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # todo :  version checks to see which c++ standards are available
  #         and warn/error if the version is incompatible
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wextra")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # todo :  version checks to see which c++ standards are available
  #         and warn/error if the version is incompatible
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -Wall -Wextra")
endif()

# add C compiler flags
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Cray")
  # todo :  set cray-specific c flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "XL")
  # todo :  set xl-specific c flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  # todo :  set intel-specific c flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  # todo :  version checks to see which c standards are available
  #         and warn/error if the version is incompatible
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # todo :  version checks to see which c++ standards are available
  #         and warn/error if the version is incompatible
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
endif()

# add Fortran compiler flags
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray")
  # todo : set cray-specific fortran flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "XL")
  # todo : set xl-specific fortran flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  # todo : set intel-specific fortran flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g noarg_temp_created")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -Wall -cpp")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -Wall -cpp")
endif()

#TODO : only include petsc/trilinos-specific files
#       if they are present (allows removal of
#       #ifdefs around entire files)
set(SOURCES
  api/m3dc1_scorec.cc
  src/BSpline.cc
  src/Expression.cc
  src/m3dc1_ls.cc
  src/m3dc1_mesh.cc
  src/m3dc1_field.cc
  src/m3dc1_numbering.cc
  src/m3dc1_sizeField.cc
  src/PolyNomial.cc
  src/ReducedQuinticImplicit.cc
  src/CMODGeoExpression.cc
  src/m3dc1_matrix.cc
  src/m3dc1_matrix_allocate.cc
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
  include/m3dc1_model.h
  src/m3dc1_numbering.h
  include/m3dc1_slnTransfer.h
  include/PolyNomial.h
  include/ReducedQuinticImplicit.h
  include/CMODGeoExpression.h
  include/m3dc1_matrix.h
  src/m3dc1_matrix_allocate.h
  include/m3dc1_sizeField.h
  include/ReducedQuinticExplicit.h
  include/slntransferUtil.h
  include/Expression.h
  include/ReducedQuintic.h)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_scorecConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/m3dc1_scorecConfig.cmake")

set(SCOREC_LIB_NAME m3dc1_scorec)

if(ENABLE_SIMMETRIX)
  set(ENABLE_PETSC "OFF")
  set(ENABLE_TRILINOS "OFF")
  find_package(Simmetrix REQUIRED)
  list(APPEND INCLUDE_DIRS ${SIMMETRIX_INCLUDE_DIRS})
  #set(LINK_LIBS ${LINK_LIBS} ${SIMMETRIX_LIBRARIES} )
  set(SCOREC_LIB_NAME ${SCOREC_LIB_NAME}_simmetrix)
  find_package(Lapack REQUIRED)
  #set(LINK_LIBS ${LINK_LIBS} ${LAPACK_LIBRARIES} )
elseif(ENABLE_TRILINOS)
  find_package(Trilinos REQUIRED)
  find_package(Boost REQUIRED)
  #include_directories(${TRILINOS_INCLUDE_DIRS})
  #set(LINK_LIBS ${LINK_LIBS} ${TRILINOS_LIBRARIES})
  #set(LINK_LIBS ${LINK_LIBS} ${BOOST_LIBRARIES})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DM3DC1_TRILINOS")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DM3DC1_TRILINOS")
  set(SCOREC_LIB_NAME ${SCOREC_LIB_NAME}_trilinos)
else()
  set(ENABLE_PETSC "ON")
endif()

if(ENABLE_PETSC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DM3DC1_PETSC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DM3DC1_PETSC")
  find_package(PETSc REQUIRED)
  find_package(Hdf5  REQUIRED)
  list(APPEND LINK_LIBS ${PETSC_LIBRARIES} ) #${HDF5_LIBRARIES} )
  list(APPEND INCLUDE_DIRS ${PETSC_INCLUDE_DIRS} ) #${HDF5_INCLUDE_DIRS}
endif()

if(ENABLE_COMPLEX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPETSC_USE_COMPLEX")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPETSC_USE_COMPLEX")
  set(SCOREC_LIB_NAME ${SCOREC_LIB_NAME}_complex)
endif()

if(ENABLE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
  set(SCOREC_LIB_NAME ${SCOREC_LIB_NAME}_openmp)
endif()

find_package(Scorec REQUIRED)
list(APPEND LINK_LIBS ${SCOREC_LIBRARIES})
list(APPEND INCLUDE_DIRS ${SCOREC_INCLUDE_DIRS})

add_library(${SCOREC_LIB_NAME} STATIC ${SOURCES})
target_include_directories(${SCOREC_LIB_NAME} PUBLIC
                                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                       $<INSTALL_INTERFACE:include>
                                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/api>
                                       $<INSTALL_INTERFACE:api>
                                       ${INCLUDE_DIRS})
target_link_libraries(${SCOREC_LIB_NAME} ${LINK_LIBS})

message(STATUS "Include directories : \n\t ${INCLUDE_DIRS}")
message(STATUS "Linking against : \n\t ${LINK_LIBS}")

INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(TARGETS ${SCOREC_LIB_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

if(ENABLE_TESTING)

  add_executable(fieldOps test/fieldOps/main.cc)
  target_link_libraries(fieldOps ${LINK_LIBS} ${SCOREC_LIB_NAME})

  if(NOT ENABLE_COMPLEX)
    add_executable(aggOrdering test/aggOrdering/main.cc)
    target_link_libraries(aggOrdering ${LINK_LIBS})
    target_link_libraries(aggOrdering ${SCOREC_LIB_NAME})
  endif()

  add_executable(field2vec2field test/f2v2f/main.cc)
  target_link_libraries(field2vec2field ${LINK_LIBS} ${SCOREC_LIB_NAME})

  add_executable(allocate test/allocate/main.cc)
  target_link_libraries(allocate ${LINK_LIBS} ${SCOREC_LIB_NAME})

  add_executable(solver test/solver/main.cc)
  target_link_libraries(solver ${LINK_LIBS})
  target_link_libraries(solver ${SCOREC_LIB_NAME})
  INSTALL(TARGETS solver RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

  if(ENABLE_SIMMETRIX)
    add_executable(m3dc1_meshgen test/meshgen/m3dc1_meshgen.cc)
    target_link_libraries(m3dc1_meshgen ${LINK_LIBS})
    target_link_libraries(m3dc1_meshgen ${SCOREC_LIB_NAME})
    INSTALL(TARGETS m3dc1_meshgen RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

    add_executable(convert_polar test/convert_polar/main.cc)
    target_link_libraries(convert_polar ${LINK_LIBS})
    target_link_libraries(convert_polar ${SCOREC_LIB_NAME})
    INSTALL(TARGETS convert_polar RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
  endif()

  add_executable(show_meshcount test/create_smb/show_meshcount.cc)
  target_link_libraries(show_meshcount ${LINK_LIBS})
  target_link_libraries(show_meshcount ${SCOREC_LIB_NAME})
  INSTALL(TARGETS show_meshcount RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endif() #enable_testing

#binary distribution package
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION "1.0.1")
include(CPack)

# Settings options for testing
enable_testing()
include(CTest)
# This will be set to ON by the CTest driver script (and only by that)
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
