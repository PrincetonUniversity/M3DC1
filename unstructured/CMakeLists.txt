cmake_minimum_required(VERSION 3.12)

# set(CMAKE_DISABLE_SOURCE_CHANGES ON)
# set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(m3dc1_unstructured LANGUAGES Fortran C CXX)


if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# enable_language(Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_Fortran_PREPROCESS "ON")

# option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

set(SOURCES
    scorec_mesh.f90
    scorec_vector.f90
    scorec_matrix.f90
    scorec_adapt.f90
    PETScInterface.cpp
    physical_mesh.f90
    read_vmec.f90
    fftw_fortran.c
    read_namelist.cpp
    gsl_wrapper.c
    # random.c
    spline.f90
    kprad.f90
    math.f90
    read_ascii.f90
    interpolate.f90
    control.f90
    iterdb.f90
    read_gyro.f90
    read_neo.f90
    region.f90
    radiation.f90
    element.f90
    field.f90
    nintegrate_mod.f90
    M3Dmodules.f90
    resistive_wall.f90
    m3dc1_nint.f90
    boundary.f90
    gyroviscosity.f90
    bootstrap.f90
    metricterms_new.f90
    two_fluid.f90
    harned_mikic.f90
    biharmonic.f90
    electric_field.f90
    pellet.f90
    parallel_heat_flux.f90
    temperature_plots.f90
    electrostatic_potential.f90
    newvar.f90
    runaway.f90
    diagnostics.f90
    read_schaffer_field.f90
    neutral_beam.f90
    coils.f90
    coil_sets.f90
    model.f90
    fit_magnetics.f90
    wall.f90
    init_common.f90
    gradshafranov.f90
    rmp.f90
    kprad_m3dc1.f90
    readgeqdsk.f90
    read_jsolver.f90
    read_dskbal.f90
    init_cyl.f90
    init_rwm.f90
    init_solovev.f90
    init_circle.f90
    init_basicj.f90
    init_tilt.f90
    init_taylor.f90
    init_force_free.f90
    init_gem.f90
    init_wave.f90
    init_gmode.f90
    init_strauss.f90
    init_mri.f90
    init_rotating_cylinder.f90
    init_eqdsk.f90
    init_dskbal.f90
    init_jsolver.f90
    init_3dwave.f90
    init_3ddiffusion.f90
    init_frs.f90
    init_ftz.f90
    init_eigen.f90
    init_intkink.f90
    init_lz.f90
    init_kstar.f90
    init_basicq.f90
    init_vmec.f90
    transport.f90
    auxiliary_fields.f90
    time_step_split.f90
    time_step_unsplit.f90
    hdf5_output.f90
    particle.f90
    time_step.f90
    output.f90
    restart_hdf5.f90
    error_estimate.f90
    adapt.f90
    input.f90
    ludef_t.f90
    init_conds.f90
    # get_pc_skip_count.f90
    hypervisc.f90
    geometry.f90
    newpar.f90
    )

set(BIN_NAME m3dc1)

# find_package(PETSc)
# if(PETSC_FOUND AND USEPETSC)
#   set(PETSC_SOURCES
#       basic_mesh.f90
#       petsc_vector.f90
#       petsc_matrix.f90)
#   set(SOURCES "${SOURCES};${PETSC_SOURCES}")
#   list(APPEND DEFS -DUSEPETSC -Dmesh_mod=basic_mesh_mod -Dvector_mod=petsc_vector_mod -Dmatrix_mod=petsc_matrix_mod)
#   set(POSTFIX "${POSTFIX}-petsc")
# elseif(PETSC_FOUND)
#   include(${CMAKE_SOURCE_DIR}/cmake/petsc.cmake)
#   petsc_get_version()
#   message(STATUS "Petsc version detected as : ${PETSC_VERSION}")
#   list(APPEND DEFS -DPETSC_FORTRAN_PETSCTRUTH_INT -DPETSC_VERSION=${PETSC_VERSION})
# endif()
find_package(PkgConfig REQUIRED)
list(APPEND PKG_CONFIG_EXECUTABLE "--static")

# pkg_check_modules (PETSC REQUIRED PETSc)
pkg_search_module(PETSC PETSc IMPORTED_TARGET REQUIRED)
message(STATUS "Found PETSc ${PETSC_VERSION}")

# foreach(ITR ${PETSC_LINK_LIBRARIES})
#     if(ITR MATCHES "liblapack\.a")
#       set(LAPACK_LIBRARY ${ITR})
#     endif()
# endforeach()
# list(FILTER PETSC_LINK_LIBRARIES EXCLUDE REGEX "liblapack\.a")
# foreach(ITR ${PETSC_LINK_LIBRARIES})
#     if(ITR MATCHES "libblas\.a")
#       list(FIND PETSC_LINK_LIBRARIES ${ITR} IDX)
#     endif()
# endforeach()
# list(INSERT PETSC_LINK_LIBRARIES ${IDX} ${LAPACK_LIBRARY})
# message(STATUS "Found PETSc ${PETSC_LINK_LIBRARIES}")

list(FILTER PETSC_LINK_LIBRARIES EXCLUDE REGEX "libpthread")
set_property(TARGET PkgConfig::PETSC PROPERTY
  INTERFACE_LINK_LIBRARIES "${PETSC_LINK_LIBRARIES}")

pkg_search_module(FFTW fftw3 IMPORTED_TARGET REQUIRED)
pkg_search_module(GSL gsl IMPORTED_TARGET REQUIRED)

#find_package(TRILINOS)
# if(TRILINOS_FOUND AND USETRILINOS)
#   list(APPEND DEFS -DM3DC1_TRILINOS)
#   set(POSTFIX "${POSTFIX}-trilinos")
# endif()

find_package(Scorec)
# if(SCOREC_FOUND AND USESCOREC)
#   if(TRILINOS_FOUND AND USETRILINOS)
#     set(SCOREC_SOURCES
#         scorec_mesh.f90
#         scorec_vector.f90
#         scorec_matrix.f90)
#   endif()
#   if(PETSC_FOUND AND USESCOREC)
#     set(SCOREC_SOURCES
#         scorec_mesh.f90
#         scorec_vector.f90
#         scorec_matrix.f90
#         PETScInterface.cpp)
#   endif()
#   set(SOURCES "${SOURCES};${SCOREC_SOURCES}")
#   list(APPEND DEFS -DUSESCOREC -Dmesh_mod=scorec_mesh_mod -Dvector_mod=scorec_vector_mod -Dmatrix_mod=scorec_matrix_mod -Dmatrix_type=scorec_matrix -Dvector_type=scorec_vector)
#   if(USERW)
#     list(APPEND DEFS -DUSERW)
#     set(POSTFIX "${POSTFIX}-rw")
#   else()
#     list(APPEND DEFS -Dglobalinsertval=insertval -Dglobalentdofs=entdofs)
#   endif()
# endif()

# if(MATRIX_DUMP)
#   list(APPEND DEFS -DMATRIX_DUMP)
# endif()

find_package(Zoltan)
#find_package(PARMETIS)
#find_package(METIS)

set(HDF5_USE_STATIC_LIBRARIES TRUE)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
if (HDF5_FOUND)
    message(STATUS "HDF5 found: ${HDF5_VERSION}")
else ()
    message(FATAL_ERROR "HDF5 not found.")
endif ()

#find_package(ZLIB)
# find_package(GSL REQUIRED)
# find_package(FFTW REQUIRED)

set(M3DC1_SCOREC_LIB_NAME m3dc1_scorec)
if(ENABLE_3D)
  set(BIN_NAME "${BIN_NAME}_3d")
  list(APPEND DEFS -DUSE3D -Dvectype=real -DMAX_PTS=125)
  # set(POSTFIX "${POSTFIX}-125")
else()
  set(BIN_NAME "${BIN_NAME}_2d")
  list(APPEND DEFS -DMAX_PTS=25)
  # set(POSTFIX "${POSTFIX}-25")
  if(ENABLE_COMPLEX)
    set(BIN_NAME "${BIN_NAME}_complex")
    list(APPEND DEFS -Dvectype=complex -DUSECOMPLEX)
    set(M3DC1_SCOREC_LIB_NAME ${M3DC1_SCOREC_LIB_NAME}_complex)
    # if(PAR)
    #   set(POSTFIX "${POSTFIX}-pix")
    # else()
    #   set(POSTFIX "${POSTFIX}-complex")
    # endif()
  else()
    list(APPEND DEFS -Dvectype=real)
  endif()
endif()

# if(USEOMP)
#   set(POSTFIX "${POSTFIX}-omp")
# endif()

execute_process(COMMAND date OUTPUT_VARIABLE DATE)
string(REGEX REPLACE "\n$" "" DATE "${DATE}")
message(STATUS "Build date : ${DATE}")
execute_process(COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/release_version OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "\n$" "" VERSION "${VERSION}")
message(STATUS "Release version : ${VERSION}")
list(APPEND DEFS -DUSESCOREC -Dmesh_mod=scorec_mesh_mod -Dvector_mod=scorec_vector_mod -Dmatrix_mod=scorec_matrix_mod -Dmatrix_type=scorec_matrix -Dvector_type=scorec_vector -DPETSC_VERSION=990 -DUSEBLAS -DDATE_BUILT='${DATE}' -DRELEASE_VERSION='${VERSION}')

# if(NOT ARCH)
#   if(M3DC1_ARCH)
#     set(ARCH ${M3DC1_ARCH})
#   else()
#     set(ARCH ${HOST})
#   endif()
# endif()
# set(POSTFIX "${ARCH}${POSTFIX}")

# include(cmake/${ARCH}.cmake)

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

#todo : build depency list as we go instead of adding all at the end (trilinos, etc)

# if(COM)
#   set(BIN_NAME "${BIN_NAME}_complex")
# endif()

set_source_files_properties(
  scorec_vector.f90
  PROPERTIES
  OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${BIN_NAME}.dir/scorec_mesh.f90.o
  )

set_source_files_properties(
  scorec_matrix.f90
  PROPERTIES
  OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${BIN_NAME}.dir/scorec_vector.f90.o
  )

set_source_files_properties(
  field.f90
  PROPERTIES
  OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${BIN_NAME}.dir/scorec_vector.f90.o
  OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${BIN_NAME}.dir/scorec_matrix.f90.o
  )

set_source_files_properties(
  M3Dmodules.f90
  PROPERTIES
  OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${BIN_NAME}.dir/scorec_mesh.f90.o
  )

set_source_files_properties(
  nintegrate_mod.f90
  PROPERTIES
  OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${BIN_NAME}.dir/scorec_mesh.f90.o
  )

add_executable(${BIN_NAME} ${SOURCES})
# set_target_properties(${BIN_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN_NAME}_mods)
set_property(TARGET ${BIN_NAME} PROPERTY LINKER_LANGUAGE Fortran)
target_include_directories(${BIN_NAME} PRIVATE ${SCOREC_INCLUDE_DIRS})
target_link_libraries(${BIN_NAME} HDF5::HDF5
  ${M3DC1_SCOREC_LIB_NAME} ${SCOREC_LIBRARIES} ${ZOLTAN_LIBRARIES}
  PkgConfig::PETSC PkgConfig::FFTW PkgConfig::GSL)
# target_link_directories(m3dc1_2d ${PETSC_STATIC_LIBRARY_DIRECTORIES})
target_compile_definitions(${BIN_NAME} PRIVATE ${DEFS})
target_compile_options(${BIN_NAME} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-r8>)
install(TARGETS ${BIN_NAME} DESTINATION ".")

#find_package(PSPLINE)
#add_executable(readgato polar.f90 readgato.f90)
#target_link_libraries(readgato ${PSPLINE_LIBRARIES})

#add_executable(read_jsolver polar.f90 read_jsolver_exec.f90)
#target_link_libraries(read_jsolver ${PSPLINE_LIBRARIES})

#add_executable(a2cc readaeqdsk.f90 a2cc.f90)
#target_link_libraries(a2cc )
#target_compile_definitions(a2cc PUBLIC ${DEFS})

#add_executable(kprad kprad.f90 kprad_test.f90)

