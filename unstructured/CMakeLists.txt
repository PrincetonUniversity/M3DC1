list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(SOURCES
    fftw_fortran.c
    read_namelist.cpp
    gsl_wrapper.c
    random.c
    spline.f90
    kprad.f90
    math.f90
    read_ascii.f90
    interpolate.f90
    control.f90
    iterdb.f90
    read_gyro.f90
    read_neo.f90
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
    restart.f90
    init_conds.f90
    get_pc_skip_count.f90
    )
# newpar.f90

set(BIN_NAME m3dc1)

find_package(PETSc)
if(PETSC_FOUND AND USEPETSC)
  set(PETSC_SOURCES
      basic_mesh.f90
      petsc_vector.f90
      petsc_matrix.f90)
  set(SOURCES "${SOURCES};${PETSC_SOURCES}")
  list(APPEND DEFS -DUSEPETSC -Dmesh_mod=basic_mesh_mod -Dvector_mod=petsc_vector_mod -Dmatrix_mod=petsc_matrix_mod)
  set(POSTFIX "${POSTFIX}-petsc")
elseif(PETSC_FOUND)
  include(${CMAKE_SOURCE_DIR}/cmake/petsc.cmake)
  petsc_get_version()
  message(STATUS "Petsc version detected as : ${PETSC_VERSION}")
  list(APPEND DEFS -DPETSC_FORTRAN_PETSCTRUTH_INT -DPETSC_VERSION=${PETSC_VERSION})
endif()

#find_package(TRILINOS)
if(TRILINOS_FOUND AND USETRILINOS)
  list(APPEND DEFS -DM3DC1_TRILINOS)
  set(POSTFIX "${POSTFIX}-trilinos")
endif()

find_package(Scorec)
if(SCOREC_FOUND AND USESCOREC)
  if(TRILINOS_FOUND AND USETRILINOS)
    set(SCOREC_SOURCES
        scorec_mesh.f90
        scorec_vector.f90
        scorec_matrix.f90)
  endif()
  if(PETSC_FOUND AND USESCOREC)
    set(SCOREC_SOURCES
        scorec_mesh.f90
        scorec_vector.f90
        scorec_matrix.f90
        PETScInterface.cpp)
  endif()
  set(SOURCES "${SOURCES};${SCOREC_SOURCES}")
  list(APPEND DEFS -DUSESCOREC -Dmesh_mod=scorec_mesh_mod -Dvector_mod=scorec_vector_mod -Dmatrix_mod=scorec_matrix_mod -Dmatrix_type=scorec_matrix -Dvector_type=scorec_vector)
  if(USERW)
    list(APPEND DEFS -DUSERW)
    set(POSTFIX "${POSTFIX}-rw")
  else()
    list(APPEND DEFS -Dglobalinsertval=insertval -Dglobalentdofs=entdofs)
  endif()
endif()

if(MATRIX_DUMP)
  list(APPEND DEFS -DMATRIX_DUMP)
endif()

#find_package(ZOLTAN)
#find_package(PARMETIS)
#find_package(METIS)
#find_package(HDF5)
#find_package(ZLIB)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

if(USE3D)
  set(BIN_NAME "${BIN_NAME}_3d")
  list(APPEND DEFS -DUSE3D -Dvectype=real -DMAX_PTS=125)
  set(POSTFIX "${POSTFIX}-125")
else()
  set(BIN_NAME "${BIN_NAME}_2d")
  list(APPEND DEFS -DMAX_PTS=25)
  set(POSTFIX "${POSTFIX}-25")
  if(COM)
    list(APPEND DEFS -Dvectype=complex -DUSECOMPLEX)
    if(PAR)
      set(POSTFIX "${POSTFIX}-pix")
    else()
      set(POSTFIX "${POSTFIX}-complex")
    endif()
  else()
    list(APPEND DEFS -Dvectype=real)
  endif()
endif()

if(USEOMP)
  set(POSTFIX "${POSTFIX}-omp")
endif()

execute_process(COMMAND date OUTPUT_VARIABLE DATE)
string(REGEX REPLACE "\n$" "" DATE "${DATE}")
message(STATUS "Build date : ${DATE}")
execute_process(COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/release_version OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "\n$" "" VERSION "${VERSION}")
message(STATUS "Release version : ${VERSION}")
list(APPEND DEFS -DNOUSE -DDATE_BUILT='${DATE}' -DRELEASE_VERSION='${VERSION}')

if(NOT ARCH)
  if(M3DC1_ARCH)
    set(ARCH ${M3DC1_ARCH})
  else()
    set(ARCH ${HOST})
  endif()
endif()
set(POSTFIX "${ARCH}${POSTFIX}")

include(cmake/${ARCH}.cmake)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

add_library(m3dc1 ${SOURCES})
target_include_directories(m3dc1 PUBLIC
                           ${SCOREC_INCLUDE_DIRS}
                           ${PETSC_INCLUDE_DIRS})
#todo detect m3cdc1_scorec name
target_link_libraries(m3dc1
                      ${DEPS}
                      m3dc1_scorec
                      ${SCOREC_LIBRARIES}
                      ${PETSC_LIBRARIES})
target_compile_definitions(m3dc1 PUBLIC ${DEFS})

#todo : build depency list as we go instead of adding all at the end (trilinos, etc)

if(USETRILINOS)
  set(BIN_NAME "${BIN_NAME}_trilinos")
endif()
if(COM)
  set(BIN_NAME "${BIN_NAME}_complex")
endif()
if(OMP)
  set(BIN_NAME "${BIN_NAME}_omp")
endif()

add_executable(${BIN_NAME} newpar.f90)
target_include_directories(${BIN_NAME} PUBLIC
                           ${CMAKE_SOURCE_DIR}/m3dc1_scorec/api
                           ${SCOREC_INCLUDE_DIRS}
                           ${PETSC_INCLUDE_DIRS})
target_link_libraries(${BIN_NAME} m3dc1 m3dc1_scorec ${SCOREC_LIBRARIES} ${PETSC_LIBRARIES} ${FFTW_LIBRARIES} ${GSL_LIBRARIES})


#find_package(PSPLINE)
#add_executable(readgato polar.f90 readgato.f90)
#target_link_libraries(readgato ${PSPLINE_LIBRARIES})

#add_executable(read_jsolver polar.f90 read_jsolver_exec.f90)
#target_link_libraries(read_jsolver ${PSPLINE_LIBRARIES})

#add_executable(a2cc readaeqdsk.f90 a2cc.f90)
#target_link_libraries(a2cc )
#target_compile_definitions(a2cc PUBLIC ${DEFS})

#add_executable(kprad kprad.f90 kprad_test.f90)
