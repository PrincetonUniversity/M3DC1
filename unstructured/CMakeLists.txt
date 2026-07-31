set(M3DC1_LIBRARY_SOURCES
  adapt.f90
  auxiliary_fields.f90
  basic_mesh.f90
  biharmonic.f90
  bootstrap.f90
  boundary.f90
  coil_sets.f90
  coils.f90
  control.f90
  diagnostics.f90
  electric_field.f90
  electrostatic_potential.f90
  element.f90
  error_estimate.f90
  field.f90
  fit_magnetics.f90
  fftw_fortran.c
  geometry.f90
  get_pc_skip_count.f90
  gsl_wrapper.c
  gradshafranov.f90
  gyroviscosity.f90
  harned_mikic.f90
  hdf5_output.f90
  hypervisc.f90
  init_3ddiffusion.f90
  init_3dwave.f90
  init_basicj.f90
  init_basicq.f90
  init_circle.f90
  init_common.f90
  init_conds.f90
  init_cyl.f90
  init_dskbal.f90
  init_eigen.f90
  init_eqdsk.f90
  init_force_free.f90
  init_frs.f90
  init_ftz.f90
  init_gem.f90
  init_gmode.f90
  init_intkink.f90
  init_jsolver.f90
  read_jsolver.f90
  readaeqdsk.f90
  init_kstar.f90
  init_lz.f90
  init_mri.f90
  init_rotating_cylinder.f90
  init_rwm.f90
  init_solovev.f90
  init_strauss.f90
  init_taylor.f90
  init_tilt.f90
  init_vmec.f90
  init_wave.f90
  input.f90
  interpolate.f90
  iterdb.f90
  kprad_m3dc1.f90
  kprad.f90
  ludef_t.f90
  m3dc1_nint.f90
  M3Dmodules.f90
  math.f90
  metricterms_new.f90
  model.f90
  neutral_beam.f90
  newvar.f90
  newpar.f90
  nintegrate_mod.f90
  output.f90
  parallel_heat_flux.f90
  parse_solver_options.c
  particle.f90
  pellet.f90
  physical_mesh.f90
  polar.f90
  radiation.f90
  random.c
  read_ascii.f90
  read_dskbal.f90
  read_gyro.f90
  read_neo.f90
  read_namelist.cpp
  read_schaffer_field.f90
  read_vmec.f90
  readgeqdsk.f90
  region.f90
  resistive_wall.f90
  restart_hdf5.f90
  rmp.f90
  runaway_advection.f90
  runaway.f90
  scorec_adapt.f90
  scorec_matrix.f90
  scorec_mesh.f90
  scorec_vector.f90
  signal_handler.f90
  spline.f90
  temperature_plots.f90
  time_step_split.f90
  time_step_unsplit.f90
  time_step.f90
  transport.f90
  two_fluid.f90
  wall.f90
)
if(M3DC1_ENABLE_GPU)
  list(REMOVE_ITEM M3DC1_SOURCES ludef_t.f90 metricterms_new.f90)
  list(APPEND M3DC1_SOURCES ludef_t_gpu.f90 metricterms_new_gpu.f90)
endif()

if(NOT DEFINED M3DC1_ADAS_DATA_DIR)
  set(M3DC1_ADAS_DATA_DIR "${CMAKE_BINARY_DIR}/adas/adf11"
    CACHE PATH "Path to OpenADAS ADF11 data directory")
endif()
set(M3DC1_OPENADAS_SOURCE_DIR "${CMAKE_BINARY_DIR}/adas/source")
set(M3DC1_OPENADAS_SOURCE_NAMES
  xxrptn
  xfelem
  i4unit
  i4fctn
  xxword
  xxcase
  xxslen
  xxdata_11
)
foreach(_openadas_source_name IN LISTS M3DC1_OPENADAS_SOURCE_NAMES)
  list(APPEND M3DC1_OPENADAS_SOURCES
    "${M3DC1_OPENADAS_SOURCE_DIR}/${_openadas_source_name}.f")
endforeach()
set_source_files_properties(${M3DC1_OPENADAS_SOURCES} PROPERTIES GENERATED TRUE)
set(M3DC1_OPENADAS_DATA_FILES
  "${M3DC1_ADAS_DATA_DIR}/scd85/scd85_ar.dat"
  "${M3DC1_ADAS_DATA_DIR}/scd89/scd89_b.dat"
  "${M3DC1_ADAS_DATA_DIR}/plt89/plt89_b.dat"
  "${M3DC1_ADAS_DATA_DIR}/plt89/plt89_ar.dat"
  "${M3DC1_ADAS_DATA_DIR}/scd96/scd96_he.dat"
  "${M3DC1_ADAS_DATA_DIR}/scd96/scd96_be.dat"
  "${M3DC1_ADAS_DATA_DIR}/scd96/scd96_c.dat"
  "${M3DC1_ADAS_DATA_DIR}/scd96/scd96_ne.dat"
  "${M3DC1_ADAS_DATA_DIR}/plt96/plt96_he.dat"
  "${M3DC1_ADAS_DATA_DIR}/plt96/plt96_be.dat"
  "${M3DC1_ADAS_DATA_DIR}/plt96/plt96_c.dat"
  "${M3DC1_ADAS_DATA_DIR}/plt96/plt96_ne.dat"
)
set(M3DC1_OPENADAS_STAMP "${CMAKE_BINARY_DIR}/adas/.openadas.stamp")
add_custom_command(
  OUTPUT
    "${M3DC1_OPENADAS_STAMP}"
    ${M3DC1_OPENADAS_SOURCES}
    ${M3DC1_OPENADAS_DATA_FILES}
  COMMAND "${CMAKE_COMMAND}"
    "-DOUT_DIR=${CMAKE_BINARY_DIR}/adas"
    "-DADF11_DIR=${M3DC1_ADAS_DATA_DIR}"
    -P "${PROJECT_SOURCE_DIR}/cmake/DownloadOpenADAS.cmake"
  COMMAND "${CMAKE_COMMAND}" -E touch "${M3DC1_OPENADAS_STAMP}"
  COMMENT "Downloading OpenADAS ADF11 data and reader sources"
  VERBATIM)
set(_download_open_adas_all)
if(M3DC1_DOWNLOAD_ADAS_DATA)
  set(_download_open_adas_all ALL)
endif()
add_custom_target(download_open_adas ${_download_open_adas_all}
  DEPENDS "${M3DC1_OPENADAS_STAMP}")

if(M3DC1_ENABLE_ADAS)
  list(APPEND M3DC1_LIBRARY_SOURCES adas_m3dc1.f90)
  add_library(openadas_reader STATIC ${M3DC1_OPENADAS_SOURCES})
  add_dependencies(openadas_reader download_open_adas)
endif()

find_package(MPI REQUIRED)
find_package(HDF5 COMPONENTS Fortran HL REQUIRED)
find_package(GSL REQUIRED)
find_package(LAPACK REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h REQUIRED)
find_library(FFTW3_LIBRARY NAMES fftw3 REQUIRED)

if(M3DC1_ENABLE_ST)
  find_library(NETCDF_C_LIBRARY NAMES netcdf
    HINTS ${NetCDF_ROOT}
    PATH_SUFFIXES lib lib64
    REQUIRED)
endif()

if(PETSC_VERSION MATCHES "^([0-9]+)\\.([0-9]+)")
  math(EXPR M3DC1_PETSC_VERSION_CODE "${CMAKE_MATCH_1} * 10 + ${CMAKE_MATCH_2}")
else()
  set(M3DC1_PETSC_VERSION_CODE 0)
endif()
string(TIMESTAMP M3DC1_DATE_BUILT "%a %b %d %H:%M:%S %Y")
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/release_version" M3DC1_RELEASE_VERSION)
string(STRIP "${M3DC1_RELEASE_VERSION}" M3DC1_RELEASE_VERSION)

add_library(m3dc1_library ${M3DC1_LIBRARY_SOURCES})

target_include_directories(m3dc1_library PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(m3dc1_library PUBLIC
  m3dc1::scorec
  MPI::MPI_Fortran
  PETSc::PETSc
  hdf5::hdf5_fortran
  hdf5::hdf5_hl_fortran
  GSL::gsl
  LAPACK::LAPACK
  ${FFTW3_LIBRARY}
)
if(M3DC1_ENABLE_ST)
  target_link_libraries(m3dc1_library PUBLIC
    ${NETCDF_C_LIBRARY}
  )
endif()
if(M3DC1_ENABLE_ADAS)
  target_link_libraries(m3dc1_library PRIVATE openadas_reader)
endif()
target_compile_definitions(m3dc1_library PRIVATE
  USESCOREC
  vector_type=scorec_vector
  matrix_type=scorec_matrix
  PETSC_VERSION=${M3DC1_PETSC_VERSION_CODE}
  DATE_BUILT="${M3DC1_DATE_BUILT}"
  RELEASE_VERSION="${M3DC1_RELEASE_VERSION}"
)
if(M3DC1_ENABLE_3D)
  target_compile_definitions(m3dc1_library PRIVATE -DUSE3D -DMAX_PTS=125)
else()
  target_compile_definitions(m3dc1_library PRIVATE -DMAX_PTS=25)
endif()
if(M3DC1_ENABLE_ST)
  target_compile_definitions(m3dc1_library PRIVATE -DUSEST)
endif()
if(M3DC1_ENABLE_COMPLEX)
  target_compile_definitions(m3dc1_library PRIVATE -Dvectype=complex -DUSECOMPLEX)
else()
  target_compile_definitions(m3dc1_library PRIVATE -Dvectype=real)
endif()
if(M3DC1_ENABLE_PARTICLE)
  target_compile_definitions(m3dc1_library PRIVATE -DUSEPARTICLES)
endif()
if(M3DC1_ENABLE_ADAS)
  target_compile_definitions(m3dc1_library PRIVATE -DUSEADAS)
endif()
set_target_properties(m3dc1_library PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_compile_options(m3dc1_library PUBLIC
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-cpp>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fdefault-real-8>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fdefault-double-8>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fallow-argument-mismatch>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
  )

# --- Build main executable name per Makefile convention ---
set(M3DC1_MAIN_NAME m3dc1)
if(M3DC1_ENABLE_3D)
  string(APPEND M3DC1_MAIN_NAME "_3d")
else()
  string(APPEND M3DC1_MAIN_NAME "_2d")
endif()
if(M3DC1_ENABLE_TRILINOS)
  string(APPEND M3DC1_MAIN_NAME "_trilinos")
endif()
if(M3DC1_ENABLE_COMPLEX)
  string(APPEND M3DC1_MAIN_NAME "_complex")
endif()
if(M3DC1_ENABLE_OPENMP)
  string(APPEND M3DC1_MAIN_NAME "_omp")
endif()
if(M3DC1_ENABLE_ST)
  string(APPEND M3DC1_MAIN_NAME "_st")
endif()
if(M3DC1_ENABLE_PARTICLE)
  string(APPEND M3DC1_MAIN_NAME "_pic")
endif()

add_executable(${M3DC1_MAIN_NAME} $<TARGET_OBJECTS:m3dc1_library>)
target_link_libraries(${M3DC1_MAIN_NAME} PRIVATE m3dc1_library)

add_custom_target(m3dc1 ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${M3DC1_MAIN_NAME} m3dc1
  DEPENDS ${M3DC1_MAIN_NAME}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Creating symlink: m3dc1 -> ${M3DC1_MAIN_NAME}")

add_executable(a2cc a2cc.f90)
target_link_libraries(a2cc PUBLIC m3dc1_library)

add_executable(mackenbach_profiles mackenbach_profiles.f90)

# add_executable(kprad kprad.f90 kprad_test.f90)

add_executable(test_region test_region.f90)
target_link_libraries(test_region PUBLIC m3dc1_library)


if(M3DC1_ENABLE_PSPLINE)
  find_package(PSPLINE REQUIRED)
  add_executable(readgato readgato.f90)
  target_link_libraries(readgato PUBLIC m3dc1_library PSPLINE::PSPLINE)

  add_executable(read_jsolver read_jsolver.f90 polar.f90)
  target_compile_definitions(read_jsolver PRIVATE READ_JSOLVER)
  target_link_libraries(read_jsolver PUBLIC PSPLINE::PSPLINE)
endif()
