# sources
set(SKELETON_BASE_SOURCES
  etimes.F90
  math.F90
  element.F90
  mesh.F90
  nintegrate.F90
  field.F90
  m3dc1_data.F90
  m3dc1_nint.F90
  matdef.F90
  main.F90
)

find_package(LAPACK REQUIRED)

add_executable(m3dc1_skeleton ${SKELETON_BASE_SOURCES})
target_compile_definitions(m3dc1_skeleton PRIVATE USE3D MAX_PTS=60)
target_link_libraries(m3dc1_skeleton PRIVATE LAPACK::LAPACK)

# --- OpenACC executable (requires PGI/NVIDIA compiler) ---
option(M3DC1_ENABLE_SKELETON_ACC "Build OpenACC skeleton" OFF)
if(M3DC1_ENABLE_SKELETON_ACC)
  set(SKELETON_ACC_SOURCES
    etimes.F90
    math.F90
    element_acc.F90
    mesh_acc.F90
    field_acc.F90
    m3dc1_data_acc.F90
    nintegrate_acc.F90
    dgemm.f
    dgemv.f
    m3dc1_nint_acc.F90
    matdef_acc.F90
    main_acc.F90
  )
  add_executable(m3dc1_skeleton_acc ${SKELETON_ACC_SOURCES})
  target_compile_definitions(m3dc1_skeleton_acc PRIVATE USE3D MAX_PTS=60)
  target_link_libraries(m3dc1_skeleton_acc PRIVATE LAPACK::LAPACK)
endif()

install(TARGETS m3dc1_skeleton
  RUNTIME DESTINATION bin
)
if(TARGET m3dc1_skeleton_acc)
  install(TARGETS m3dc1_skeleton_acc RUNTIME DESTINATION bin)
endif()