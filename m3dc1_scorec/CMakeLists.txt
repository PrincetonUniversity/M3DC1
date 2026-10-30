# Doxygen documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in"
    "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
endif()

# Generate config file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_scorecConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/m3dc1_scorecConfig.cmake")

#Source and header files
set(SOURCES
  api/m3dc1_scorec.cc
  src/BSpline.cc
  src/Expression.cc
  src/m3dc1_ls.cc
  src/m3dc1_field.cc
  src/m3dc1_mesh.cc
  src/m3dc1_adapt.cc
  src/m3dc1_entity.cc
  src/m3dc1_sizeField.cc
  src/PolyNomial.cc
  src/ReducedQuinticImplicit.cc
  src/CMODGeoExpression.cc
  src/m3dc1_matrix.cc
  src/m3dc1_model.cc
  src/m3dc1_slnTransfer.cc
  #src/ReducedQuinticExplicit.cc
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
  include/m3dc1_slnTransfer.h
  include/PolyNomial.h
  include/ReducedQuinticImplicit.h
  include/CMODGeoExpression.h  
  include/m3dc1_matrix.h
  include/m3dc1_sizeField.h
  #include/ReducedQuinticExplicit.h
  include/slntransferUtil.h
  include/Expression.h
  include/ReducedQuintic.h)


find_package(SCOREC CONFIG REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED COMPONENTS C)

add_library(m3dc1_scorec ${SOURCES})
add_library(m3dc1::scorec ALIAS m3dc1_scorec)
target_link_libraries(m3dc1_scorec PUBLIC SCOREC::core PETSc::PETSc LAPACK::LAPACK MPI::MPI_C)
target_include_directories(m3dc1_scorec PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/api"
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_compile_definitions(m3dc1_scorec PRIVATE M3DC1_PETSC)

# --- Utility executables (SCOREC tools) ---
# Utilities that only depend on SCOREC core
set(M3DC1_UTILITY_CORE_SOURCES
  utility/check_smb.cc
  utility/create_smb.cc
  utility/split_smb.cc
  utility/vtk_order.cc
)

foreach(_util_src IN LISTS M3DC1_UTILITY_CORE_SOURCES)
  get_filename_component(_util_name ${_util_src} NAME_WE)
  add_executable(${_util_name} ${_util_src})
  target_link_libraries(${_util_name} PRIVATE SCOREC::core MPI::MPI_C)
endforeach()

# show_meshcount depends on the m3dc1_scorec library (m3dc1_mesh, m3dc1_model, etc.)
add_executable(show_meshcount utility/show_meshcount.cc)
target_link_libraries(show_meshcount PRIVATE m3dc1_scorec)

# Install rules
install(TARGETS m3dc1_scorec
  EXPORT m3dc1_scorecTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include)

install(FILES ${HEADERS} DESTINATION include)

# Install utility executables
set(M3DC1_UTILITY_NAMES check_smb create_smb split_smb vtk_order show_meshcount)
install(TARGETS ${M3DC1_UTILITY_NAMES} DESTINATION bin)

# CPack configuration for binary distribution
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION "1.0.1")
include(CPack)
