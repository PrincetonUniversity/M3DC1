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
find_package(PETSc REQUIRED)

add_library(m3dc1_scorec ${SOURCES})
add_library(m3dc1::scorec ALIAS m3dc1_scorec)
target_link_libraries(m3dc1_scorec PRIVATE SCOREC::core PETSc::PETSc)
#target_include_directories(m3dc1_scorec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/api ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(m3dc1_scorec PUBLIC api include)
target_compile_definitions(m3dc1_scorec PRIVATE M3DC1_PETSC)

