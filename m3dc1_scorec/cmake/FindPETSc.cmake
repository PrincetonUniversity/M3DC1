find_package(PkgConfig REQUIRED QUIET)
if(DEFINED PETSC_DIR)
  if(DEFINED PETSC_ARCH)
    set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  else()
    set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  endif()
endif()

# we give an internal name _petsc
# so we can fill up the PETSC_VARIABLE based
# on static or not
pkg_check_modules(_petsc PETSc QUIET)

if(_petsc_FOUND AND _petsc_VERSION)
  set(PETSC_VERSION ${_petsc_VERSION})
endif()

# note, there are a number of additional properties that
# can be extracted / set on a target. We just do the basic
# set for now. See: https://cmake.org/cmake/help/latest/module/FindPkgConfig.html
if(PETSC_LINK_STATIC)
  set(PETSC_LIBRARIES ${_petsc_STATIC_LIBRARIES})
  set(PETSC_INCLUDE_DIRS ${_petsc_STATIC_INCLUDE_DIRS})
  set(PETSC_LIBRARY_DIRS ${_petsc_STATIC_LIBRARY_DIRS})
  set(PETSC_LDFLAGS ${_petsc_STATIC_LDFLAGS})
  set(PETSC_LDFLAGS_OTHER ${_petsc_STATIC_LDFLAGS_OTHER})
elseif(_petsc_FOUND)
  set(PETSC_LIBRARIES ${_petsc_LIBRARIES})
  set(PETSC_INCLUDE_DIRS ${_petsc_INCLUDE_DIRS})
  set(PETSC_LIBRARY_DIRS ${_petsc_LIBRARY_DIRS})
  set(PETSC_LDFLAGS ${_petsc_LDFLAGS})
  set(PETSC_LDFLAGS_OTHER ${_petsc_LDFLAGS_OTHER})
endif()

if(NOT TARGET PETSc::PETSc)
  add_library(PETSc::PETSc INTERFACE IMPORTED GLOBAL)
  set_target_properties(PETSc::PETSc PROPERTIES INTERFACE_LINK_LIBRARIES "${PETSC_LIBRARIES}"
                                     INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIRS}"
                                     INTERFACE_LINK_DIRECTORIES "${PETSC_LIBRARY_DIRS}"
                                     INTERFACE_LINK_OPTIONS "${PETSC_LDFLAGS_OTHER}")
endif()



include(FindPackageHandleStandardArgs)
# TODO consider adding version check logic
find_package_handle_standard_args(PETSc
  REQUIRED_VARS PETSC_LIBRARIES PETSC_INCLUDE_DIRS
  VERSION_VAR PETSC_VERSION)

#message(STATUS "STATIC_LIBRARIES=${_petsc_STATIC_LIBRARIES}")
#message(STATUS "STATIC_LDFLAGS=${_petsc_STATIC_LDFLAGS}")
#message(STATUS "STATIC_LDFLAGS_OTHER=${_petsc_STATIC_LDFLAGS_OTHER}")