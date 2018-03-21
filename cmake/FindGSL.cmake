find_library(GSL_LIB gsl)
find_library(GSL_CBLAS_LIB gslcblas)

set(GSL_LIBRARIES ${GSL_LIB} ${GSL_CBLAS_LIB})
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_LIBRARIES)
mark_as_advanced(GSL_LIBRARIES)