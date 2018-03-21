find_library(FFTW_LIB fftw3)
find_library(FFTW_THREADS_LIB fftw3_threads)

set(FFTW_LIBRARIES ${FFTW_LIB} ${FFTW_THREADS_LIB})
message(STATUS "FFTW : ${FFTW_LIBRARIES}")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES)
mark_as_advanced(FFTW_LIBRARIES)