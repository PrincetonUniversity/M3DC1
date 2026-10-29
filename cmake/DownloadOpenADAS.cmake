# Download OpenADAS ADF11 data and reader sources
# Called from unstructured/CMakeLists.txt with:
#   -DOUT_DIR=<path>     : output directory (CMAKE_BINARY_DIR/adas)
#   -DADF11_DIR=<path>   : where ADF11 data files are placed

set(OA_DATA_URL "https://open.adas.ac.uk/download/adf11")
set(OA_SOURCE_DIR "${OUT_DIR}/source")

# ---- Download and extract reader source code ----
message(STATUS "Downloading OpenADAS reader sources...")
file(DOWNLOAD
  "https://open.adas.ac.uk/code/xxdata_11.tar.gz"
  "${OUT_DIR}/xxdata_11.tar.gz"
  STATUS _dl_status
  SHOW_PROGRESS
)
list(GET _dl_status 0 _dl_code)
if(NOT _dl_code EQUAL 0)
  message(FATAL_ERROR "Failed to download xxdata_11.tar.gz")
endif()

# Extract into a temporary directory, then move .for → .f to source/
set(_extract_dir "${OUT_DIR}/_extract_tmp")
file(REMOVE_RECURSE "${_extract_dir}")
file(MAKE_DIRECTORY "${_extract_dir}")
file(MAKE_DIRECTORY "${OA_SOURCE_DIR}")

execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xzf "${OUT_DIR}/xxdata_11.tar.gz"
  WORKING_DIRECTORY "${_extract_dir}"
  RESULT_VARIABLE _tar_res
)
if(NOT _tar_res EQUAL 0)
  message(FATAL_ERROR "Failed to extract xxdata_11.tar.gz")
endif()

# Find all .for files (recursively, handles possible subdirectories)
# and copy/rename them to OA_SOURCE_DIR as .f
file(GLOB_RECURSE _for_files "${_extract_dir}/*.for")
foreach(_f IN LISTS _for_files)
  get_filename_component(_name "${_f}" NAME_WE)
  file(COPY "${_f}" DESTINATION "${OA_SOURCE_DIR}")
  file(RENAME "${OA_SOURCE_DIR}/${_name}.for" "${OA_SOURCE_DIR}/${_name}.f")
endforeach()

# Cleanup
file(REMOVE_RECURSE "${_extract_dir}")
file(REMOVE "${OUT_DIR}/xxdata_11.tar.gz")
message(STATUS "OpenADAS reader sources extracted to ${OA_SOURCE_DIR}")

# ---- Download ADF11 data files ----
set(_adf11_files
  "scd85/scd85_ar.dat"
  "scd89/scd89_b.dat"
  "plt89/plt89_b.dat"
  "plt89/plt89_ar.dat"
  "scd96/scd96_he.dat"
  "scd96/scd96_be.dat"
  "scd96/scd96_c.dat"
  "scd96/scd96_ne.dat"
  "plt96/plt96_he.dat"
  "plt96/plt96_be.dat"
  "plt96/plt96_c.dat"
  "plt96/plt96_ne.dat"
)

foreach(_f IN LISTS _adf11_files)
  get_filename_component(_subdir "${_f}" DIRECTORY)
  file(MAKE_DIRECTORY "${ADF11_DIR}/${_subdir}")
  message(STATUS "Downloading ${_f}...")
  file(DOWNLOAD
    "${OA_DATA_URL}/${_f}"
    "${ADF11_DIR}/${_f}"
    STATUS _dl_status
  )
  list(GET _dl_status 0 _dl_code)
  if(NOT _dl_code EQUAL 0)
    message(FATAL_ERROR "Failed to download ${_f}")
  endif()
endforeach()

message(STATUS "OpenADAS download complete")