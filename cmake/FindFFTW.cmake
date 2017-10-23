# - Find FFTW
# Find FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h PATHS /usr /usr/local)

find_library (FFTW_LIB NAMES fftw3 PATHS /usr/lib /usr/local/lib )
find_library (FFTWF_LIB NAMES fftw3f PATHS /usr/lib /usr/local/lib )
find_library (FFTWT_LIB NAMES fftw3_threads PATHS /usr/lib /usr/local/lib )
find_library (FFTWTF_LIB NAMES fftw3f_threads PATHS /usr/lib /usr/local/lib )

if (FFTWT_LIB)
  set (FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWT_LIB})
endif()
set (FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_LIB})

if (FFTWTF_LIB)
  set (FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWTF_LIB})
endif()
set (FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTWF_LIB})


include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
