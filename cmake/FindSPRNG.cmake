#
# this module look for SPRNG (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# SPRNG_INCLUDE_PATH = where sprng.h can be found
# SPRNG_LIBRARY      = the library to link against (sprng etc)
#

SET(TRIAL_LIBRARY_PATHS
  /usr/local/lib
  /sw/lib
  ${CMAKE_SOURCE_DIR}/lib
  $ENV{SPRNG_HOME}/lib
  $ENV{HOME}/libs/sprng/lib
)

SET(TRIAL_INCLUDE_PATHS
  /usr/local/include
  /usr/include/sprng
  /sw/include
  ${CMAKE_SOURCE_DIR}/include
  $ENV{SPRNG_HOME}/include
  $ENV{HOME}/libs/sprng/include
)

FIND_LIBRARY(SPRNG_LIBRARY sprng ${TRIAL_LIBRARY_PATHS})
FIND_PATH(SPRNG_INCLUDE_DIR sprng.h ${TRIAL_INCLUDE_PATHS})

IF(SPRNG_INCLUDE_DIR AND SPRNG_LIBRARY)
  SET(SPRNG_FOUND 1 CACHE BOOL "Found sprng library")
  MESSAGE(STATUS "SPRNG_INCLUDE_DIR=${SPRNG_INCLUDE_DIR}")
  MESSAGE(STATUS "SPRNG_LIBRARY=${SPRNG_LIBRARY}")
ELSE(SPRNG_INCLUDE_DIR AND SPRNG_LIBRARY)
  SET(SPRNG_FOUND 0 CACHE BOOL "Not found sprng library")
ENDIF(SPRNG_INCLUDE_DIR AND SPRNG_LIBRARY)

MARK_AS_ADVANCED(
  SPRNG_INCLUDE_DIR
  SPRNG_LIBRARY
  SPRNG_FOUND
)
