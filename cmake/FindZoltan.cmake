# - Try to find zoltan
# Once done this will define
#  ZOLTAN_FOUND - System has ZOLTAN
#  ZOLTAN_INCLUDE_DIRS - The ZOLTAN include directories
#  ZOLTAN_LIBRARIES - The libraries needed to use ZOLTAN
#  ZOLTAN_DEFINITIONS - Compiler switches required for using ZOLTAN


SET(TRIAL2_LIBRARY_PATHS
        /usr/local/lib
        /sw/lib
        ${CMAKE_SOURCE_DIR}/lib
        $ENV{HOME}/libs/zoltan/3.82/lib
        )

SET(TRIAL2_INCLUDE_PATHS
        /usr/local/include
        /usr/include/sprng
        /sw/include
        ${CMAKE_SOURCE_DIR}/include
        $ENV{home}/libs
        $ENV{HOME}/libs/zoltan/3.82/include
        )

FIND_LIBRARY(ZOLTAN_LIBRARY zoltan ${TRIAL2_LIBRARY_PATHS})
FIND_PATH(ZOLTAN_INCLUDE_DIR zoltan.h ${TRIAL2_INCLUDE_PATHS})


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ZOLTAN_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    ZOLTAN  
    DEFAULT_MSG
    ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR
)

mark_as_advanced(ZOLTAN_INCLUDE_DIR ZOLTAN_LIBRARY )
