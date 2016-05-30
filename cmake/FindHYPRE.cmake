# - Try to find HYPRE
#

SET(TRIAL2_LIBRARY_PATHS
        /usr/local/lib
        /sw/lib
        ${CMAKE_SOURCE_DIR}/lib
        /home/czaki/libs/hypre/2.10/lib
        /Users/grzegorzbokota/libs/hypre/2.10.1/lib
        )

SET(TRIAL2_INCLUDE_PATHS
        /usr/local/include
        /usr/include/sprng
        /sw/include
        ${CMAKE_SOURCE_DIR}/include
        $ENV{home}/libs
        /home/czaki/libs/hypre/2.10/include
        /Users/grzegorzbokota/libs/hypre/2.10.1/include
        )

FIND_LIBRARY(HYPRE_LIBRARY HYPRE ${TRIAL2_LIBRARY_PATHS})
FIND_PATH(HYPRE_INCLUDE_DIR HYPRE.h ${TRIAL2_INCLUDE_PATHS})