cmake_minimum_required(VERSION 3.2)
project(timothy)
#set(CMAKE_PREFIX_PATH /home/czaki/libs/Zoltan/3.82)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})


find_package(Zoltan)
find_package (HYPRE REQUIRED)
find_package(SPRNG REQUIRED)
#include_directories(${ZOLTAN_INCLUDE_DIR})
#include_directories(${HYPRE_INCLUDE_DIR})
#include_directories(${SPRNG_INCLUDE_DIR})
include_directories("/usr/include")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -lgmp")


set(SOURCE_FILES
    src/cells.c
    src/cells.h
    src/chemf.c
    src/chemf.h
    src/comm.c
    src/comm.h
    src/compute.c
    src/compute.h
    src/domdec.c
    src/domdec.h
    src/fields.c
    src/fields.h
    src/global.h
    src/grid.c
    src/grid.h
    src/init.c
    src/init.h
    src/inline.h
    src/interp.c
    src/interp.h
    #src/io.c
    src/io.h
    src/main.c
    src/octree.c
    src/octree.h
    src/potential.c
    src/potential.h
    src/random.c
    src/random.h
    src/stats.c
    src/stats.h
    #src/tempf.c
    #src/tempf.h
    src/utils.c
    src/utils.h
    src/dicts.h
    src/ini_parser/ini_parser.c
    src/ini_parser/ini_parser.h
    #src/functions.cpp src/functions.h
    src/local_endian.h
    src/fake_io.c)

set(VALIDATOR_SOURCE src/ini_parser/ini_parser.c src/ini_parser/ini_parser.h src/validator.cpp src/validator/validator.cpp src/validator/validator.h src/ini_manipulator/ini_manipulator.cpp src/ini_manipulator/ini_manipulator.h src/functions.cpp src/functions.h src/local_endian.h)

add_executable(timothy ${SOURCE_FILES})
target_include_directories(timothy PUBLIC ${ZOLTAN_INCLUDE_DIR})
target_include_directories(timothy PUBLIC ${HYPRE_INCLUDE_DIR})
target_include_directories(timothy PUBLIC ${SPRNG_INCLUDE_DIR})
target_link_libraries (timothy ${HYPRE_LIBRARY})
target_link_libraries (timothy ${ZOLTAN_LIBRARY})
target_link_libraries (timothy ${SPRNG_LIBRARY})
target_link_libraries (timothy m)
target_link_libraries(timothy ${MPI_C_LIBRARIES})

add_executable(validator ${VALIDATOR_SOURCE})
target_link_libraries (validator ${ZOLTAN_LIBRARY})
target_include_directories(validator PUBLIC ${ZOLTAN_INCLUDE_DIR})
target_link_libraries(validator ${MPI_C_LIBRARIES})
target_link_libraries(validator ${MPI_CXX_LIBRARIES})


if(MPI_COMPILE_FLAGS)
    set_target_properties(timothy PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(validator PROPERTIES
        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(timothy PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
    set_target_properties(validator PROPERTIES
        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/Doxyfile @ONLY)
    add_custom_target(doc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif(DOXYGEN_FOUND)
set( CMAKE_VERBOSE_MAKEFILE on )