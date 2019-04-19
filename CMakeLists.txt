cmake_minimum_required (VERSION 2.6)

# Project title
project (toptor)

# Manual assignation
#SET(TOPGEN_PBGL "yes")

# Set the compiler
if(TOPGEN_PBGL)
    message(STATUS "Compiling with PBGL")
    SET(CMAKE_CXX_COMPILER "mpic++")
else()
    message(STATUS "Compiling with BGL. Use -DTOPGEN_PBGL to enable the PBGL support")
    SET(CMAKE_CXX_COMPILER "g++")
endif()

SET(CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -Wall -Wno-deprecated -Wwrite-strings -no-pie")

# Installation by default
IF (NOT CMAKE_INSTALL_PREFIX)
    SET(CMAKE_INSTALL_PREFIX /usr/local/topgen)
ENDIF()
message(STATUS "Install to ${CMAKE_INSTALL_PREFIX}")

# Include topgen headers
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Include MPI libraries
IF(TOPGEN_PBGL)

    # Include boost libraries
    FIND_PACKAGE( Boost 1.65 COMPONENTS program_options serialization system graph graph_parallel mpi thread regex REQUIRED )
    IF(Boost_FOUND)
        message(STATUS "Found Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
        message(STATUS "Found Boost_LIBRARIES: ${Boost_LIBRARIES}")
        INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
    else()
        message(WARNING "Boost NOT found")
    ENDIF()

    #FIND_PACKAGE(MPI REQUIRED)
    #IF (MPI_FOUND)
    #    message(STATUS "Found MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
    #    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
    #ENDIF()
ELSE()
    # Include boost libraries
    FIND_PACKAGE( Boost 1.65 COMPONENTS program_options serialization system graph REQUIRED )
    IF(Boost_FOUND)
        message(STATUS "Found Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
        message(STATUS "Found Boost_LIBRARIES: ${Boost_LIBRARIES}")
        INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
    else()
        message(WARNING "Boost NOT found")
    ENDIF()

ENDIF()

# Enable GTest tests
#enable_testing()
#find_package( GTest REQUIRED )

#if(GTEST_FOUND)
#    message(STATUS "Found GTEST_INCLUDE_DIRS: ${GTEST_INCLUDE_DIRS}")
#    message(STATUS "Found GTEST_LIBRARIES: ${GTEST_LIBRARIES}")
#
#    include_directories( ${GTEST_INCLUDE_DIRS})
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
#endif()

if (TOPGEN)
    find_path(TOPGEN_INCLUDE_DIR NAMES topgen.hpp HINTS ${TOPGEN}/include HINTS ${TOPGEN}/include/topgen)
    if (NOT TOPGEN_INCLUDE_DIR)
        message(FATAL_ERROR "I cannot find the file topgen.hpp. Set -DTOPGEN in CMakeLists.txt")
    endif()
    message(STATUS "topgen.hpp found at ${TOPGEN_INCLUDE_DIR}")
    include_directories(${TOPGEN_INCLUDE_DIR})

    find_library(TOPGEN_LIBRARY NAMES topgen HINTS ${TOPGEN}/lib)
    if (NOT TOPGEN_LIBRARY)
        message(FATAL_ERROR "I cannot find the file libtopgen")
    endif()
    message(STATUS "libtopgen fount at ${TOPGEN_LIBRARY}")
    link_libraries(${TOPGEN_LIBRARY})
endif()

# Include source files
file(GLOB TOPTOR_SOURCES "src/toptor.cpp")

# Build the tool
add_executable(toptor ${TOPTOR_SOURCES})

IF(TOPGEN_PBGL)
    #target_include_directories(toptor PRIVATE ${MPI_CXX_INCLUDE_PATH} )
    #target_compile_options(toptor PRIVATE ${MPI_CXX_COMPILE_FLAGS} )
    target_link_libraries(toptor ${TOPGEN_LIBRARY} m  ${Boost_LIBRARIES})
ELSE()
    target_link_libraries(toptor ${TOPGEN_LIBRARY} m  ${Boost_LIBRARIES})
ENDIF()

install(TARGETS toptor DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

