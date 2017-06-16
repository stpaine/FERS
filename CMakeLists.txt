CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(fers)

SET(CMAKE_VERBOSE_MAKEFILE ON)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_BINARY_DIR}/config.h )

#Check for build configuration 
IF(NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "No build configuration specified, defaulting to Release")
    SET(CMAKE_BUILD_TYPE "Release")
ENDIF()


#Set compiler flags
IF (CMAKE_COMPILER_IS_GNUCXX)
   MESSAGE(STATUS "Setting general compiler flags for detected compiler: gnu-g++")	
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -pthread -ffast-math")

    IF(CMAKE_BUILD_TYPE MATCHES "Release")
    
        MESSAGE(STATUS "Setting g++ flags for Release configuration")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")   ## Optimize
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")  ## Strip binary
    
        ELSEIF(CMAKE_BUILD_TYPE MATCHES "Debug")
            MESSAGE(STATUS "Setting g++ flags for Debug configuration")
    
    ENDIF()

ENDIF ()

#FERS is multithreaded, boost needs to be informed
ADD_DEFINITIONS(-D_REENTRANT)

ADD_SUBDIRECTORY( fftwcpp )
ADD_SUBDIRECTORY( src )

ADD_SUBDIRECTORY( utilities )
ADD_SUBDIRECTORY( test )

