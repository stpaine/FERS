INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${fers_SOURCE_DIR}/src )
INCLUDE_DIRECTORIES( ${fers_SOURCE_DIR}/fftwcpp )

########### noisegen_test ###############

SET(noisegen_test_SRCS
   noisegen_test.cpp
   ${fers_SOURCE_DIR}/src/rsnoise.cpp
   ${fers_SOURCE_DIR}/src/rsnoise.h
   ${fers_SOURCE_DIR}/src/rsdsp.cpp
   ${fers_SOURCE_DIR}/src/rsdsp.h
   ${fers_SOURCE_DIR}/src/rspython.cpp
   ${fers_SOURCE_DIR}/src/rspython.h
   ${fers_SOURCE_DIR}/src/rsparameters.cpp
   ${fers_SOURCE_DIR}/src/rsparameters.h
   ${fers_SOURCE_DIR}/src/rsgeometry.cpp
   ${fers_SOURCE_DIR}/src/rsgeometry.h
   ${fers_SOURCE_DIR}/src/rsdebug.cpp
   ${fers_SOURCE_DIR}/src/rsdebug.h
   )

ADD_EXECUTABLE(noisegen_test ${noisegen_test_SRCS})

#Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(noisegen_test fftwcpp)

#See if the FFTW3 libraries are available
FIND_LIBRARY(FERS_LIB_FFTW3 NAMES fftw3)

IF (FERS_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(noisegen_test ${FERS_LIB_FFTW3})
ENDIF (FERS_LIB_FFTW3) 

#Check if we are building on a multithreaded FFTW
FIND_LIBRARY(FERS_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (FERS_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(noisegen_test ${FERS_LIB_FFTW3_THREADS})
ENDIF (FERS_LIB_FFTW3_THREADS)

#Include code to check for python libs

INCLUDE(FindPythonLibs)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(noisegen_test ${PYTHON_LIBRARIES})
   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF (PYTHONLIBS_FOUND)

#Link to Boost Thread and System

IF (WIN32)
   SET(Boost_USE_STATIC_LIBS ON) ##Link Boost statically
ELSE ()
   SET(Boost_USE_STATIC_LIBS OFF) ##Link Boost with shared libs
ENDIF ()
FIND_PACKAGE(Boost COMPONENTS system thread REQUIRED)

MESSAGE(STATUS "Linking libraries for Boost: " ${Boost_LIBRARIES})

IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(noisegen_test ${Boost_LIBRARIES} )
ENDIF ()


########### timing_test ###############

SET(timing_test_SRCS
   timing_test.cpp
   ${fers_SOURCE_DIR}/src/rsnoise.cpp
   ${fers_SOURCE_DIR}/src/rsnoise.h
   ${fers_SOURCE_DIR}/src/rsdsp.cpp
   ${fers_SOURCE_DIR}/src/rsdsp.h
   ${fers_SOURCE_DIR}/src/rspython.cpp
   ${fers_SOURCE_DIR}/src/rspython.h
   ${fers_SOURCE_DIR}/src/rsparameters.cpp
   ${fers_SOURCE_DIR}/src/rsparameters.h
   ${fers_SOURCE_DIR}/src/rsgeometry.cpp
   ${fers_SOURCE_DIR}/src/rsgeometry.h
   ${fers_SOURCE_DIR}/src/rsdebug.cpp
   ${fers_SOURCE_DIR}/src/rsdebug.h
   ${fers_SOURCE_DIR}/src/rstiming.cpp
   ${fers_SOURCE_DIR}/src/rstiming.h
   )

ADD_EXECUTABLE(timing_test ${timing_test_SRCS})

#Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(timing_test fftwcpp)

#See if the FFTW3 libraries are available
FIND_LIBRARY(FERS_LIB_FFTW3 NAMES fftw3)

IF (FERS_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(timing_test ${FERS_LIB_FFTW3})
ENDIF (FERS_LIB_FFTW3) 

#Check if we are building on a multithreaded FFTW
FIND_LIBRARY(FERS_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (FERS_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(timing_test ${FERS_LIB_FFTW3_THREADS})
ENDIF (FERS_LIB_FFTW3_THREADS)

#Include code to check for python libs

INCLUDE(FindPythonLibs)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(timing_test ${PYTHON_LIBRARIES})
   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF (PYTHONLIBS_FOUND)

#Link to Boost Thread and System
IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(timing_test ${Boost_LIBRARIES} )
ENDIF ()

########### clock_test ###############

SET(clock_test_SRCS
   clock_test.cpp
   ${fers_SOURCE_DIR}/src/rsnoise.cpp
   ${fers_SOURCE_DIR}/src/rsnoise.h
   ${fers_SOURCE_DIR}/src/rsdsp.cpp
   ${fers_SOURCE_DIR}/src/rsdsp.h
   ${fers_SOURCE_DIR}/src/rspython.cpp
   ${fers_SOURCE_DIR}/src/rspython.h
   ${fers_SOURCE_DIR}/src/rsparameters.cpp
   ${fers_SOURCE_DIR}/src/rsparameters.h
   ${fers_SOURCE_DIR}/src/rsgeometry.cpp
   ${fers_SOURCE_DIR}/src/rsgeometry.h
   ${fers_SOURCE_DIR}/src/rsdebug.cpp
   ${fers_SOURCE_DIR}/src/rsdebug.h
   ${fers_SOURCE_DIR}/src/rstiming.cpp
   ${fers_SOURCE_DIR}/src/rstiming.h
   )

ADD_EXECUTABLE(clock_test ${clock_test_SRCS})

#Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(clock_test fftwcpp)

#See if the FFTW3 libraries are available
IF (FERS_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(clock_test ${FERS_LIB_FFTW3})
ENDIF (FERS_LIB_FFTW3) 

#Check if we are building on a multithreaded FFTW
IF (FERS_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(clock_test ${FERS_LIB_FFTW3_THREADS})
ENDIF (FERS_LIB_FFTW3_THREADS)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(clock_test ${PYTHON_LIBRARIES})
ENDIF (PYTHONLIBS_FOUND)

#Link to Boost thread and system
IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(clock_test ${Boost_LIBRARIES} )
ENDIF ()


########### filter_test ###############

SET(filter_test_SRCS
   filter_test.cpp
   ${fers_SOURCE_DIR}/src/rsnoise.cpp
   ${fers_SOURCE_DIR}/src/rsnoise.h
   ${fers_SOURCE_DIR}/src/rsdsp.cpp
   ${fers_SOURCE_DIR}/src/rsdsp.h
   ${fers_SOURCE_DIR}/src/rspython.cpp
   ${fers_SOURCE_DIR}/src/rspython.h
   ${fers_SOURCE_DIR}/src/rsparameters.cpp
   ${fers_SOURCE_DIR}/src/rsparameters.h
   ${fers_SOURCE_DIR}/src/rsgeometry.cpp
   ${fers_SOURCE_DIR}/src/rsgeometry.h
   ${fers_SOURCE_DIR}/src/rsdebug.cpp
   ${fers_SOURCE_DIR}/src/rsdebug.h
   )

ADD_EXECUTABLE(filter_test ${filter_test_SRCS})

#Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(filter_test fftwcpp)

#See if the FFTW3 libraries are available
FIND_LIBRARY(FERS_LIB_FFTW3 NAMES fftw3)

IF (FERS_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(filter_test ${FERS_LIB_FFTW3})
ENDIF (FERS_LIB_FFTW3) 

#Check if we are building on a multithreaded FFTW
FIND_LIBRARY(FERS_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (FERS_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(filter_test ${FERS_LIB_FFTW3_THREADS})
ENDIF (FERS_LIB_FFTW3_THREADS)

#Include code to check for python libs

INCLUDE(FindPythonLibs)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(filter_test ${PYTHON_LIBRARIES})
   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF (PYTHONLIBS_FOUND)

#Link to Boost Thread and System

IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(filter_test ${Boost_LIBRARIES} )
ENDIF ()

########### install files ###############

