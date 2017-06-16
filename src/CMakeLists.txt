INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${fers_SOURCE_DIR}/fftwcpp )

########### next target ###############

SET(fers_SRCS
   rsdebug.cpp
   rsmain.cpp
   rsplatform.cpp
   rsradar.cpp
   rstarget.cpp
   rsworld.cpp
   xmlimport.cpp
   rsgeometry.cpp
   rspath.cpp
   rssim.cpp
   rsthreadedsim.cpp
   rspulserender.cpp
   rsradarwaveform.cpp
   rsparameters.cpp
   rsobject.cpp
   rsantenna.cpp
   rsportable.cpp
   rsresponse.cpp
   rsnoise.cpp
   rssignal.cpp
   rstiming.cpp
   rshdf5.cpp
   rspython.cpp
   rsinterp.cpp
   rsdsp.cpp
   rsmultipath.cpp
   rspattern.cpp
   rspolarize.cpp
   rsdebug.h
   rspath.h
   rsportable.h
   rsradar.h
   rssim.h
   rsworld.h
   rsgeometry.h
   rsplatform.h
   rstarget.h
   xmlimport.h
   rsradarwaveform.h
   rsparameters.h
   rssignal.h
   rsthreadedsim.h
   rsobject.h
   rsantenna.h
   rsresponse.h
   rsnoise.h
   rstiming.h
   rshdf5.h
   rspython.h
   rsinterp.h
   rsdsp.h
   rsmultipath.h
   rspattern.h
   rspolarize.h
   )

ADD_EXECUTABLE(fers ${fers_SRCS})

#Link to the HDF5 libraries
FIND_LIBRARY(FERS_LIB_HDF5 NAMES hdf5 )

IF (FERS_LIB_HDF5)
   TARGET_LINK_LIBRARIES(fers ${FERS_LIB_HDF5})
   ADD_DEFINITIONS(-DHAVE_LIBHDF5)
ENDIF (FERS_LIB_HDF5)

FIND_LIBRARY(FERS_LIB_HDF5_HL NAMES hdf5_hl )

IF (FERS_LIB_HDF5_HL)
   TARGET_LINK_LIBRARIES(fers ${FERS_LIB_HDF5_HL})
ENDIF (FERS_LIB_HDF5_HL)

#Link to tinyxml libraries
FIND_LIBRARY(FERS_LIB_TINYXML NAMES tinyxml)

IF (FERS_LIB_TINYXML)
   TARGET_LINK_LIBRARIES(fers ${FERS_LIB_TINYXML})
ENDIF (FERS_LIB_TINYXML) 

#Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(fers fftwcpp)

#See if the FFTW3 libraries are available
FIND_LIBRARY(FERS_LIB_FFTW3 NAMES fftw3)

IF (FERS_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(fers ${FERS_LIB_FFTW3})
ENDIF (FERS_LIB_FFTW3) 

#Check if we are building on a multithreaded FFTW
FIND_LIBRARY(FERS_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (FERS_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(fers ${FERS_LIB_FFTW3_THREADS})
ENDIF (FERS_LIB_FFTW3_THREADS)

#Include code to check for python libs
#We must use Python2
SET(Python_ADDITIONAL_VERSIONS 2.7 2.6 2.5)
INCLUDE(FindPythonLibs)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(fers ${PYTHON_LIBRARIES})
   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF (PYTHONLIBS_FOUND)

#Link to Boost Thread and System
IF(WIN32)
   SET(Boost_USE_STATIC_LIBS ON) ##Link Boost statically
ELSE ()
   SET(Boost_USE_STATIC_LIBS OFF) ##Link Boost with shared libs
ENDIF ()
FIND_PACKAGE(Boost COMPONENTS system thread REQUIRED)

MESSAGE(STATUS "Linking libraries for Boost: " ${Boost_LIBRARIES})

IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(fers ${Boost_LIBRARIES} )
ENDIF ()

########### install files ###############

INSTALL( TARGETS fers DESTINATION bin)
