# $Id$

# find the Arpack library 
#
# Joerg F. Unger, BAM 2013
#
# variables used by this module (can be also defined as environment variables):
#   ARPACK_ROOT - preferred installation prefix for searching for ARPACK
#   ARPACK_FIND_STATIC_LIBRARY - searches for static libraries (UNIX only)
#   ARPACK_DEBUG - print debug messages
#
# variables defined by this module
#   ARPACK_FOUND - defines whether metis was found or not
#   ARPACK_INCLUDE_DIR - ARPACK include directory
#   ARPACK_LIBRARIES   - ARPACK libraries


# initialize variables
MESSAGE(STATUS "Checking for ARPACK Library ...")
# check if ARPACK_ROOT is set
IF(NOT ARPACK_ROOT AND NOT $ENV{ARPACK_ROOT} STREQUAL "")
  SET(ARPACK_ROOT $ENV{ARPACK_ROOT})
ENDIF(NOT ARPACK_ROOT AND NOT $ENV{ARPACK_ROOT} STREQUAL "")

# convert path to unix style path and set search path
IF(ARPACK_ROOT)
  FILE(TO_CMAKE_PATH ${ARPACK_ROOT} ARPACK_ROOT)
  SET(_ARPACK_LIBRARIES_SEARCH_DIRS ${ARPACK_ROOT}/lib ${ARPACK_ROOT} ${_ARPACK_LIBRARIES_SEARCH_DIRS})
ENDIF(ARPACK_ROOT)

# search for ARPACK library
IF(UNIX AND ARPACK_FIND_STATIC_LIBRARY)
  SET(ARPACK_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
ENDIF(UNIX AND ARPACK_FIND_STATIC_LIBRARY)
FIND_LIBRARY(ARPACK_LIBRARIES NAMES arpack HINTS ${_ARPACK_LIBRARIES_SEARCH_DIRS})
IF(UNIX AND ARPACK_FIND_STATIC_LIBRARY)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ${ARPACK_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
ENDIF(UNIX AND ARPACK_FIND_STATIC_LIBRARY)

# handle the QUIETLY and REQUIRED arguments
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ARPACK DEFAULT_MSG ARPACK_LIBRARIES)

if(ARPACK_DEBUG)
  MESSAGE(STATUS "ARPACK_FOUND=${ARPACK_FOUND}")
  MESSAGE(STATUS "ARPACK_LIBRARIES=${ARPACK_LIBRARIES}")
endif(ARPACK_DEBUG)
