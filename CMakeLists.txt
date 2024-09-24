################################################################################
# Project: RooUnfold
################################################################################

if(${NO_DEFAULT_ROOUNFOLD})
  return()
endif()

cmake_minimum_required(VERSION 3.0)

# The project version number.
project(RooUnfold)
set(VERSION_MAJOR   3   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR)

find_package( ROOT COMPONENTS Tree  Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser)

set(CMAKE_BUILD_TYPE Debug)

option (RooUnfoldGenerateCMakeConfig "generate a RooUnfoldConfig.cmake file for dependency tracking with other packages" OFF)
# See https://cmake.org/cmake/help/latest/policy/CMP0077.html
# For compatability with old CMake versions, this variable is unset, and then set again.
unset(RooUnfoldTests)
option (RooUnfoldTests "enable testing" OFF)
option (RooUnfoldEnableProfiling "enable profiling" OFF)

find_package( ROOT COMPONENTS Tree Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser REQUIRED)

set(ROOUNFOLD_SOURCE_DIR src)
set(ROOUNFOLD_HEADER_DIR src)
include(cmake/version.cmake)
include(cmake/source.cmake)

if(ROOT_VERSION VERSION_LESS 6.20)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS} -g -DNO_WRAPPERPDF")
else()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS} -g")
endif()

execute_process( COMMAND ln -sf ${RooUnfoldHeaders} -t ${CMAKE_CURRENT_BINARY_DIR} )

include(cmake/PlainROOT.cmake)

if(RooUnfoldTests)  
  include(cmake/testing.cmake)
endif()

include(cmake/doxygen.cmake)

