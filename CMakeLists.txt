################################################################################
# Project: RooUnfold
################################################################################

cmake_minimum_required(VERSION 2.8)

# The project version number.
project(RooUnfold)
set(VERSION_MAJOR   3   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR)

find_package( ROOT COMPONENTS Tree Unfold Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser)

include(cmake/source.cmake)

if(ROOT_VERSION VERSION_LESS 6.20)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS} -g -DNO_WRAPPERPDF")
endif()

execute_process( COMMAND ln -sf ${RooUnfoldHeaders} -t ${CMAKE_CURRENT_BINARY_DIR} )

if(NOT PlainROOT_BUILD)
  include(cmake/ATLAS.cmake)
endif()

if(NOT ATLAS_BUILD)
  include(cmake/PlainROOT.cmake)
endif()

include(cmake/testing.cmake)

include(cmake/doxygen.cmake)

