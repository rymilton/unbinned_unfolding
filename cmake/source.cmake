file(GLOB RooUnfoldLinkDef ${ROOUNFOLD_HEADER_DIR}/*_LinkDef.h)
file(GLOB RooUnfoldSources ${ROOUNFOLD_SOURCE_DIR}/*.cxx)
file(GLOB RooUnfoldHeaders ${ROOUNFOLD_HEADER_DIR}/R*.h ${ROOUNFOLD_HEADER_DIR}/RooUnfold/TUnfold/*.h src/*.tpp)
file(GLOB RooUnfoldHeadersForLinkDef RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${ROOUNFOLD_HEADER_DIR}/*.h
  ${ROOUNFOLD_HEADER_DIR}/RooUnfold/TUnfold/*.h
  src/*.tpp)
file(GLOB RooUnfoldLinkDef2 RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src ${ROOUNFOLD_HEADER_DIR}/*_LinkDef.h)
list(REMOVE_ITEM RooUnfoldHeadersForLinkDef ${RooUnfoldLinkDef2})

file(GLOB RooUnfoldExecSources test/src/RooUnfoldTest.cxx test/src/RooUnfoldTest2D.cxx test/src/RooUnfoldTest3D.cxx)
file(GLOB RooUnfoldUnitTests test/*.cxx)

# for installing header files in separate (sub)directories
file(GLOB RooUnfoldHeadersOnly ${ROOUNFOLD_HEADER_DIR}/*.h)
file(GLOB TUnfoldHeaders ${ROOUNFOLD_HEADER_DIR}/RooUnfold/TUnfold/*.h)

if(RooUnfoldEnableProfiling)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS} -g -ftest-coverage -fprofile-arcs ")
else()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS} -g ")
endif()

set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)
