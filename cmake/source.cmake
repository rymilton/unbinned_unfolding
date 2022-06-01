file(GLOB RooUnfoldLinkDef RELATIVE ${CMAKE_SOURCE_DIR}/src ${ROOUNFOLD_HEADER_DIR}/*_LinkDef.h)
file(GLOB RooUnfoldSources RELATIVE ${CMAKE_SOURCE_DIR} ${ROOUNFOLD_SOURCE_DIR}/*.cxx)
file(GLOB RooUnfoldHeaders RELATIVE ${CMAKE_SOURCE_DIR}/src ${ROOUNFOLD_HEADER_DIR}/*.h src/*.tpp)
list(REMOVE_ITEM RooUnfoldHeaders ${RooUnfoldLinkDef})

file(GLOB RooUnfoldExecSources test/src/RooUnfoldTest.cxx test/src/RooUnfoldTest2D.cxx test/src/RooUnfoldTest3D.cxx)
file(GLOB RooUnfoldUnitTests test/*.cxx)

# -fprofile-arcs -ftest-coverage
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS} -g -ftest-coverage  -fprofile-arcs ")
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)
