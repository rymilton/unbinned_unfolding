
foreach(UnitTest ${RooUnfoldUnitTests})
get_filename_component(Test ${UnitTest} NAME_WE)
add_executable(${Test} ${UnitTest})
target_link_libraries(${Test} PUBLIC tests gcov)
endforeach()

add_test(Bayes test_bayes)
#add_test(SVD test_svd)
#add_test(Invert test_invert)
#add_test(BinByBin test_bbb)

file(GLOB test_methods "test/test_methods.py")
add_test(
NAME test_methods
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
COMMAND python ${test_methods}
)

file(GLOB test_fakes "test/test_fakes.py")
add_test(
NAME test_fakes
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
COMMAND python ${test_fakes}
)

file(GLOB test_bin_correlation "test/test_bin_correlation.py")
add_test(
NAME test_bin_correlation
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
COMMAND python ${test_bin_correlation}
)

file(GLOB test_uncertainty "test/test_uncertainty.py")
add_test(
NAME test_uncertainty
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
COMMAND python ${test_uncertainty}
)

file(GLOB test_overflow "test/test_overflow.py")
add_test(
NAME test_overflow
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
COMMAND python ${test_overflow}
)

file(GLOB test_2D "test/test_2D.py")
add_test(
NAME test_2D
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
COMMAND python ${test_2D}
)

file(GLOB test_3D "test/test_3D.py")
add_test(
NAME test_3D
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
COMMAND python ${test_3D}
)