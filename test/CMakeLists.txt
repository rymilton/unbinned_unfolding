if((NOT ${RooUnfoldTests}) OR (TARGET RooUnfoldTests))
  return()
endif()

include_directories ("${ROOT_INCLUDE_DIRS}")
add_library(RooUnfoldTests
        src/generate_test.cxx
        src/response_test.cxx
        src/bayes_test.cxx
        src/RooUnfold_unit_tests.cxx)
target_link_libraries(RooUnfoldTests PUBLIC ROOT::Matrix ROOT::Core ROOT::Hist RooUnfold)
