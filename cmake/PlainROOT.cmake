find_package( ROOT COMPONENTS Tree Unfold Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser)

if(${ROOT_FOUND})
  message("Setup using plain ROOT")    
  set(PlainROOT_BUILD 1)
  
  execute_process( COMMAND ln -sf ${RooUnfoldHeaders} -t ${CMAKE_CURRENT_BINARY_DIR} )
  set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
  file(WRITE ${SETUP} "#!/bin/bash\n")
  file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )
  
  
  # register all the files and directories
  include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/src")
  include_directories ("${ROOT_INCLUDE_DIRS}")
  
  file(APPEND ${SETUP} "export PATH=\${PATH}:${CMAKE_CURRENT_BINARY_DIR}\n")  
  file(APPEND ${SETUP} "export PYTHONPATH=\${PYTHONPATH}:${CMAKE_CURRENT_BINARY_DIR}\n")
  file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_CURRENT_BINARY_DIR}\n")  
  
  # generate the dictionary source code
  if(DEFINED ROOT_USE_FILE)
    include(${ROOT_USE_FILE})
  else()
    include("${ROOTSYS}/cmake/RootMacros.cmake")
  endif()
  ROOT_GENERATE_DICTIONARY(G__RooUnfold ${RooUnfoldHeaders} LINKDEF ${RooUnfoldLinkDef} OPTIONS ${EXTRA_FLAGS} -Wno-inconsistent-missing-override)
  
  # register the shared object to include both sources and dictionaries
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${RooUnfold_FLAGS}")  
  add_library( RooUnfold SHARED ${RooUnfoldSources}  G__RooUnfold.cxx)
  
  # link everything together at the end
  target_link_libraries( RooUnfold ${ROOT_LIBRARIES} )
  
  # Add all targets to the build-tree export set
  export(TARGETS RooUnfold FILE "${PROJECT_BINARY_DIR}/RooUnfoldTargets.cmake")
  
  # Export the package for use from the build-tree
  # (this registers the build-tree with a global CMake-registry)
  export(PACKAGE RooUnfold)
  
  set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
  set(CONF_LIBRARY_DIRS "${PROJECT_BINARY_DIR}")
  set(CONF_LIBRARIES    RooUnfold)
  configure_file(RooUnfoldConfig.cmake.in
    "${PROJECT_BINARY_DIR}/RooUnfoldConfig.cmake" @ONLY)
  
  # Install the RooUnfoldConfig.cmake and RooUnfoldConfigVersion.cmake
  install(FILES
    "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/RooUnfoldConfig.cmake"
    "${PROJECT_BINARY_DIR}/RooUnfoldConfigVersion.cmake"
    DESTINATION "${PROJECT_SOURCE_DIR}" COMPONENT dev)

  if(${RooUnfoldTests})    
    include(CTest)
    enable_testing()
    
    foreach(ExecSource ${RooUnfoldExecSources})
      get_filename_component(ExecName ${ExecSource} NAME_WE)    
      add_executable( ${ExecName} ${ExecSource} )
      target_link_libraries ( ${ExecName} RooUnfold ${ROOT_LIBRARIES} gcov)
    endforeach()
    
    add_subdirectory(test)
  endif()
  file(GLOB pyfiles "python/*.py")
  execute_process(
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/RooUnfold
    )
  foreach(pyfile ${pyfiles})
    execute_process(
      COMMAND ln -sf ${pyfile} ${CMAKE_CURRENT_BINARY_DIR}/RooUnfold
      )
  endforeach()
  
  
endif()
