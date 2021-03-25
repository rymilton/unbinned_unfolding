#
# This is a template for a CMakeLists.txt file that can be used in a client
# project (work area) to set up building ATLAS packages against the configured
# release.
#

find_package(AnalysisBase QUIET)
find_package(AthAnalysis QUIET)

if (${AnalysisBase_FOUND} OR ${AthAnalysis_FOUND})
  message ("Configuring for build within analysis release")
  set (ATLAS_BUILD 1)
else() 
  set (ATLAS_BUILD 0)
endif()

if(${ATLAS_BUILD})

# Set the minimum required CMake version:
cmake_minimum_required( VERSION 3.0 FATAL_ERROR )

#register the package_filter file as a dependent
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS $ENV{ACM_PACKAGE_FILTER_FILE} )
set( ATLAS_PACKAGE_FILTER_FILE $ENV{ACM_PACKAGE_FILTER_FILE} )

# If there's a directory called AtlasCMake in the project,
# and the user didn't specify AtlasCMake_DIR yet, then let's
# give it a default value.
if( IS_DIRECTORY ${CMAKE_SOURCE_DIR}/Build/AtlasCMake AND
      NOT AtlasCMake_DIR AND NOT ENV{AtlasCMake_DIR} )
   set( AtlasCMake_DIR ${CMAKE_SOURCE_DIR}/Build/AtlasCMake )
endif()

# If there's a directory called AtlasLCG in the project,
# and the user didn't specify LCG_DIR yet, then let's
# give it a default value.
if( IS_DIRECTORY ${CMAKE_SOURCE_DIR}/Build/AtlasLCG AND
      NOT LCG_DIR AND NOT ENV{LCG_DIR} )
   set( LCG_DIR ${CMAKE_SOURCE_DIR}/Build/AtlasLCG )
endif()

# Pick up a local version of the AtlasCMake code if it exists:
find_package( AtlasCMake QUIET )

# Find the project that we depend on:
find_package( $ENV{AtlasProject} )

# Set up CTest:
atlas_ctest_setup()

# These next lines are a temporary fix for issue: ATLINFR-2388
if( DEFINED ENV{WorkDir_PLATFORM} )
    set( ATLAS_PLATFORM $ENV{WorkDir_PLATFORM} )
endif()

# Set up a work directory project:
atlas_project( WorkDir 1.0.0
   USE $ENV{AtlasProject} $ENV{AtlasVersion}
   FORTRAN $ENV{ATLAS_PROJECT_EXTRA})

# Set up the runtime environment setup script(s):
lcg_generate_env( SH_FILE ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh )
install( FILES ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh
   DESTINATION . )

  atlas_subdir( RooUnfold )

  atlas_add_root_dictionary( RooUnfold _dictSource
    ROOT_HEADERS
    ${RooUnfoldHeaders}
    ${RooUnfoldLinkDef})

  # Ensure that the 'RooUnfold' folder expected by atlas_add_library exists at source level
  # by linking the appropriate files from src/ 
  set( _RooUnfold_header_dir ${CMAKE_CURRENT_SOURCE_DIR}/RooUnfold )
  file(MAKE_DIRECTORY ${_RooUnfold_header_dir})
  execute_process( COMMAND ln -sf ${RooUnfoldHeaders} -t ${_RooUnfold_header_dir} )

  
  atlas_add_library( RooUnfold
    ${RooUnfoldHeaders} ${RooUnfoldSources} ${_dictSource}
    PUBLIC_HEADERS RooUnfold
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES})

  foreach(ExecSource ${RooUnfoldExecSources})
    get_filename_component(ExecName ${ExecSource} NAME_WE)    
    atlas_add_executable( ${ExecName} ${ExecSource}
      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} RooUnfold ${CMAKE_CURRENT_SOURCE_DIR}/examples ${CMAKE_CURRENT_SOURCE_DIR}/src
      LINK_LIBRARIES ${ROOT_LIBRARIES} RooUnfold)      
  endforeach()

 
# Set up CPack:
atlas_cpack_setup()

endif()
