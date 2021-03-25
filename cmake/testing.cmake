file(GLOB Tests "test/*.sh")
foreach(TestScript ${Tests})
  get_filename_component(TestName ${TestScript} NAME)
  add_test(
    NAME ${TestName}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND bash ${TestScript}
    )
endforeach()

