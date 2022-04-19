if (${ROOT_VERSION_MAJOR} VERSION_LESS 6.24)
    set(RooUnfold_FLAGS "${RooUnfold_FLAGS} -DNO_WRAPPERPDF")
endif()

