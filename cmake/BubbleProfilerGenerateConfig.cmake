function(generate_bubbleprofiler_config_files)

  include(CMakePackageConfigHelpers)

  configure_package_config_file(
    "${CMAKE_MODULE_PATH}/BubbleProfilerConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/${BUBBLEPROFILER_INSTALL_CMAKE_DIR}/BubbleProfilerConfig.cmake"
    INSTALL_DESTINATION "${BUBBLEPROFILER_INSTALL_CMAKE_DIR}"
    PATH_VARS
    BUBBLEPROFILER_INSTALL_ARCHIVE_DIR
    BUBBLEPROFILER_INSTALL_CMAKE_DIR
    BUBBLEPROFILER_INSTALL_INCLUDE_DIR
    BUBBLEPROFILER_INSTALL_LIBRARY_DIR
    BUBBLEPROFILER_INSTALL_RUNTIME_DIR
    )

  install(FILES "${PROJECT_BINARY_DIR}/${BUBBLEPROFILER_INSTALL_CMAKE_DIR}/BubbleProfilerConfig.cmake"
    DESTINATION ${BUBBLEPROFILER_INSTALL_CMAKE_DIR})

  install(EXPORT BubbleProfilerTargets
    NAMESPACE BubbleProfiler::
    DESTINATION ${BUBBLEPROFILER_INSTALL_CMAKE_DIR})

  write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/${BUBBLEPROFILER_INSTALL_CMAKE_DIR}/BubbleProfilerConfigVersion.cmake"
    VERSION "${BUBBLEPROFILER_VERSION}"
    COMPATIBILITY SameMajorVersion)

  install(FILES "${PROJECT_BINARY_DIR}/${BUBBLEPROFILER_INSTALL_CMAKE_DIR}/BubbleProfilerConfigVersion.cmake"
    DESTINATION ${BUBBLEPROFILER_INSTALL_CMAKE_DIR})

endfunction()
