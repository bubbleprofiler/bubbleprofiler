find_package(Doxygen)

set(BUBBLEPROFILER_DOC_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/../README.md"
  "${CMAKE_CURRENT_SOURCE_DIR}/../include"
  "${CMAKE_CURRENT_SOURCE_DIR}/../src"
  "${CMAKE_CURRENT_SOURCE_DIR}/../examples"
)

if(DOXYGEN_FOUND)
  set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  set(DOXYGEN_INPUT "")
  foreach(s ${BUBBLEPROFILER_DOC_SOURCES})
    set(DOXYGEN_INPUT "${DOXYGEN_INPUT} \"${s}\"")
  endforeach()

  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

  add_custom_target(docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
endif()
