# add compiler warnings to compiler flags

set(GNU_CXX_WARNING_FLAGS
  "-Wall"
  "-pedantic"
  "-Wextra"
  "-Wcast-align"
  "-Woverloaded-virtual"
  "-Wnon-virtual-dtor"
  "-Wint-in-bool-context"
)

set(CLANG_CXX_WARNING_FLAGS
  "-Wall"
  "-pedantic"
  "-Wextra"
  "-Wcast-align"
  "-Woverloaded-virtual"
  "-Wnon-virtual-dtor"
)

set(INTEL_CXX_WARNING_FLAGS
  "-w3"
  "-Wnon-virtual-dtor"
)

set(MSVC_WARNING_FLAGS
  "/W3"
)

set(WARNING_FLAGS_TO_ADD "")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  list(APPEND WARNING_FLAGS_TO_ADD ${GNU_CXX_WARNING_FLAGS})
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  list(APPEND WARNINGS_FLAGS_TO_ADD ${INTEL_CXX_WARNING_FLAGS})
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # remove default setting
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  endif()
  list(APPEND WARNING_FLAGS_TO_ADD ${MSVC_WARNING_FLAGS})
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  list(APPEND WARNING_FLAGS_TO_ADD ${CLANG_CXX_WARNING_FLAGS})
endif()

if(${WARNINGS_AS_ERRORS})
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    list(APPEND WARNING_FLAGS_TO_ADD "-Werror")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    list(APPEND WARNING_FLAGS_TO_ADD "-Werror")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    list(APPEND WARNING_FLAGS_TO_ADD "/WX")
  elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    list(APPEND WARNING_FLAGS_TO_ADD "-Werror")
  endif()
endif()

set(CXX_WARNING_FLAGS "")
foreach(flag ${WARNING_FLAGS_TO_ADD})
  add_cxx_flag_if_supported("${flag}" CXX_WARNING_FLAGS)
endforeach()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if(CXX_WARNING_FLAGS MATCHES "-Wint-in-bool-context")
    string(REGEX REPLACE "-Wint-in-bool-context" "-Wno-int-in-bool-context"
      CXX_WARNING_FLAGS  "${CXX_WARNING_FLAGS}")
  endif()
endif()

separate_arguments(CXX_WARNING_FLAGS)
