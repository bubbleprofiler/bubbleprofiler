# adds a C++ compiler flag if it is supported by the current compiler

include(CheckCXXCompilerFlag)
include(CMakeParseArguments)

function(add_cxx_flag_if_supported flag output_flags)
  CMAKE_PARSE_ARGUMENTS(add_flag_if_supported "REQUIRED" "" "" ${ARGN})
  string(REGEX REPLACE "[ -/=:]" "_" TEST_FLAG_NAME "${flag}")
  unset(HAS_${TEST_FLAG_NAME})
  CHECK_CXX_COMPILER_FLAG("${flag}" HAS_${TEST_FLAG_NAME})
  if(add_flag_if_supported_REQUIRED AND NOT HAS_${TEST_FLAG_NAME})
    message(FATAL_ERROR
      "The C++ compiler does not support the required flag \"${flag}\"")
  endif()
  if(HAS_${TEST_FLAG_NAME})
    set(${output_flags} "${${output_flags}} ${flag}" PARENT_SCOPE)
  else()
    message(STATUS "C++ compiler does not support ${flag}")
  endif()
endfunction(add_cxx_flag_if_supported)
