# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#.rst:
# FindCLN
# ----------
#
# Find CLN include dirs and libraries.
#
# This module sets the following variables:
#
# ::
#
#   CLN_FOUND - set to true if the library is found
#   CLN_INCLUDE_DIRS - list of required include directories
#   CLN_LIBRARIES - list of libraries to be linked
#

function(_cln_get_version _cln_version_h _found_major _found_minor _found_patch)
  file(READ "${_cln_version_h}" _cln_version_header)

  string(REGEX MATCH "define[ \t]+CL_VERSION_MAJOR[ \t]+([0-9]+)" _cln_major_version_match "${_cln_version_header}")
  set(${_found_major} "${CMAKE_MATCH_1}" PARENT_SCOPE)
  string(REGEX MATCH "define[ \t]+CL_VERSION_MINOR[ \t]+([0-9]+)" _cln_minor_version_match "${_cln_version_header}")
  set(${_found_minor} "${CMAKE_MATCH_1}" PARENT_SCOPE)
  string(REGEX MATCH "define[ \t]+CL_VERSION_PATCHLEVEL[ \t]+([0-9]+)" _cln_patch_version_match "${_cln_version_header}")
  set(${_found_patch} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction(_cln_get_version)

find_package(PkgConfig)
pkg_check_modules(PC_CLN QUIET CLN)

find_path(CLN_INCLUDE_DIR
  NAMES cln/cln.h
  PATHS ${PC_CLN_INCLUDE_DIRS}
)

find_library(CLN_LIBRARY
  NAMES libcln cln
  PATHS ${PC_CLN_LIBRARY_DIRS}
)

if(CLN_INCLUDE_DIR)
  _cln_get_version(
    ${CLN_INCLUDE_DIR}/cln/version.h
    CLN_VERSION_MAJOR
    CLN_VERSION_MINOR
    CLN_VERSION_PATCHLEVEL)
  set(CLN_VERSION ${CLN_VERSION_MAJOR}.${CLN_VERSION_MINOR}.${CLN_VERSION_PATCHLEVEL})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLN
  FOUND_VAR CLN_FOUND
  REQUIRED_VARS
  CLN_LIBRARY
  CLN_INCLUDE_DIR
  VERSION_VAR CLN_VERSION
)

if(CLN_FOUND)
  set(CLN_LIBRARIES ${CLN_LIBRARY})
  set(CLN_INCLUDE_DIRS ${CLN_INCLUDE_DIR})
  set(CLN_DEFINITIONS ${PC_CLN_CFLAGS_OTHER})
endif()

mark_as_advanced(
  CLN_INCLUDE_DIR
  CLN_LIBRARY
)
