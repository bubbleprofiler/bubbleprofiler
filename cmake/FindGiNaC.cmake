# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#.rst:
# FindGiNaC
# ----------
#
# Find GiNaC include dirs and libraries.
#
# This module sets the following variables:
#
# ::
#
#   GINAC_FOUND - set to true if the library is found
#   GINAC_INCLUDE_DIRS - list of required include directories
#   GINAC_LIBRARIES - list of libraries to be linked
#

find_package(CLN 1.2.2 REQUIRED)

find_package(PkgConfig)
pkg_check_modules(PC_GINAC QUIET GINAC)

find_path(GINAC_INCLUDE_DIR
  NAMES ginac/ginac.h
  PATHS ${PC_GINAC_INCLUDE_DIRS}
)

find_library(GINAC_LIBRARY
  NAMES libginac ginac
  PATHS ${PC_GINAC_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GINAC
  FOUND_VAR GINAC_FOUND
  REQUIRED_VARS
  GINAC_LIBRARY
  GINAC_INCLUDE_DIR
  VERSION_VAR GINAC_VERSION
)

if(GINAC_FOUND)
  set(GINAC_LIBRARIES ${GINAC_LIBRARY})
  set(GINAC_INCLUDE_DIRS ${GINAC_INCLUDE_DIR})
  set(GINAC_DEFINITIONS ${PC_GINAC_CFLAGS_OTHER})
endif()

mark_as_advanced(
  GINAC_INCLUDE_DIR
  GINAC_LIBRARY
)

