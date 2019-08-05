cmake_minimum_required(VERSION 2.8.12)

project(BubbleProfiler)
enable_language(CXX)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(BUBBLEPROFILER_MAJOR_VERSION "1")
set(BUBBLEPROFILER_MINOR_VERSION "0")
set(BUBBLEPROFILER_PATCH_VERSION "0")
set(BUBBLEPROFILER_VERSION "${BUBBLEPROFILER_MAJOR_VERSION}.${BUBBLEPROFILER_MINOR_VERSION}.${BUBBLEPROFILER_PATCH_VERSION}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

include(set_cxx11)
set_cxx11()

# destinations for build products
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
endif()
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
endif()

include(BubbleProfilerBuildOptions)

# destinations for installed files
if(NOT BUBBLEPROFILER_INSTALL_ARCHIVE_DIR)
  set(BUBBLEPROFILER_INSTALL_ARCHIVE_DIR lib
    CACHE PATH
    "The directory relative to CMAKE_INSTALL_PREFIX where BubbleProfiler archives are installed")
endif()

if(NOT BUBBLEPROFILER_INSTALL_CMAKE_DIR)
  set(BUBBLEPROFILER_INSTALL_CMAKE_DIR share/BubbleProfiler
    CACHE PATH
    "The directory relative to CMAKE_INSTALL_PREFIX where BubbleProfiler CMake files are installed")
endif()

if(NOT BUBBLEPROFILER_INSTALL_INCLUDE_DIR)
  set(BUBBLEPROFILER_INSTALL_INCLUDE_DIR include
    CACHE PATH
    "The directory relative to CMAKE_INSTALL_PREFIX where BubbleProfiler headers are installed")
endif()

if(NOT BUBBLEPROFILER_INSTALL_LIBRARY_DIR)
  set(BUBBLEPROFILER_INSTALL_LIBRARY_DIR lib
    CACHE PATH
    "The directory relative to CMAKE_INSTALL_PREFIX where BubbleProfiler libraries are installed")
endif()

if(NOT BUBBLEPROFILER_INSTALL_RUNTIME_DIR)
  set(BUBBLEPROFILER_INSTALL_RUNTIME_DIR bin
    CACHE PATH
    "The directory relative to CMAKE_INSTALL_PREFIX where BubbleProfiler executables are installed")
endif()

# if not using the default executable suffix for the platform,
# set it to ".x"
if(NOT USE_NATIVE_EXE_SUFFIX)
  if("${CMAKE_EXECUTABLE_SUFFIX}" STREQUAL "")
    set(CMAKE_EXECUTABLE_SUFFIX .x)
  endif()
endif()

find_package(Boost 1.53.0 REQUIRED
  COMPONENTS
  filesystem
  program_options
  system
  )
find_package(CLN 1.2.2 REQUIRED)
find_package(Eigen3 3.1 REQUIRED)
find_package(GiNaC REQUIRED)
find_package(GSL 1.15 REQUIRED)
find_package(NLopt REQUIRED)
find_package(Threads REQUIRED)
find_package(effpotential REQUIRED)

# add compiler warnings if requested
if(ENABLE_COMPILER_WARNINGS)
  include(add_cxx_flag_if_supported)
  include(CompilerWarnings)
endif(ENABLE_COMPILER_WARNINGS)

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC BUBBLEPROFILER_GIT_COMMIT)
if(BUBBLEPROFILER_GIT_COMMIT MATCHES "GITDIR-NOTFOUND")
  set(BUBBLEPROFILER_GIT_COMMIT "")
endif()

configure_file(
  include/bubbleprofiler_version.hpp.in
  ${PROJECT_BINARY_DIR}/include/bubbleprofiler_version.hpp
  @ONLY
)

add_library(bubbler
  src/action.cpp
  src/algebraic_potential.cpp
  src/thermal_potential.cpp
  src/basic_logger.cpp
  src/euclidean_action.cpp
  src/field_profiles.cpp
  src/gaussian_potential.cpp
  src/gsl_interpolator.cpp
  src/instream_profile_guesser.cpp
  src/integration_utils.cpp
  src/kink_profile_guesser.cpp
  src/log_message.cpp
  src/logging_manager.cpp
  src/nlopt_optimizer.cpp
  src/numeric.cpp
  src/observers.cpp
  src/perturbations_ode_system.cpp
  src/relative_convergence_tester.cpp
  src/restricted_quartic_potential.cpp
  src/rotation.cpp
  src/shooting.cpp
  src/shooting_profile_guesser.cpp
)

set(bubbler_headers_to_install
  include/action.hpp
  include/algebraic_potential.hpp
  include/thermal_potential.hpp
  include/basic_logger.hpp
  ${PROJECT_BINARY_DIR}/include/bubbleprofiler_version.hpp
  include/default_integration_policy.hpp
  include/eigen_state_utils.hpp
  include/error.hpp
  include/euclidean_action.hpp
  include/field_profiles.hpp
  include/gaussian_potential.hpp
  include/generic_perturbative_profiler.hpp
  include/gsl_interpolator.hpp
  include/gsl_root_finder.hpp
  include/instream_profile_guesser.hpp
  include/integration_policy.hpp
  include/integration_utils.hpp
  include/kink_profile_guesser.hpp
  include/logging_manager.hpp
  include/log_message.hpp
  include/math_wrappers.hpp
  include/nlopt_optimizer.hpp
  include/numeric.hpp
  include/observers.hpp
  include/perturbations_ode_system.hpp
  include/perturbative_profiler.hpp
  include/potential.hpp
  include/profile_convergence_tester.hpp
  include/profile_guesser.hpp
  include/raii_guard.hpp
  include/relative_convergence_tester.hpp
  include/restricted_quartic_potential.hpp
  include/root_finder.hpp
  include/rotation.hpp
  include/shooting.hpp
  include/shooting_profile_guesser.hpp
  include/univariate_interpolation.hpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(bubbler
    PRIVATE ${CXX_WARNING_FLAGS})
endif()

target_include_directories(bubbler
  PUBLIC
  ${Boost_INCLUDE_DIRS}
  ${CLN_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${GINAC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${NLopt_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${BUBBLEPROFILER_INSTALL_INCLUDE_DIR}>
)

target_link_libraries(bubbler
  INTERFACE
  ${NLopt_LIBRARIES}
  ${GINAC_LIBRARIES}
  ${CLN_LIBRARIES}
  ${Boost_LIBRARIES}
  ${GSL_LIBRARIES}
  ${CMAKE_DL_LIBS}
)

add_executable(bubbler_run
  src/run_cmd_line_potential.cpp
)

set_target_properties(bubbler_run
  PROPERTIES OUTPUT_NAME run_cmd_line_potential
)

target_include_directories(bubbler_run PRIVATE ${Boost_INCLUDE_DIRS})

if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(bubbler_run "-pthread")
endif()

target_link_libraries(bubbler_run
    bubbler
    ${Boost_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(bubbler_run
    PRIVATE ${CXX_WARNING_FLAGS})
endif()

install(TARGETS
  bubbler bubbler_run
  EXPORT BubbleProfilerTargets
  RUNTIME DESTINATION ${BUBBLEPROFILER_INSTALL_RUNTIME_DIR}
  LIBRARY DESTINATION ${BUBBLEPROFILER_INSTALL_LIBRARY_DIR}
  ARCHIVE DESTINATION ${BUBBLEPROFILER_INSTALL_ARCHIVE_DIR}
)

install(FILES
  ${bubbler_headers_to_install}
  DESTINATION ${BUBBLEPROFILER_INSTALL_INCLUDE_DIR}
)

add_library(BubbleProfiler::bubbler ALIAS bubbler)
add_executable(BubbleProfiler::bubbler_run ALIAS bubbler_run)

export(TARGETS bubbler bubbler_run
  NAMESPACE BubbleProfiler::
  FILE "${PROJECT_BINARY_DIR}/${BUBBLEPROFILER_INSTALL_CMAKE_DIR}/BubbleProfilerTargets.cmake"
)
export(PACKAGE BubbleProfiler)

add_subdirectory(docs EXCLUDE_FROM_ALL)
add_subdirectory(examples EXCLUDE_FROM_ALL)

if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(test EXCLUDE_FROM_ALL)
endif(ENABLE_TESTING)

# Custom targets for generating archived source distributions
if(BUBBLEPROFILER_GIT_COMMIT)
  add_custom_target(archive-version
    COMMAND "${CMAKE_COMMAND}"
    -DGIT_ARCHIVE_REPO="${PROJECT_SOURCE_DIR}"
    -DGIT_ARCHIVE_REF="v${BUBBLEPROFILER_VERSION}"
    -DGIT_ARCHIVE_PREFIX_ROOT="${PROJECT_NAME}"
    -DGIT_ARCHIVE_OUTPUT_ROOT="${PROJECT_BINARY_DIR}/${PROJECT_NAME}"
    -DGIT_ARCHIVE_EXTENSION="tar.gz"
    -P "${CMAKE_MODULE_PATH}/git_archive.cmake")

  add_custom_target(archive-head
    COMMAND "${CMAKE_COMMAND}"
    -DGIT_ARCHIVE_REPO="${PROJECT_SOURCE_DIR}"
    -DGIT_ARCHIVE_REF="HEAD"
    -DGIT_ARCHIVE_PREFIX_ROOT="${PROJECT_NAME}"
    -DGIT_ARCHIVE_OUTPUT_ROOT="${PROJECT_BINARY_DIR}/${PROJECT_NAME}"
    -DGIT_ARCHIVE_EXTENSION="tar.gz"
    -P "${CMAKE_MODULE_PATH}/git_archive.cmake")
endif()

include(BubbleProfilerGenerateConfig)
generate_bubbleprofiler_config_files()
