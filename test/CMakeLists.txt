set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/bin)

add_executable(bubble_unit
  catch_main.cpp
  mock_potentials.cpp
  test_euclidean_action.cpp
  test_field_profiles.cpp
  test_gaussian_potential.cpp
  test_instream_profile_guesser.cpp
  test_kink_profile_guesser.cpp
  test_perturbative_profiler_kink_ansatz.cpp
  test_root_finders.cpp
  test_rotation.cpp
  test_spline_interpolation.cpp
  test_undershoot_overshoot_profiler.cpp
  test_univariate_interpolation.cpp
)

target_include_directories(bubble_unit
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/catch
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(bubble_unit PRIVATE ${CXX_WARNING_FLAGS})
endif()

if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(bubble_unit "-pthread")
endif()

target_link_libraries(bubble_unit
  BubbleProfiler::bubbler
  ${CMAKE_THREAD_LIBS_INIT})

include(ParseAndAddCatchTests)
ParseAndAddCatchTests(bubble_unit)

add_executable(BubbleProfiler::bubble_unit ALIAS bubble_unit)

add_custom_target(check COMMAND "ctest")
add_dependencies(check bubble_unit)
