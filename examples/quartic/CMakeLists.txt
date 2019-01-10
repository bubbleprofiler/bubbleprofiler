add_executable(quartic
  action.cpp
)

target_link_libraries(quartic
  BubbleProfiler::bubbler
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(quartic PRIVATE ${CXX_WARNING_FLAGS})
endif()


add_executable(quartic_tabulate
  tabulate.cpp
)

target_link_libraries(quartic_tabulate
  BubbleProfiler::bubbler
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(quartic_tabulate PRIVATE ${CXX_WARNING_FLAGS})
endif()

add_executable(BubbleProfiler::quartic ALIAS quartic)
add_executable(BubbleProfiler::quartic_tabulate ALIAS quartic_tabulate)

add_dependencies(examples quartic quartic_tabulate)
