add_executable(fubini
  general_fubini.cpp
  generalized_fubini_observer.cpp
  generalized_fubini_potential.cpp
)

target_link_libraries(fubini
  BubbleProfiler::bubbler
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(fubini
    PRIVATE ${CXX_WARNING_FLAGS}
    )
endif()

add_executable(BubbleProfiler::fubini ALIAS fubini)

add_dependencies(examples fubini)
