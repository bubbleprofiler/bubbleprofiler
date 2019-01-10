add_executable(logarithmic
  logarithmic.cpp
  logarithmic_observer.cpp
  logarithmic_potential.cpp
)

target_link_libraries(logarithmic
  BubbleProfiler::bubbler
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(logarithmic
    PRIVATE ${CXX_WARNING_FLAGS}
    )
endif()

add_executable(BubbleProfiler::logarithmic ALIAS logarithmic)

add_dependencies(examples logarithmic)
