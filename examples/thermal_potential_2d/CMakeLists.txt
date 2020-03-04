add_executable(thermal_potential_2d
  thermal_potential.cpp simple_2d_potential.cpp
)

target_link_libraries(thermal_potential_2d
  BubbleProfiler::bubbler
  effpotential
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(thermal_potential_2d
    PRIVATE ${CXX_WARNING_FLAGS}
    )
endif()

add_executable(BubbleProfiler::thermal_potential_2d ALIAS thermal_potential_2d)

add_dependencies(examples thermal_potential_2d)
