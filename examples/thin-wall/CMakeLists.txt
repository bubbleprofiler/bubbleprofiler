add_executable(thin
  thin_wall.cpp
  thin_wall_observer.cpp
  thin_wall_potential.cpp
)

target_link_libraries(thin
  BubbleProfiler::bubbler
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(thin
    PRIVATE ${CXX_WARNING_FLAGS}
    )
endif()

add_executable(BubbleProfiler::thin ALIAS thin)

add_dependencies(examples thin)
