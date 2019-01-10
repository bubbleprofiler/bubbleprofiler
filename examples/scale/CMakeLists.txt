add_executable(scale
  scale.cpp
)

target_link_libraries(scale
  BubbleProfiler::bubbler
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(scale
    PRIVATE ${CXX_WARNING_FLAGS}
    )
endif()

add_executable(BubbleProfiler::scale ALIAS scale)

add_dependencies(examples scale)
