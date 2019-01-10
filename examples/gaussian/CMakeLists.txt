add_executable(gaussian
  gaussian_potential_test.cpp
)

target_link_libraries(gaussian
  BubbleProfiler::bubbler
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(gaussian
    PRIVATE ${CXX_WARNING_FLAGS}
    )
endif()

add_executable(BubbleProfiler::gaussian ALIAS gaussian)

add_dependencies(examples gaussian)
