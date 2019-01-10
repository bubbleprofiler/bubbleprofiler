configure_file(
  example_2field.sh.in
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/example_2field.sh
  @ONLY
)

configure_file(
  example_3field.sh.in
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/example_3field.sh
  @ONLY
)

add_custom_target(perturbative)
add_dependencies(perturbative BubbleProfiler::bubbler_run)

add_custom_command(TARGET perturbative
  COMMAND
  ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/make_plots.py
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/make_plots.py
)

add_dependencies(examples perturbative)
