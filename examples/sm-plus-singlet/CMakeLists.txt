configure_file(
  sm_plus_singlet.py.in
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sm_plus_singlet.py
  @ONLY
)

add_custom_target(sm_plus_singlet)
add_dependencies(sm_plus_singlet BubbleProfiler::bubbler_run)

add_custom_command(TARGET sm_plus_singlet
  COMMAND
  ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/sm_plus_singlet_tn_plot.py
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sm_plus_singlet_tn_plot.py
)

add_dependencies(examples sm_plus_singlet)
