configure_file(
  textfile_ansatz.sh.in
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/textfile_ansatz.sh
  @ONLY
  )

add_custom_target(textfile_ansatz)
add_dependencies(textfile_ansatz BubbleProfiler::bubbler_run)

add_custom_command(TARGET textfile_ansatz
  COMMAND
  ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/AnsatzTab.csv
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/AnsatzTab.csv
)

add_dependencies(examples textfile_ansatz)
