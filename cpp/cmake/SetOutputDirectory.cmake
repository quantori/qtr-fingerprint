function(set_output_directory target)

    set(${target}_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../build/bin/)

    set_target_properties(${target} PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${${target}_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${${target}_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${${target}_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${${target}_RUNTIME_OUTPUT_DIRECTORY})

endfunction()