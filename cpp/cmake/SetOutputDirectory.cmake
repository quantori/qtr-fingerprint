function(set_output_directory target)
    
    # COPY AND PASTE FROM INDIGO setup.cmake
    ###############################################################################
    # If cross-compiling, we are able to set another OS
    if (NOT CMAKE_SYSTEM_NAME_LOWER)
        string(TOLOWER ${CMAKE_SYSTEM_NAME} CMAKE_SYSTEM_NAME_LOWER)
    endif()
    # If cross-compiling, we are able to set another Arch, like "aarch64"
    if (NOT CMAKE_SYSTEM_PROCESSOR_LOWER)
        string(TOLOWER ${CMAKE_SYSTEM_PROCESSOR} CMAKE_SYSTEM_PROCESSOR_LOWER)
    endif()
    if (CMAKE_SYSTEM_PROCESSOR_LOWER STREQUAL "amd64")
        set(CMAKE_SYSTEM_PROCESSOR_LOWER "x86_64")
    endif()
    if (CMAKE_SYSTEM_PROCESSOR_LOWER STREQUAL "arm64")
        set(CMAKE_SYSTEM_PROCESSOR_LOWER "aarch64")
    endif()
    if (CMAKE_SYSTEM_NAME_LOWER STREQUAL "msys")
        set(CMAKE_SYSTEM_NAME_LOWER "windows")
    endif()
    ###############################################################################

    set(${target}_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../dist/lib/${CMAKE_SYSTEM_NAME_LOWER}-${CMAKE_SYSTEM_PROCESSOR_LOWER})

    set_target_properties(${target} PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${${target}_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${${target}_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${${target}_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${${target}_RUNTIME_OUTPUT_DIRECTORY})

endfunction()