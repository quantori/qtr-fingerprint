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

file(GLOB Headers
        "*.h")

file(GLOB Sources
        "*.cpp")

message("Headers = ${Headers}")
message("Sources = ${Sources}")

add_executable(app main.cpp ${Headers} ${Sources})

target_include_directories(app
    PUBLIC .
    PUBLIC ../third_party/Indigo/api/cpp/src
)

target_link_directories(app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(app indigo-cpp)

set(APP_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../dist/lib/${CMAKE_SYSTEM_NAME_LOWER}-${CMAKE_SYSTEM_PROCESSOR_LOWER})

set_target_properties(app PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${APP_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${APP_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${APP_RUNTIME_OUTPUT_DIRECTORY}
                      RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${APP_RUNTIME_OUTPUT_DIRECTORY}
)