cmake_minimum_required(VERSION 3.12)

set(This qtrlib)

file(GLOB_RECURSE ${This}-sources
        "*.h"
        "*.c"
        "*.hpp"
        "*.cpp"
        "../third_party/Indigo/api/c/bingo-nosql/src/*.cpp"
        "../third_party/Indigo/api/c/indigo/src/*.cpp")


add_library(${This} STATIC ${${This}-sources})

include(../third_party/Indigo/api/indigo-version.cmake)
configure_file(../third_party/Indigo/api/c/indigo/src/indigo_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/indigo_version.h)

target_link_libraries(${This} glog::glog)
target_link_libraries(${This} indigo indigo-cpp bingo-nosql indigo-core)
target_link_libraries(${This} absl::flags absl::flags_parse)

# Include the directory containing indigo_version.h to the include directories of the target
target_include_directories(${This} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/src)

include(../third_party/Indigo/api/indigo-version.cmake)
configure_file(../third_party/Indigo/api/c/indigo/src/indigo_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/indigo_version.h)

target_include_directories(${This}
        PUBLIC .
        PUBLIC ../third_party/abseil-cpp
        PUBLIC profiler
        PUBLIC search_engine
        PUBLIC utils
        PUBLIC fingerprints
        PUBLIC data_io
        PUBLIC io
        PUBLIC ball_tree
        PUBLIC smiles_table
        PUBLIC answer_filtering
        PUBLIC query_data
        )