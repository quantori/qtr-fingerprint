cmake_minimum_required(VERSION 3.12)

set(This qtrlib)

option(QTR_PROFILING "" OFF)

file(GLOB_RECURSE ${This}-sources
        "*.h"
        "*.c"
        "*.hpp"
        "*.cpp")


add_library(${This} STATIC ${${This}-sources})

if (QTR_PROFILING)
    target_compile_definitions(${This} PUBLIC QTR_PROFILING)
    message(STATUS "Enabling QTR_PROFILING because QTR_PROFILING=ON")
endif (QTR_PROFILING)

target_link_libraries(${This} glog::glog)
target_link_libraries(${This} indigo-core indigo-cpp bingo-nosql bingo-nosql-object)
target_link_libraries(${This} absl::flags absl::flags_parse)

target_include_directories(${This}
        PUBLIC .
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
        PUBLIC mmf
        )