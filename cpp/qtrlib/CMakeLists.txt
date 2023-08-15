cmake_minimum_required(VERSION 3.12)

set(This qtrlib)

file(GLOB_RECURSE ${This}-sources
        "*.h"
        "*.c"
        "*.hpp"
        "*.cpp")


add_library(${This} STATIC ${${This}-sources})

target_link_libraries(${This} glog::glog)
target_link_libraries(${This} indigo-core indigo-cpp bingo-nosql)
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
        )