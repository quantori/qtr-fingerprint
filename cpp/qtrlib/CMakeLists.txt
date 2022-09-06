set(This qtrlib)

file(GLOB_RECURSE ${This}-sources
        "*.h"
        "*.c"
        "*.hpp"
        "*.cpp"
        )

add_library(${This} STATIC ${${This}-sources})

target_include_directories(${This}
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/Indigo/api/c/bingo-nosql
        PUBLIC ../third_party/Indigo/api/c/indigo

        PUBLIC profiler
        PUBLIC search_engine
        PUBLIC utils
        PUBLIC fingerprints
        PUBLIC splitter_tree
        PUBLIC columns_selection
        PUBLIC data_io
        PUBLIC io
        PUBLIC ball_tree
        PUBLIC ball_tree/split_bit_selectors
        )

target_link_libraries(${This} glog::glog)