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
        PUBLIC ../third_party/rdkit/Code

        PUBLIC profiler
        PUBLIC search_engine
        PUBLIC utils
        PUBLIC fingerprints
        PUBLIC data_io
        PUBLIC io
        PUBLIC ball_tree
        PUBLIC smiles_table
        PUBLIC answer_filtering
        )

target_link_libraries(${This} glog::glog)
target_link_libraries(parse_sdf indigo indigo-cpp)
