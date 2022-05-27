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
        )

target_link_libraries(${This} glog::glog)