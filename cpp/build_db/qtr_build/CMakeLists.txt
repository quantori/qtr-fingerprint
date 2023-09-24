cmake_minimum_required(VERSION 3.12)

set(This qtr_build)

file(GLOB_RECURSE ${This}-sources
        "*.h"
        "*.c"
        "*.hpp"
        "*.cpp")

add_library(${This} STATIC ${${This}-sources})

target_link_libraries(${This} qtrlib)

target_include_directories(${This}
        PUBLIC .
        PUBLIC database_builders
        )