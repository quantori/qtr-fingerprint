cmake_minimum_required(VERSION 3.12)

set(This qtr_run)

file(GLOB_RECURSE ${This}-sources
        "*.h"
        "*.c"
        "*.hpp"
        "*.cpp")

add_library(${This} STATIC ${${This}-sources})

target_link_libraries(${This} PRIVATE qtrlib)
target_link_libraries(${This} PRIVATE Crow::Crow)
target_link_libraries(${This} PRIVATE bingo-nosql-object)

target_include_directories(${This}
        PUBLIC .
        PUBLIC search_data
        PUBLIC modes
        )