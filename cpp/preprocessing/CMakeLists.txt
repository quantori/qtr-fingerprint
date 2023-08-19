set(CMAKE_CXX_STANDARD 20)
set(This preprocessing)

file(GLOB_RECURSE sources *.cpp *.h)

add_executable(${This} main.cpp ${sources})

target_include_directories(${This}
        PUBLIC .
        PUBLIC preprocessors
        )

target_link_directories(${This} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(${This} bingo-nosql-object)
target_link_libraries(${This} absl::flags absl::flags_parse)
target_link_libraries(${This} qtrlib glog::glog)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(${This})
