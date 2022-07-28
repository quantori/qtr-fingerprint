set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)

add_executable(splitter_tree main.cpp)

target_include_directories(splitter_tree
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/c/indigo
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/abseil-cpp
        )

target_link_directories(splitter_tree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(splitter_tree indigo indigo-cpp OpenMP::OpenMP_CXX)
target_link_libraries(splitter_tree absl::flags absl::flags_parse)
target_link_libraries(splitter_tree qtrlib glog::glog)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(splitter_tree)
