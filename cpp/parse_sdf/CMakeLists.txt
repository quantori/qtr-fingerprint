set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

add_executable(parse_sdf main.cpp)

target_include_directories(parse_sdf
        PUBLIC .
        )

target_link_directories(parse_sdf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(parse_sdf bingo-nosql-object OpenMP::OpenMP_CXX)
target_link_libraries(parse_sdf absl::flags absl::flags_parse)
target_link_libraries(parse_sdf qtrlib glog::glog)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(parse_sdf)
