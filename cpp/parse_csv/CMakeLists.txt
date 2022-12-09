set(CMAKE_CXX_STANDARD 20)
set(This parse_csv)
find_package(OpenMP REQUIRED)

add_executable(${This} main.cpp)

target_include_directories(${This}
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/c/indigo
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/abseil-cpp
        )

target_link_directories(${This} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(${This} indigo indigo-cpp OpenMP::OpenMP_CXX)
target_link_libraries(${This} absl::flags absl::flags_parse)
target_link_libraries(${This} qtrlib glog::glog)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(${This})
