set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

add_executable(playground main.cpp)

target_include_directories(playground
        PUBLIC .
        )

target_link_directories(playground PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(playground OpenMP::OpenMP_CXX)
target_link_libraries(playground absl::flags absl::flags_parse)
target_link_libraries(playground qtrlib glog::glog)
target_link_libraries(playground bingo-nosql-object)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(playground)
