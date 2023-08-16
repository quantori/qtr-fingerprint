set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

file(GLOB_RECURSE sources *.cpp *.h)

add_executable(build_db ${sources})

target_include_directories(build_db
        PUBLIC .
        )

target_link_directories(build_db PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(build_db OpenMP::OpenMP_CXX)
target_link_libraries(build_db absl::flags absl::flags_parse)
target_link_libraries(build_db qtrlib glog::glog)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(build_db)
