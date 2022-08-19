set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

add_executable(fingerprints main.cpp)

target_include_directories(fingerprints
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/c/indigo
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/abseil-cpp
        )

target_link_directories(fingerprints PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(fingerprints indigo indigo-cpp OpenMP::OpenMP_CXX)
target_link_libraries(fingerprints absl::flags absl::flags_parse)
target_link_libraries(fingerprints qtrlib glog::glog)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(fingerprints)
