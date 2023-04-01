set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

file(GLOB RdKitSrc "../third_party/rdkit/Code/*.cpp")
add_executable(playground main.cpp
        ../third_party/Indigo/core/indigo-core/common/gzip/gzip_output.cpp
        ${RdKitSrc})

target_include_directories(playground
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/c/indigo
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/Indigo/core/indigo-core
        PUBLIC ../third_party/Indigo/core/indigo-core/common
        PUBLIC ../third_party/Indigo/third_party/object_threadsafe
        PUBLIC ../third_party/abseil-cpp
        PUBLIC ../third_party/rdkit/Code
        )


target_link_directories(playground PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(playground indigo indigo-cpp OpenMP::OpenMP_CXX)
target_link_libraries(playground absl::flags absl::flags_parse)
target_link_libraries(playground qtrlib glog::glog)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(playground)
