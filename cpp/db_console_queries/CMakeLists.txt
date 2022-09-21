set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

add_executable(db_console_queries main.cpp)

target_include_directories(db_console_queries
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/c/indigo
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/abseil-cpp
        )

target_link_directories(db_console_queries PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(db_console_queries indigo indigo-cpp OpenMP::OpenMP_CXX)
target_link_libraries(db_console_queries absl::flags absl::flags_parse)
target_link_libraries(db_console_queries qtrlib glog::glog)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(db_console_queries)
