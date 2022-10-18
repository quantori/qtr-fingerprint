add_executable(app main.cpp)

target_include_directories(app
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/Indigo/api/c/indigo
        PUBLIC ../third_party/Indigo/api/c/bingo-nosql
        PUBLIC ../profiler
        PUBLIC ../third_party/abseil-cpp
        PUBLIC ../third_party/Crow
        )

target_link_directories(app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(app indigo-cpp qtrlib glog::glog)
target_link_libraries(app absl::flags absl::flags_parse)
target_link_libraries(app Crow::Crow)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(app)