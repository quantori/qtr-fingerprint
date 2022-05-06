add_executable(app main.cpp)

target_include_directories(app
    PUBLIC .
    PUBLIC ../third_party/Indigo/api/cpp/src
    PUBLIC ../profiler
    PUBLIC ../search_engine
)

target_link_directories(app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(app profiler search_engine indigo-cpp)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(app)