add_library(profiler STATIC
    CompleteSearchEngineProfiler.cpp
    SearchEngineProfiler.cpp)

target_include_directories(profiler
    PUBLIC .
    PUBLIC ../search_engine
    PUBLIC ../third_party/Indigo/api/cpp/src
)

#target_link_directories(search_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
#target_link_libraries(search_engine indigo-cpp)
