add_library(profiler STATIC
    CompleteSearchEngineProfiler.cpp
    SearchEngineProfiler.cpp)

target_include_directories(profiler
    PUBLIC .
    PUBLIC ../../third_party/Indigo/api/cpp/src
    PUBLIC ../search_engine
)
