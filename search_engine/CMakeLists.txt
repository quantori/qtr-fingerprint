add_library(search_engine STATIC MockSearchEngine.cpp)

target_include_directories(search_engine
    PUBLIC .
    PUBLIC ../third_party/Indigo/api/cpp/src
)

#target_link_directories(search_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
#target_link_libraries(search_engine indigo-cpp)

