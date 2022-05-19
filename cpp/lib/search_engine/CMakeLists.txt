add_library(search_engine STATIC
        MockSearchEngine.cpp
        BingoSearchEngine.cpp)

target_include_directories(search_engine
    PUBLIC .
    PUBLIC ../../third_party/Indigo/api/cpp/src
    PUBLIC ../../third_party/Indigo/api/c/bingo-nosql
    PUBLIC ../../third_party/Indigo/api/c/indigo
    PUBLIC ../utils
)

target_link_libraries(search_engine utils)

