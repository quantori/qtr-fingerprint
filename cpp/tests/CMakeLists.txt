enable_testing()

add_executable(tests 
    "main.cpp"
    "common.cpp"
    "BingoSearchEngineTests.cpp"
)

target_include_directories(tests
    PUBLIC ../lib/search_engine 
    PUBLIC ../lib/utils 
    PUBLIC ../third_party/Indigo/api/c/indigo-inchi
    PUBLIC ../third_party/Indigo/third_party/googletest/include
)

target_link_libraries(tests search_engine utils indigo-cpp gtest)

include(GoogleTest)
gtest_discover_tests(tests)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(tests)