enable_testing()

add_executable(tests 
    "main.cpp"
    "common.cpp"
    "BingoSearchEngineTests.cpp"
)

target_include_directories(tests
    PUBLIC ../third_party/Indigo/api/c/indigo-inchi
    PUBLIC ../third_party/Indigo/third_party/googletest/include
)

target_link_libraries(tests qtrlib indigo-cpp gtest  glog::glog)

include(GoogleTest)
gtest_discover_tests(tests)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(tests)