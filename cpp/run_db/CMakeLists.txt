set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

file(GLOB_RECURSE sources *.cpp *.h)

add_executable(run_db main.cpp ${sources})

target_include_directories(run_db
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/c/indigo
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/abseil-cpp
        PUBLIC ../third_party/Crow
        )

target_link_directories(run_db PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(run_db indigo indigo-cpp bingo-nosql indigo-core OpenMP::OpenMP_CXX)
target_link_libraries(run_db absl::flags absl::flags_parse)
target_link_libraries(run_db qtrlib glog::glog)
target_link_libraries(run_db Crow::Crow)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(run_db)
