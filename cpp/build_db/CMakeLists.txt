set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE sources *.cpp *.h)

include(../cmake/AddRDKit.cmake)

add_executable(build_db ${sources})

target_include_directories(build_db
        PUBLIC .
        )

target_link_directories(build_db PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)

target_link_libraries(build_db qtr_build)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(build_db)
