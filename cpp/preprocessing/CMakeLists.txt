set(CMAKE_CXX_STANDARD 20)
set(This preprocessing)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/AddRDKit.cmake)

file(GLOB_RECURSE sources *.cpp *.h)

add_executable(${This} main.cpp ${sources})

target_include_directories(${This}
        PUBLIC .
        PUBLIC qtr_preprocessing/preprocessors
        )

target_link_directories(${This} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(${This} qtr_preprocessing)
target_link_libraries(${This} ${LIBS} ${RDKit_LIBS})

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)

set_output_directory(${This})
