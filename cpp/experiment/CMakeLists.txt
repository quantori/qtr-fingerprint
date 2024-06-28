set(CMAKE_CXX_STANDARD 20)

set(This experiment)

find_package(glog REQUIRED)
find_package(absl CONFIG REQUIRED)

file(GLOB_RECURSE sources *.cpp *.h)

add_executable(${This} main.cpp ${sources})

target_link_libraries(${This} PRIVATE glog::glog)
target_link_libraries(${This} PRIVATE absl::flags absl::flags_parse)
target_link_libraries(${This} PRIVATE ${LIBS} ${RDKit_LIBS})

target_include_directories(${This}
        PUBLIC fingerprints
        PUBLIC search_engines
        PUBLIC parsers
)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(${This})