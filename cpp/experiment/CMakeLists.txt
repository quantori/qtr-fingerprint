set(CMAKE_CXX_STANDARD 20)

set(This experiment)
find_package(glog REQUIRED)

file(GLOB_RECURSE sources *.cpp *.h)

include(../cmake/AddRDKit.cmake)

add_executable(${This} main.cpp ${sources})

target_link_libraries(${This} PRIVATE glog::glog)

target_include_directories(${This}
        PUBLIC Fingerprints
        PUBLIC SearchEngines
)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(${This})