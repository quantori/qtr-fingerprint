set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(ServiceExe qtr-service)

# Source files
file(GLOB_RECURSE SERVICE_SOURCES ServiceManager.cpp main.cpp)

add_executable(${ServiceExe} ${SERVICE_SOURCES})

target_link_libraries(${ServiceExe} PRIVATE qtr-core)
target_link_libraries(${ServiceExe} PRIVATE Crow::Crow)
target_link_libraries(${ServiceExe} PRIVATE glog::glog)
target_link_libraries(${ServiceExe} PRIVATE absl::flags absl::flags_parse)
target_link_libraries(${ServiceExe} PRIVATE indigo-core indigo-cpp bingo-nosql bingo-nosql-object)

target_include_directories(${ServiceExe}
        PRIVATE .
        PRIVATE ../core
        PRIVATE ../third_party/Crow/include
        PRIVATE ../third_party/Indigo/api/c/bingo-nosql/
)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(${ServiceExe}) 