set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(This experiment)

#find_package(glog REQUIRED)
#find_package(absl CONFIG REQUIRED)
find_package(TBB REQUIRED)

file(GLOB_RECURSE EXPERIMENT_SOURCES *.cpp *.h)

add_executable(${This} ${EXPERIMENT_SOURCES})

target_link_libraries(${This} PRIVATE qtr-core)
target_link_libraries(${This} PRIVATE glog::glog)
target_link_libraries(${This} PUBLIC absl::flags absl::flags_parse)
target_link_libraries(${This} PRIVATE ${LIBS} ${RDKit_LIBS})
target_link_libraries(${This} PRIVATE indigo-core indigo-cpp bingo-nosql bingo-nosql-object)
target_link_libraries(${This} PRIVATE TBB::tbb)

target_include_directories(${This}
        PUBLIC .
        PUBLIC ../core
        PUBLIC ../core/dataset
        PUBLIC ../third_party/Indigo/api/c/bingo-nosql/src/
)

if(QTR_PROFILING)
    target_compile_definitions(${This} PRIVATE QTR_PROFILING=1)
endif()

target_compile_options(experiment PRIVATE -mavx2 -mfma -mbmi2)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(${This})