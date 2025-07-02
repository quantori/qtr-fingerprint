set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(This qtr-core)

find_package(TBB REQUIRED)

file(GLOB_RECURSE CORE_SOURCES *.cpp *.h)

add_library(${This} STATIC ${CORE_SOURCES})

target_link_libraries(${This} PRIVATE glog::glog)
target_link_libraries(${This} PUBLIC absl::flags absl::flags_parse)
target_link_libraries(${This} PRIVATE ${LIBS} ${RDKit_LIBS})
target_link_libraries(${This} PRIVATE indigo-core indigo-cpp bingo-nosql bingo-nosql-object)
target_link_libraries(${This} PRIVATE TBB::tbb)

target_include_directories(${This}
        PUBLIC .
        PUBLIC frameworks
        PUBLIC search
        PUBLIC search/engines
        PUBLIC search/algorithms
        PUBLIC search/utils
        PUBLIC io
        PUBLIC utils
        PUBLIC stats
        PUBLIC benchmarking
        PUBLIC dataset
        PUBLIC ../third_party/Indigo/api/c/bingo-nosql/src/
)

if(QTR_PROFILING)
    target_compile_definitions(${This} PRIVATE QTR_PROFILING=1)
endif()

target_compile_options(${This} PRIVATE -mavx2 -mfma -mbmi2)
