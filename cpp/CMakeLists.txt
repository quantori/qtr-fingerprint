cmake_minimum_required(VERSION 3.13)

project(QtrFingerprint LANGUAGES CXX)

## Sanitizers and debug-friendly flags
#set(SANITIZER_FLAGS
#        -fsanitize=address
#        -fsanitize=undefined
#        -fsanitize=leak
#        -fsanitize=null
#        -fsanitize=bounds
#        -fsanitize=vptr
#        -fno-sanitize=alignment
#        -fno-sanitize-recover
#        -fno-omit-frame-pointer
#        -fno-common
#        -fstack-protector
#)
#
#message("Added Sanitizer flags for debug")
#add_compile_options(${SANITIZER_FLAGS})
#add_link_options(${SANITIZER_FLAGS})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -march=native")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "9.4")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized") # Disable -Werror=maybe-uninitialized flag,
    # needed to compile gtest on new gcc compilers
endif ()

set(QTR_PROFILING OFF CACHE BOOL "Enable profiling")

SET(BUILD_BINGO OFF CACHE BOOL "BUILD_BINGO Indigo option")
SET(BUILD_BINGO_ELASTIC OFF CACHE BOOL "BUILD_BINGO_ELASTIC Indigo option")
SET(BUILD_INDIGO_WRAPPERS OFF CACHE BOOL "BUILD_INDIGO_WRAPPERS Indigo option")
SET(BUILD_INDIGO_UTILS OFF CACHE BOOL "BUILD_INDIGO_UTILS Indigo option")
SET(WITH_STATIC ON CACHE BOOL "WITH_STATIC Indigo option")
SET(ABSL_PROPAGATE_CXX_STD ON)

include(cmake/AddRDKit.cmake)
SET(WITH_GTEST OFF)

add_subdirectory(third_party/Indigo)
add_subdirectory(third_party/glog)
add_subdirectory(third_party/Crow)
add_subdirectory(third_party/abseil-cpp)

add_subdirectory(core)
add_subdirectory(experiment)
add_subdirectory(service)
