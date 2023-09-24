cmake_minimum_required(VERSION 3.13)

project(QtrFingerprint LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "9.4")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized") # Disable -Werror=maybe-uninitialized flag,
    # needed to compile gtest on new gcc compilers
endif ()

SET(BUILD_BINGO OFF CACHE BOOL "BUILD_BINGO Indigo option")
SET(BUILD_BINGO_ELASTIC OFF CACHE BOOL "BUILD_BINGO_ELASTIC Indigo option")
SET(BUILD_INDIGO_WRAPPERS OFF CACHE BOOL "BUILD_INDIGO_WRAPPERS Indigo option")
SET(BUILD_INDIGO_UTILS OFF CACHE BOOL "BUILD_INDIGO_UTILS Indigo option")
SET(WITH_STATIC ON CACHE BOOL "WITH_STATIC Indigo option")
SET(ABSL_PROPAGATE_CXX_STD ON)


add_subdirectory(run_db)
add_subdirectory(run_db/qtr_run)
add_subdirectory(build_db)
add_subdirectory(build_db/qtr_build)
add_subdirectory(preprocessing)
add_subdirectory(preprocessing/qtr_preprocessing)
add_subdirectory(playground)
add_subdirectory(qtrlib)
add_subdirectory(tests)

add_subdirectory(third_party/Indigo)
add_subdirectory(third_party/glog)
add_subdirectory(third_party/Crow)
add_subdirectory(third_party/abseil-cpp)
