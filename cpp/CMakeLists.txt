cmake_minimum_required(VERSION 3.13)

project(QtrFingerprint LANGUAGES CXX)

find_package (glog 0.6.0 REQUIRED)

SET(BUILD_BINGO OFF CACHE BOOL "BUILD_BINGO Indigo option")
SET(BUILD_BINGO_ELASTIC OFF CACHE BOOL "BUILD_BINGO_ELASTIC Indigo option")
SET(BUILD_INDIGO_WRAPPERS OFF CACHE BOOL "BUILD_INDIGO_WRAPPERS Indigo option")
SET(BUILD_INDIGO_UTILS OFF CACHE BOOL "BUILD_INDIGO_UTILS Indigo option")
SET(WITH_STATIC ON CACHE BOOL "WITH_STATIC Indigo option")
SET(ABSL_PROPAGATE_CXX_STD ON)



add_subdirectory(app)
add_subdirectory(fingerprints)
add_subdirectory(qtrlib)
add_subdirectory(third_party/Indigo)
add_subdirectory(third_party/abseil-cpp)
