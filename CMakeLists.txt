cmake_minimum_required(VERSION 3.13)

project(QtrFingerprint LANGUAGES CXX)


SET(BUILD_BINGO OFF CACHE BOOL "BUILD_BINGO Indigo option")
SET(BUILD_BINGO_ELASTIC OFF CACHE BOOL "BUILD_BINGO_ELASTIC Indigo option")
SET(BUILD_INDIGO_WRAPPERS OFF CACHE BOOL "BUILD_INDIGO_WRAPPERS Indigo option")
SET(BUILD_INDIGO_UTILS OFF CACHE BOOL "BUILD_INDIGO_UTILS Indigo option")
SET(WITH_STATIC ON CACHE BOOL "WITH_STATIC Indigo option")

add_subdirectory(third_party/Indigo)

add_subdirectory(profiler)
add_subdirectory(search_engine)
add_subdirectory(app)
add_subdirectory(fingerprints)
add_subdirectory(engine-on-bingo-nosql)