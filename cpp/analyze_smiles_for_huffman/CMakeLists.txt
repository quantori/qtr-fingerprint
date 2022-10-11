add_executable(analyze_smiles_for_huffman main.cpp)

target_include_directories(analyze_smiles_for_huffman
        PUBLIC .
        PUBLIC ../third_party/Indigo/api/cpp/src
        PUBLIC ../third_party/Indigo/api/c/indigo
        PUBLIC ../third_party/Indigo/api/c/bingo-nosql
        PUBLIC ../profiler
        PUBLIC ../third_party/abseil-cpp
        )

target_link_directories(analyze_smiles_for_huffman PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(analyze_smiles_for_huffman indigo-cpp qtrlib glog::glog)
target_link_libraries(analyze_smiles_for_huffman absl::flags absl::flags_parse)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)
set_output_directory(analyze_smiles_for_huffman)