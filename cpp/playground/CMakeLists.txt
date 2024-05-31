set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

#----------- rdkit example part
set(Boost_NO_BOOST_CMAKE ON)

cmake_minimum_required( VERSION 3.5 )

set(RDBASE $ENV{RDBASE})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${RDBASE}/Code/cmake/Modules")

# note that if you haven't installed/built the toolkit with CoordGen, you'll
# have problems with this.
add_definitions("-DRDK_BUILD_COORDGEN_SUPPORT=ON")
add_compile_options(-Wall)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package( Boost COMPONENTS iostreams filesystem system)
find_package( Cairo REQUIRED )

# specify where CMake can find the RDKit libraries
include_directories ( ${RDBASE}/Code ${CAIRO_INCLUDE_DIRS} )
link_directories ( ${RDBASE}/lib )

set(RDKit_LIBS RDKitChemReactions RDKitFileParsers RDKitSmilesParse RDKitDepictor
	       RDKitRDGeometryLib RDKitRDGeneral RDKitSubstructMatch RDKitSubgraphs
	       RDKitMolDraw2D RDKitGraphMol RDKitDistGeometry RDKitDistGeomHelpers
  	       RDKitMolAlign RDKitOptimizer RDKitForceField RDKitForceFieldHelpers
	       RDKitAlignment RDKitForceField  RDKitMolTransforms RDKitEigenSolvers )

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)

set( LIBS ${RDKIT_LIBRARIES} Boost::iostreams ${RDKit_THREAD_LIBS}
     ${CAIRO_LIBRARIES} z  )

include_directories(${RDKIT_INCLUDE_DIR})
# end of rdkir exmaple

add_executable(playground main.cpp)
target_include_directories(playground
        PUBLIC .
        )

target_link_libraries(playground ${LIBS} ${RDKit_LIBS})

target_link_directories(playground PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)
target_link_libraries(playground OpenMP::OpenMP_CXX)
target_link_libraries(playground absl::flags absl::flags_parse)
target_link_libraries(playground qtrlib glog::glog)
target_link_libraries(playground bingo-nosql-object)

#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../third_party/rdkit/Docs/Book/C++Examples ${CMAKE_CURRENT_BINARY_DIR}/../third_party/rdkit/Docs/Book/C++Examples)

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SetOutputDirectory.cmake)


set_output_directory(playground)
