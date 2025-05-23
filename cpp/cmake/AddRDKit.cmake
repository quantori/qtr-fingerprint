set(Boost_NO_BOOST_CMAKE ON)

cmake_minimum_required( VERSION 3.5 )

set(RDBASE "${CMAKE_SOURCE_DIR}/third_party/rdkit/")
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
link_directories ( ${RDBASE}/build/lib )

set(RDKit_LIBS RDKitChemReactions RDKitFileParsers RDKitSmilesParse RDKitDepictor
        RDKitRDGeometryLib RDKitRDGeneral RDKitSubstructMatch RDKitSubgraphs
        RDKitMolDraw2D RDKitGraphMol RDKitDistGeometry RDKitDistGeomHelpers
        RDKitMolAlign RDKitOptimizer RDKitForceField RDKitForceFieldHelpers
        RDKitAlignment RDKitForceField  RDKitMolTransforms RDKitEigenSolvers
        RDKitFingerprints RDKitDataStructs RDKitSubstructLibrary)

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)

set( LIBS ${RDKIT_LIBRARIES} Boost::iostreams ${RDKit_THREAD_LIBS}
        ${CAIRO_LIBRARIES} z  )

include_directories(SYSTEM ${RDKIT_INCLUDE_DIR})