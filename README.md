# qtr-fingerprint

A molecular substructure search engine that provides fast searching capability for chemical compounds. The project includes both a benchmarking framework and a web service API for molecular searches.

## Quick Start with Docker 

```bash
# Build and run the service
docker-compose up -d
```

The service will be available at `http://localhost:8080`.

## API Usage

The service now supports **collection-based** molecular searches:

1. **Create Collection**: `POST /collections`
2. **Upload CSV**: `POST /collections/{id}/upload` 
3. **Build Index**: `POST /collections/{id}/build`
4. **Search**: `POST /fastquery` with `collectionId`

## Project Overview

This project is a proof of concept implementation demonstrating the application of BallTree data structures to the problem of chemical fingerprint indexing, which enables efficient molecular substructure searching. The codebase focuses on benchmarking the performance of this approach against traditional search methods and provides a REST API service for real-time molecular searches.

## Project Structure

The project is organized as follows:

* `cpp/core` - Core C++ library containing search algorithms and frameworks
  * `search/engines` - Search engine implementations
  * `search/algorithms` - Search algorithms (BallTree, etc.)
  * `frameworks` - Molecular frameworks adapters (RDKit, Indigo)
  * `dataset` - Dataset handling and storage
  * `io` - Input/output utilities and parsers
  * `benchmarking` - Benchmarking infrastructure
  * `stats` - Statistics collection
  * `utils` - Utility functions

* `cpp/experiment` - Benchmarking application
* `cpp/service` - Web service implementation
* `cpp/build` - Build directory for compiled binaries

## Search Engines

The project supports multiple search engine types:

* `BallTreeRDKit` - BallTree search engine using RDKit framework
* `BallTreeIndigo` - BallTree search engine using Indigo framework  
* `RDKit` - Direct RDKit SubstructLibrary search
* `Indigo` - Direct Indigo Bingo NoSQL search

## Experiment

The experiment application accepts the following command line arguments:

* `--SearchEngineType` - Type of search engine to be tested (BallTreeRDKit, BallTreeIndigo, RDKit, Indigo)
* `--MaxResults` - Maximum number of results to retrieve for each query
* `--TimeLimit` - Time limit in seconds for each query
* `--QueriesFile` - File containing the queries to be tested in the experiment
* `--DatasetDir` - Directory containing the dataset (CSV files)
* `--QueriesStatisticFile` - File where query statistics will be written
* `--SearchEngineStatisticFile` - File where search engine statistics will be written

## Requirements

* CMake 3.13 or higher
* C++20 compatible compiler (GCC 9.4 or higher recommended)
* Required libraries:
  * libfreetype6-dev
  * libfontconfig1-dev
  * libasio-dev
  * libgflags-dev
  * libtbb-dev
  * Boost libraries

Install the required libraries on Ubuntu with:
```bash
apt-get install libfreetype6-dev libfontconfig1-dev libasio-dev libgflags-dev libtbb-dev libboost-all-dev
```

## Build Instructions

1. Clone the repository:
   ```bash
   git clone https://github.com/quantori/qtr-fingerprint.git
   cd qtr-fingerprint
   git submodule update --init --recursive
   ```

2. Build RDKit (see [RDKIT_BUILD.md](RDKIT_BUILD.md) for detailed instructions):
   ```bash
   cd cpp/third_party/rdkit
   mkdir build && cd build
   cmake -DPy_ENABLE_SHARED=1 -DRDK_INSTALL_INTREE=ON -DRDK_INSTALL_STATIC_LIBS=OFF -DRDK_BUILD_CPP_TESTS=ON -DRDK_BUILD_INCHI_SUPPORT=ON -DRDKIT_RDINCHILIB_BUILD=ON ..
   make -j
   ```

   Note: 
   - You may need to specify the numpy location: `-DPYTHON_NUMPY_INCLUDE_PATH="$(python -c 'import numpy ; print(numpy.get_include())')"`
   - You may need to specify the boost location: `-DBOOST_ROOT="/path/to/boost"`

3. Build the qtr-fingerprint code:
   ```bash
   cd ../../../  # Return to cpp directory
   cmake -DCMAKE_BUILD_TYPE=Release -S ./ -B ./cmake-build-release
   cmake --build ./cmake-build-release --target experiment -j
   ```

4. The compiled executable will be located in `cpp/cmake-build-release/bin/`

## Docker Build

You can also use the provided Dockerfile to build the project:

```bash
./build_docker.sh
```

This will create a Docker image with all dependencies installed and the project built.

## Example Usage

Run the experiment with:

```bash
./cpp/cmake-build-release/bin/experiment --SearchEngineType=BallTreeRDKit --MaxResults=100 --TimeLimit=60 --QueriesFile=path/to/queries.txt --DatasetDir=path/to/dataset --QueriesStatisticFile=queries_stats.csv --SearchEngineStatisticFile=engine_stats.csv
```

## Benchmarking and Research

The results described in [this article](https://arxiv.org/abs/2310.02022) were obtained using [this dataset](https://www.dropbox.com/scl/fi/5je47quy4naarr1svflgt/compound_libraries.tar.gz?rlkey=6n05eexhfatcjpvgyvehya8ks&dl=0). The set of queries can be found in [this file](data/queries_3488_good.txt).
