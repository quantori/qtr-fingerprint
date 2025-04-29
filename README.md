# qtr-fingerprint

A molecular substructure search engine that provides fast searching capability for chemical compounds.

## Project Overview

This project is a proof of concept implementation demonstrating the application of BallTree data structures to the problem of chemical fingerprint indexing, which enables efficient molecular substructure searching. The codebase focuses on benchmarking the performance of this approach against traditional search methods.

## Project Structure

The project is organized as follows:

* `cpp/experiment` - Main C++ codebase containing the experiment framework
  * `search/engines` - Search engine implementations
  * `search/algorithms` - Search algorithms
  * `frameworks` - Molecular frameworks adapters (RDKit, Indigo)
  * `dataset` - Dataset handling
  * `io` - Input/output utilities
  * `benchmarking` - Benchmarking infrastructure
  * `stats` - Statistics collection
  * `utils` - Utility functions

## Search Engines

The project supports multiple search engine types:

* `BallTreeRDKit` - BallTree search engine using RDKit framework
* `BallTreeIndigo` - BallTree search engine using Indigo framework
* `RDKit` - Direct RDKit search
* `Indigo` - Direct Indigo search

## Command Line Arguments

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

The original benchmarks compared `QtrRam` and `BingoNoSQL` implementations. With the current code structure, the equivalent comparisons would be between `BallTreeRDKit`/`BallTreeIndigo` and `RDKit`/`Indigo` engines.

Benchmarks can be run with parameters similar to:

* `--SearchEngineType=BallTreeRDKit` (or other supported engine type)
* `--MaxResults=10000`  
* `--TimeLimit=60`

## Logging Configuration

Configure the logger using environment variables:

1. `GLOG_log_dir=` - Directory for log files
2. `GLOG_alsologtostderr=true/false` - Whether to also log to stderr
3. `GLOG_logtostdout=true/false` - Whether to log to stdout
4. `GLOG_minloglevel=0` - Minimum log level (0=INFO, 1=WARNING, 2=ERROR, 3=FATAL)
