# Use the official Ubuntu 22.04 as a base image
FROM ubuntu:22.04

# Set environment variables for non-interactive apt installs
ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages and dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    python3 \
    python3-pip \
    libfreetype6-dev \
    libfontconfig1-dev \
    libasio-dev \
    libgflags-dev \
    libtbb-dev \
    g++-9 \
    && apt-get clean

# Set GCC-9 as the default compiler
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 100 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 100

# Install Boost and numpy
RUN pip3 install numpy

# Clone the qtr-fingerprint repository
COPY ./ /qtr-fingerprint

# Move to qtr-fingerprint directory
WORKDIR /qtr-fingerprint

# Build RDKit
WORKDIR /qtr-fingerprint/cpp/third_party/rdkit

RUN apt-get install -y  \
    catch2 \
    libboost-all-dev

# Configure and build RDKit (adjust paths if necessary)
RUN mkdir build \
    && cd build  \
    && cmake -DPy_ENABLE_SHARED=1  \
    -DRDK_INSTALL_INTREE=ON  \
    -DRDK_INSTALL_STATIC_LIBS=OFF  \
    -DRDK_BUILD_CPP_TESTS=ON  \
    -DRDK_BUILD_INCHI_SUPPORT=ON  \
    -DRDKIT_RDINCHILIB_BUILD=ON  \
    -DBOOST_ROOT=/usr/include/boost .. \
    && make -j15

## Build qtr-fingerprint code
WORKDIR /qtr-fingerprint/cpp

RUN cmake -DCMAKE_BUILD_TYPE=Release -S ./ -B ./cmake-build-release \
    && cmake --build ./cmake-build-release --target preprocessing -j 15

# Executables will be located in qtr-fingerprint/cpp/cmake-build-release/bin
WORKDIR /qtr-fingerprint/cpp/cmake-build-release/bin

CMD ["/bin/bash"]
