FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

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
    catch2 \
    libboost-all-dev \
    && apt-get clean

RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 100 \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 100

COPY ./ /qtr-fingerprint
WORKDIR /qtr-fingerprint

RUN cd cpp/third_party/rdkit  \
    && mkdir build && cd build \
    && cmake \
      -DPy_ENABLE_SHARED=1 \
      -DRDK_INSTALL_INTREE=ON \
      -DRDK_INSTALL_STATIC_LIBS=OFF \
      -DRDK_BUILD_CPP_TESTS=ON  \
      -DRDK_BUILD_INCHI_SUPPORT=ON  \
      -DRDKIT_RDINCHILIB_BUILD=ON  \
      .. \
    && make -j

WORKDIR /qtr-fingerprint

RUN cd cpp \
    && cmake -DCMAKE_BUILD_TYPE=Release -S ./ -B ./build \
    && cmake --build ./cmake-build-release --target preprocessing -j

ENV PATH="/qtr-fingerprint/cpp/cmake-build-release/bin:${PATH}"

CMD ["/bin/bash"]
