FROM nvidia/cuda:13.0.1-devel-ubuntu24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        ninja-build \
        libtinyxml2-dev \
        zlib1g-dev \
        libomp-dev \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY CMakeLists.txt ./
COPY include ./include
COPY src ./src
COPY tests ./tests

ARG CUDA_ARCH=86
RUN cmake -S . -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_TESTING=ON \
        -DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCH} \
    && cmake --build build -j \
    && ctest --test-dir build --output-on-failure

ENV PATH=/src/build:${PATH}
WORKDIR /work
ENTRYPOINT ["cuda-demux"]
