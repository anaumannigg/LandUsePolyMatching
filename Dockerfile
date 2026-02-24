# ---------- BUILD STAGE ----------
FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libgdal-dev \
    libmpfr-dev \
    libgmp-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /tmp

# ----------------------------
# Install Boost 1.83.0
# ----------------------------
RUN wget https://archives.boost.io/release/1.83.0/source/boost_1_83_0.tar.gz \
    && tar -xzf boost_1_83_0.tar.gz \
    && cd boost_1_83_0 \
    && ./bootstrap.sh \
    && ./b2 install --prefix=/usr/local -j$(nproc)

# ----------------------------
# Install CGAL 6.0.1
# ----------------------------
RUN wget https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1.tar.xz \
    && tar -xf CGAL-6.0.1.tar.xz \
    && mkdir CGAL-6.0.1/build \
    && cd CGAL-6.0.1/build \
    && cmake .. -DCMAKE_BUILD_TYPE=Release \
    && cmake --build . -j$(nproc) \
    && cmake --install .

# ----------------------------
# Clone HiGHS into external/
# ----------------------------
RUN mkdir -p /app/external \
    && git clone --depth=1 https://github.com/ERGO-Code/HiGHS.git /app/external/HiGHS


WORKDIR /app

# Copy necessary files into container
COPY CMakeLists.txt .
COPY include/ include/
COPY src/ src/

# Force HiGHS backend
RUN cmake -B build -S . -DFORCE_HIGHS=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/usr/local -DBUILD_SHARED_LIBS=OFF
RUN cmake --build build -j$(nproc)

# ---------- RUNTIME STAGE ----------
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    libgdal-dev \
    libmpfr-dev \
    libgmp-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app/build
COPY --from=builder /app/build/TCPolygonMatching /app/build/TCPolygonMatching

# Copy exact Boost + CGAL installations
COPY --from=builder /usr/local /usr/local

# Create expected directories
RUN mkdir -p /app/input /app/output

ENTRYPOINT ["/app/build/TCPolygonMatching"]