FROM bubbleprofiler/bubbleprofiler-base-ubuntu-gcc:latest

RUN apt-get update -y \
   && apt-get install -y \
   cmake \
   cmake-extras \
   libboost-all-dev \
   libeigen3-dev \
   libginac-dev \
   libgsl-dev \
   libnlopt-dev \
   && apt-get clean -y \
   && rm -rf /var/lib/apt/lists/*
