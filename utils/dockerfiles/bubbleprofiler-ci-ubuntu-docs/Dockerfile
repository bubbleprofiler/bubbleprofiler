FROM bubbleprofiler/bubbleprofiler-ci-ubuntu-gcc:latest

RUN apt-get update -y \
   && apt-get install -y \
   doxygen \
   && apt-get clean -y \
   && rm -rf /var/lib/apt/lists/*
