FROM ubuntu:latest

RUN apt-get update -y \
   && apt-get install -y \
   gcc g++ gfortran binutils make \
   && apt-get clean -y \
   && rm -rf /var/lib/apt/lists/*
