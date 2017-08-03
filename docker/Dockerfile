FROM ubuntu:16.04
RUN apt-get -qq update
RUN apt-get install -y software-properties-common
RUN apt-get -qq update
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y
RUN apt-get -qq update
RUN apt-get install -y g++-6 clang libeigen3-dev libiomp-dev swig3.0 doxygen python3-dev python3-numpy libboost-dev libopenblas-dev libmetis-dev libmumps-seq-dev libann-dev libarpack2-dev
RUN apt-get install -y lcov curl texlive-font-utils libglu1-mesa-dev gmsh python3-numpy
WORKDIR /project
