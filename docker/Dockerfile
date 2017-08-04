FROM ubuntu:16.04
RUN apt-get -qq update
RUN apt-get install -y software-properties-common
RUN apt-get -qq update
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y
RUN apt-get -qq update

# Toolchain?
RUN apt-get install -y cmake g++-6 clang 
# Nuto deps
RUN apt-get install -y libeigen3-dev libiomp-dev swig3.0 doxygen python3-dev python3-numpy libopenblas-dev libmetis-dev libmumps-seq-dev libann-dev libarpack2-dev gmsh 
# boost
RUN apt-get install -y libboost-filesystem-dev libboost-system-dev libboost-mpi-dev libboost-test-dev
# coverage and documentation
RUN apt-get install -y lcov curl texlive-font-utils

WORKDIR /project
