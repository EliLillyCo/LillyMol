FROM gcc:10.2

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux

# Install Eigen
# (No need to install zlib since zlib1g and zlib1g-dev are alreday in the gcc base image)
RUN apt-get update && \ 
    apt-get upgrade -y && \
    apt-get install -y libeigen3-dev && \
    cd /usr/include && \
        ln -sf eigen3/Eigen Eigen

RUN apt-get install npm -y && \
    npm install -g @bazel/bazelisk && \
    apt-get install cmake -y

RUN apt-get install python3-pip -y && \
    pip3 install pandas && \
    apt-get install libblas-dev -y && \
    apt-get install liblapack-dev -y && \
    pip3 install scipy

COPY . ./LillyMol

WORKDIR /LillyMol

RUN mkdir bin && \
    mkdir bin/Linux && \
    cd src && \
    ./build_third_party.sh

WORKDIR /LillyMol/src

RUN cp /tmp/WORKSPACE .

WORKDIR /LillyMol/src/build_deps

RUN sed -i 's/\/workspaces\/LillyMolPrivate\/bin\/Linux\//\/LillyMol\/bin\/Linux\//g' install.bzl

WORKDIR /LillyMol

ENV LD_LIBRARY_PATH=/LillyMol/third_party/lib

RUN mkdir /node && \
    mkdir /node/scratch

WORKDIR /LillyMol/src

RUN ./build_from_src.sh

WORKDIR /LillyMol

RUN rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

# RUN apt-get install python3-pip -Y