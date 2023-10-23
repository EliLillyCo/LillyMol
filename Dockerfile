FROM gcc:13.2

RUN apt-get update && \ 
    apt-get upgrade -y

RUN apt-get install npm -y && \
    npm install -g @bazel/bazelisk

RUN apt-get install libblas-dev -y && \
    apt-get install liblapack-dev -y 

RUN apt-get install python3-minimal -y && \
    apt-get install python3-pandas python3-scipy python3-absl python3-pybind11 python3-protobuf -y

RUN rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

COPY . ./LillyMol

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux

# This step probably not necessary since now, third party
# dependencies are all linked static. Maybe that will change.
ENV LD_LIBRARY_PATH=/LillyMol/third_party/lib

WORKDIR /LillyMol/src

ENV BUILD_BDB=1
ENV BUILD_PYTHON=1

RUN ./build_third_party.sh
RUN ./build_from_src.sh

WORKDIR /LillyMol
