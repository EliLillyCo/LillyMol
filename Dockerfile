# builder stage
FROM gcc:13.2 AS build

RUN apt-get update && \ 
    apt-get upgrade -y

RUN apt-get install npm -y && \
    npm install -g @bazel/bazelisk && \
    apt-get install libblas-dev -y && \
    apt-get install liblapack-dev -y 

RUN apt-get install python3-minimal -y && \
    apt-get install python3-pandas python3-scipy python3-absl python3-pybind11 python3-protobuf -y && \
    rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

COPY . ./LillyMol

WORKDIR /LillyMol/src

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux \
    BUILD_BDB=1 \
    BUILD_PYTHON=1

RUN ./update_bazel_configs.sh && ./build_third_party.sh && ./build_from_src.sh

# final stage
FROM ubuntu:mantic AS final

RUN apt-get update && \ 
    apt-get upgrade -y && \
    apt-get install python3-minimal -y && \
    apt-get install python3-pandas python3-scipy python3-absl python3-pybind11 python3-protobuf -y && \
    rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

COPY --from=build /LillyMol /LillyMol

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux

WORKDIR /LillyMol
