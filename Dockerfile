# builder stage
FROM python:3.11.8 AS build

RUN apt-get update && \ 
    apt-get upgrade -y

RUN apt-get install npm -y && \
    npm install -g @bazel/bazelisk && \
    apt-get install libblas-dev liblapack-dev libzmq3-dev -y

RUN pip install pandas scipy absl-py pybind11 protobuf

COPY . ./LillyMol

WORKDIR /LillyMol/src

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux \
    BUILD_BDB=1 \
    BUILD_PYTHON=1

RUN ./build_linux.sh

# final stage
FROM python:3.11.8-slim AS final

RUN apt-get update && \ 
    apt-get upgrade -y && \
    apt-get install libgomp1 -y

RUN pip install pandas scipy absl-py pybind11 protobuf

COPY --from=build /LillyMol /LillyMol

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux

WORKDIR /LillyMol
