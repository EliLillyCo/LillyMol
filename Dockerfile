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
    BUILD_PYTHON=1 \
    BUILD_GO=1

RUN apt-get install -y golang

RUN ./build_linux.sh

# Remove executables currently not being used.
RUN ./uninstall.sh

# final stage
FROM python:3.11.8-slim AS final

RUN apt-get update && \ 
    apt-get upgrade -y && \
    apt-get install build-essential libgomp1 ruby-dev protobuf-compiler -y && \
    rm -rf /var/lib/apt/lists/*

RUN gem install google-protobuf -v 3.21.12

# Note to Xuyan - this line might be a duplicate of line 11
RUN pip install pandas scipy absl-py pybind11 protobuf

COPY --from=build /LillyMol /LillyMol

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux

WORKDIR /LillyMol

RUN protoc -I=. --ruby_out=. test/lillymol_tests.proto

