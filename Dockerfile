# builder stage
ARG PYTHON_VERSION='3.12'
FROM python:${PYTHON_VERSION}-trixie AS build

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
        cmake \
        golang \
        libblas-dev \
        liblapack-dev \
        npm \
        protobuf-compiler \
        xz-utils && \
    rm -rf /var/lib/apt/lists/*

RUN npm install -g @bazel/bazelisk && \
    npm cache clean --force

RUN pip install pandas scipy absl-py nanobind protobuf

COPY . ./LillyMol

WORKDIR /LillyMol/src

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux \
    BUILD_BDB=1 \
    BUILD_PYTHON=1 \
    BUILD_GO=1 \
    PYTHON_VERSION=${PYTHON_VERSION}

RUN ./build_linux.sh

# Remove executables currently not being used.
RUN ./uninstall.sh

# final stage
ARG PYTHON_VERSION='3.12'
FROM python:${PYTHON_VERSION}-slim-trixie AS final

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install build-essential libgomp1 ruby-dev protobuf-compiler -y && \
    rm -rf /var/lib/apt/lists/*

RUN gem install google-protobuf -v 3.21.12

RUN pip install pandas scipy absl-py nanobind protobuf

COPY --from=build /LillyMol /LillyMol

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux

WORKDIR /LillyMol
