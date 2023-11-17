FROM mcr.microsoft.com/devcontainers/cpp:1.0.4-debian-12
# Install the xz-utils package
RUN apt-get update && apt-get install -y xz-utils

RUN apt-get install npm -y && \
    npm install -g @bazel/bazelisk && \
    apt-get install libblas-dev -y && \
    apt-get install liblapack-dev -y 

RUN apt-get install python3-minimal -y && \
    apt-get install python3-pandas python3-scipy python3-absl python3-pybind11 python3-protobuf awscli -y && \
    rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python
