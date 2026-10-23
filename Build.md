# Building LillyMol

## Table of Contents

* [Quick Start (Ubuntu)](#quick-start-ubuntu)
* [System Requirements](#system-requirements)

  * [macOS](#macos)
  * [Python](#python)
  * [Ruby](#ruby)
* [Detailed Instructions](#detailed-instructions)

  * [Configuring for Bazel](#configuring-for-bazel)
  * [C++ Dependencies](#c-dependencies)
* [Build](#build)
* [Troubleshooting](#troubleshooting)

  * [Bazel](#bazel)
  * [Protocol Buffers](#protocol-buffers)
  * [cmake](#cmake)


# Quick Start (Ubuntu)

```bash
sudo apt install git gcc-12 g++-12 wget unzip xz-utils zlib1g-dev npm
sudo npm install -g @bazel/bazelisk

git clone https://github.com/EliLillyCo/LillyMol
cd LillyMol

make
```
A successful build should result in executables appearing in `LillyMol/bin/Linux`.

A typical build may take 15 minutes, more if optional functionality is enabled.
The first build may take longer because Bazel downloads and compiles several
third-party dependencies.
A full build with third-party dependencies may require 5-10GB of free disk space.

> WARNING
> Bazel performs poorly or may fail entirely when its cache is located on NFS.
> If needed, configure `--output_user_root` to use local disk storage. See below

The following platforms are supported

| Platform | Status | Notes |
|---|---|---|
| Ubuntu 24+ | Fully supported | Primary development platform |
| RHEL/Rocky 8+ | Fully supported | Internal development platform |
| macOS | Experimental | First successful builds Sep 2025 |
| Windows | Unsupported | WSL may work |

In addition LillyMol can build in [Docker](/Dockerfile).

## System Requirements.

The primary build system used for LillyMol is [Bazel](https://bazel.build/).
You might also choose to use [Bazelisk](https://github.com/bazelbuild/bazelisk)
which makes keeping Bazel up to date easier. That is strongly encouraged.

Within a GitHub CodeSpace, this worked to install Bazelisk.
```
sudo apt install npm
sudo npm install -g @bazel/bazelisk
```
If you use the module system
```
module load bazelisk
```
If you are NOT building the python bindings, Bazel or Bazelisk is equivalent.

The software requires a gcc version of at least version 12. This version of LillyMol
uses some fairly recent c++ features, which require a recent compiler. The software
has been tested with gcc15.

If you use the module system
```
module load gcc12
module load bazelisk
module load git
```
Other system components that are needed

* wget
* unzip
* libz-dev
* xz

## macOS
LillyMol runs under [macOS](docs/macOS.md).

### Python
If you build LillyMol with 'make all' or by setting shell variable BUILD_PYTHON,
python bindings will be built. See [python](docs/python/Build.md).

### Ruby
Many scripts in the contrib directory are written in ruby. When protocol
buffers are involved there may be extra complications. See [ruby](docs/ruby/Build.md).

# Detailed Instructions

If you have Bazelisk and gcc installed, there is a reasonable possibility that
issuing `make` in the top level directory will work (but see note above
about NFS filesystems).

```
# Inside Lilly use the private repo
git clone https://github.com/EliLillyCo/LillyMol

cd /path/to/LillyMol
make
```
Executables will be in `bin/$(uname)` and libraries in `lib`. More details
below. There is no concept of installation prefix, everything remains
in the repo, although the 'install.sh' script will copy binaries to
${LILLYMOL_HOME}/bin/Linux.

*Note* by default neither Python bindings nor Berkeley DB dependencies
are built. If you wish to build either of those 
```
make python
make berkeleydb
```
or
```
make all
```

If you look at [Makefile](Makefile) you will see that all it is doing
is invoking src/build_linux.sh with different shell variables set.

### Configuring for Bazel
Within the src directory, the file `MODULE.bazel` configures the build environment
for `Bazel`. If you are building python bindings, see [BUILD](docs/python/Build.md),
[MODULE.bazel](src/MODULE.bazel) and 
[.bazelrc](src/.bazelrc) for information about configuring different versions of python.

### C++ Dependencies.
The preferred way of using third party software is via the Bazel
Module system. Most of the external dependencies needed are handled
via that mechanism. Today that includes

- **absl**: Google's c++ library - we use crc32c and some data structures.
- **benchmark**: Google's c++ benchmarking tool suite.
- **eigen**: matrix operations
- **googletest**: Google's c++ unit tests
- **protobuf**: Google's Protocol Buffers
- **highwayhash**: Google's fast fingerprint hash.
- **onetbb**: Threaded Building Blocks for multi-threading
- **pybind11**: c++ interface to python.
- **re2**: Google's regular expression library.
- **tomlplusplus**: a TOML library.
- **zlib**: compression
- **zstd**: compression - in development.

The complete listing is in the file [MODULE.bazel](src/MODULE.bazel).

There are several other dependencies which could be installed on the system,
which would considerably simplify the build configuration, but during
development we have frequently encountered machines that either lacked
needed software, or could
not be updated to the versions needed: lack of privileges, compatibility...
Several external dependencies are therefore downloaded and managed explicitly.

These third party dependencies are downloaded and built by the
[build_linux.sh](src/build_linux.sh) script, which will populate the `third_party`
directory (next to src) and then download, build and install the following dependencies

- **BerkeleyDb**: used for key/value database - if configured.
- **f2c/libf2c**: there is some fortran in LillyMol.
- **xgboost**: used for XGBoost models - if configured.
- **InChi**: if InChi bindings - if configured
- **nlopt**: an optimisation library - if configured.

Note that BerkeleyDB and Python bindings are only built if requested. 
In [Makefile](/Makefile) you will see use of the shell variables
'BUILD_PYTHON' and 'BUILD_BDB' which if set, enables building of
these optional features. These can be set any time.

Several utilities use 'omp' for parallel processing. That is only
supported as a system install.

Running `build_linux.sh` may be a lengthy process. 

# Build
Once third party dependencies are built, the install_linux.sh script will
begin compiling LillyMol with Bazel. The first phase usually involves building
the dependencies specified in MODULE.bazel.

Bazel needs to be able to store its cache on a local disk, *not* NFS. When building
inside Lilly, I have used `--output_user_root=/node/scratch/${USER}` to
use local scratch storage for Bazel's cache. Note that if there is a
recycling policy in place for the cache, you may see unexpected outcomes.
Purge the cache completely to start afresh. 

If outside Lilly, the 'build_linux.sh' script will check to
see if your HOME directory is on an NFS mounted file system, and if so, will
specify /tmp for Bazel's cache. This is almost certainly not what you want,
so edit 'build_linux.sh' to specify a local directory for
`--output_user_root`. Again, only needed if you are on an NFS file system.
You can also enter this value in Bazel's configuration file `.bazelrc`.

By default, Bazel will use all cores available on the local machine.
If needed, limit the number of cores with the `--jobs` option inside
'build_linux.sh'

Optionally set shell variables BUILD_BDB, BUILD_XGBOOST and BUILD_PYTHON to enable
building of optional features, or just use `make all`.

Once the Bazel preconditions are set, do the build, test and installs
```
cd src                 # you might already be here
./build_linux.sh       # takes a while
```

The script will

1. run the C++ unit tests,
2. build all executables
3. build the python bindings
4. install executables into the `/path/to/LillyMol/bin/$(uname)` directory
5. copy python related shared libraries to /path/to/LillyMol/lib (if BUILD_PYTHON)
6. run python unit tests (if BUILD_PYTHON)

Step 5 is done via the [copy_shared_libraries.sh](/src/copy_shared_libraries.sh)
script. It also copies some python compiled protos. Adjust as needed.

For anyone interested in doing their own development, a typical build
inside Lilly might be (change the path for test_env)

```
bazelisk --output_user_root=/node/scratch/${USER}
         build
         --jobs=8
         -c opt
         --cxxopt=-DGIT_HASH=\"$(git rev-parse --short --verify HEAD)\"
         --cxxopt=-DTODAY=\"$(date +%Y-%b-%d)\"
         Molecule_Tools:all     <-  or some other target
```

Most will want to put this in a small shell script, and/or add to .bazelrc where
possible.

When building for release, it is convenient to include the git hash and
the date of the build in the executables. That is not necessary, omit those if not needed.
Note that because the date is included with cxxopt, this _will_ cause a daily
recompile. While this is hardly desirable, the benefits are many.

## Troubleshooting

### Bazel
Bazel is updated regularly. LillyMol attempts to follow the Google philosophy
of "live at head", staying as close as possible to new versions. Bazel
introduces breaking changes which can cause problems for LillyMol. You can
freeze the Bazel version by adding a file '.bazelversion' to the 'LillyMol/src'
directory with the version of Bazel to use.

### Protocol Buffers
While protos are now foundational to LillyMol, they can introduce build
difficulties.

If only C++ tools are built, Bazel usually takes care of compilation of
protocol buffers. It does this by downloading and building the protobuf version specified
in [MODULE.bazel](src/MODULE.bazel). This means that if there is an older
/usr/bin/protoc on the system it will be ignored and not interfere.

If building the python bindings there is the possibility of incompatibilities
between the python version used by Bazel and the installed python protobuf package.
This will show up at run-time.
Experiment with changing the 'protobuf' version specified in [MODULE.bazel](src/MODULE.bazel).
See also [protobuf](https://protobuf.dev/support/version-support/) for more information.

### cmake
The distribution contains `cmake` infrastructure, that is currently
not functional.  We have not been able to make it work,
usually as a result of conflicting protocol buffer versions on the
system (see previous).  Work is ongoing to get cmake working for the public release.

The 'cmake' branch also contains a non-working cmake attempt, this time
guided by ChatGPT. The main problem again seems to be systems with obsolete
protoc system installs and the difficulty of having cmake use an alternate.
