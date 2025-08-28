# Building LillyMol

## REQUIREMENTS:

**Note**: This is significantly different from prior versions.

LillyMol is primarily developed on RedHat and Ubuntu systems. It has
been built on a variety of Cloud environments. It does not yet build on
a Mac - that is a work in progress.

The primary build system used for LillyMol is [bazel](https://bazel.build/).
You might also choose to use [bazelisk](https://github.com/bazelbuild/bazelisk)
which makes keeping bazel up to date easier. That is strongly encouraged.

Within a GitHub CodeSpace, this worked to install bazelisk.
```
sudo apt install npm
sudo npm install -g @bazel/bazelisk
```
If you use the module system
```
module load bazelisk
```
If you are NOT building the python bindings, bazel or bazelisk is equivalent.

The software requires a gcc version of at least version 10. This version of LillyMol
uses some fairly recent c++ features, which require a recent compiler. The software
has been tested with gcc13.

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

### Python
If you wish to build the python bindings, you will need a recent version of
python. Development was done with python3.11 and has not been
tested on any other version, although we have no reason to believe
it will not work with other versions. You will need to install
```
pip install pybind11 absl-py protobuf
apt install python-dev
```

If you wish to use the xgboost QSAR model building tools in LillyMol,
also pip install xgboost, scikit-learn, matplotlib and pandas.

Make sure that python-dev and libblas-dev are installed at system level.

```
sudo apt install python-dev libblas-dev
```

Things seem to work seamlessly in virtualenv.

Note that with the default build (below) Python bindings are not built,
but 'make all' will.

# TLDR
If you have bazelisk and gcc installed, there is a reasonable possibility that
issuing `make` in the top level directory will work (but see note below
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
LILLYMOL_HOME/bin/Linux.

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
is sequentially invoking the three scripts discussed below, with
different shell variables set.

### Configuring for bazel
Within the src directory, the file `MODULE.bazel` configures the build environment
for `bazel`. If you are building python bindings, this file needs to be updated
to reflect the location of your local python. The script `update_python_in_module_bazel.py`
does this automatically from the Makefile.

### C++ Dependencies.
There are several dependencies which could be installed on the system,
which would considerably simplify the build configuration, but during
development we have frequently found ourselves on machines that could
not be updated to the versions we needed, or where we lacked privileges, or...
So external dependencies are downloaded and managed explicitly.
The preferred way of using third party software is via the Bazel
Module system. Most of the external dependencies needed are handled
via that mechanism. Today that includes

- **absl**: Google's c++ library - we use crc32c and some data structures.
- **benchmark**: Google's c++ benchmarking tool suite.
- **eigen**: matrix operations
- **googletest**: Google's c++ unit tests
- **protobuf**: Google's Protocol Buffers
- **highwayhash**: Google's fast fingerprint hash.
- **re2**: Google's regular expression library.
- **tomlplusplus**: a TOML library.
- **tbb**: Threaded Building Blocks for multi-threading
- **zlib**: compression

The complete listing is in the file [MODULE.bazel](src/MODULE.bazel).

Other third party dependencies are downloaded and built by the
[build_linus.sh](src/build_linus.sh) script, which will create a `third_party`
directory (next to src) and then download, build and install the following dependencies

- **BerkeleyDb**: used for key/value databases
- **f2c/libf2c**: there is some fortran in LillyMol.
- **xgboost**: used for XGBoost models.
- **inchii**: if InChi bindings are needed.
- **nlopt**: an optimisation library.

Note that BerkeleyDB and Python bindings are only built if requested. 
In [Makefile](/Makefile) you will see use of the shell variables
'BUILD_PYTHON' and 'BUILD_BDB' which if set, enables building of
these optional features. These can be set any time.

Running `build_linux.sh` may be a lengthy process. 

Note too that installing these external dependencies and running bazel may require
considerable amounts of disk space. For example at the time of writing my
'third_party' directory contains 1.2GB and my bazel temporary area contains 4.2GB.

### Python Bindings
During building of external dependencies (with build_third_party.sh
and if BUILD_PYTHON is set)
the script `update_python_in_module_bazel.py` will examine your python
installation and get information about the include path. With that
info it will update MODULE.bazel with new values for the 'path'
attributes of the python related features.

Note that if it does not find a pybind11 installation, the build
will continue, but the python related parts of the build will subsequently fail.

You can of course manually update MODULE.bazel to point to your
python installation. See the 'new_local_repository' sections for 'python'
and 'pybind11'

# Build
Once the third party dependencies have been built, and MODULE.bazel and
install.bzl configured, LillyMol building can begin. 

Bazel needs to be able to store its cache on a local disk, *not* NFS. When building
inside Lilly, I have used `--output_user_root=/node/scratch/${USER}` to
use local scratch storage for bazel's cache. Note that if there is a
recycling policy in place for the cache, you may see unexpected outcomes.
Purge the cache completely to start afresh. 

If outside Lilly, the 'build_linux.sh' script will check to
see if your HOME directory is on an NFS mounted file system, and if so, will
specify /tmp for bazel's cache. This is almost certainly not what you want,
so edit 'build_linux.sh' to specify a local directory for
`--output_user_root`. Again, only needed if you are on an NFS file system.
You can also enter this value in bazel's configuration file `.bazelrc`.

By default, bazel will use all cores available on the local machine.
If needed, limit the number of cores with the `--jobs` option inside
'build_linux.sh' (sorry no command line options here).

Optionally set shell variables BUILD_BDB and BUILD_PYTHON to enable
building of optional features.

Once the bazel preconditions are set, do the build, test and installs
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

For anyone interesting in doing their own development, a typical build
inside Lilly might be (change the path for test_env)

```
bazelisk --output_user_root=/node/scratch/${USER}
         build
         --jobs=8
         -c opt
         --cxxopt=-DGIT_HASH=\"$(git rev-parse --short --verify HEAD)\"
         --cxxopt=-DTODAY=\"$(date +%Y-%b-%d)\"
         --test_env=${C3TK_DATA_PERSISTENT}=/full/path/to/LillyMolQueries
         Molecule_Tools:all     <-  or some other target
```

Most will want to put this in a small shell script, and/or add to .bazelrc where
possible.

When building for release, it is convenient to include the git hash and
the date of the build in the executables. That is not necessary, omit those if not needed.
Note that because the date is included with cxxopt, this _will_ cause a daily
recompile. While this is hardly desirable, the benefits are many.

## cmake
The distribution contains `cmake` infrastructure, that is currently
not functional.  Within Lilly we have not been able to make it work,
usually as a result of conflicting protcol buffer versions on the
system.  Work is ongoing to get cmake working for the public release.
