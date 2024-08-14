# Welcome to the Eli Lilly LillyMol implementation.

## Background
LillyMol is a C++ library for Cheminformatics. This repo also contains a variety of useful
command line tools that have been built with LillyMol.

LillyMol does only a subset of Cheminformatics tasks, but tries to do those tasks
efficiently and correctly.

LillyMol has some novel approaches to substructure searching, reaction enumeration and
chemical similarity. These have been developed over many years, driven by the needs
of Computational and Medicinal Chemists at Lilly and elsewhere.

Recent work has focussed on making *de-novo* molecule construction and there are
several tools desiged to either support or complement A/I driven molecule
generation.

LillyMol is fast and scalable, with modest memory requirements.

This release includes a number of C++ unit tests. All
tests can be run with address sanitizer, with no problems reported.

The file [Molecule_Tools/introduction.cc](src/Molecule_Tools/introduction.cc) provides
an introduction to LillyMol for anyone wishing to develop with C++.

This release includes first steps towards more extensive documentation of LillyMol,
see the [docs](/docs) directory. More work is needed on this front. Most parts
of LillyMol have been feature stable for a long time.

## Python
This release include a python interface to LillyMol via pybind11. This first release
includes most Molecule related functionality, substructure searching and reaction
enumeration. In the pybind directory there are some *_test.py files that exemplify
much of the current functionality. Documentation is in [docs](/docs/python/LillyMolPython.md).
This should be sufficient support for a great many tasks involving querying
or manipulation of molecules at the connection table level.

The current roadmap for the python interface primarily involves two directions

* Enabling gfp fingerprints for similarity calculations.
* Making existing LillyMol applications available.

There is already a Julia interface to an earlier version of LillyMol, and this
release will soon be adapted to support Julia.

## REQUIREMENTS:

**Note**: This is significantly different from prior versions.

LillyMol is primarily developed on RedHat and Ubuntu systems.

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
module load gcc10
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
Note that with the default build (below) Python bindings are not built.


Make sure that python-dev and libblas-dev are installed.

```
sudo apt install python-dev libblas-dev
```
Things seem to work seamlessly in virtualenv.

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
in the repo.

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
Within the src directory, the file `WORKSPACE` configures the build environment
for `bazel`. If you are building python bindings, this file needs to be updated
to reflect the location of your local python. The script `update_bazel_configs.sh`
does this automatically from the Makefile.

### Installation Directory
There is an 'install' target in the BUILD files, and defined in
[build_deps/install.bzl](src/build_deps/install.bzl).
This is where LillyMol executables will be installed when
the 'install' run target is run. Again `update_bazel_configs.sh`
will update this to '/path/to/LillyMol/bin/$(uname)'.
Check to see that the update has
been done correctly and adjust if not, or to set another location.

```
tail build_deps/install.bzl
```
That file contains other mechanisms for specifying the install directory.
But remember, every time `make` is run, that file will be automatically
updated again. Remove the calls to `update_bazel_configs.sh` from
the Makefile if needed.

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
- **eigen**: matrix operations
- **googletest**: Google's c++ unit tests
- **protobuf**: Google's Protocol Buffers
- **re2**: Google's regular expression library
- **tbb**: Threaded Building Blocks for multi-threading
- **zlib**: compression

The complete listing is in the file [MODULE.bazel](src/MODULE.bazel).

Other third party dependencies are downloaded and built by the
[build_third_party.sh](src/build_third_party.sh) script, which will create a `third_party`
directory (next to src) and then download, build and install the following dependencies

- **BerkeleyDb**: used for key/value databases
- **f2c/libf2c**: there is some fortran in LillyMol.

Running 'build_third_party.sh' needs to be done once.

Note that BerkeleyDB and Python bindings are only built if requested. 
In [Makefile](/Makefile) you will see use of the shell variables
'BUILD_PYTHON' and 'BUILD_BDB' which if set, enables building of
these optional features. These can be set any time.

Running `build_third_party.sh` may be a lengthy process. It can be re-run at
any time thereafter. For those repos that are cloned GitHub repos, it will
pull a new version and build. Remove the entire `third_party` directory and
re-run the script and all dependencies are downloaded and rebuilt. If there is
an individual dependency that you would like to rebuild, just remove it from
the `third_party` directory, run the script again and it will be rebuilt.

Note too that installing these external dependencies and running bazel may require
considerable amounts of disk space. For example at the time of writing my
'third_party' directory contains 1.2GB and my bazel temporary area contains 2.2GB.


Note that [.bazelrc](/src/.bazelrc) contains a hardware restriction to quite old
Intel hardware. You should update update `--cxxopt` to reflect
your hardware. Using `--cxxopt=-march=native --cxxopt=-mtune=native` is likely
what you want. Build for the local hardware.

### Python Bindings
During building of external dependencies (with build_third_party.sh
and if BUILD_PYTHON is set)
the script `update_python_in_workspace.py` will examine your python
installation and get information about the include path. With that
info it will update WORKSPACE with new values for the 'path'
attributes of the python related features.

Note that if it does not find a pybind11 installation, the build
will continue, but the python related parts of the build will subsequently fail.

You can of course manually update WORKSPACE to point to your
python installation. See the 'new_local_repository' sections for 'python'
and 'pybind11'

Note that we recently observed a change in how shared
libraries are handled by bazel. For now, there is a .bazelversion file that 
freezes the bazel version until we figure out how to handle shared libraries
with newer versions of bazel. Having a .bazelversion file makes use of bazelisk
superfluous, but once we figure out the new shared library stuff, the bazel
version will again be allowed to float via bazelisk. The current way
shared libraries are handled is not ideal, and causes some undesirable
behaviour in LillyMol/python.

# Build
Once the third party dependencies have been built, and WORKSPACE and
install.bzl configured, LillyMol building can begin. 

Bazel needs to be able to store its cache on a local disk, *not* NFS. When building
inside Lilly, I have used `--output_user_root=/node/scratch/${USER}` to
use local scratch storage for bazel's cache. Note that if there is a
recycling policy in place for the cache, you may see unexpected outcomes.
Purge the cache completely to start afresh. 

If outside Lilly, the 'build_from_src.sh' script (below) will check to
see if your HOME directory is on an NFS mounted file system, and if so, will
specify /tmp for bazel's cache. This is almost certainly not what you want,
so edit 'build_from_src.sh' to specify a local directory for
`--output_user_root`. Again, only needed if you are on an NFS file system.
You can also enter this value in bazel's configuration file `.bazelrc`.

By default, bazel will use all cores available on the local machine.
If needed, limit the number of cores with the `--jobs` option inside
'build_from_src.sh' (sorry no command line options here).

Optionally set shell variables BUILD_BDB and BUILD_PYTHON to enable
building of optional features.

Once the bazel preconditions are set, do the build, test and installs
```
cd src                    # you might already be here
./build_from_src.sh       # takes a while
```

The script will

1. run the C++ unit tests,
2. build all executables
3. build the python bindings
4. install executables into the `/path/to/LillyMol/bin/$(uname)` directory (build_deps/install.bzl)
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
