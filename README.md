# Welcome to the Eli Lilly LillyMol implementation.

## Background
LillyMol is a C++ library for Cheminformatics. This repo also contains a variety of useful
command line tools that have been built with LillyMol.

Lillymol does only a subset of Cheminformatics tasks, but tries to do those tasks
efficiently and correctly.

Lillymol has some novel approaches to substructure searching, reaction enumeration and
chemical similarity. These have been developed over many years, driven by the needs
of Computational and Medicinal Chemists at Lilly and elsewhere.

LillyMol is fast. The public release now includes a number of C++ unit tests. All
tests can be run with address sanitizer, with no problems reported.

The file [Molecule_Tools/introduction.cc](src/Molecule_Tools/introduction.cc) provides
an introduction to LillyMol for anyone wishing to develop with C++.

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
which makes keeping bazel up to date easier. 

Within a GitHub CodeSpace, this worked to install bazelisk.
```
sudo apt install npm
sudo npm install -g @bazel/bazelisk
```
If you use the module system
```
module load bazelisk
```

You **WILL** need to edit bazel's configuration file [WORKSPACE](/src/WORKSPACE) in order to
make things work on your file system. More below.

Support for `cmake` is rudimentary. It sometimes works. Improvements would
be welcome. Some third party dependencies may depend on recent versions of `cmake`.

If you use the module system
```
module load cmake
```

## Dependencies.

There are several dependencies which could be installed on the system,
which would considerably simplify the build configuration, but during
development we have frequently found ourselves on machines that could
not be updated to the versions we needed, or where we lacked priveleges.
So external dependencies are managed explicitly.

The software requires a gcc version of at least version 10. This version of LillyMol
uses some fairly recent c++ features, which require a recent compiler. 

If you use the module system
```
module load gcc10
module load bazelisk
```
The preferred way of using third party software is via the Bazel
Module system. Most of the external dependencies needed are handled
via that mechanism. Today that includes

- **absl**: Google's c++ library - we use crc32c and some data structures.
- **googletest**: c++ unit tests
- **protobuf**: handling of protocol buffers
- **re2**: Google's regular expression library
- **zlib**: compression

The complete listing is in the file [MODULE.bazel](src/MODULE.bazel).


Other third party dependencies are downloaded and built by the
[build_third_party.sh](src/build_third_party.sh) script, which will create a `third_party`
directory and build the following dependencies.

- **BerkeleyDb**: used for key/value databases
- **f2c/libf2c**: there is some fortran in LillyMol.

The `build_third_party.sh` script will download, configure, build and install
the requisite third party tools. 
See the `WORKSPACE` file if you would prefer to do a local install of
those tools instead. The difference is that if you have bazel silently
handle those, they are not readily available to you if you decide you
need them for other purposes. Doing a local install will also avoid
possibly needless recompilations.

By default the third_party tools are installed next to the `src` directory.

Note that installing these external dependencies may require considerable
amounts of disk space. For example at the time of writing my bazel
temporary area contains 847MB.

Note too that there are some third party dependencies that are made
part of the LillyMol repo, xmlParser, jama, tnt, and various code
snippets.

`WORKSPACE` configures the build environment for `bazel`. The absolute path of
the dependencies must be configured in that file. Wherever you see a `path=`
directive in WORKSPACE, that just be updated to reflect where the repo
is on your system. The `build_third_party` script does attempt to generate
an updated WORKSPACE file. At the end of that script, it places a possibly
useful new WORKSPACE file in `/tmp`. It is up to you to review that file, and
if it looks OK, copy it to the `src` directory, overwriting the existing
WORKSPACE file. The WORKSPACE file that is in the repo will **not** work.

If some third party tools do not build, you may be able to skip building
them during the build by use of --build_tag_filters=. For example to *not*
build any tool that involves BerkeleyDb `--build_tag_filters=-berkeleydb`
should do that. Tags are specified in the respective BUILD files.

```
cd src
./build_third_party.sh   # takes a while
ls -l /tmp/WORKSPACE  # should have just been created
vi /tmp/WORKSPACE   # inspect and update if needed
cp /tmp/WORKSPACE WORKSPACE   # press the new version into service.
```
If you are lucky, the WORKSPACE file in /tmp will require no changes.

There is an install target in the BUILD files, and defined in [build_deps/install.bzl](src/build_deps/install.bzl).
This is used to copy built executables to a pre-defined location. There are
two ways of specifying this.

1. hard code the directory into the `default = ""` specification at the bottom of
   `install.bzl`
2. use the `--action_env=` directive below to specify BINDIR.

Generally it is easier to do a one-time change to install.bzl, although you then
need to remember to not commit your install path to the repo. Running `build_third_party.sh`
will again create what is likely to be an ok version of this file in `/tmp/install.bzl`. If
it looks OK, copy it to the `build_deps` directory.

Running `build_third_party.sh` may be a lengthy process. It can be re-run at
any time thereafter. For those repos that are cloned GitHub repos, it will
pull a new version and build. Remove the entire `third_party` directory and
re-run the script and all dependencies are downloaded and rebuilt. If there is
an individual dependency that you would like to rebuild, just remove it from
the `third_party` directory, run the script again and it will be rebuilt.

Note that [.bazelrc](src/.bazelrc) contains a hardware restriction to quite old
Intel hardware. You should update update `--cxxopt` to reflect
your hardware. Using `--cxxopt=-march=native --cxxopt=-mtune=native` is likely
what you want. Build for the local hardware.

## Building
When building for release, it is convenient to include the git hash and
the date of the build in
the executables. Make sure you are building from a clean repo.

Bazel needs to be able to store its cache on a local disk. When building
inside Lilly, I have used `--output_user_root=/node/scratch/${USER}` to
use local scratch storage for bazel's cache. Note that if there is a
recycling policy in place for the cache, you may see unexpected outcomes.
Remove the cache completely to start afresh.

By default, bazel will use all cores available on the local machine.
If needed, limit the number of cores with the `--jobs` option.

There is an install run target built into the BUILD files, and you must
decide upon the directory to which that target should copy executables.
See previous section.

While there is a script that should do the build for you, underneath these
are the commands being issued. Skip to the next section unless interested.

A typical build command might be (change the path for test_env)

```
bazelisk --output_user_root=/node/scratch/${USER}
         build
         --jobs=8
         --enable_bzlmod
         --experimental_cc_shared_library
         -c opt
         --cxxopt=-DGIT_HASH=\"$(git rev-parse --short --verify HEAD)\"
         --cxxopt=-DTODAY=\"$(date +%Y-%b-%d)\"
         //Molecule_Tools:all     <-  targets listed below
         --action_env=BINDIR=/path/for/binaries
         --test_env=C3TK_DATA_PERSISTENT=/full/path/to/LillyMolQueries
```

Again, the binary installation directory can be either specified here, or hard coded in
`install.bzl`. If hard coded in `install.bzl` omit the `--action_env` option - that is
what I do.

Note that with the advent of python bindings, we observed a change in how shared
libraries are handled by bazel. For now, there is a .bazelversion file that 
freezes the bazel version until we figure out how to handle shared libraries
with newer versions of bazel. Having a .bazelversion file makes use of bazelisk
superfluous, but once we figure out the new shared library stuff, the bazel
version will be allowed to float via bazelisk.

## Targets to Test and Build

First run `build_third_party.sh` to ensure dependencies are built.

Actual tests, builds and installs can be handled by the `build_from_src.sh` script.

The script `build_from_src.sh` will

1. run the C++ unit tests,
2. build all executables
3. install executables into the `bin/Linux` directory (build_deps/install.bzl)

## cmake
There is `cmake` infrastructure
included, but that may, or may not, work - within Lilly we have not
been able to make it work, usually as a result of conflicting protcol buffer
versions on the system, but externally it often works, although it has not
been tested recently.

## Recipe to build (inside Lilly).

```bash
module load gcc10
module load bazelisk
cd ./src

# Do these steps once
./build_third_party.sh
# Examine /tmp/WORKSPACE
# cp /tmp/WORKSPACE WORKSPACE
# Examine /tmp/install.bzl
# cp /tmp/install.bzl build_deps
./build_from_src.sh

# Rebuild Molecule_Tools any time
bazelisk --output_user_root=/node/scratch/$USER build --jobs=8 --enable_bzlmod --experimental_cc_shared_library -c opt --cxxopt=-DGIT_HASH=\"$(git rev-parse --short --verify HEAD)\" --cxxopt=-DTODAY=\"$(date +%Y-%b-%d)\" //Molecule_Tools:all --test_env=C3TK_DATA_PERSISTENT=/full/path/to/LillyMol/queries
```
Note that because the date is included with cxxopt, this _will_ cause a daily
recompile. While this is hardly desirable, knowing when an executable was compiled
can be useful.