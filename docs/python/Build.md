# Building LillyMol Python

The public LillyMol Python modules are built from the nanobind bindings. During
local development, build them with Bazel and stage the shared objects into the
repository `lib` directory used by `run_python.sh`.

## Required Python Packages

Use a recent Python 3. Development has been exercised with Python 3.13, and the
bindings are expected to work with other supported Python 3 versions.

```shell
python -m pip install nanobind absl-py protobuf numpy
```

At system level, make sure the Python headers and BLAS development libraries are
installed. Package names vary by distribution; on Debian/Ubuntu systems this is
typically:

```shell
sudo apt install python3-dev libblas-dev
```

`numpy` is needed for the array-returning fingerprint, descriptor, and coordinate
helpers. `protobuf` is needed for proto-backed configuration objects and tests.

## Optional Python Packages

If you use the xgboost QSAR model building tools in LillyMol:

```shell
python -m pip install xgboost scikit-learn matplotlib pandas optuna
```

If you use the molecular property profile tool:

```shell
python -m pip install scipy numpy matplotlib
```

Virtual environments and `uv` both work well. For example:

```shell
uv venv /path/to/lillymol_venv
uv pip install --python /path/to/lillymol_venv/bin/python numpy protobuf
```

## Local Development

From the source checkout:

```shell
cd ${LILLYMOL_HOME}/src
bazel build --config=native -c opt //nanobind:all
./copy_shared_libraries.sh ../lib
```

Then run Python through the top-level wrapper, which sets `PYTHONPATH` and
`LD_LIBRARY_PATH` for the staged extension modules and LillyMol shared libraries:

```shell
${LILLYMOL_HOME}/run_python.sh my_script.py
```

Use `-c opt` for normal testing and timing. Bazel's default `fastbuild` is much
slower and can produce tiny floating point differences from the optimised command
line tools.

### Optional 2D depiction module

The `lillymol_depict` module depends on Schrodinger's coordgen library. Build it
through the normal build script by setting both Python and depiction support:

```shell
BUILD_PYTHON=1 BUILD_SCHRODINGER_2D=1 ./build_linux.sh
```

For an existing coordgen installation in `third_party`, the focused Bazel build
and test are:

```shell
bazel build -c opt //nanobind:lillymol_depict
bazel test -c opt //nanobind:lillymol_nb_depict_test
./copy_shared_libraries.sh ../lib
```

The main `lillymol` extension does not link coordgen and remains available when
depiction support is not built.

## Wheel Packaging

The wheel infrastructure lives in `${LILLYMOL_HOME}/python`. It packages
Bazel-built extension modules and private LillyMol shared libraries; setuptools
does not compile the C++ code.

```shell
cd ${LILLYMOL_HOME}/src
bazel build -c opt //nanobind:all
./copy_shared_libraries.sh ../lib

cd ${LILLYMOL_HOME}/python
./scripts/stage_wheel_files.sh
python -m pip install build wheel setuptools
python -m build --wheel
```

The wheel is written to `${LILLYMOL_HOME}/python/dist`. The generated
`python/prebuilt`, `python/build`, and `python/dist` directories are build
artifacts. Commit only the wheel-building infrastructure, not staged binaries or
generated wheels.

The wheel should install the nanobind-backed `lillymol` module, the
`lillymol_bdb` BerkeleyDB helper module when BDB support is built, and small
compatibility modules such as `lillymol_io` and `lillymol_tools` that re-export
from `lillymol`.
