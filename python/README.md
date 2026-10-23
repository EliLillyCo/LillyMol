# LillyMol Python wheel staging

This directory contains the minimal infrastructure needed to build a Python
wheel from Bazel-built LillyMol nanobind modules. Bazel still does the C++ build;
setuptools only packages already-built binaries and Python helper files.

The normal non-wheel runtime workflow is:

```shell
cd ${LILLYMOL_HOME}/src
bazel build -c opt --config=native \
  //nanobind:lillymol \
  //nanobind:lillymol_gfp_server \
  //nanobind:lillymol_bdb
./copy_shared_libraries.sh ${LILLYMOL_HOME}/lib
```

For portable/container builds, omit `--config=native`.

`run_python.sh` uses `${LILLYMOL_HOME}/lib` via `PYTHONPATH` and
`LD_LIBRARY_PATH`. The wheel path is different: staged shared objects are copied
into the wheel and, on Linux when `patchelf` is available, patched to find
private LillyMol shared libraries via `$ORIGIN`.

## Build A Wheel

From `${LILLYMOL_HOME}`:

```shell
cd src
bazel build -c opt --config=native \
  //nanobind:lillymol \
  //nanobind:lillymol_gfp_server \
  //nanobind:lillymol_bdb
./copy_shared_libraries.sh ../lib

cd ../python
./scripts/stage_wheel_files.sh
python -m pip install build wheel setuptools
python -m build --wheel
```

With `uv`, the last build step can be:

```shell
uv build --wheel
```

The wheel is written to `python/dist/`.

## What Gets Staged

`scripts/stage_wheel_files.sh` creates `python/prebuilt/` containing:

- nanobind extension modules: `lillymol.so`, `lillymol_bdb.so`, and
  `lillymol_gfp_server.so`;
- pure-Python compatibility modules such as `lillymol_tools.py` and
  `lillymol_io.py`, which re-export symbols from the canonical `lillymol`
  module;
- private LillyMol shared libraries from `${LILLYMOL_HOME}/lib`;
- selected generated protobuf Python modules and GFP HTTP helper modules.

`python/prebuilt/`, `python/build/`, `python/dist/`, and egg-info directories are
generated artifacts and are ignored by git.

## Runtime Dependencies

The base wheel declares runtime dependencies on `numpy` and `protobuf`.
The GFP HTTP server helpers additionally need:

```shell
pip install 'lillymol[http]'
```

which installs `fastapi` and `uvicorn`.

## Notes

If a wheel installs but importing a module fails with a missing private shared
library, inspect the installed extension with `ldd`. Linux wheels may need a
repair step such as `auditwheel repair`, or staged extension modules may need an
`$ORIGIN` runtime path so they can find private shared libraries shipped in the
same wheel. The staging script applies that rpath with `patchelf` when available.

The `src/pybind` directory is retained for historical reference during the
changeover, but new wheels are built from `src/nanobind` targets.
