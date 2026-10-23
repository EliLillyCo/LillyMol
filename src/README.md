# LillyMol Source Tree

This directory contains the LillyMol C++ source, command-line tools, tests, and
Python bindings. Most development work starts here.

For repository-level background, see [../README.md](../README.md). For build
setup details, see [../Build.md](../Build.md). For compact agent guidance, see
[../AGENTS.md](../AGENTS.md).

## Directory Map

- `Foundational/` contains low-level utilities used throughout LillyMol,
  including strings, command-line parsing, I/O helpers, and general data
  structures.
- `Molecule_Lib/` contains the core molecule representation, SMILES and structure
  file I/O, aromaticity, chirality, substructure searching, reactions,
  standardisation, and atom typing.
- `Molecule_Tools/` contains most command-line molecular tools and reusable tool
  libraries.
- `Molecule_Tools_Bdb/` contains BerkeleyDB-backed tools and libraries.
- `Utilities/GFP_Tools/` contains fingerprint, similarity, nearest-neighbour,
  and distance-matrix tools.
- `nanobind/` contains the active Python bindings. New Python functionality
  should generally be added here.
- `pybind/` contains the former pybind11 bindings. It remains for reference
  during the transition, but is no longer the active binding target.
- `Utilities/`, `General/`, and other subdirectories contain older or more
  specialised utilities. Inspect local BUILD files and nearby code before making
  changes.

## Building

Prefer Bazel builds from this directory. A portable optimised build looks like:

```shell
bazel build -c opt //Molecule_Tools:fileconv
```

For local development or benchmarking on the current machine, use the native
configuration:

```shell
bazel build -c opt --config=native //Molecule_Tools:fileconv
```

`--config=native` enables local CPU-specific compiler options. Do not use it for
portable or container build artifacts unless native CPU targeting is intended.

Bazel outputs are written under `bazel-bin`. A built command-line tool can be run
from its package path, for example:

```shell
bazel-bin/Molecule_Tools/fileconv -v input.smi
```

Many LillyMol tools infer input type from filename suffixes. Many also support
`-` for standard input or standard output, but check the specific tool before
relying on that behavior.

## Testing

Use focused Bazel tests while developing C++ libraries and tools:

```shell
bazel test -c opt --config=native //Molecule_Tools:ring_system_shape_test
```

Executable regression tests live under `../test` and are usually run with the
Ruby test driver:

```shell
ruby ../test/run_all_test.rb -bindir bazel-bin/Molecule_Tools fileconv
```

Build with `-c opt` before running executable regression tests. Some expected
outputs depend on optimised floating-point behavior and can differ in the last
printed digit from `fastbuild` binaries.

Small, focused tests are usually preferable to large end-to-end fixtures when
changing an older tool with many options.

## Python Bindings

The active Python bindings are built with nanobind. The installed module name is
`lillymol`; BerkeleyDB-dependent bindings are in `lillymol_bdb`.

Build the bindings and stage the shared libraries into the repository `lib`
directory:

```shell
bazel build -c opt --config=native //nanobind:lillymol
./copy_shared_libraries.sh
```

Run Python scripts through the repository wrapper:

```shell
../run_python.sh script.py
```

The wrapper sets `PYTHONPATH` and `LD_LIBRARY_PATH` to use the staged libraries.
It also checks `../lib/lillymol.so.soabi` against the active Python interpreter
ABI, which helps catch mismatches between the Python used to build the extension
and the Python used at runtime.

Focused nanobind tests can be run with:

```shell
bazel test -c opt --config=native //nanobind:lillymol_test --test_output=errors
```

Many Python examples and API notes are in [../docs/python](../docs/python).

## Documentation

User documentation is under [../docs](../docs). Common starting points are:

- [../docs/Molecule_Tools](../docs/Molecule_Tools) for command-line tools.
- [../docs/Molecule_Lib](../docs/Molecule_Lib) for core molecule-library
  behavior.
- [../docs/GFP](../docs/GFP) for fingerprints, similarity, and distance tools.
- [../docs/python](../docs/python) for Python bindings.
- [../docs/Workflows](../docs/Workflows) for higher-level task recipes.

When updating tool behavior, update the corresponding documentation if the
user-visible semantics change. Prefer concrete command examples and keep option
descriptions tied to observed code behavior.

## Development Notes

Follow the style of nearby code. Much of LillyMol predates modern C++; avoid
large modernization passes while making targeted behavioral changes.

Older code commonly uses LillyMol container and string types such as `IWString`,
`const_IWSubstring`, `IWString_and_File_Descriptor`, `resizable_array`,
`resizable_array_p`, and `Set_of_Atoms`. Prefer existing local types and helper
APIs when working in those areas.

Many command-line tools have broad option surfaces and some file-scope global
state. Preserve command-line compatibility unless a change explicitly requires a
break. For reusable behavior, prefer moving logic into a library with a thin CLI
wrapper so it can also be tested directly and exposed to Python when needed.

Be careful with output files in tools that make more than one input pass. Opening
an output file too early can truncate an input file that has not yet been read.
For molecule-writing tools, check for input/output clobbering where nearby code
uses helpers such as `would_overwrite_input_files`.

Use `rg` for source searches. Inspect the relevant source, tests, and docs before
changing behavior.
