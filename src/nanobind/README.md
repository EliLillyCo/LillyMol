# Nanobind bindings

This directory contains the nanobind implementation of the LillyMol Python
bindings. Nanobind is now the production Python binding path. The older
`src/pybind` tree remains in the repository for reference, but it is no longer
the active build/test path. Python code should not mix pybind11-bound and
nanobind-bound LillyMol objects in one process. Avoid adding dependencies on
pybind-bound classes.

`Atom`, `Bond`, and `Ring` objects returned from a `Molecule` are borrowed views
into that molecule. They are valid while the parent molecule remains alive and
structurally unchanged. Structural mutation of the parent molecule can invalidate
previously returned child objects; the bindings do not add extra ownership or
copying solely to guard against that case. APIs that return independent objects
should do so explicitly, as `rings()` does by returning copied rings.

Current module:

```shell
bazel build -c opt //nanobind:lillymol
bazel test -c opt //nanobind:lillymol_test
```

The binding implementation is split by functional area:

- `lillymol_nb.cc` contains only the `NB_MODULE` entry point;
- `lillymol_nb_internal.h` contains shared declarations for the bindings;
- `lillymol_nb_common.cc` contains shared helper implementations;
- `lillymol_nb_io.cc`, `lillymol_nb_molecule.cc`,
  `lillymol_nb_atom_bond.cc`, `lillymol_nb_set_of_atoms.cc`,
  `lillymol_nb_substructure.cc`, `lillymol_nb_descriptors.cc`,
  `lillymol_nb_standardise.cc`, and `lillymol_nb_fingerprint.cc` contain the
  Python-visible bindings. The filenames retain `_nb` for now to keep the
  changeover diff focused; the installed module name is `lillymol`.

To stage the Python modules:

```shell
./copy_shared_libraries.sh ../lib
../run_python.sh
```

The default destination is the repository-level `lib` directory when run from
`src`. Pass an explicit directory to override it.

The Python 3.13 interpreter installed by `uv` is externally managed, so install
extra packages into a virtual environment rather than modifying that base
interpreter. One working setup is:

```shell
uv venv --python /home/ian/.local/bin/python3.13 /home/ian/lillymol_py313_venv
uv pip install --python /home/ian/lillymol_py313_venv/bin/python numpy protobuf
```

The nanobind tests currently validate:

- constructing, copying, indexing, and iterating over a `Molecule`;
- common `Molecule` accessors, RDKit-style aliases, mutators, bond lookup,
  formal-charge, isotope, Hydrogen, ring, aromaticity, symmetry, path, SMILES,
  atom-renumbering, organic-element, and chirality helpers;
- `MolFromSmiles` returning an owned `Molecule` or `None`;
- basic string, scalar, tuple, optional, and `list[int]` conversions;
- molecule-owned `Atom` and `Bond` references from `atom(i)`, `mol[i]`,
  `for atom in mol`, `atom[i]`, `for bond in atom`, `bond(i)`, and `bonds()`;
- copied `Chiral_Centre` values, `ChiralType`, `is_actually_chiral`, and
  `tetrahedral_chirality`;
- borrowed `Atom` and `Bond` views reflecting parent molecule property mutations;
- sequence-like `Set_of_Atoms` and `Ring` atom-number iteration;
- `Ring` bindings, including both molecule-owned references from `ring(i)` and
  copied rings from `rings()`;
- `Molecule.__contains__` overloads for atomic numbers, element symbols, and
  `SubstructureQuery`;
- `SubstructureQuery` and `SubstructureResults` bindings with copied embeddings
  and borrowed-result iteration via `for embedding in sresults`;
- `TSubstructure` query-list workflows, including batch search, match counts,
  query modifiers, matched-atom reporting, and isotope labelling;
- transient SMARTS convenience helpers (`HasSubstructMatch`,
  `CountSubstructMatches`, `GetSubstructMatches`);
- GIL release around substructure search calls;
- simple descriptor helper dependencies (`alogp`, `xlogp`, `tpsa`);
- basic `Chemical_Standardisation` access via `Standardise`;
- list-returning and NumPy-returning linear, ECFP, and atom-pair fingerprint
  helpers plus `tanimoto`;
- basic molecule file I/O via `Reader`, `Writer`, `MolReaderContext`, and
  `MolWriterContext`, with `ReaderContext` and `ContextWriter` retained as
  compatibility aliases;
- `slurp` and the small global MDL/SDF reader option helpers exposed by the
  Python I/O module.

See `changeover_checklist.md` for remaining changeover packaging tasks and
intentional compatibility decisions.

Remaining cleanup work is mostly deleting obsolete pybind build/test targets and
optionally renaming the implementation files to drop the `_nb` suffix.
