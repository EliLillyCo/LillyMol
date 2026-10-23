# LillyMol Python Binding Benchmarks

These benchmarks time small Python-facing LillyMol workflows. They are timing
tools, not pass/fail regression tests.

The production Python binding is the nanobind-backed `lillymol` module. Build it
first from `src`:

```shell
bazel build -c opt //nanobind:lillymol
```

Run the current binding benchmark:

```shell
PYTHONPATH="$(pwd)/bazel-bin/nanobind:$(pwd)/bazel-bin" \
  python3 nanobind/benchmarks/benchmark_lillymol_bindings.py --binding nanobind
```

The script still has a `--binding pybind` mode for historical comparisons in
worktrees where the old pybind targets are intentionally built, but the hard
changeover no longer builds or tests those targets.

For a useful ChEMBL-derived input set, a few thousand molecules is enough for
quick iteration and 50k-100k molecules is useful for steadier whole-workflow
measurements. Use a normal LillyMol SMILES file, one molecule per line:

```shell
python3 nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding nanobind \
  --smiles-file /path/to/random_chembl.smi \
  --max-mols 10000
```

Descriptor benchmarks require `LILLYMOL_HOME` to point at a tree containing the
standard query data. NumPy-specific benchmarks are skipped if NumPy is not
installed or if the selected binding does not expose the needed method.

Checked-in baseline runs:

- `results_builtin.md`: 10 built-in molecules, useful as a smoke benchmark.
- `results_rand10k.md`: 10,000 random molecules from `/home/ian/rand10k.smi`,
  useful as the current larger comparison baseline.

Raw JSON is stored beside each summary.
