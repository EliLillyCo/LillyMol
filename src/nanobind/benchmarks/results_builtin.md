# Built-in LillyMol Binding Benchmark Results

These are small smoke-benchmark results for the pybind11 and nanobind Python
bindings. They use the 10 built-in molecules in
`benchmark_lillymol_bindings.py`, so the numbers are useful for quick local
comparison rather than final performance claims.

## Run Context

- Date: 2026-08-21
- Source checkout HEAD at run time: `24e9f3da`
- Python: `3.13.7`, `/home/ian/lillymol_py313_venv/bin/python`
- NumPy: `2.5.2`
- Platform: `Linux-6.8.0-137-generic-x86_64-with-glibc2.39`
- Molecules: 10 built-in SMILES
- Loops per benchmark: 20
- Raw pybind results: `results_builtin_pybind.json`
- Raw nanobind results: `results_builtin_nanobind.json`

Commands:

```shell
PYTHONPATH="$PWD/bazel-bin/pybind:$PWD/bazel-bin" \
LILLYMOL_HOME=/home/ian/LillyMol_IWFORK \
/home/ian/lillymol_py313_venv/bin/python \
  nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding pybind --loops 20 --json \
  > nanobind/benchmarks/results_builtin_pybind.json

PYTHONPATH="$PWD/bazel-bin/nanobind:$PWD/bazel-bin" \
LILLYMOL_HOME=/home/ian/LillyMol_IWFORK \
/home/ian/lillymol_py313_venv/bin/python \
  nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding nanobind --loops 20 --json \
  > nanobind/benchmarks/results_builtin_nanobind.json
```

## Common Benchmarks

Lower `best_s` is faster. The ratio is `nanobind best_s / pybind best_s`.

| Benchmark | pybind best s | nanobind best s | Ratio |
| --- | ---: | ---: | ---: |
| atom_index_access | 0.000094 | 0.000024 | 0.26 |
| atom_iteration | 0.000196 | 0.000053 | 0.27 |
| coordinates_list | 0.000115 | 0.000035 | 0.30 |
| coordinates_numpy | 0.000116 | 0.000065 | 0.56 |
| ecfp | 0.000155 | 0.000117 | 0.76 |
| ecfp_numpy | 0.000216 | 0.000141 | 0.65 |
| iwdescr | 0.003110 | 0.002964 | 0.95 |
| iwdescr_list | 0.003148 | 0.002918 | 0.93 |
| linear_fingerprint | 0.000386 | 0.000357 | 0.92 |
| linear_fingerprint_numpy | 0.000445 | 0.000360 | 0.81 |
| parse | 0.000203 | 0.000059 | 0.29 |
| scalar_properties | 0.000014 | 0.000003 | 0.18 |
| substructure_count | 0.000121 | 0.000106 | 0.88 |
| unique_smiles | 0.000009 | 0.000002 | 0.26 |

## Skips

The only skipped benchmark in this run was pybind
`substructure_embeddings`, because the pybind `SubstructureQuery` module does
not expose `substructure_search_match_lists`. Nanobind completed all benchmark
cases.
