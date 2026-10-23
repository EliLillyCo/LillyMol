# rand10k LillyMol Binding Benchmark Results

These results use `/home/ian/rand10k.smi`, a 10,000 molecule random sample. This
is a more useful Python binding benchmark than the built-in smoke set because it
runs enough molecules for per-call overhead to show up consistently.

## Run Context

- Date: 2026-08-21
- Source checkout HEAD at run time: `24e9f3da`
- Python: `3.13.7`, `/home/ian/lillymol_py313_venv/bin/python`
- NumPy: `2.5.2`
- Platform: `Linux-6.8.0-137-generic-x86_64-with-glibc2.39`
- Input: `/home/ian/rand10k.smi`
- Molecules: 10,000
- Loops per benchmark: 5
- Raw pybind results: `results_rand10k_pybind.json`
- Raw nanobind results: `results_rand10k_nanobind.json`

Commands:

```shell
PYTHONPATH="$PWD/bazel-bin/pybind:$PWD/bazel-bin" \
LILLYMOL_HOME=/home/ian/LillyMol_IWFORK \
/home/ian/lillymol_py313_venv/bin/python \
  nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding pybind --smiles-file /home/ian/rand10k.smi \
  --loops 5 --json \
  > nanobind/benchmarks/results_rand10k_pybind.json

PYTHONPATH="$PWD/bazel-bin/nanobind:$PWD/bazel-bin" \
LILLYMOL_HOME=/home/ian/LillyMol_IWFORK \
/home/ian/lillymol_py313_venv/bin/python \
  nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding nanobind --smiles-file /home/ian/rand10k.smi \
  --loops 5 --json \
  > nanobind/benchmarks/results_rand10k_nanobind.json
```

## Common Benchmarks

Lower `best_s` is faster. The ratio is `nanobind best_s / pybind best_s`.

| Benchmark | pybind best s | nanobind best s | Ratio |
| --- | ---: | ---: | ---: |
| atom_index_access | 0.173033 | 0.063464 | 0.37 |
| atom_iteration | 0.440423 | 0.105229 | 0.24 |
| coordinates_list | 0.140071 | 0.068161 | 0.49 |
| coordinates_numpy | 0.156259 | 0.098843 | 0.63 |
| ecfp | 0.242689 | 0.246199 | 1.01 |
| ecfp_numpy | 0.308462 | 0.278444 | 0.90 |
| iwdescr | 2.768553 | 2.655178 | 0.96 |
| iwdescr_list | 2.277313 | 2.273640 | 1.00 |
| linear_fingerprint | 1.027003 | 1.041146 | 1.01 |
| linear_fingerprint_numpy | 1.063631 | 1.037910 | 0.98 |
| parse | 0.305769 | 0.049755 | 0.16 |
| scalar_properties | 0.009010 | 0.002534 | 0.28 |
| substructure_count | 0.186422 | 0.170045 | 0.91 |
| unique_smiles | 0.004524 | 0.002784 | 0.62 |

## Skips

The only skipped benchmark in this run was pybind
`substructure_embeddings`, because the pybind `SubstructureQuery` module does
not expose `substructure_search_match_lists`. Nanobind completed all benchmark
cases.
