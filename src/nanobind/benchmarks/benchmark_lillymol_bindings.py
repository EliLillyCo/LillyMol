#!/usr/bin/env python3
"""Microbenchmarks for LillyMol pybind11 vs nanobind bindings.

The benchmark names are intentionally small and binding-oriented. For each
benchmark, `items` is the number of molecules or operation groups processed, and
`loops` is the number of full benchmark repetitions performed by the harness.
"""

from __future__ import annotations

import argparse
import gc
import importlib
import json
import os
import platform
import statistics
import sys
import time
from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

BUILTIN_SMILES = [
    "CCO ethanol",
    "c1ccccc1 benzene",
    "CC(=O)Oc1ccccc1C(=O)O aspirin",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine",
    "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4 atorvastatin_like",
    "O=C(NCc1ccccc1)c1ccc(O)cc1 benzamide",
    "CCN(CC)C tertiary_amine",
    "C1CCCCC1 cyclohexane",
    "COc1cc(C=O)ccc1O vanillin",
    "N[C@@H](CC(=O)O)C(=O)O aspartic_acid",
]


@dataclass
class Binding:
    name: str
    core: object
    fingerprint: object
    tools: object
    query: object


def _import_optional(name: str):
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError:
        return None


def load_binding(name: str) -> Binding:
    if name == "pybind":
        core = importlib.import_module("lillymol")
        fingerprint = _import_optional("lillymol_fingerprint") or core
        tools = _import_optional("lillymol_tools") or core
        query = _import_optional("lillymol_query") or core
        return Binding(name, core, fingerprint, tools, query)
    if name == "nanobind":
        core = importlib.import_module("lillymol")
        return Binding(name, core, core, core, core)
    raise ValueError(f"Unrecognised binding '{name}'")


def read_smiles(path: str | None, max_mols: int) -> list[str]:
    if path is None:
        smiles = list(BUILTIN_SMILES)
    else:
        smiles = []
        with open(path, "r", encoding="utf-8") as reader:
            for line in reader:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                smiles.append(line)
                if max_mols and len(smiles) >= max_mols:
                    break
    if not smiles:
        raise ValueError("No molecules available for benchmarking")
    return smiles


def mol_from_smiles(binding: Binding, smiles: str):
    mol = binding.core.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse {smiles!r}")
    return mol


def build_molecules(binding: Binding, smiles: Sequence[str]):
    return [mol_from_smiles(binding, smi) for smi in smiles]


def maybe_numpy():
    try:
        import numpy as np  # pylint: disable=import-outside-toplevel
        return np
    except ModuleNotFoundError:
        return None


def coordinate_vector(mol) -> list[float]:
    coords: list[float] = []
    for i in range(mol.natoms()):
        coords.extend([float(i), float(i % 3), float((i * 7) % 5)])
    return coords


@dataclass
class Benchmark:
    name: str
    func: Callable[[Binding, list[str], list[object]], int]
    min_mols: int = 1


class SkipBenchmark(Exception):
    pass


def bench_parse(binding: Binding, smiles: list[str], _mols: list[object]) -> int:
    total = 0
    for smi in smiles:
        total += mol_from_smiles(binding, smi).natoms()
    return total


def bench_scalar_properties(_binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    total = 0
    for mol in mols:
        total += mol.natoms()
        total += mol.nrings()
        total += mol.nedges()
    return total


def bench_unique_smiles(_binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    total = 0
    for mol in mols:
        total += len(mol.unique_smiles())
    return total


def bench_atom_index_access(_binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    total = 0
    for mol in mols:
        for i in range(mol.natoms()):
            total += mol.atomic_number(i)
            total += mol.ncon(i)
    return total


def bench_atom_iteration(_binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    total = 0
    for mol in mols:
        for atom in mol:
            total += atom.atomic_number()
            total += atom.ncon()
    return total


def make_query(binding: Binding, smarts: str):
    query = binding.query.SubstructureQuery()
    ok = query.build_from_smarts(smarts)
    if not ok:
        raise ValueError(f"Cannot build SMARTS {smarts!r}")
    return query


def bench_substructure_count(binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    query = make_query(binding, "c1ccccc1")
    total = 0
    for mol in mols:
        total += int(query.substructure_search(mol))
    return total


def bench_substructure_embeddings(binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    query = make_query(binding, "[#6]")
    if not hasattr(query, "substructure_search_match_lists"):
        raise SkipBenchmark("SubstructureQuery.substructure_search_match_lists unavailable")
    total = 0
    for mol in mols:
        total += len(query.substructure_search_match_lists(mol))
    return total


def bench_linear_fingerprint(binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    if not hasattr(binding.fingerprint, "linear_fingerprint"):
        raise SkipBenchmark("linear_fingerprint unavailable")
    total = 0
    for mol in mols:
        fp = binding.fingerprint.linear_fingerprint(mol)
        total += len(fp)
    return total


def bench_linear_fingerprint_numpy(binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    np = maybe_numpy()
    if np is None:
        raise SkipBenchmark("NumPy unavailable")
    if hasattr(binding.fingerprint, "linear_fingerprint_numpy"):
        func = binding.fingerprint.linear_fingerprint_numpy
    elif hasattr(binding.fingerprint, "linear_fingerprint"):
        func = binding.fingerprint.linear_fingerprint
    else:
        raise SkipBenchmark("linear fingerprint unavailable")
    total = 0
    for mol in mols:
        fp = func(mol)
        total += int(np.asarray(fp).sum())
    return total


def bench_ecfp(binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    if not hasattr(binding.fingerprint, "ECFingerprintCreator"):
        raise SkipBenchmark("ECFingerprintCreator unavailable")
    creator = binding.fingerprint.ECFingerprintCreator(1024)
    total = 0
    for mol in mols:
        total += len(creator.fingerprint(mol))
    return total


def bench_ecfp_numpy(binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    np = maybe_numpy()
    if np is None:
        raise SkipBenchmark("NumPy unavailable")
    if not hasattr(binding.fingerprint, "ECFingerprintCreator"):
        raise SkipBenchmark("ECFingerprintCreator unavailable")
    creator = binding.fingerprint.ECFingerprintCreator(1024)
    method = getattr(creator, "fingerprint_numpy", creator.fingerprint)
    total = 0
    for mol in mols:
        total += int(np.asarray(method(mol)).sum())
    return total


def bench_coordinates_list(_binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    total = 0
    for mol in mols:
        coords = coordinate_vector(mol)
        mol.set_coordinates(coords)
        total += len(mol.get_coordinates())
    return total


def bench_coordinates_numpy(_binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    np = maybe_numpy()
    if np is None:
        raise SkipBenchmark("NumPy unavailable")
    total = 0
    for mol in mols:
        coords = np.asarray(coordinate_vector(mol), dtype=np.float32)
        if hasattr(mol, "set_coordinates_numpy"):
            mol.set_coordinates_numpy(coords)
            total += mol.get_coordinates_numpy().shape[0]
        else:
            mol.set_coordinates(coords)
            total += np.asarray(mol.get_coordinates(), dtype=np.float32).shape[0]
    return total


def bench_iwdescr(binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    np = maybe_numpy()
    if np is None:
        raise SkipBenchmark("NumPy unavailable")
    if "LILLYMOL_HOME" not in os.environ:
        raise SkipBenchmark("LILLYMOL_HOME not set")
    if not hasattr(binding.tools, "IWDescr"):
        raise SkipBenchmark("IWDescr unavailable")
    calc = binding.tools.IWDescr()
    total = 0
    for mol in mols:
        total += np.asarray(calc.process(mol)).shape[0]
    return total


def bench_iwdescr_list(binding: Binding, _smiles: list[str], mols: list[object]) -> int:
    np = maybe_numpy()
    if np is None:
        raise SkipBenchmark("NumPy unavailable")
    if "LILLYMOL_HOME" not in os.environ:
        raise SkipBenchmark("LILLYMOL_HOME not set")
    if not hasattr(binding.tools, "IWDescr"):
        raise SkipBenchmark("IWDescr unavailable")
    calc = binding.tools.IWDescr()
    return int(np.asarray(calc.process_list(mols)).size)


BENCHMARKS = [
    Benchmark("parse", bench_parse),
    Benchmark("scalar_properties", bench_scalar_properties),
    Benchmark("unique_smiles", bench_unique_smiles),
    Benchmark("atom_index_access", bench_atom_index_access),
    Benchmark("atom_iteration", bench_atom_iteration),
    Benchmark("substructure_count", bench_substructure_count),
    Benchmark("substructure_embeddings", bench_substructure_embeddings),
    Benchmark("linear_fingerprint", bench_linear_fingerprint),
    Benchmark("linear_fingerprint_numpy", bench_linear_fingerprint_numpy),
    Benchmark("ecfp", bench_ecfp),
    Benchmark("ecfp_numpy", bench_ecfp_numpy),
    Benchmark("coordinates_list", bench_coordinates_list),
    Benchmark("coordinates_numpy", bench_coordinates_numpy),
    Benchmark("iwdescr", bench_iwdescr),
    Benchmark("iwdescr_list", bench_iwdescr_list),
]


def run_one(benchmark: Benchmark, binding: Binding, smiles: list[str], mols: list[object], loops: int):
    times = []
    checksum = None
    for _ in range(loops):
        gc.collect()
        start = time.perf_counter()
        value = benchmark.func(binding, smiles, mols)
        elapsed = time.perf_counter() - start
        if checksum is None:
            checksum = value
        elif checksum != value:
            raise RuntimeError(
                f"Benchmark {benchmark.name} checksum changed: {checksum} vs {value}")
        times.append(elapsed)
    return {
        "benchmark": benchmark.name,
        "items": len(smiles),
        "loops": loops,
        "best_s": min(times),
        "median_s": statistics.median(times),
        "mean_s": statistics.fmean(times),
        "checksum": checksum,
    }


def print_table(results: list[dict]) -> None:
    print(f"{'benchmark':28s} {'items':>7s} {'loops':>5s} {'best_s':>11s} {'median_s':>11s} {'items/s(best)':>14s} checksum")
    for row in results:
        rate = row["items"] / row["best_s"] if row["best_s"] else 0.0
        print(
            f"{row['benchmark']:28s} {row['items']:7d} {row['loops']:5d} "
            f"{row['best_s']:11.6f} {row['median_s']:11.6f} {rate:14.1f} {row['checksum']}"
        )


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binding", choices=["pybind", "nanobind"], required=True)
    parser.add_argument("--smiles-file", help="Optional LillyMol SMILES input file")
    parser.add_argument("--max-mols", type=int, default=0,
                        help="Maximum molecules to read from --smiles-file; 0 means all")
    parser.add_argument("--loops", type=int, default=5,
                        help="Number of timing loops for each benchmark")
    parser.add_argument("--benchmark", action="append",
                        help="Benchmark name to run. May be repeated. Default: all")
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of a text table")
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    binding = load_binding(args.binding)
    smiles = read_smiles(args.smiles_file, args.max_mols)
    mols = build_molecules(binding, smiles)

    selected = BENCHMARKS
    if args.benchmark:
        wanted = set(args.benchmark)
        selected = [b for b in BENCHMARKS if b.name in wanted]
        missing = wanted - {b.name for b in selected}
        if missing:
            raise ValueError(f"Unknown benchmarks: {sorted(missing)}")

    results = []
    skipped = []
    for benchmark in selected:
        try:
            results.append(run_one(benchmark, binding, smiles, mols, args.loops))
        except SkipBenchmark as err:
            skipped.append({"benchmark": benchmark.name, "reason": str(err)})
        except ModuleNotFoundError as err:
            skipped.append({"benchmark": benchmark.name, "reason": str(err)})

    if args.json:
        payload = {
            "binding": args.binding,
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "smiles_file": args.smiles_file,
            "molecules": len(smiles),
            "loops": args.loops,
            "results": results,
            "skipped": skipped,
        }
        print(json.dumps(payload, indent=2))
    else:
        print(f"binding: {args.binding}")
        print_table(results)
        if skipped:
            print("\nskipped:")
            for row in skipped:
                print(f"  {row['benchmark']}: {row['reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
