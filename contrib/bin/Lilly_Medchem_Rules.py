#!/usr/bin/env python3
"""Python driver for the Lilly Medchem Rules pipeline.

This is intentionally a conservative translation of Lilly_Medchem_Rules.rb.
It preserves the overall shell pipeline construction and calls the existing
LillyMol executables rather than reimplementing chemistry functionality.

Important compatibility choice:
  * Wrapper-owned paths/files are shell-quoted.
  * Pass-through fragments from -tp and -iwd are appended unquoted, so a user can
    write -iwd 'v1 v2' and have v1 and v2 arrive as separate downstream tokens.
"""

from __future__ import annotations

import argparse
import os
import platform
import re
import shlex
import subprocess
import sys
from pathlib import Path


DEFAULT_LOWER_ATOM_COUNT_CUTOFF = 7
DEFAULT_SOFT_UPPER_ATOM_COUNT_CUTOFF = 25
DEFAULT_HARD_UPPER_ATOM_COUNT_CUTOFF = 40


class Odm:
    """Matches demerit names either by exact name or by one RX= regex."""

    def __init__(self, omit_specs: list[str]):
        self._names: set[str] = set()
        self._rx: re.Pattern[str] | None = None

        for spec in omit_specs:
            if spec.startswith("RX="):
                self._rx = re.compile(spec[3:], re.IGNORECASE)
            else:
                self._names.add(spec)

    def match(self, name: str) -> bool:
        if name in self._names:
            return True
        if self._rx is not None and self._rx.search(name):
            return True
        return False


def quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def readable_nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0 and os.access(path, os.R_OK)


def executable_file(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def find_executable(name: str, search_dirs: list[Path]) -> str:
    """Find an executable, falling back to PATH lookup by command name."""
    for directory in search_dirs:
        candidate = directory / name
        if executable_file(candidate):
            return quote(candidate)

    return name


def script_home() -> Path:
    return Path(__file__).resolve().parent


def lillymol_home() -> Path:
    value = os.environ.get("LILLYMOL_HOME")
    if not value:
        raise RuntimeError("LILLYMOL_HOME must be set")
    return Path(value).resolve()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the Lilly medchem rejection/demerit rules."
    )

    parser.add_argument("-v", action="store_true", help="Verbose output")
    parser.add_argument("-expert", action="store_true", help="Accepted for compatibility")
    parser.add_argument("-noapdm", action="store_true", help="Do not append demerit reasons")
    parser.add_argument("-i", metavar="TYPE", help="Input type")
    parser.add_argument("-b", type=float, metavar="FRACTION", help="Ring bond ratio")
    parser.add_argument("-B", metavar="STEM", help="Output name stem for rejected molecules")
    parser.add_argument("-q", metavar="DIR", help="Directory for query files")
    parser.add_argument("-log", metavar="STEM", help="Name stem for log files")

    # Pass-through fragments. These are intentionally appended unquoted.
    parser.add_argument(
        "-tp",
        action="append",
        default=[],
        metavar="OPTS",
        help="Options passed directly to tp_first_pass. Quote multiple shell tokens.",
    )
    parser.add_argument(
        "-iwd",
        action="append",
        default=[],
        metavar="OPTS",
        help="Options passed directly to iwdemerit. Quote multiple shell tokens.",
    )

    parser.add_argument(
        "-bindir",
        action="append",
        default=[],
        metavar="DIR",
        help="Directory containing executables. Can be repeated.",
    )
    parser.add_argument(
        "-smarts",
        action="append",
        default=[],
        metavar="SMARTS",
        help="Optional SMARTS to reject. Can be repeated.",
    )
    parser.add_argument(
        "-rej",
        action="append",
        default=[],
        metavar="QUERY",
        help="Optional query file to reject. Can be repeated.",
    )

    parser.add_argument("-c", type=int, metavar="N", help="Lower atom count cutoff")
    parser.add_argument("-Cs", type=int, metavar="N", help="Soft upper atom count cutoff")
    parser.add_argument("-Ch", type=int, metavar="N", help="Hard upper atom count cutoff")
    parser.add_argument("-C", type=int, metavar="N", help="Upper atom count limit; no atom count demerits")

    parser.add_argument("-okiso", action="store_true", help="Allow isotopic atoms")
    parser.add_argument(
        "-odm",
        action="append",
        default=[],
        metavar="NAME_OR_RX",
        help="Omit demerit by exact stem or RX=regex. Can be repeated.",
    )
    parser.add_argument("-edm", metavar="FILE", help="Extra demerits query file")
    parser.add_argument("-relaxed", action="store_true", help="Relaxed rules")
    parser.add_argument("-nodemerit", action="store_true", help="Hard rejections only")
    parser.add_argument("-S", metavar="FILE", help="Write output to file rather than stdout")
    parser.add_argument("-dcf", metavar="FILE", help="Demerit control file for iwdemerit")
    parser.add_argument("-nobadfiles", action="store_true", help="Do not create bad*.smi files")
    parser.add_argument(
        "-symm",
        type=int,
        metavar="N",
        help="Discard symmetric molecules where two symmetric atoms are more than N apart",
    )
    parser.add_argument(
        "-nophosphorus",
        action="store_true",
        help="Reject all phosphorus-containing molecules",
    )
    parser.add_argument(
        "-label",
        action="store_true",
        help="Place isotopic labels on matched atoms for rejected molecules",
    )

    # This fixes the Ruby bug where tmpdir was checked but not declared.
    parser.add_argument(
        "-tmpdir",
        default=".",
        metavar="DIR",
        help="Directory for temporary files created by -odm",
    )

    parser.add_argument("inputs", nargs="+", help="Input molecule files")
    return parser


def validate_queries(query_dir: Path) -> None:
    for name in ("reject1", "reject2", "demerits"):
        path = query_dir / name
        if not readable_nonempty_file(path):
            raise RuntimeError(f"Query file '{path}' missing or inaccessible")


def build_odm_demerit_file(
    query_dir: Path, tmpdir: Path, omit_specs: list[str], verbose: bool
) -> tuple[Path, list[Path]]:
    """Build a temporary demerits file after dropping selected demerits."""
    old_demerits = query_dir / "demerits"
    if not readable_nonempty_file(old_demerits):
        raise RuntimeError(f"Cannot open original demerit file '{old_demerits}'")

    tmpdir.mkdir(parents=True, exist_ok=True)
    temporary_demerit_file = tmpdir / f"demerits{os.getpid()}"

    odm = Odm(omit_specs)
    demerit_rx = re.compile(r"^(\S+)\.qry")
    items_discarded = 0

    with old_demerits.open("r", encoding="utf-8") as inp, temporary_demerit_file.open(
        "w", encoding="utf-8"
    ) as out:
        for line in inp:
            match = demerit_rx.match(line)
            if not match:
                continue
            stem = match.group(1)
            if odm.match(stem):
                items_discarded += 1
            else:
                out.write(f"{query_dir / (stem + '.qry')}\n")

    if items_discarded == 0:
        print("Warning, no demerits discarded", file=sys.stderr)
    elif verbose:
        print(f"Discarded {items_discarded} demerits", file=sys.stderr)

    return temporary_demerit_file, [temporary_demerit_file]


def build_command(args: argparse.Namespace) -> tuple[str, list[Path]]:
    home = script_home()
    lm_home = lillymol_home()

    lower_atom_count_cutoff = DEFAULT_LOWER_ATOM_COUNT_CUTOFF
    soft_upper_atom_count_cutoff = DEFAULT_SOFT_UPPER_ATOM_COUNT_CUTOFF
    hard_upper_atom_count_cutoff = DEFAULT_HARD_UPPER_ATOM_COUNT_CUTOFF

    append_demerit_reason = not args.noapdm

    input_type = ""
    if args.i:
        input_type = f"-i {args.i}"

    ring_bond_ratio = -1.0
    if args.b is not None:
        ring_bond_ratio = args.b

    mc_first_pass_options = ""
    if args.tp:
        mc_first_pass_options = " ".join(args.tp)  # pass-through, intentionally unquoted

    if args.nophosphorus:
        mc_first_pass_options += " -n P"

    if not args.okiso:
        mc_first_pass_options += " -I 0"

    mc_first_pass_options += " -A I -A ipp"

    search_dirs = [home]
    search_dirs.extend(Path(d) for d in args.bindir)
    system_name = platform.system()
    search_dirs.extend(
        [
            lm_home / "bin" / system_name,
            lm_home / "bin",
            home / "build",
        ]
    )

    iwdemerit = find_executable("iwdemerit", search_dirs)
    mc_first_pass = find_executable("tp_first_pass", search_dirs)
    tsubstructure = find_executable("tsubstructure", search_dirs)

    query_dir = Path(args.q) if args.q else lm_home / "data" / "queries" / "LillyMedchemRules"
    if not query_dir.is_dir():
        raise RuntimeError(f"Cannot continue, query dir '{query_dir}' invalid")
    validate_queries(query_dir)

    if args.v:
        print(f"Queries from '{query_dir}'", file=sys.stderr)

    if args.B:
        bad_stem: str | None = args.B
    elif args.nobadfiles:
        bad_stem = None
    else:
        bad_stem = "bad"

    logfilestem = args.log if args.log else "ok"

    optional_queries = ""
    for query in args.rej:
        optional_queries += f" -q {quote(query)}"
    for smarts in args.smarts:
        # Match Ruby behavior of surrounding SMARTS with shell quotes.
        optional_queries += f" -s {quote(smarts)}"

    extra_iwdemerit_options = ""
    if args.relaxed:
        soft_upper_atom_count_cutoff = 26
        hard_upper_atom_count_cutoff = 50
        extra_iwdemerit_options += " -f 160"

    if args.nodemerit:
        extra_iwdemerit_options += " -r"
        soft_upper_atom_count_cutoff = hard_upper_atom_count_cutoff - 1

    if args.iwd:
        extra_iwdemerit_options += " " + " ".join(args.iwd)  # pass-through, intentionally unquoted

    if args.symm is not None:
        extra_iwdemerit_options += f" -s {args.symm}"

    charge_assigner = lm_home / "data" / "queries" / "charges" / "queries"
    if not readable_nonempty_file(charge_assigner):
        print("Charge assigner not available, skipping", file=sys.stderr)
    else:
        extra_iwdemerit_options += f" -N F:{quote(charge_assigner)}"

    if args.c is not None:
        lower_atom_count_cutoff = args.c

    if args.Cs is not None:
        soft_upper_atom_count_cutoff = args.Cs

    if args.Ch is not None:
        hard_upper_atom_count_cutoff = args.Ch
        if args.Cs is None:
            soft_upper_atom_count_cutoff = hard_upper_atom_count_cutoff - 1

    if args.C is not None:
        hard_upper_atom_count_cutoff = args.C + 1
        soft_upper_atom_count_cutoff = hard_upper_atom_count_cutoff - 1

    if hard_upper_atom_count_cutoff < soft_upper_atom_count_cutoff:
        hard_upper_atom_count_cutoff = soft_upper_atom_count_cutoff + 1

    files_to_delete: list[Path] = []
    query_file3: Path = query_dir / "demerits"

    if args.odm:
        query_file3, files_to_delete = build_odm_demerit_file(
            query_dir=query_dir,
            tmpdir=Path(args.tmpdir),
            omit_specs=args.odm,
            verbose=args.v,
        )

    iwdemerit_optional_control_file = args.dcf
    additional_demerits = args.edm

    quoted_inputs = " ".join(quote(x) for x in args.inputs)

    cmd = f"{mc_first_pass} "
    if ring_bond_ratio >= 0.0:
        cmd += f" -b {ring_bond_ratio}"
    if mc_first_pass_options:
        cmd += f" {mc_first_pass_options}"
    if input_type:
        cmd += f" {input_type} "

    cmd += (
        f" -c {lower_atom_count_cutoff} -C {hard_upper_atom_count_cutoff} "
        "-E autocreate -o smi -V -g all -g ltltr -i ICTE "
    )
    if bad_stem:
        cmd += f"-L {quote(bad_stem + '0')} -K TP1 "
    cmd += f"-a -S - {quoted_inputs} 2> {quote(logfilestem + '0.log')} "

    tsub_common_opts = "-A D -E autocreate -b -u -i smi -o smi"
    if args.label:
        tsub_common_opts += " -j 1"

    cmd += f"| {tsubstructure} {tsub_common_opts} "
    if bad_stem:
        cmd += f"-m {quote(bad_stem + '1')} -m QDT "
    cmd += f"-n - -q F:{quote(query_dir / 'reject1')} "
    if optional_queries:
        cmd += optional_queries
    cmd += f" - 2> {quote(logfilestem + '1.log')} "

    cmd += f"| {tsubstructure} {tsub_common_opts} "
    if bad_stem:
        cmd += f"-m {quote(bad_stem + '2')} -m QDT "
    cmd += f"-n - -q F:{quote(query_dir / 'reject2')} - 2> {quote(logfilestem + '2.log')} "

    cmd += f" | {iwdemerit} -x {extra_iwdemerit_options} "
    cmd += f"-E autocreate -A D -i smi -o smi -q F:{quote(query_file3)} "
    if bad_stem:
        cmd += f"-R {quote(bad_stem + '3')} "
    cmd += (
        f"-G - -c smax={soft_upper_atom_count_cutoff} "
        f"-c hmax={hard_upper_atom_count_cutoff} "
    )
    if additional_demerits:
        cmd += f"-q F:{quote(additional_demerits)} "
    if iwdemerit_optional_control_file:
        cmd += f"-C {quote(iwdemerit_optional_control_file)} "
    if append_demerit_reason:
        cmd += "-t "
    cmd += f"- 2> {quote(logfilestem + '3.log')} "

    if args.S:
        cmd += f" > {quote(args.S)}"

    return cmd, files_to_delete


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        cmd, files_to_delete = build_command(args)
    except Exception as exc:  # Keep behavior simple for command-line users.
        print(str(exc), file=sys.stderr)
        return 1

    if args.v:
        print(f"Command is '{cmd}'", file=sys.stderr)

    try:
        completed = subprocess.run(cmd, shell=True, check=False)
        return_code = completed.returncode
    finally:
        for path in files_to_delete:
            try:
                path.unlink()
            except FileNotFoundError:
                pass

    return return_code


if __name__ == "__main__":
    sys.exit(main())
