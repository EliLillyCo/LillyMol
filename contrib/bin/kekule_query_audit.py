#!/usr/bin/env python3
"""Find substructure queries whose answers depend on the Kekule matching mode.

By default LillyMol keeps the Kekule single and double bonds of an aromatic ring,
and a query bond can match an aromatic bond through its Kekule form. With
`tsubstructure -M nokekule` aromatic bonds match only aromatic query bonds. A query
that gives different answers in the two modes depends on the setting, either
deliberately or by accident, and it changes behaviour if the setting changes.

This script runs every query in a set of query files over a set of molecules in the
default mode and with -M nokekule (and optionally -M fkekule), and reports the queries
that differ. See docs/Molecule_Tools/tsubstructure.md, "Aromatic bonds and Kekule forms".

Examples

  # counts only, a few thousand molecules, every query under data/queries
  kekule_query_audit.sh --head 5000 molecules.smi

  # a systematic sample of a big file, and compare the labelled atoms, not just the counts
  kekule_query_audit.sh --every 100 --labels --jobs 16 all.smi > audit.tsv

  # only some queries, and say which query lists use each one
  kekule_query_audit.sh --used-in molecules.smi data/queries/hbonds data/queries/charges

Two levels of comparison

  counts   (the default) the number of molecules that match each query.
  --labels also puts isotopes on the matched atoms and compares the labelled output of
           every molecule. This is stricter, and slower. It can find a query that
           matches the same molecules but different atoms, which matters for queries
           whose matched atoms are used, such as the hydrogen bond acceptor queries.

A query is reported as

  same        identical in every mode tested.
  COUNT_DIFF  a different number of molecules match.
  LABEL_DIFF  the same number match, but the labelled atoms differ (with --labels).
  ERROR       tsubstructure could not use the query, in any mode. These are reported
              because a query that cannot be read cannot be checked.

Differences are not necessarily bugs. Some queries depend on the default deliberately,
for example hbonds/imine.qry originally found aromatic ring nitrogens through their
Kekule double bonds. A query that you have written to be independent of the mode should
come out as 'same', and --fail-on-change makes the exit status 1 if any query does not.

Only the standard library is used, and tsubstructure is run as a subprocess, so no
LillyMol python bindings are needed.
"""

import argparse
import concurrent.futures
import hashlib
import itertools
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

MATCH_RE = re.compile(r"(\d+) molecules read, (\d+) molecules match")

# Mode name and the tsubstructure arguments that select it. The default mode must be first.
MODES = [("default", []), ("nokekule", ["-M", "nokekule"]), ("fkekule", ["-M", "fkekule"])]

# How a query file is given to tsubstructure, by suffix.
QUERY_SUFFIXES = {".qry": "", ".txtproto": "PROTO:", ".textproto": "PROTO:"}

# Suffixes looked for when a directory is searched. Other .textproto files, such as
# data/queries/pharmacophore/pharmacophore.textproto, are not substructure queries.
# A .textproto file is used if it is named on the command line.
DIRECTORY_SUFFIXES = {".qry", ".txtproto"}


def find_tsubstructure(explicit):
  if explicit:
    return explicit
  found = shutil.which("tsubstructure")
  if found:
    return found
  home = os.environ.get("LILLYMOL_HOME")
  if home:
    candidate = Path(home) / "bin" / platform.system() / "tsubstructure"
    if candidate.exists():
      return str(candidate)
  return None


def find_queries(paths):
  """All query files in the files and directories given, in a stable order.

  A directory is searched for .qry and .txtproto files. A file named on the command
  line is used whatever its suffix.
  """
  queries = []
  for p in paths:
    p = Path(p)
    if p.is_dir():
      for root, dirs, files in os.walk(p):
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        for f in sorted(files):
          if Path(f).suffix in DIRECTORY_SUFFIXES and not f.startswith("."):
            queries.append(Path(root) / f)
    elif p.is_file():
      queries.append(p)
    else:
      sys.exit(f"{p}: no such file or directory")
  return queries


def sample_molecules(source, dest, head, every):
  """Write a subset of the molecules in source to dest, returning how many."""
  n = 0
  with open(source, errors="replace") as src, open(dest, "w") as out:
    lines = itertools.islice(src, head) if head else src
    for i, line in enumerate(lines):
      if i % every == 0:
        out.write(line)
        n += 1
  return n


def run_one(tsubstructure, query, mode_args, molecules, labels, timeout):
  """(matches, digest, error). The digest is of the labelled output, or None."""
  qarg = QUERY_SUFFIXES.get(query.suffix, "") + str(query.resolve())
  cmd = [tsubstructure] + mode_args + ["-q", qarg]
  if labels:
    cmd += ["-j", "1", "-j", "same", "-m", "-", "-m", "NONMX"]
  cmd.append(str(molecules))
  # stderr goes to a file. A query that cannot be read can write a great deal there, and
  # a full pipe would stop tsubstructure while we are still waiting on stdout.
  with tempfile.TemporaryFile() as errfile:
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errfile,
                            cwd=query.resolve().parent)
    digest = hashlib.md5()
    try:
      for chunk in iter(lambda: proc.stdout.read(1 << 20), b""):
        digest.update(chunk)
      proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
      proc.kill()
      proc.wait()
      return None, None, "timeout"
    errfile.seek(0)
    err = errfile.read().decode(errors="replace")
  m = MATCH_RE.search(err)
  if not m:
    # Exit status is not informative for tsubstructure, so the report line decides.
    first = next((ln for ln in err.splitlines() if ln.strip()), "no output")
    return None, None, first[:120]
  return int(m.group(2)), (digest.hexdigest() if labels else None), None


def audit_query(args, tsubstructure, molecules, query):
  modes = MODES if args.fkekule else MODES[:2]
  results = {}
  for name, margs in modes:
    results[name] = run_one(tsubstructure, query, margs, molecules, args.labels, args.timeout)

  errors = [r[2] for r in results.values() if r[2]]
  row = {"query": query, "matches": {n: r[0] for n, r in results.items()}}
  if errors:
    row["status"] = "ERROR"
    row["note"] = errors[0]
    return row

  base_matches, base_digest, _ = results["default"]
  status = "same"
  for name, (matches, digest, _) in results.items():
    if name == "default":
      continue
    if matches != base_matches:
      status = "COUNT_DIFF"
    elif args.labels and digest != base_digest and status == "same":
      status = "LABEL_DIFF"
  row["status"] = status
  if args.labels:
    row["labels_same"] = {n: (r[1] == base_digest) for n, r in results.items() if n != "default"}
  return row


def used_in(query):
  """Names of the list files, files with no extension, that mention the query.

  The queries that LillyMol tools read are usually listed in such files, for example
  data/queries/charges/positive. Only the query's own directory is searched.
  """
  names = {query.name, query.name.split(".")[0] + ".qry"}
  found = []
  for f in sorted(query.parent.iterdir()):
    if f.is_file() and f.suffix == "" and not f.name.startswith("."):
      try:
        tokens = f.read_text(errors="replace").split()
      except OSError:
        continue
      if any(Path(t).name in names for t in tokens):
        found.append(f.name)
  return found


def main():
  parser = argparse.ArgumentParser(
      description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("molecules", help="smiles file to search")
  parser.add_argument("queries", nargs="*",
                      help="query files (.qry, .txtproto) or directories. Default: ${LILLYMOL_HOME}/data/queries")
  parser.add_argument("--head", type=int, default=0, help="only use the first N molecules")
  parser.add_argument("--every", type=int, default=1,
                      help="use every Nth molecule, a systematic sample (default 1, all)")
  parser.add_argument("--labels", action="store_true",
                      help="also compare the labelled atoms, not just the number of matches")
  parser.add_argument("--fkekule", action="store_true", help="also test -M fkekule")
  parser.add_argument("--all", action="store_true", help="list every query, not just those that differ")
  parser.add_argument("--used-in", action="store_true", help="add the query lists that name each query")
  parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="parallel tsubstructure runs")
  parser.add_argument("--timeout", type=int, default=3600, help="seconds allowed for each run")
  parser.add_argument("--tsubstructure", help="path to tsubstructure (default: PATH, then LILLYMOL_HOME)")
  parser.add_argument("--fail-on-change", action="store_true",
                      help="exit with status 1 if any query depends on the mode")
  args = parser.parse_args()

  if args.every < 1:
    sys.exit("--every must be at least 1")

  tsubstructure = find_tsubstructure(args.tsubstructure)
  if not tsubstructure:
    sys.exit("Cannot find tsubstructure. Put it on PATH, set LILLYMOL_HOME, or use --tsubstructure")

  query_paths = args.queries
  if not query_paths:
    home = os.environ.get("LILLYMOL_HOME")
    if not home:
      sys.exit("No queries given and LILLYMOL_HOME is not set")
    query_paths = [Path(home) / "data" / "queries"]
  queries = find_queries(query_paths)
  if not queries:
    sys.exit("No query files found")

  with tempfile.TemporaryDirectory(prefix="kekule_audit_") as tmp:
    sample = Path(tmp) / "molecules.smi"
    n = sample_molecules(args.molecules, sample, args.head, args.every)
    if n == 0:
      sys.exit("No molecules to search")
    modes = [m[0] for m in (MODES if args.fkekule else MODES[:2])]
    print(f"{len(queries)} queries, {n} molecules, modes {', '.join(modes)}, "
          f"{'labelled atoms compared' if args.labels else 'counts only'}", file=sys.stderr)

    rows = []
    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
      futures = [pool.submit(audit_query, args, tsubstructure, sample, q) for q in queries]
      for fut in concurrent.futures.as_completed(futures):
        rows.append(fut.result())
        done += 1
        if done % 25 == 0 or done == len(queries):
          print(f"  {done} of {len(queries)} queries done", file=sys.stderr)

  rows.sort(key=lambda r: str(r["query"]))
  differing = [r for r in rows if r["status"] in ("COUNT_DIFF", "LABEL_DIFF")]
  errors = [r for r in rows if r["status"] == "ERROR"]

  header = ["query", "status"] + [f"matches_{m}" for m in modes]
  header += ["change_nokekule"]
  if args.labels:
    header += ["labels_same_" + m for m in modes[1:]]
  if args.used_in:
    header += ["used_in"]
  header += ["note"]
  print("\t".join(header))
  for r in rows:
    if r["status"] == "same" and not args.all:
      continue
    m = r["matches"]
    cells = [str(r["query"]), r["status"]] + [str(m.get(x, "")) if m.get(x) is not None else "" for x in modes]
    change = (m["nokekule"] - m["default"]) if m.get("nokekule") is not None and m.get("default") is not None else ""
    cells.append(str(change))
    if args.labels:
      cells += [str(r.get("labels_same", {}).get(x, "")) for x in modes[1:]]
    if args.used_in:
      cells.append(",".join(used_in(Path(r["query"]))))
    cells.append(r.get("note", ""))
    print("\t".join(cells))

  print(f"\n{len(rows)} queries: {len(rows) - len(differing) - len(errors)} same, "
        f"{len(differing)} depend on the mode, {len(errors)} could not be run", file=sys.stderr)
  return 1 if (args.fail_on_change and differing) else 0


if __name__ == "__main__":
  sys.exit(main())
