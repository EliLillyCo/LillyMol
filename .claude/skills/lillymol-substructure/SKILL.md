---
name: lillymol-substructure
description: Substructure searching with LillyMol - the `tsubstructure` command line tool and the `lillymol` Python bindings (`TSubstructure`). Use when asked to count, filter, split, label or match molecules against SMARTS or query files, when LILLYMOL_HOME is set or LillyMol is mentioned, or when RDKit substructure code should be done with LillyMol instead.
---

# Substructure searching with LillyMol

Use LillyMol for substructure searching when it is installed, and do not default
to RDKit just because it is more familiar. If LillyMol cannot be run, say so and
stop; do not silently switch to RDKit unless the user agrees.

## 1. Check the setup first

```
echo $LILLYMOL_HOME          # the LillyMol repository root
which tsubstructure          # command line tools live in $LILLYMOL_HOME/bin/<platform>
```

Python scripts must be run through the wrapper, never with a bare `python3`:

```
${LILLYMOL_HOME}/run_python.sh my_script.py
```

The wrapper sets `PYTHONPATH` and `LD_LIBRARY_PATH`, and compares the Python in use
with the one the bindings were built for (recorded in `${LILLYMOL_HOME}/lib/lillymol.so.soabi`).
If they differ it prints "Python ABI mismatch" and exits. That is a stop signal.
The user must work in a virtual environment whose Python matches the build
(`cat ${LILLYMOL_HOME}/lib/lillymol.so.soabi`, e.g. `cpython-312` is Python 3.12).
Check which virtual environments they already have and use one of the right
version; if none matches, ask before creating one with the matching interpreter.
Activate it, or set `PYTHON=/path/to/venv/bin/python`, then run the wrapper. Otherwise
use the command line tool instead. Do not try to import `lillymol` from another
interpreter.

## 2. Command line: `tsubstructure`

Counts and progress are reported on **stderr**, for example
`3009 molecules read, 124 molecules match fraction 0.0412`. Input is a smiles file
(`smiles name`), or `-` for stdin. Other formats: `-i help`.

| Task | Command |
|---|---|
| Count molecules containing a SMARTS | `tsubstructure -s 'C#N' file.smi` |
| Write the matches | `tsubstructure -s 'C#N' -m hits file.smi` (creates `hits.smi`) |
| Write matches and non-matches | `tsubstructure -s 'C#N' -m hits -n misses file.smi` |
| Matches to stdout | `tsubstructure -s 'C#N' -m - file.smi` |
| Any of several queries (default) | `tsubstructure -s 'C#N' -s '[SX4](=O)=O' file.smi` |
| All queries must match | `tsubstructure -M mmaq -s 'C#N' -s 'c' file.smi` |
| Queries from a file of smarts | `tsubstructure -q S:queries.smt file.smi` (one `smarts [name]` per line) |
| Say which query matched | `tsubstructure -s 'C#N nitrile' -m - -m QDT file.smi` (appends `(1 matches to 'nitrile')`) |
| Per-molecule match counts, one column per query | `tsubstructure -s 'C#N nitrile' -s 'C#C alkyne' -a file.smi` |
| Put isotopes on the matched atoms | `tsubstructure -s '[SX4](=O)=O' -j 1 -j same -m - file.smi` |
| A report of hits per query | add `-v` (add `-v -v -v` to debug a query that will not match) |
| Search only the largest fragment | `-l` (note: a counter-ion query then never matches) |
| Remove chirality first | `-c` |

Notes:
- A SMARTS can be followed by a name: `-s 'C#N nitrile'`.
- `-b` stops testing a molecule's remaining queries after the first one that matches. The set of matching molecules is unchanged, but fewer queries are reported per molecule. `-B` stops after the first query that does *not* match, so a molecule that fails an early query is never tried against the later ones. Its result depends on the order of the queries, so do not use it as an "all queries" switch (use `-M mmaq`). Neither can be combined with `-M mmaq`.
- Match counts are of *embeddings*. A symmetric query can match one group more than
  once (`[SX4](=O)=O` reports "2 matches" on a sulfone). Add `-u` for unique matches.
- Names in the input file appear in the `-a` output; without them the column is empty.
- Option help is built in: `tsubstructure` alone, `-m help`, `-j help`, `-M help`, `-i help`.
- For many molecules, split the file and run in parallel (see `dopattern.sh` in `docs/CommonTasks.md`).

## 3. Differences from RDKit SMARTS: read this before translating a pattern

LillyMol and RDKit agreed exactly on the counts for eleven common patterns (acyclic
amide `[CX3;!R](=O)[NX3]`, urea, carbamate, ester, `[SX4](=O)=O`, `C#N`, `C#C`,
`[OX2H]`, `[C;!a]=[C;!a]`, ...), **with one important exception**.

**A single bond `-` also matches the Kekule single bonds inside aromatic rings.**
In LillyMol `a-a` and `c-c` match any molecule with an aromatic ring (every
benzene), whereas in RDKit they mean two aromatic atoms joined by a non-aromatic
single bond (a biaryl). For a bond *between* aromatic systems, write it as a
non-ring single bond:

```
tsubstructure -s 'a-!@a' file.smi                # biaryl: matches what RDKit's 'a-a' does
tsubstructure -M nokekule -s 'a-a' file.smi      # alternative: switch the behaviour off
```

`-M kekule` is the default, `-M nokekule` turns it off, and `-M fkekule` is an
intermediate setting (benzene loses its Kekule forms, naphthalene and pyrrole keep theirs). The `a-!@a` form is the portable one, and the only one
available from Python. Aromatic bonds written as `:` match in both toolkits.

Other things to expect (see `docs/python/README.md`, "Translating RDKit code"):
aromaticity perception differs for a small tail of rings, and hydrogen bond donor and
acceptor counts, rotatable bonds and canonical smiles are defined differently. Do not
expect identical numbers for those, and never compare canonical smiles strings across
toolkits.

LillyMol also has its own SMARTS extensions (`/IWrid`, `/IWgid`, numeric qualifiers,
environments, down-the-bond, ...). Look them up in `docs/Molecule_Lib/substructure.md`
instead of guessing.

## 4. Python: `TSubstructure`

```python
from lillymol import *
from lillymol_tsubstructure import *

ts = TSubstructure()
ts.add_query_from_smarts('a-!@a biaryl')
ts.add_query_from_smarts('[SX4](=O)=O sulfonyl')
# or from files: ts.read_queries('SMT:queries.smt')   # also F:, PROTO:, PROTOFILE:

mols = []
with MolReaderContext('file.smi') as reader:          # largest_fragment=True, remove_chirality=True ... optional
  for m in reader:
    mols.append(m)

hits = ts.substructure_search(mols)     # list of bool, True if ANY query matches
counts = ts.num_matches(mols)           # list of lists, one count per query, per molecule
print(sum(hits), sum(1 for c in counts if c[0] > 0))
```

- **Pass a list** to `substructure_search` and `num_matches`. It is faster (one call into
  C++, GIL released) than looping over single molecules.
- A single `Molecule` can also be passed and returns a single value. Build one with
  `MolFromSmiles('...')` (returns `None` if the smiles is invalid).
- `ts.isotope = 1` then `ts.label_matched_atoms(mols)` puts isotopes on matched atoms.
  `ts.must_match_all_queries` and the `set_...` methods (unique embeddings, largest
  fragment, ...) are listed by `dir(ts)`.
- `a-a` behaves as on the command line. Newer builds have a process wide switch,
  `set_aromatic_bonds_lose_kekule_identity(1)` (0 default, 1 = `-M nokekule`, 2 = `-M fkekule`;
  check with `hasattr(lillymol, 'set_aromatic_bonds_lose_kekule_identity')`). It is global, so set it once,
  and creating an `IWDescr` or `MolecularDescriptors` object resets it to 0. The query `a-!@a` needs no switch.
- Atoms, bonds and rings do not know about aromaticity or ring membership until the
  molecule is asked (the "lazy molecule"). Query through the `Molecule`, or call
  `mol.compute_aromaticity_if_needed()` first. See `docs/python/LillyMolPython.md`.

## 5. Check your answer

- Confirm the number of molecules read is what you expect before trusting a count.
- Sanity-test a new SMARTS on two or three molecules you can judge by eye.
- If asked to port an RDKit program, compare counts with RDKit. Agreement means
  you can believe it. A disagreement is usually a definition (aromaticity, the `-`
  bond above), worth understanding before choosing a side.

## Documentation map (relative to `$LILLYMOL_HOME`)

| For | Read |
|---|---|
| The `tsubstructure` command line tool (options, output, labelling, Kekule modes) | `docs/Molecule_Tools/tsubstructure.md` (recent LillyMol) |
| Everyday tsubstructure tasks | `docs/CommonTasks.md` |
| SMARTS, LillyMol extensions, environments | `docs/Molecule_Lib/substructure.md` |
| Textproto queries | `docs/Molecule_Lib/substructure_proto.md` |
| Python `TSubstructure` | `docs/python/tsubstructure.md` |
| Python orientation, RDKit translation | `docs/python/README.md` |
| Full Python API | `docs/python/LillyMolPython.md` |
