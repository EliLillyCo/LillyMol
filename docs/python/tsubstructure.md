# tsubstructure (python)

## tsubstructure
The c++ substructure searching tool `tsubstructure` is one of the most commonly
used LillyMol utilities. Due to the accumulated needs of many years of use it
is quite comprehensive and complex.

The python implementation of `tsubstructure` is designed to capture some
of the most common use cases into a python interface.

### tsubstructure
`tsubstructure` takes two sets of inputs

1. One of more queries
2. A set of molecules to be searched

The most common task is to separate the molecules into those that match any
of the queries from those that do not match any of the queries, but many other
uses are possible.

A typical python usage might be
```
from lillymol import *
from lillymol_tsubstructure import *

ts = TSubstructure()
ts.add_query_from_smarts("C carbon")

mol = MolFromSmiles("C1CC1")
ts.substructure_search(mol)
```
returns `True`, since the query matches the molecule.

A `TSubstructure` object can hold any number of queries
```
from lillymol import *
from lillymol_tsubstructure import *

ts = TSubstructure()
ts.add_query_from_smarts("C carbon")
ts.add_query_from_smarts("c aromatic carbon")
ts.add_query_from_smarts("N nitrogen")
ts.add_query_from_smarts("n aromatic nitrogen")
ts.add_query_from_smarts("F Fluorine")

mol = MolFromSmiles("Cc1ncc(N)cc1")
ts.substructure_search(mol)
```
again returns `True` since at least one of the queries matches the input.

If you need more information about which of the queries (if any) have matched
the input, try
```
from lillymol import *
from lillymol_tsubstructure import *

ts = TSubstructure()
ts.add_query_from_smarts("C carbon")
ts.add_query_from_smarts("c aromatic carbon")
ts.add_query_from_smarts("N nitrogen")
ts.add_query_from_smarts("n aromatic nitrogen")
ts.add_query_from_smarts("F Fluorine")

mol = MolFromSmiles("Cc1ncc(N)cc1")
ts.num_matches(mol)
```
returns `[1, 5, 1, 1, 0]`. The first query, aliphatic Carbon, matched one atom in the structure,
the second query, aromatic Carbon, matched 5 atoms...

If you have multiple molecules, `substructure_search` will return a list of Booleans.
```
from lillymol import *
from lillymol_tsubstructure import *

ts = TSubstructure()
ts.add_query_from_smarts("C carbon")
ts.add_query_from_smarts("N nitrogen")
ts.add_query_from_smarts("O oxygen")

smiles = ["C", "N", "O", "[U]", "CN"]
mols = [MolFromSmiles(smi) for smi in smiles]
ts.substructure_search(mols)
```
returns `[True, True, True, False, True]`, one list entry per molecule.

### Query Files
All the examples shown here show queries being loaded as smarts. The method `read_queries`
will read a file of queries just the same mechanisms as as the `-q` option to `tsubstructure`.
For example
```
ts = TSubstructure()
ts.read_query("file.qry")          # read an old style query file
ts.read_query("F:file")            # read a file containing names of old style query files
ts.read_query("SMT:file")          # read a file of smarts - one per line
ts.read_query("PROTO:file.qry")    # read a single textproto query
ts.read_query("PROTOFILE:file")    # read a file containing names of textproto query files
```
These may be the most convenient way of loading queries into a `TSubstructure` object.

## num_matches
Sometimes knowing whether or not a molecule matches any of the queries is good enough.
But sometimes it may be important to know which of the queries have matched a given molecule.
```
from lillymol import *
from lillymol_tsubstructure import *

ts = TSubstructure()
ts.add_queries_from_smarts(["C carbon", "N nitrogen", "O oxygen"])

smiles = ["CCC(N)N"]
mols = [MolFromSmiles(smi) for smi in smiles]
print(ts.num_matches(mols))
```
returns `[[3, 2, 0]]`. We see that the input molecule had 3 matches
to Carbon, 2 matches to the Nitrogen query and no matches to the Oxygen
query. There is also a signature that accepts a single molecule as argument,
in which case the value returned is just a list, rather than a list of lists.

Prefer the list form when you have more than a handful of molecules. It does the
whole loop inside one c++ call, and releases the GIL while it does, so python
threads searching concurrently actually overlap. See [Performance](#performance).

## Label Matched Atoms
A common use case for `tsubstructure` is to place isotopic labels on the matched atoms.
The most common use cases are enabled via python.

To place an isotopic label on all atoms matched by any of the queries
```
from lillymol import *
from lillymol_tsubstructure import *

ts = TSubstructure()
ts.add_query_from_smarts("C carbon")
ts.add_query_from_smarts("F-c Fluorine")
ts.add_query_from_smarts("[U] Uranium")
ts.isotope = 1

mol = MolFromSmiles("Cc1ncc(F)cc1")
nhits = ts.label_matched_atoms(mol)
```
`nhits` will be 2, since two of the queries match the input. `mol` will have
isotopes applied and the unique smiles will be `[1F][1c]1c[n]c([1CH3])cc1`.

A list of molecules can be labelled in one call, which returns the number of
molecules where something matched
```
ts = TSubstructure()
ts.add_query_from_smarts("c1ccccc1 benzene")
ts.isotope = 5

mols = [MolFromSmiles(smi) for smi in ["c1ccccc1C", "CCO", "c1ccccc1N"]]
print(ts.label_matched_atoms(mols))          # 2
print([mol.smiles() for mol in mols])
```
prints `2`, then
`['[5CH]1=[5CH][5CH]=[5CH][5CH]=[5C]1C', 'CCO', '[5CH]1=[5CH][5CH]=[5CH][5CH]=[5C]1N']`.
The labels are applied to the molecules you passed in - the ethanol, which does not
match, is untouched.

This did not used to work, and the reason is worth knowing. The binding took a
`std::vector<Molecule>`, which makes the binding layer copy every molecule in
the list at the boundary, so the labels went onto the copies and were discarded -
there was no way to get the result back. Taking a `std::vector<Molecule*>`
instead copies nothing, and the same change made the batch `substructure_search`
and `num_matches`
substantially faster. See [Performance](#performance).

The other kind of matched atom is where the atoms are labelled by the number
of the query that matches that atom.
```
from lillymol import *
from lillymol_tsubstructure import *

ts = TSubstructure()
ts.add_query_from_smarts("C(=O)-[OH] acid")
ts.add_query_from_smarts("O=N=O Nitro")
ts.add_query_from_smarts("[U] Uranium")
ts.add_query_from_smarts("[OD2]-[CD1] methoxy")
ts.set_label_by_query_number(True)

mol = MolFromSmiles("OC(=O)c1ccc(N(=O)=O)cc1")
ts.label_matched_atoms(mol)
print(mol.unique_smiles())
```
prints `[1OH][1C](=[1O])c1ccc([2N](=[2O])=[2O])cc1` which can be a useful
way of identifying functional groups. Note however that `TSubstructure` does
not support any concept of exclusive matching - use `substituent_model` for
that.

## Aromatic Bonds and Kekule Forms
Inside an aromatic ring, LillyMol keeps the underlying Kekule bonds as well as the
aromatic perception. By default a bond in a query can therefore match an aromatic
bond in either way. The consequence that surprises people is that a single bond
in a query matches the Kekule single bonds *inside* aromatic rings.
```
from lillymol import *

ts = TSubstructure()
ts.add_query_from_smarts("a-a")
benzene = MolFromSmiles("c1ccccc1 benzene")
biphenyl = MolFromSmiles("c1ccccc1-c1ccccc1 biphenyl")

print(ts.substructure_search([benzene, biphenyl]))     # [True, True]
```
`a-a` matches benzene, because its ring bonds have a Kekule single form. In many
other toolkits, RDKit among them, `a-a` would only match two aromatic atoms joined by a
single bond that is not part of an aromatic ring, which is what is usually meant. That
meaning is available by turning the Kekule matching off
```
set_aromatic_bonds_lose_kekule_identity(1)
print(ts.substructure_search([benzene, biphenyl]))     # [False, True]
```
The function takes one of three values

| Mode | Behaviour |
| ---- | --------- |
| `0` | The default. An aromatic bond matches aromatic query bonds, and also single and double query bonds via its Kekule form |
| `1` | An aromatic bond matches only aromatic query bonds. `a-a` is two aromatic atoms joined by a non aromatic bond |
| `2` | Like `1`, but some aromatic rings keep their Kekule forms (see below) |

and `aromatic_bonds_lose_kekule_identity()` returns the current value. Any other
value raises `ValueError`. These are the same three behaviours as `-M kekule`,
`-M nokekule` and `-M fkekule` in [tsubstructure](../Molecule_Tools/tsubstructure.md).

Mode `2` is more intricate. The rings it treats as having no Kekule form are decided
by the library and not by the query. Searching `a-a` against a set of molecules
```
mols = [MolFromSmiles(smi) for smi in ["c1ccccc1 benzene", "c1ccc2ccccc2c1 naphthalene",
        "c1cc[nH]c1 pyrrole", "c1ccoc1 furan", "O=c1cccc[nH]1 pyridone",
        "c1ccccc1-c1ccccc1 biphenyl", "c1ccc(cc1)-c1ccc[nH]1 phenylpyrrole"]]
ts = TSubstructure()
ts.add_query_from_smarts("a-a")
for mode in (0, 1, 2):
  set_aromatic_bonds_lose_kekule_identity(mode)
  hits = ts.substructure_search(mols)
  print(mode, [mol.name() for mol, hit in zip(mols, hits) if hit])
```
prints
```
0 ['benzene', 'naphthalene', 'pyrrole', 'furan', 'pyridone', 'biphenyl', 'phenylpyrrole']
1 ['biphenyl', 'phenylpyrrole']
2 ['naphthalene', 'pyrrole', 'furan', 'pyridone', 'biphenyl', 'phenylpyrrole']
```
so in mode `2` only benzene loses its Kekule single bonds. If the goal is to find a bond
between two aromatic systems, use mode `1`, or write the query so that it does not depend
on the mode: `a-!@a` (a single bond that is not in a ring) gives the same answer in
all three modes.

Three things to know, because the setting is not an attribute of a `TSubstructure`

* **It is process wide.** Every `TSubstructure` in the process, in every thread, uses
  the current value. Set it once, before starting any threads that search.
* **It is applied when a search is done**, and when queries are built from molecules, so
  queries that were added before the call are affected by it. Set it before searching.
* **Creating an `IWDescr` or `MolecularDescriptors` object resets it to `0`.** If
  you use descriptors and substructure searching in the same program, set the mode
  after the descriptor objects have been created.

Queries that depend on the mode are not unusual in older query files. If you change the
mode, check the queries you rely on, because they may have been written with the
default in mind.

## Performance
As a test of performance, load 70 queries from the Lilly Medchem Rules
```
ts = TSubstructure()
ts.read_queries("F:/home/you/path/to/LillyMol/data/queries/LillyMedchemRules/reject1")
logging.info("Read %d queries", ts.number_queries())
```
and then run these queries against 2000 random Chembl molecules
```
matches = 0
for smi in smiles:
  mol = MolFromSmiles(smi)
  if ts.substructure_search(mol):
    matches += 1

print(f"Matched {matches}")
```
There are four ways to run that, and they do not perform alike. Measured with 74
queries against 2000 molecules, best of five runs on a modern 32 core machine,
built `-c opt`:

| | | |
| --- | --- | --- |
| loop over single molecules, as above | 0.136 s | 1.0x |
| **a list of Molecule** | **0.072 s** | **1.9x** |
| a list of smiles | 0.099 s | 1.4x |
| bulk `MolFromSmiles` then a list of Molecule | 0.116 s | 1.2x |

```
# the fastest of the four
mols = [MolFromSmiles(smi) for smi in smiles]
result = ts.substructure_search(mols)
matches = op.countOf(result, True)
```
Passing the list does the whole loop inside one c++ call, crossing the boundary
once instead of once per molecule, and releases the GIL while it runs so that
python threads searching concurrently genuinely overlap.

Passing a list of smiles builds the molecules inside c++, which avoids creating
python Molecule objects at all, at the cost of never having one available. It used
to be the fastest of the four; it no longer is, and it re-parses each smiles.

An older version of this document reported these as 0.44, 0.40, 0.30 and 0.36
seconds on 2012 SandyBridge hardware, when the list forms took a
`std::vector<Molecule>` and the binding layer therefore copy constructed every
molecule at the boundary. A Molecule copy is about 5.4 microseconds,
comparable to a whole search, so the copying dominated and the ordering came out
differently. The
bindings now take `std::vector<Molecule*>`, which copies nothing. Do not restore
the by value form.

Whether to use the python implementation or the command line version depends on
the size of the dataset being processed.


## Summary
The python implementation of `tsubstructure` functionality provides many of the
most commonly used features of command line `tsubstructure`, with modest
performance trade-offs.
