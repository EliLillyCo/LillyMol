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

Unfortunately, multiple molecules cannot be processed this way, due to some
complexities in how lists and vectors are passed between python and c++, but
obviously molecules can be processed in a loop or list comprehension.

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

## Performance
As a test of performance, load 70 queries from the Lilly Medchem Rules
```
ts = TSubstructure()
ts.read_queries("F:/home/you/path/to/LillyMol/data/queries/medchemrules/reject1")
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
runs in about 0.44 seconds on SandyBridge hardware released in 2012. Running the same
queries through command line `tsubstructure` takes about 0.26 seconds.
Running the vectorized form of the query
```
mols = [MolFromSmiles(smi) for smi in smiles]
result = ts.substructure_search(mols)
matches = op.countOf(result, True)
```
runs in about 0.40 seconds, indicating that minimising the number of times the c++/python barrier
is crossed may be beneficial. Further gains are possible if smiles interpretation is done
entirely within c++ by passing a list of smiles to `TSubstructure`.
```
smiles = [ 2000 smiles ]
result = ts.substructure_search(smiles)
matched = op.countOf(result, True)
print(f"matched {matched} queries")
```
which runs in 0.30 seconds. But this has the disadvantage of never having a Molecule
object available within python. Which might be fine.

Yet another variant which relies on a bulk conversion from smiles to Molecule
```
smiles = [ 2000 smiles ]
mols = MolFromSmiles(smiles)
result = ts.substructure_search(mols)
matched = op.countOf(result, True)
print(f"matched {matched} queries")
```
runs in 0.36 seconds.

As you contemplate whether to use the python `tsubstructure`
implementation or the command line version will depend on the size of the dataset being
processed.


## Summary
The python implementation of `tsubstructure` functionality provides many of the
most commonly used features of command line `tsubstructure`, with modest
performance trade-offs.
