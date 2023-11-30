# Synthetic Precedent

## Background
The idea of using fingerprints computed across collections as an indicator
of the likely synthetic feasibility of a new molecule was developed independently
many years by both Ertl (Novartis) and Watson (Lilly). The LillyMol implementation
uses a BerkeleyDB database to store bits. While the current implementation is for
Extended Connectivity (Morgan) type fingerprints, there is no requirement that
this kind of fingerprint be used. [smi2rings](/docs/Molecule_Tools/smi2rings.md)
is a tool that implements this idea for rings, rather than fingerprints.

In practice we find that using both fingerprints and rings to assess precedent
for a new molecule are very effective at identifying molecules likely to
be discarded by chemists as being unrealistic or too hard to make.

## TLDR
```
from lillymol import *
from lillymol_tools import *

# Instantiate an empty set of synthetic precedent databases
sp = SyntheticPrecedentDatabases()

# Add some databases - see below for how to build.
sp.add_database("/path/to/chembl.bdb")
sp.add_database("/path/to/corporate.bdb")

mol = MolFromSmiles("C1CC(=O)NC(=O)C1N2CC3=C(C2=O)C=CC=C3N Lenalidomide")
print(sp.per_shell_data(mol))
```
might print something like


