# Structure Database

This describes the python bindings for the structure database functionality
described in [structure database](/docs/Molecule_Tools/structure_database.md).


## TLDR
```
from lillymol import *
from lillymol_tools import *

db = StructureDatabase()
# Open the corporate database for matching proprietary structures
db.open_read("/path/to/corporate/database/structures")

# Optionally also open other databases
db.open_read("/path/to/chembl/structures")

mol = MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O aspirin")
db.lookup(mol)
```
might yield
```
CHEMBL25:F%CHEMBL1697753:F%CHEMBL2296002:F%CHEMBL3833325:F%CHEMBL3833404
```
The exact lookup has yielded 5 matches. Note that if more than one database
was open, the result might look like
```
CHEMBL25...|PROPRIETARY1:PROPRIETARY2:...
```
where a vertical bar appears between data retrieved from different databases.

The first 'CHEMBL25' is an exact match to the starting molecule. The others
all match with an 'F%' qualifier, which means they matched because those
starting molecules have lost one or more smaller fragments, leaving someting
that matches our smiles.

Longer term, a protocol buffer will be returned which avoids awkward, and
possibly expensive string processing.

## Details
The database used must have been build by `structure_database_load`, and the
same naming conventions apply.

There can be a second parameter to lookup(). It is an or'd set of values
from the LookupParams enumeration
```
EXACT
STRIP
NOCHIRAL
GRAPH
```
These describe changes that are applied to the input molecule before the lookup
is done. Using `EXACT` means no change at all. That is the default. Note that
these changes are NOT applied to your molecule - internally a copy is made and
the changes are applied to that copy, so your input molecule will not be altered
by this function.

So functionally, starting with a multi-fragment molecule
```
mol = MolFromSmiles("CCC.C")

db.lookup(mol, LookupParams.STRIP)

# should yield the same result as

mol.reduce_to_largest_fragment()
db.lookup(mol)
```
but in the second case, your molecule is altered. If you are stripping to the
largest fragment anyway, the second way would be more efficient.

In the future we hope to be able to avoid the Molecule copy in most cases where
an exact lookup is requested.
