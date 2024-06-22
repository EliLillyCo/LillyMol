# Synthetic Precedent

## Background
The idea of using EC type fingerprints computed across collections as an indicator
of the likely synthetic feasibility of a new molecule was developed independently
many years by both Ertl (Novartis) and Watson (Lilly). The LillyMol implementation
uses a BerkeleyDB database to store bits. While the current implementation is for
Extended Connectivity (Morgan) type fingerprints, there is no requirement that
this kind of fingerprint be used. [smi2rings](/docs/Molecule_Tools/smi2rings.md)
is a tool that implements this idea for rings, rather than fingerprints. Work
is underway to implement this concept for dicer fragments. Any plausible
molecular subset can work.

In practice we find that using both fingerprints and rings to assess precedent
for a new molecule are very effective at identifying molecules likely to
be discarded by chemists as being unrealistic or too hard to make.

Fundamental to use of this tool is to understand what an EC, Morgan, type
fingerprint is, [Rogers and Hahn](https://pubs.acs.org/doi/10.1021/ci100050t)
for example.

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
```
[212404, 365, 3, 1]
```
This indicates that among the zero radius shells in the molecule, which
are just a single atom, the rarest such atom has 212k instances in the
databases examined.

At radius 1, which will be a central atom together with the atoms directly
attached to it, the rarest such grouping of atoms has 365 instances in the
databases. As the radius increases, the number of examples in the databases
should decrease monotonically - although there is a low probability of bit
collisions that might upset this.

And for a shell of radius 3, diameter 6, the rarest grouping of atoms has
just one examplar - likely the molecule itself which is in Chembl.

The fundamental idea behind the method is that if there are arrangements of
atoms in the molecule that are minimally precedented, that raises the risk
that this molecule could be hard to make, and be rejected by Chemists.

Clearly at radius 3, we are encountering shells that can be very specific
to a particular molecule, so most use of the tool tends to focus on
radius 2. Where again, there is no "right" answer for how to interpret
results. Candidate molecules can be sorted by their synthetic precedent
score, and the lowest ranked candidates de-prioritised.

## Other considerations
This technique does not work for small molecules. For example if we
lookup propane
```
mol = MolFromSmiles("CCC")
sp.per_shell_data(mol)
```
we might get `[3350136, 126, 0, 0]`. This shows a radius 2 miss, which
would indicate a never before seen arrangement of atoms in this molecule!
But instead, this molecule is so small that a radius 2 shell cannot be
formed, and so there are zero precedents for radius 2 shells.

And as expected, if we evaluate a molecule containing an atom not in
the database at all
```
mol = MolFromSmiles("[Au]")
sp.per_shell_data(mol)
```
we get the expected `[0, 0, 0, 0]`.

It is important to note that this tool comes with some important
limitations. There might be molecules that are very difficult to make,
but which are economically important, and so are in fact common.
Those molecules, and their
fingerprints, will show up in the database as being well precedented,
but they are nevertheless hard to make.

## Databases
This tool depends on the prior construction of BerkeleyDB databases
containing the fingerprints computed over collections that represent
molecules that have been made. Usually collections like your
corporate collection, chembl, and catalogues of molecules that are
available for immediate shipment are all good candidates.

We do not recommend using virtual compounds for databases, since
these may contain molecules that have been hypothesised, but which
in fact may not react as planned, leading to erroneous conclusions.

The EC type fingerprints stored in the database are stored with a user
defined atom typing. Early work with this tool showed that simplistic
atom typing, atomic number, for example did not adequately capture
atomic context. We recommend using a custom atom typeing `-P sfx` when
building a database. This atom type consists of the following atomic
features

1. atomic number
2. ring status - aromatic, aliphatic, fused
3. formal charge
4. number of connections
5. Hydrogens attached.

This combination of atomic properties seems to do a good job of
capturing the relevant information about an atom in this context. Again
there is no 'right' answer for which atom typing to use.

Building a database might be done via
```
iwecfp_database_load -R 3 -d chembl.spc.bdb -M noex -P sfx -g all -l -v chembl.smi
```
creates `chembl.spc.bdb` which can be used by the `add_database` method above.

The tool can create databases where the identity of the first molecule
exemplifying a feature is stored, but in practice that has not proven
useful. Therefore `-M noex` is specified.

Feel free to experiment with different atom typing schemes to use with the
`-P` option [atom typing](/docs/Molecule_Lib/atom_typing.md). The atom
typing is stored in the database with the key `_ATYPE`. Examine via
```
iwbdb_fetch -d /path/to/collection/collection.spc.bdb -K _ATYPE
```

## De-Novo Generated Molecules.
One of the primary use cases for this tool is prioritising virtual and
de-novo generated molecules.  While use of filters such as
[Lilly Medchem Rules](https://github.com/IanAWatson/Lilly-Medchem-Rules) is
very important in removing likely reactive or unstable molecules, these
rule sets often cannot deal with absurd molecules coming from de-novo
generation.

For example the Lilly Medchem Rules pass a bridged benzene ring
`c12ccc(C1)cc2` without complaint. The molecule is valence OK and
does not contain any obviously reactive or unstable functionality.
But most filtering lists were generated based on input from Chemists
and nobody ever thought of blocking such absurd molecules.

Rather than burden filtering tools like the Lilly Medchem Rules
with what would likely be a limitless stream of outrageous, but
benign structures, tools like this one, ring rarity, and other precedent
based tools work well to stamp out these molecules. 

But even then, results can be unpredictable.
For example if we score the bridged benzene molecule above in this tool
we get `[636703, 19935, 1721, 3]` which would indicate a well
precedented molecule - all individual components of this molecule
are common. But looking up this ring system in a database of precedented
rings will quickly flag this as unprecedented. It is also possible that
a different atom typing in building the database could capture enough
information to identify this molecule as unprecedented.

## Performance
Note, this is all for the command line version, `iwecfp_database_lookup`
but the python results should be similar.

The throughput of this technique can be relatively high. Obviously
the process can be I/O bound since a database is involved. Currently
the Chembl database is 216MB, which is small enough to be read
into RAM. Running a 10k molecule query against a 'cold' database
takes about 7 seconds. Repeating that query soon after, while the
operating system is likely caching the data, drops run-time to
about 5 seconds.

The option combination `-M slurp=n` allows caching the database
contents into RAM. This can significantly speed up processing. The
number given to slurp is the minimum number of examples to load
into RAM. Bits are stored in the database together with the number
of examples of that bit. Specifying `-M slurp=1` will result in
the entire database loading to RAM. There will be an initial
hit to performance, but after initialisation, processing can be
fast. For example, slurping all of Chembl, and then processing
100k molecules is done in 10 seconds.

If you specify a value other than 1 for the slurp value, then
bits with fewer than that number of examples will not be considered,
they will show up as missing bits.
This might be fine if you are running with a radius of 2 and
want to discard molecules with fewer than (say) 5 examples,
`-R 2 -M slurp=5` will accomplish that.

Within python, the slurp functionality is available via the
slurp method,
```
sp.add_database("/path/to/chembl.bdb")
sp.add_database("/path/to/corporate.bdb")
sp.slurp(1)
```
slurps the entire contents of all databases to RAM.

Note that if you are never going to query radius 3 bits, the
size of the database, and size of any cache, will be substantially
reduced.

## Summary
Assessing the synethic viability of molecules is hard. There are
many retrosynthetic tools, both old and newer AI driven approaches,
that can be used for assessing synethic feasibility. Most of those
suffer from poor performance, and/or very low throughput.

The approaches outlined here can be a useful component of a 
high performance virtual and de-novo molecule generation
scheme.

