# Tools

The python bindings come with a small number of tools, mostly derived
from tools otherwise avaialble as c++ executables.

## Selimsteg
Selimsteg is an anadrome if getsmiles. At Lilly a variety of selimsteg
tools are used to fetch the smiles for an identifier, from BerkeleyDB
databases which store key/value relationships between identifier and
smiles. For example
```
selimsteg_chembl -K CHEMBL45466
```
returns `C(=O)NCCCC CHEMBL45466'. Giving a file containing a list of id's
will result in a list of smiles being returned.

There is a python binding that can be used with such databases.

### Usage
```
from lillymol import *
from lillymol_tools import *

selimsteg = Selimsteg()
selimsteg.open_database("/path/to/chembl/selimsteg.bdb")
selimsteg.get_smiles("CHEMBL45466")
```
returns "C(=O)NCCCC CHEMBL45466".

Requesting an invalid id will result in None being returned.

Since the first thing you are likely to want to do with the smiles is
to form a Molecule, there is also
```
mol = selimsteg.get_molecule("CHEMBL45466")
```
and there is also a version that given a list of identifiers, returns a list of
Molecules.
```
ids = ["id1", "id2", "id3"....]
mols = selimsteg.get_molecules(ids)
```
If any invalid ids are presented, that molecule will be empty - no atoms. The
length of `mols` will be the same as the length of `ids`.

A suitable BerkeleyDB database can be built from a smiles via
```
iwbdb_load -c 2 -C 1 -d collection.bdb collection.smi
```
where we have specified that the key is in column 2 and the
value in column 1. Note that for larger collections, > 1M,
it can be very advantageous to build the database in /dev/shm
and then move the resulting database to its final location.
Generally BerkeleyDB databases seem to exhibit poor scalability
when loading larger datasets - perhaps it is just a matter of
tuning parameters...

Performance tends to be very good, but obviously is limited by I/O
capabilities.

## Position3D
This is also covered in the [trxn](/docs/Molecule_Tools/trxn.md) documentation.

A common operation in 3d is to react an arbitrary fragment with a core that
is fixed. Exploring different substituents at given position(s) on a scaffold.
In order to do this, the fragment must be translated and rotated so
that the atoms across which the new bond will be formed are facing each other.

It is necessary to identify the atom in the scaffold and the atom in the
fragment that will form the bond. There are `a1` and `a2` in the pseudo code
below.

The important thing to note is that `a2` must be converted from an atom number
in `fragment` to an atom number in the merged molecule, so it is added
to the number of atoms in the scaffold - which could also have been
derived from `scaffold.natoms()`.

We make a copy of `scaffold` because it will be changed and we want to
re-use it across other fragments.
```
a1 = 4
a2 = len(scaffold) + 5

mcopy = copy.copy(scaffold)
mcopy += fragment
Position3D(mcopy, a1, 1.35, a2)
mcopy.add_bond(a1, a2, BondType.SINGLE_BOND)
```
At this stage, we should have a 3d plausible arrangement of atoms `a1` and
`a2`, separated by 1.35 Angstroms.

The next task might be to perform a conformational search around the newly
formed bond. This could of course be done by successively setting the
dihedral angle, but that would be inefficient. Instead we can do a one-time
generation of all coordinates implied by a conformational scan. 
```
angle = 30
bump_check = 2.5
coords = mcopy.dihedral_scan(a1, a2, angle, bump_check)

for conf in coords:
  mcopy.set_coordindates(conf)
```
which privides a very cheap means of generating conformers. The number of
conformers returned will be influenced by the value assigned to `bump_check`
which specifies a minimum distance separating any moving or non moving atom.
Now, of course these conformers are entirely naieve, they know nothing about
their energetic implications beyond passing the bump check. More expensive
methods will be needed in order to evaluate each conformer.

## Unique Molecules
A common task in Cheminformatics is to de-duplicate a set of molecules.
This task is complicated by the problems inherent in chemical representations.
For example, if one person draws an acid in charged form, and the other in neutral
form, should those be considered identical? While that may be obvious, what
about tautomeric variants? Or the Sodium salt and the Potassium salt? Should
enantiomers be considered the same? What do you do if a chiral centre is not
marked? How should that be compared with a version that has marked chirality?

There is no 'right' answer to any of these questions, it all depends on the
task at hand. It also depends on how reliable you believe information about
chirality to be. For example, in a random sample of Chembl molecules, about
17% have atoms that are chiral, but there is no chirality marked in the
smiles. Does that mean that they have been determined to be racemic, or
has there never been a determination? Or is it a pure form, just that nobody
knows which form it is? 

The same kinds of arguments can be applied to cis-trans bonding.

With this in mind, LillyMol provides a variety of means for determining
molecular equivalence.

### TLDR
```
from lillymol import *
from lillymol_tools import UniqueMolecules

methane = MolFromSmiles("C methane")

unique_molecules = UniqueMolecules()

unique_molecules.is_unique(methane)   # returns True
unique_molecules.is_unique(methane)   # returns False
```
This works by building an internal hash of unique smiles.
If the unique smiles has not been encountered before, it will
be added to the hash, and the function returns True. If the
unique smiles has been seen before, it returns False.

Obvously the hash can grow in size considerably.

Note that by default, incoming molecules undergo chemical
standardisation, a potentially expensive process, and not necessary
if your molecules have already been standardised. See below for how to
turn off this feature.
### Variations
By default, chirality is considered
```
m1 = MolFromSmiles("C(O)[C@@H](N)C")
m2 = MolFromSmiles("C(O)[C@H](N)C")

u = UniqueMolecules()
u.is_unique(m1)
u.is_unique(m2)
```
both return true. But if you specify  `u.set_include_chiral_info(False)`
then chirality will be stripped before the unique smiles is computed
and hash lookup done.

Note that even though chirality is dropped during the comparison, your
starting molecule is unchanged - a copy of the molecule is created
internally.

Comparisons involving fragments are done in a similar way, and the
setting is controlled by setting
```
    u.set_strip_to_largest_fragment(True)
```

Comparisons involving isotopes are done in a similar way, and the
setting is controlled by setting
```
    u.set_consider_isotopes(False)
```
There is an option for considering isotopes, but converting all
isotopic atoms to a constant form.

As noted previously, by default molecules are standardised before their
unique smiles is computed. If your molecules are already standardised,
this costly process can be shut off
```
    u.set_standardize_molecules(False)
```

### Performance
Uniqueness determination is unavoidably an expensive operation,
involving chemical standardisation, possible fragment removal and
other changes, followed by a unique smiles determination.
Determining uniqueness for 50k random molecules from
Chembl on the command line
```
unique_molecules -g all -l -c -v rand50k.smi
```
takes 3.1 seconds. Doing the same thing via the python bindings
```
u = UniqueMolecules()
u.set_include_chiral_info(False)
u.set_strip_to_largest_fragment(True)

for mol in mols:
  u.is_unique(mol)
```
takes 5.5 seconds.

Interestingly, during development I experimented with using a
python Set to hold the hashed unique smiles, and that seemed to
perform better than the same data structure in c++. That surprised
me greatly, so perhaps there is something not performant in this
implementation.

## Dicer
[dicer](docs/Molecule_Tools/dicer.md) is one of the more interesting and useful
tools in LillyMol, able to efficiently implement a recursive fragment generation
that has proven useful across a variety of contexts. Python bindings for some
rudimentary dicer-like functionality is available via python bindings. While
the performance is slower than the dicer executable, the convenience may make
it worthwhile.

```
from lillymol import *
from lillymol_io import *
from lillymol_tools import *

dicer = Dicer()
dicer.set_label_join_points(1)
dicer.set_max_bonds_to_break(3)
dicer.set_   other dicer options...

mol = MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O aspirin")
frags = dicer.dice(mol)
```
yields
```
{'O=C(Oc1[1cH]cccc1)C': 1,
 '[1OH]C(=O)C': 1,
 'O[1CH]=O': 1,
 '[1cH]1[1cH]cccc1': 1,
 'OC(=O)c1[1cH]cccc1': 1,
 '[1OH][1CH]=O': 1,
 'OC(=O)c1c(O[1CH]=O)cccc1': 1,
 'O=[1CH]Oc1[1cH]cccc1': 1,
 '[1CH4]': 1
}
```
Where the result is a dictionary from unique smiles of the fragment to the number of
times that fragment is encountered in the starting molecule - in this case all fragments
are found only once.

There is a very important performance tweak that needs to be considered. During
development, dicing 2000 random Chembl molecules with dicer
```
dicer -M 16 -c -B brcb -B bscb -B nbamide -I 1 -k 3 rand.smi
```
takes about 3.2 seconds. For the default options of the python tool, this same calculation
takes 29 seconds. `dicer` goes to considerable lengths to optimise fragment identification,
and not all those methods are implemented in the python bindings - perhaps that will change
over time.

The default python implementation forms all plausible fragments, and then checks their
unique smiles against a hash. This turns out to be very expensive. But it also
enables computation of the number of times each fragment occurs in each starting
molecule. If that is not important, one of the duplication detection concepts from
dicer is applied, and the run time drops to 4.4 seconds, which is not significantly
different from `dicer` itself. While count data is still returned, it will be inaccurate,
although potentially more useful.

The optimisation is to form a bit vector of the atoms in each proposed fragment.
That bit vector is checked against bit vectors for all previously found fragments,
and if it has alrady been formed, the fragment is ignored. So what this means is
that where the count is > 1 in this case, that means the same unique smiles is
generated in different parts of the molecule. Again, this may be more useful.

Enable this feature by setting
```
  dicer.set_determine_fragment_counts(False)
```
again noting that the resulting fragment counts are not exhaustive,
but unique.

For the python version, it is recommended that chirality be removed before
processing. The reason for this is that as bonds are removed and replaced,
chirality will likely be destroyed in unpredictable ways
```
  mol.remove_all_chiral_centres()
```
before calling `dicer.dice(mol)`.

The following options to the Dicer object are implemented
### set_max_bonds_to_break(nbonds)
The maximum number of bonds to simultaneously break. By default this is set to 1, which
will generate substituents only.

### set_min_fragment_size(natoms) set_max_fragment_size(natoms)
Set constraints on the size of the fragments included in the results. Note that during
computation larger fragments may be generated, but they will be recursively broken,
potentially generating fragments satisfying the size constraints.

The default maximum fragment size is 16 atoms. Increasing this will likely
have significant performance implications.

### set_break_cc_bonds(True)
By default C-C bonds are not broken by Dicer. Set this to `True` and C-C single
bonds will be broken when the default bond breakage rules are used - no external
queries.

### set_break_amide_bonds(True)
By default, the C-N and S-N bonds in amides and sulfonamides are not broken.
That can be changed by this setting.

TODO: perhaps acids too should be prevented from fragmenting under the default rules.

### set_label_join_points(isotope)
An isotopic label applied to the join points.

### set_determine_fragment_counts(False)
Set to `False` to enable faster computation and unique fragment counts. See above.

### set_accumulate_global_fragment_count(True)
Set to True to enable accumulation of a count of all features encountered across all
molecules processed by the Dicer object. Use 

### get_global_fragment_count()
Retrieve the global fragment count Dictionary - again a mapping from fragment
unique smiles to number of molecules containing that fragment.

## Bonds to Break
The python Dicer object has very similar bond breaking rules as does `dicer`. But
for complete flexibility, you can specify your own rules by specifying one or
more substructure queries. The assumption is that the first two matched atoms
define a bond to be broken.
### add_bond_break_smarts(smarts)
Add a bond breaking substructure query. Again, first two matched atoms define
the bond.

### add_bond_break_query(fname)
Add a bond breaking substructure query via a textproto query. Actually this
is parsed by the same function as the -q option to `tsubstructure` so
there should be considerable flexibility.

### set_perceive_symmetry_equivalent_matches(bool)
If you have specified external queries for the bonds to be broken, are
symmetry equivalent matches included or not? This can be difficult because
something that might be symmetric in the starting molecule may no longer
be symmetric once fragmentation has started. But it may offer some efficiency
advantages. Generally recommend running with and without this turned on to
gauge the impact.

### Atom Typing
`dicer` offers considerable complexity wrt atom typing. That is not yet
implemented in the python version.

### Performance
The run-time of dicer is a strong function of the number of bonds that
are simultaneously broken - which corresponds to the depth of the
recursion. Here is the timing (seconds) and number of fragments generated when 
dicing 2k random Chembl molecules
for different values of the maximum number of bonds broken.

| max break | time | generated |
| --------  | ---- | --------- |
| 1 | 0.76 | 17648 |
| 2 | 2.0 | 47113 |
| 3 | 4.4 | 58238 |
| 4 | 10.6 | 79186 |
| 5 | 30 | 83826 |
| 6 | 68 | 85298 |
| 7 | 116 | 85656 |
| 8 | 156 | 85745 |

We see that as the number of bonds simultaneously split increases, the
run-time increases substantially, while the number of new fragments identified
decreases rapidly. This reflects the inherent graph simplicity of most
drug-like molecules in Chembl.

Most practical uses of dicer set a max bonds broken of 3, which we see
has a pleasing combination of good run time and identifying a good fraction
of the discoverable fragments. Fragments with 4 or more connections are
also harder to work with in any subsequent manipulation.

The c++ dicer tool takes 144 seconds to process the 8 breakable bonds task above.

### Recap
Recap fragmentation does not do any recursive/overlapping breaking. Breakable bonds
are identified and all are simultaneously broken and the resulting fragments written.
Enable this with the `set_work_like_recap` method. For example
```
O=C(N(C)C1=C(N=C2N1C=C(C=C2)C(=O)NCCOC1=CC=C(OC)C=C1)CC)CC1=CC=CC=C1 CHEMBL2354634
```
given bond breaking query `[R]!@*`, and isotope 1, generates
```
[1cH]1ccccc1 1
[1CH3]C 1
[1OH]C 1
[1cH]1cc[1cH]cc1 1
[1OH]CCN[1CH]=O 1
[n]1c2[n](c[1cH]cc2)[1cH][1cH]1 1
O=C([1NH]C)[1CH3] 1
```
These are non-overlapping fragments comprising all the atoms in the molecule.
Of course if there are size constraints on the fragments produced, some fragments
may not be generated.

## Development
It is unclear that the Dictionary mapping smiles to count currently returned is optimal. Better might
be a `DicedMolecule` proto, which should be more efficient. Or perhaps the actual Molecule
forms of the fragments should be returned in some way - perhaps just as a vector of Molecule's.
