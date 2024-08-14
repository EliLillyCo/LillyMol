# Molecule Filter

## Motivation
This tool is designed to quickly filter large collections of molecules as
efficiently as possible.

Internally tests are performed roughly in order of expense of computation.
For example, if reading a smiles, the number of atoms in the smiles can
be very quickly determined by examing the smiles string, without performing a much
more expensive Molecule interpretation. Note however that this will only
be beneficial to the extent that molecules get rejected for size. If hardly
any molecules are discarded for size, using this approach may be worse.
Some tests are very expensive, xlogp for example, and those will be
performed after all cheaper tests.

Note that such a task is 'pleasingly parallel' and generally this task is
used in combination with batch processing large collections in sharded parallel
form.


## HOWTO
Built a textproto configuration file based on the proto
[proto](/src/Molecule_Tools/molecule_filter.proto).
```
molecule_filter -F config.textproto:w


## Implementation

The tool is drive by a [proto](/src/Molecule_Tools/molecule_filter.proto) which
provides the configuration information for how the filtering is to be done. An
example that sets all available properties might be
```
min_natoms: 10
max_natoms: 40
min_heteroatom_count: 1
min_heteroatom_fraction: 0.05
max_heteroatom_fraction: 0.90
min_nrings: 1
max_nrings: 6
min_aromatic_ring_count: 1
min_aliphatic_ring_count: 0
max_aliphatic_ring_count: 4
min_rotatable_bonds: 2
max_rotatable_bonds: 10
max_ring_system_size: 3
max_aromatic_rings_in_system: 3
min_tpsa: 20
max_tpsa: 180
min_alogp: -2
max_alogp: 6
min_xlogp: -2
max_xlogp: 6
min_hba: 2
max_hba: 10
min_hbd: 0
max_hbd: 5
largest_ring_size: 7
exclude_non_organic: true
exclude_isotopes: true
max_halogen_count: 7
max_distance: 30
min_sp3_carbon: 1
max_aromatic_density: 0.8
max_chiral: 2
max_number_fragments: 2
```
Note that as written this makes no sense. It would be silly to filter
on both alogp and xlogp. Choose one. That said, the xlogp implementation
is more expensive than the alogp implementation, so if a logP filter is
required, and speed is of the utmost concern, alogp would be preferred.
Which one might provide better predictions is beyond the scope of this 
document!

A typical invocation might be
```
molecule_filter -F config.textproto -v large.smi > passed.smi
```
Using the -v option results in a summary of failure modes being written
to stderr at the end of processing. If you care about the molecules that
have been discarded, add the `-B rejected` option, to specify a file where
rejected molecules are written.

The tool takes 62 seconds to process 2.4M molecules in a recent Chembl,
reporting (because of the -v option)
```
Read 2409270 molecules, passed 1488785 0.61794
7905 non organic 
1901 isotope 
6582 too few atoms 10
253993 too many atoms 40
18088 too few rings 1
10255 too many rings 6
282 too few heteroatoms 1
2170 min heteroatom fraction 0.05
13 max heteroatom fraction 0.9
41203 too few aromatic rings 1
64474 too many aromatic rings 4
0 too few aliphatic rings 0
901 too many aliphatic rings 4
22132 ring systems too large 3
0 too many aromatic rings in system 3
7964 ring too large 7
84062 too few rotatable bonds 2
60752 too many rotatable bonds 10
14284 low TPSA 20
3559 high TPSA 180
294 low ALOGP -2
56821 high ALOGP 6
13937 too few HBA 2
34722 too many HBA 10
0 too few HBD 0
13070 too many HBD 5
1761 too many halogens 7
0 molecules too long 30
625 too few CSP3 1
47640 aromatic density too high 0.8
151095 too many chiral centres 2
18800 too many fragments 2
```
Note that some rules may show zero matches because they might overlap
with previously run rules. Note that the order above is the order of
the properties in the proto definition. The actual order of computation
is as in the source file, and is generally from cheapest to most expensive.
As we learn more, the order of execution may change.
