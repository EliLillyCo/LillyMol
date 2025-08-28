# Molecule Tools

This directory contains command line tools useful for various cheminformatics tasks.

Some tools are more extensively described in the [docs](docs/Molecule_Tools) directory.

A listing of the tools with a brief description

## _fragment_filter
Part of the [zof.sh](../contrib/bin/zof.sh) fragment filtering script.
This tool implements fragment filtering rules useful for identifying fragment
like molecules.

## abraham
Computes various descriptors from Michael Abraham. Used for both descriptor
computation and fingerprint generation.

## activity_consistency
Tool for reconciling multiple experimental values
[activity_consistency](../docs/Molecule_Tools/activity_consistency.md)

## align_molecule
Aligns molecule(s) so as to align them via a substructure match with a template
molecule.

## alogp
Implementation of alogp. Used for both descriptor generation (make_descriptors)[../docs/Molecule_Tools/make_descriptors.md]
and fingerprint generation (gfp_make)[../docs/Utilities/gfp_make.md]

## atom_triples
Experimental tool for generating atom triple fingerprints - as an extension
of atom pairs. This needs to be finished and tested before widespread use.

## chirality_fingerprint
Generate a fingerprint that sees different isomers get different fingerprints.
Based on an EC type expansion around chiral atoms.
## common_names
Tool for grouping duplicate molecules
(common_names)[../docs/Molecule_Tools/common_names.md]

## dbf
Descriptor generation. Stands for Distance Between Features. This can generate
either 2D (bonds between) or 3D (Angstroms) descriptors, based on features
defined by either a substructure query or by a charge assigner or donor
acceptor assigner.
## dicer
Tool for recursively splitting a molecule into overlapping fragments
[dicer](../docs/Molecule_Tools/dicer.md).

## dicer_to_topological_types
Given a set of dicer fragments, generated with
```
dicer ... -B fragstat=fragments.textproto -B fragstatproto file.smi
```
convert `fragstatproto.textproto` into individual files based
on the number of attachment points, the number of atoms in
single attachment fragments (substituents) or the number
of bonds between attachment points (linkers).

Useful for supporting de-novo molecule generation.

## echoqry
Debugging tool for helping to debug query files. Reads and writes
a query file. Useful for converting from old form to textproto form.

## echorxn
Debugging tool for helping to debug reaction files. Reads and writes
a query file. Useful for converting from old form to textproto form.

## ec_fingerprint
A new generation EC fingerprint generator. Quite similar to the older
tool `iwecfp` but with some useful extra functionality.

## enumeration
Tool to systematically add fragments to a molecule.

## extended_atom_pairs
Generates atom pairs, mostly used as a fingerprint generator, supporting
the `-MAP` option to [gfp_make](../docs/Utilities.gfp_make.md).

## ez_descriptor, ez_fingerprint, ez_fingerprint_v2
Generates a fingerprint based on E/Z bonds present in a molecule.
While old, it has not been used much just because such information is
frequently unavailable, or unreliable.

Unclear which version is really the best choice.

## fileconv
Tool for generic reading and writing connection table files
[fileconv](../docs/Molecule_Tools/fileconv.md). Very useful.

## fingerprint_substructure
Tool that generates a fingerprint for a substructure of a molecule,
and attached atoms (fingerprint_substructure)[../docs/Molecule_Tools/fingerprint_substructure.md].

## firstatom
Given a smarts, it writes the molecule with the matched atom as the first
atom in the connection table.

## get_coordinates
Given a smarts, extract the coordinates of the matched atoms from a file
of molecules. Probably depends on the molecules having been aligned if
actual geometric positions are to be examined.

## get_substituents
Tool for collecting substituents that match a given criterion from
a set of molecules. Very useful for de-novo molecule creation.
[get_substituents](../docs/Molecule_Tools/get_substituents.md).
## gfp_erg
Implements the ErG reduced graph fingerprints, part of [gfp_make](../docs/Utilities/gfp_make.md)

## ghose_crippen
An implementaton of Ghose Crippen descriptors, the -ghose option to
[make_descriptors](../docs/Molecule_Tools/make_descriptors.md)

## grease
Once upon a time this was part of the Lilly Medchem Rules. A step that tried to identify
molecules that had large contiguous collections of likely hydrophobic atoms.

## grep_molecule
Find exact matches to a molecule(s) in a file of molecules.

## grid_fingerprint
Generates a fingerprint of a molecule based on its proximity to points in
a grid [grid_fingerprint](../docs/Molecule_Tools/grid_fingerprint.md).

## hydrophobic_sections
Descriptor generator [make_descriptors](../docs/Molecule_Tools/make_descriptors.md). This
is the -hpo option. Identifies likely hydrophobic sections of molecules.

## id_chirality
Identify unmarked chiral centres, and invalid chiral centres. Can also
enumerate all chiral possibilities.

## iwdemerit
Part of Lilly Medchem Rules. Applies demerits to molecules.
## iwdescr
Descriptor generation. This is the -w option to 
[make_descriptors](../docs/Molecule_Tools/make_descriptors.md).

## iwecfp
Fingerprint generator. This is the -EC option to
[gfp_make](../docs/Utilities/gfp_make.md). One of the best
fingerprints in LillyMol.

## iwecfp_intermolecular
A version of EC fingerprinting that only fingerprints regions of
a molecule that might be spatially close to another molecule - likely
a protein.
## iwfp
Primary linear fingerprint generator, the -IW option to
[gfp_make](../docs/Utilities/gfp_make.md).

## iwpathd
A failed idea around t-shaped fingerprints. Works but is too slow
to be useful. Should be re-implemented, and tested to see if it is
useful.

## jwcats
The -CATS option to [gfp_make](../docs/Utilities/gfp_make.md). Often
works very well in svmfp models.

## jwdip
The -jwdip option to [make_descriptors](../docs/Molecule_Tools/make_descriptors.md].

## jwdist
The -jwdist option to [make_descriptors](../docs/Molecule_Tools/make_descriptors.md].

## jwestate
The -estate option to [make_descriptors](../docs/Molecule_Tools/make_descriptors.md].

## jwmedv
The -medv option to [make_descriptors](../docs/Molecule_Tools/make_descriptors.md].

## jwmolconn
The -jwmc option to [make_descriptors](../docs/Molecule_Tools/make_descriptors.md].
Molecular connectivity descriptors.

## jwmorse
The -morse option to [make_descriptors](../docs/Molecule_Tools/make_descriptors.md].

## jwsadb
The -jurs option to [make_descriptors](../docs/Molecule_Tools/make_descriptors.md].
3D surface based descriptors.

## linear_fingerprint
A new generation linear fingerprint generator, with extra flexibility.

## long_molecules
Identify molecules that are "long". Not really needed, the same effect can
be discerned by examining descriptors coming from `iwdescr`.

## maccskeys
Implementation of MACCS keys. This implementation has diverged considerably
from the original definitions. This is the -MK option to
[gfp_make](../docs/Utilities/gfp_make.md).

## make_these_molecules
Makes a subset of the molecules defined by a large virtual library, without
enumerating the library.

## medchemwizard
Implementation of Medchem Wizard.

## minor_changes
Make small changes to a molecule. Useful in de-novo molecule
construction.

## mkfrag
Convert the individual fragments in a molecule to separate molecules.
So CC.C becomes
```
CC
C
```

## mol2qry
Convert a molecule to a query
[mol2qry](../docs/Molecule_Tools/mol2qry.md).

## mol2SAFE
Implementation of the [SAFE](https://arxiv.org/pdf/2310.10773.pdf).

## molecular_abstraction
Apply various molecular abstractions
[molecular_abstraction](../docs/Molecule_Tools//molecular_abstraction.md)

## molecular_grid
Generate a grid that can surround a molecule.
[molecular_grid](../docs/Molecule_Tools/molecular_grid.md). Need
to complete that documentation.

## molecular_merge
Concatenate the molecules in an input stream to one, multi-fragment molecule.
The opposite of `mkfrag`.

## molecular_scaffold
Generate the scaffold of a molecule. Prefer using 
```
molecular_abstraction -a 'scaffold(WRITE)' file.smi
```

## molecular_transformations
Apply multiple reactions to a set of molecules.

## molecular_variants
Tool that suppresses molecules that are minor variants of a molecule
previously encountered in a file. The variations are defined by a
set of reactions.

## molecule_filter
Tool for qickly filtering a set of molecules based on properties.
[molecule_filter](../docs/Molecule_Tools/molecule_filter.md).

## molecule_subset
Similar to `fingerprint_substructure`. Generates molecules that
are derived from a substructure match.
[molecule_subset](../docs/Molecule_Tools/molecule_subset.md).

## molecules_from_reagents

## msort msort_parallel
Sort molecules by various properties
[msort](../docs/Molecule_Tools/msort.md).

## numbered_smiles
Write a smiles with atom numbers as isotopes, or atom map numbers.
Most useful for debugging.

## overlapping_fragment_model
A Fragment Additivity model implementation where it is possible for
fragments to verlap.

## parsimonious_set
Given a set of diced molecules, produce a small subset that exemplifies
as many of the frgments identified by dicing.

## pharmacophore_2d
Convert substructure matches to queries
[pharmacophore_2d](../docs/Molecule_Tools/pharmacophore2d.md).

## preferred_smiles
Mostly used internally at Lilly to generate smiles in specific forms.

## psafp
The -PSA fingerprint in [gfp_make](../Utilities/gfp_make.md).

## pubchem_fingerprints
Implementation of the Pubchem fingerprint definition. Available
from within [gfp_make](../Utilities/gfp_make.md).

## r1r2etc
Used for identifying substituents.

## random_geometric_changes
Mostly useful for testing. Apply random geometric changes.

## random_molecular_permutations
Apply changes from a pre-defined set of known molecular permutations
[random_molecular_permutations](../docs/Molecule_Tools/random_molecular_permutations.md).

## random_smiles
For each smiles, generate random versions that describe the same
molecule.

## reduced_graph
A reduced graph implementation based on [Harper](https://doi.org/10.1021/ci049860f).

## remove_and_label
Another tool for identifying substituents and cores.

## remove_matched_atoms
Given substructure matches, remove them. Prefer to use
```
molecular_abstraction -a 'rmat(smarts WRITE)'
```
but `remove_matched_atoms` has more options and flexiblity.

## representations
Write a tabular file of standardised forms. Useful for uniqueness
and database loads.

## retrosynthesis
Part of the retrosynthesis tool set. Needs documentation.

## rgroup
Another tool for identifying R groups. Likely better to use
get_substituents.

## ring_extraction
Part of the ring_replacement suite (../docs/Molecule_Tools/ring_replacement.md).
This tool extracts rings from known collections.

## ring_fingerprint
Generates the -RING fingerprints within [gfp_make](../docs/Utilities/gfp_make.md).

## ring_replacement
The actual ring replacement part of the ring_replacement tool suite
[ring_replacement](../docs/Molecule_Tools/ring_replacement.md)

## ring_replacement_collate
Given multiple sets of rings extracted by ring_extraction from different
collections, aggregate those separate files to be as if they came from a
single collection.

## ring_size_fingerprint
The -RZE option to [gfp_make](../docs/Utilities/gfp_make.md).

## ring_substitution
The -RS option to [gfp_make](../docs/Utilities/gfp_make.md). This is
a surprisingly powerful fingerprint, frequently showing up in svmfp
models.

## ring_trimming
Removes various kinds of rings from a molecule. Built for a specific
project, maybe not generally useful.

## rotatable_bond_fingerprint
Generate a fingerprint based on rotatable bonds in a molecule. Not
widely used.

## rotatable_bonds
Compute the rotatable bonds in a molecule. A great deal of flexibility
in how rotatable bonds are identified.

## rule_of_five
Implementation of Lipinski's rule of five.

## rxn_fingerprint
Fingerprint a reaction.

## rxn_reverse
Reverse a reaction.

## rxn_signature
Generate a unique representation for a reaction - defined by a context around
the reaction core.

## rxn_standardize
Perform various standardisations on a reaction.

## rxn_substructure_search
Substructure search across reactions.
## same_structures
Examine two files for differences. Where tokens differ, see if they can be
interpreted as structures and compare the unique smiles. Useful in tests
across versions where the unique smiles might change.

## smiles_mutation
A "wild" de-novo structure generator that works by making string mutations
to smiles. Almost all structures generated are rejected, but can produce interesting
molecules. Needs to be re-implemented to increase the rate of success.

## sp3_filter
Filter molecules by fraction of sp3 atoms.

## substituent_model
Implementation of group additive model.

## substitutions
Another R group identification tool.

## substructure_match_fraction
Compares molecules by working out the fraction of matched atoms from a query.
Needs documentation, but possibly obscure.

## substructure_mcs
Looks for molecules within a set that are subsets of other molecules, possibly via
application of one or more reactions.

## superimpose_by_matched_atoms
Superimpose molecules based on the matches to a substructure query. See also
`align_molecule`.

## tautomer_generation
Attempt at tautomer generation. Not robust, instead LillyMol tools depend more
on chemical standardisation which forces tautomeric forms into common forms.

## temperature
Generates the -MPR (Molecular PRoperty) fingerprint in
[gfp_make](../docs/Utilities/gfp_make.md)

## tnass
Optimised fragment substituent model that has dependencies between
fragments.

## topotorsion
Create topological torsions - 4 atom paths.

## topotorsion_fingerprints
The -TT (Topological Torsion) fingerprints to [gfp_make](../docs/Utilities/gfp_make.md)

## tp1_summarise
Part of Lilly Medchem Rules. Summarise the results from running the Medchem Rules.

## tp_first_pass
Part of Lilly Medchem Rules. This is the first component of the pipeline, which performs
simple filtering based on counts and other easy/fast to compute properties.

## trxn
General purpose reaction implemenation [trxn](../docs/Molecule_Tools/trxn.md).

## tshadow
Generates 3 dimensional "shadow" type descriptors. Molecules are aligned along
the X axis and various area based features computed. Often works very
well in models. This is the -shadow option to [make_descriptors](../docs/Molecule_Tools/make_descriptors.md).

## tsmiles
Tester for unique smiles - not for general use.

## tstandardise
Tester for chemical standardisation - not for general use.

## tsubstructure
Substructure search tool [tsubstructure](../docs/Molecule_Tools/tsubstructure.md).

## tsubstructure_summarise_hits
Postprocess tool once `tsubstructure -m QDT` has been run.

## tsymmetry
Identifies molecules containing high degrees of symmetry. These
were considered undesirable for acquisition.

## unique_molecules
Remove duplicate molecules [unique_molecules](../docs/Molecule_Tools/unique_molecules.md)
## verloop

## xlogp
Implementation of xlogp - also generates the -XLOGP fingerprint
in [gfp_make](../docs/Utilities/gfp_make.md).

## xray_structure_compare
3D structure comparison tool.
