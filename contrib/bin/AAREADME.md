# contrib/bin

This directory contains scripts that are used with LillyMol tools. 

A quick explantion of each is described here, some may be more fully described
with their own entry in the [docs](/docs) directory. Almost all will report usage information
by executing the script with no arguments or options.

Many of the tools in this directory will have two forms. A .sh form and a .rb or .py form.
In all cases, the .sh version is a wrapper than will call the underlying Ruby or Python
script, sometimes assembling a set of default arguments. Gererally prefer to use the
shell wrapper.

## Lilly_Medchem_Rules.sh
Run Lilly Medchem Rules using the executables built with LillyMol rather than the original repo.
We try to keep consistency between implementations, but there will always be differences due
to things like changes in standardisation, changes in aromaticity definitions, bugs.

## QED.sh
An implementation of QED [docs](/docs/Molecule_Tools/QED.md)

## abraham.sh
An implementation of some of Michael Abraham's parameters. These features often perform
very well in tree models, and fingerprint models. References at
[abraham.cc](/src/Molecule_Tools/abraham.cc)

## ap.sh
Another group additive model from Michael Abraham. Constants are stored at
[abraham](/data/abraham_hbond_constants)

## assign_donor_acceptors.sh
Hydrogen bond donor and acceptor queries are stored at [bruns](/data/queries/hbonds/AAREADME.md).

## assign_formal_charges.sh
Formal charge assignments. The queries that drive the assignments are stored in
[bruns](/data/queries/charges/AAREADME.md)

## chgfp
Generates linear fingerprints on molecules. The charge assigner is applied to the molecules before
fingerprint generation. TODO:ianwatson Needs to be made more aware of the charges applied.

## concat_files.sh
The LillyMol file joiner. Shares feature with the Linux 'join' command, but more flexibility and
taylor made for many of the kinds of files generated in LillyMol.

## confusion_matrix
 confusion matrix generator.

## consecutive_reactions
Runs consecutive reactions on molecules. Useful for positional analogue scanning
type work, where a specific number of changes are to be applied.

## dbf.sh
Distance Between Features. The wrapper script adds the LillyMol charge and donor/acceptor queries.
The resulting data consists of statistics about the 2D and/or 3D distances between those assigned
features. Note that 3D features only make sense if your input molecules are 3D. Any set of queries
can be used to describe the features.

## dopattern
Quite possibly the most useful tool I have ever written [dopattern](/docs/General/dopattern.md]. If you
ever have to perform an identical operation across multiple files or directories, dopattern
can probably make things easier.

## fileconv
LillyMol's tool for reading, optionally changing and filtering, and writing streams of molecules
[fileconv](/docs/Molecule_Tools/fileconv.md). Extremely useful.

## fragment_filter
This ruby script filters a set of molecules to those that might rightly be considered to be "fragments".
Note that there are no hard and fast rules for "what is a fragment".

## get_coordinates.py
A LillyMol python tool for extracting the coordinates of matched atoms from a molecule. To some extent
it mirros the functionality of the c++ tool of the same name.

## gfp_erg
This generates the -ERG fingerprint. Mostly this is only called from gfp_make. The initial implemenation
was from Nick Stiefl, back in the late 1990's.

## gfp_make
LillyMol's fingerprint generator [gfp_make](/docs/GFP/gfp_make.md).

## ghose_crippen.sh
Implements a fragment additivity model from Ghose and Crippen. Mostly used by [make_descriptors](/docs/Molecule_Tools/make_descriptors.md)

## graph_edit_changes
A wrapper to random_molecular_permutations which provides a limited set of transformations
which closely approximate graph edit changes only. Possibly useful for evolving molecules, but
more likely of theoretical interest.

## ha.sh
Heteratom descriptors. Describes the relative spacing between heteroatoms.
This is the -ha descriptor in [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).
Not intended for general use.

## hb.sh
Hydrogen Bonding descriptors. Describes the relative spacing between Hydrogeon bonding features.
This is the -hb descriptor in [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).
Not intended for general use.

## hydrophobic_sections
Describes regions of molecules that likely consist of hydrophobic atoms. This iw the
-hpo option to  [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).

## iwdescr
2D molecular descriptors [iwdescr](/docs/Molecule_Tools/iwdescr.md). These are the -w
descriptors in  [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).

## jwcats.sh
The -cats descriptor in  [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).

## jwdist.sh
The -jwdist descriptor in  [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).
Generally these features describe the disposition of atoms along the longest spatial
path in the 3D molecule.

## jwestate
The -estate descriptor in  [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).

## jwmedv.sh
The -jwmedv descriptor in  [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).

## jwmolconn.sh
The -jwmc descriptor in  [make_descriptors](/docs/Molecule_Tools/make_descriptors.md). These
are molecular connectivity and graph theory features.

## maccskeys.sh
An evolution of the original MACCS keys. Many of the original definitions have been changed
and new queries added. Currently there are 6*32 = 192 features defined. This is the -mk
option to  [make_descriptors](/docs/Molecule_Tools/make_descriptors.md) and the -MK and -MK2
option to [gfp_make](/docs/GFP/gfp_make.md).

## make_descriptors
Molecular descriptor generation [make_descriptors](/docs/Molecule_Tools/make_descriptors.md).

## medchem_wizard
Makes reaction based transformations to molecules. Uses a set of isostere like reactions from
a paper from Abbot to generate new molecules, likely strongly relaed to starting molecules.
This tool reliably generates interesting new molecules.

## minor_changes
This tool exhaustively generates minor changes to a starting molecule. [minor_changes](/docs/Molecule_Tools/minor_changes.md)

## mybenzene.py
A python script that can generate an arbitrary sized recursive smarts - used for testingtesting.
For fun, just execute it and you will receive a list of recursive smarts that all match
benzene. We have tested LillyMol with a 2MB smarts generated by this tool and they
work as expected.

## pubchem_fingerprints
An implementation of molecular fingerprints published by the Pubchem project
[Pubchem](https://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.pdf). Also
the -PUBCHEM option to gfp_make.

## random_molecular_permutations
Randomly permutes starting molecules [random_molecular_permutations](../docs/Molecule_Tools/random_molecular_permutations.md)

## rdkit2gfp
Converts an RDKit fingerprint to gfp form.

## reduced_graph
A reduced graph implemenation. While interesting, reduced graphs only occasionally seem to be useful.

## replace_sidechain
Uses the rings extracted by [get_substituents](/docs/Molecule_Tools/get_substituents.md} to
elaborate molecules with known sidechains from Chembl. Can be a nice way of enumerating new molecules.

## rf_make rf_evaluate
Build and evaluate Random Forest models. Underlying building and scoring are dong with SKLearn. These
scripts know about descriptor computation, and make_descriptors, so a smiles file can be directly
scored. The script will figure out which descriptors are needed, and compute the sets
necessary.

## smiles_mutation
A 'wild and crazy' molecule generator. Randomly permuts smiles strings as text, and then attempts
molecular interpretation on the results. Mostly of theoretical interest.

## speadplot.jl
Generate plots showing the distances encountere during any of the gfp_spread* variants. Generate 
fingerprints, run gfp_spread*, post-process with nplotnn -S <fname> and send <fname> to spreadplot.

## svmfp_summarise_results
Given multiple outputs from iwstats, summarise and aggregate those results.

## topotorshion
Tolological Torsion descriptors and fingerprints. Used by gfp_make and make_descriptors.

## xgbd_make xgbd_evaluate
Build and score XGBoost models. These scripts know about the descriptors produced by
make_descriptors, so a smiles file can be scored directly.

## xgboost_variable_studies.rb xgboost_variable_studies_analyse.rb 
Shuffles descriptor files in order to examine the effect of feature ordering when
building XGBoost models. The first script generates some number of shuffled files,
then builds a model on each shuffling. The *_analyse.rb script then aggregates
the variable importance data stored with each of those models.
