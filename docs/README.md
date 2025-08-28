# LillyMol Tools
There are a great number of LillyMol command line tools. Over the years
some have proven more useful than others, but all have played useful roles
on one or more projects. Here is an alphabetic listing of tools with a very
brief description. Some of the tools have their own documentation.

Some tools are seldom used, and/or have been superseded by others. So you
will see things that seem duplicative. Work is ongoing to clean deprecate
obsolete tools and replace with their more modern counterparts.

## abraham

## activity_consistency
[activity_consistency](/docs/Molecule_Tools/activity_consistency.md)

## align_molecule
Crude 3D alignment tool. Uses a substructure query to identify the
atoms that will be placed at the origin, along the X axis and on 
the Y axis. Mostly used to aid visualisation.

## atom_triples
Fingerprint generator. Logical extension of the idea of atom pairs.

## average
Find the average of a column of numeric values, but has a lot of
flexibility - including weighted averages.
```
average -d w_natoms -d w_rotbond rand.w
```

## bucketise
Convert continuous values into discrete ranges. Note that the
column(s) being processed are replaced by their bucketised values
```
bucketise -v -b 10 -j -P xlogp rand.w
```
## buildsmidb_bdb
Builds A BerkeleyDB database with the smiles, or smiles variants
as the key, and the name as the value. See
[in_database](Molecule_Tools/in_database.md).

## byte_offset_index
Obscure. For a large file, generate the byte offsets where various
fractions of the file start. Useful for having multiple processes
process a large file without splitting it.

## chirality_fingerprint
Generates a fingerprint based on chirality. Initial results are disappointing.

## class_label_translation
Infrastructure. Maintain mappings between text and numeric class labels. 

## common_names
Identify groups of structures that are identical under various smiles assumptions.
[common_names](Molecule_Tools/common_names.md).

## concat_files
Join descriptor files. Any number of files can be processed and there is
no assumption about sorting. 

## correlate
Enforce a correlation structure between columns of data. Unrelated to
Cheminformatics, was developed to support Clinical Trial simulations, where
for example, we might have a distribution of patient weights and patient
heights, and we need to enforce a correlation structure beween those values.

## dbf
Distance betwen features fingerprints. Can do both topological and 3D
distances. The default script uses LillyMol pharmacaphoric features
and computes bond distances, and/or geometric distances between feastures.

## descriptor_file_cat
Concatenate descriptor files. Knows about header records, so if the
second+ files are in different order, it can properly handle them.

## descriptor_file_filter
Filter a descriptor file based on values in the columns. Specifically set
up to work with descriptor files. Any number of filtering conditions can
be used, and it can recognise columns by name.
```
descriptor_file_filter -e 'w_natoms<40 && w_xlogp<6' -v /rand.w
```
## descriptor_file_select_rows
Select a subset of the rows from a descriptor file.

## descriptor_file_sort
Sorts descriptor files - helpful because it knows about header records and
column names.

## descriptor_file_to_01_fingerprints
Convert numeric value(s) in a descriptor file to fingerprints. The input
file must consist entirely of positive numbers. If a binary fingerprint
is being produced, all non zero values set the bit. If a sparse counted
fingerprint is being generated, counts are used.

## descriptor_file_to_svm_lite
Convert a descriptor file to svm-lite form. Needs to be combined with
an activity file when writing a training set. When writing a test set,
it checks that the features present are consistent with what was present
when the training set was processed.

## descriptor_file_transpose
Useful for debugging - columnar data can be hard to see outside a
proper data analysis environment.

## descriptors_to_fingerprint
Complex. Converts a descriptor file to fingerprint form.

## dicer
Very useful but very complex molecular fragmentation.
[dicer](/docs/Molecule_Tools/dicer.md)

## difference_sort
Obscure. Sort a file based on the difference between two columns. Mostly
useful for dealing with observed vs predicted situations.

## distribution
For a given feature, generate a distribution of the values encountered.
```
distribution -d w_natoms rand.w
```
note similarities with bucketise.

## ec_fingerprint
A new generation ec fingerprint generator. Mostly works but needs more
development.

## echoqry
Read a substructure query and write it. Mostly useful for debugging,
but can also be used for converting old style query files to new
textproto forms.

## echorxn
Read a reaction file and write it. Mostly useful for debugging, but can
also be used for converting old style reaction files to new textproto forms.

## enough_inventory_bdb
Special purpose BerkeleyDB tool for reading Lilly inventory data.

## exon_correlated
Tool for computing correlations under conditions or quite large data,
where reading into R/python might be prohibitive.

## extended_atom_pairs
Robust atom pair fingerprint generator.

## ez_descriptor
Generate molecular descriptors based on E/Z bond specifications. Generally
not used, because the reliability/presence of such information can be low.

## ez_fingerprint
Fingerprint based on E/Z bond specifications. See previous comment.
## ez_fingerprint_v2
Newer generation of E/Z fingerprinting. Still unsatisfactory. E/Z
data is not robustly handled in LillyMol.

## feature_scaling
Infrastructure. Scales and unscales numeric data - for use in models.

## fetch_sdf
## fetch_sdf_quick
Extremely useful utilities for fetching a subset of molecules from an .sdf
file. The difference is that `fetch_sdf` retrieves the records in the same
order as they appear in the file containing the list of items requested,
whereas `fetch_sdf_quick`
returns molecules in the order they appear in the structure file.

## fetch_smiles
## fetch_smiles_quick
Fetching utilities for smiles. By combining with the `-X` and `-Y` options
these become generally useful differencing tools.

## fetch_tdt
## fetch_tdt_quick
Fetching utilities for TDT files.

## fileconv
Generally useful structure processing tool [fileconv](/docs/Molecule_Tools/fileconv.md)

## fingerprint_substructure
Fingerprint just a subset of the atoms in a molecule. The subset is defined by
a substructure query, and an optional radius. Useful for identifying unique
reagents - to a given radius. [fingerprint_substructure](fingerprint_substructure.md)

## firstatom
Re-write a smiles with a given atom as the first atom in the smiles.

## _fragment_filter
Part of the zof.sh fragment filtering script.

## get_coordinates
For a given set of matched atoms, extract the coordinates from molecules.

## get_substituents
For a core defined by a substructure, identify the substituents attached.
[get_substituents](get_substituents.md)

## gfp_add_descriptors
Obscure. Adds information from a descriptor file to a fingerprint file.

## gfp_distance_filter
Filter one set of fingerprints according to their distance to another set.
Hard to use, but useful in certain circumstances.

## gfp_distance_matrix
Generate a distance matrix from a set of fingerprints.

## gfp_distance_matrix_iwdm (deprecated).
Generate a binary distance matrix file from fingerprints.

## gfp_erg
Fingerprints from Erg reduced graph representations. Only used within gfp_make
where Erg fingerprints sometimes do well in models.

## gfp_flatten_counted
Reduce all counted fingerprints to max count 1.

## gfp_histogram
Generate histogram of distances in a collection.

## gfp_incremental_diversity
Report extra diversity added as new molecules are added to a set. Useful for
tracking the progress of a discovery/optimisation program.

## gfp_iterative_expansion

## gfp_leader
## gfp_leader_standard
## gfp_leader_tbb
Leader, or sphere exclusion. Extremely useful tools that implement sphere
exclusion. This is generally the preferred method for dealing with the problem
of having a set or ranked molecules from which a small subset must be selected.
Sort by desirability and run leader.
[gfp_leader](/docs/GFP/gfp_leader.md)

## gfp_lnearneighbours
## gfp_lnearneighbours_standard
Tool designed for comparing a small set of molecules, needles, against a larger 
set, haystack.

```
gfp_lnearneighbours_standard -p needles.gfp -n 100 haystack.gfp.gz > needles.nn
nn2csv needles.nn > needles.nn.csv
```
[gfp_nearneighbours](/docs/GFP/gfp_nearneighbours.md).

## gfp_mcs
An attempt at using the gfp framework for an approximate MCS solution. Needs work,
preferred means of tackling the approximate MCS problem is with [dicer](dicer.md).

## gfp_naive_bayesian
An often powerful modeling method. Make sure that more sophisticated models
actually perform better than this.

## gfp_nearneighbours
Generally use gfp_lnearneighbours instead. This version loads the haystack
into memory and then reads the needles one at a time. Turns out this is the
slow way of doing things.
[gfp_nearneighbours](/docs/GFP/gfp_nearneighbours.md).

## gfp_nearneighbours_single_file
## gfp_nearneighbours_single_file_tbb
Nearest neighbours within a single file - no needles/haystack distinction.

## gfp_pairwise_distances
For given pairs of molecules, return the distances. Useful if you only need
information about specific pairs from a larger set, although you may be
better off extracting the fingerprints of interest from the larger set
first. Or extract the subset of interest and run that alone.

## gfp_profile_activity_by_bits
Extremely useful tool for examining the tendency of bits to differntially appear
across classes, or with different activity values.

## gfp_single_linkage
A simple clustering tool.

## gfp_sparse_to_fixed
Convert sparse fingerprints to fixed form.

## gfp_spread
## gfp_spread_omp
## gfp_spread_standard
Max/min picker. Useful for maximum diversity selections.
[gfp_spread](/docs/GFP/gfp_spread.md)

## gfp_spread_buckets
A bucket constrained version of spread. Molecules are assigned to buckets and
as molecules are selected, selections are forced to uniformly sample
from the buckets.
[gfp_spread](/docs/GFP/gfp_spread.md)

## gfp_standalone
Can do similarity calculations using smiles. Computes only the default
set of fingerprints internally. 

## gfp_svmfp_score
## gfp_svmfp_score_tbb
Infrastructure. Used for scoring svmfp models.

## gfp_to_descriptors
## gfp_to_descriptors_multiple
Convert gfp fingerprint forms to descriptor file format. Always
use `gfp_to_descriptors_multiple`. [gfp_to_descriptors_multiple](/docs/GFP/gfp_to_descriptors_multiple.md).

## gfp_to_svm_lite.v2
Convert gfp to svm-lite format.

## ghose_crippen
Implements the Ghose-Crippen fragment additivity models.

## grease
Identifies large likely hydrophobic sections in molecules.

## grep_by_column
Look for text in certain column(s) only.

## grid_fingerprint
Compute a fingerprint of how a molecule interacts with a fixed grid.

## grid_overlap
Compute a tanomito overlap between two molecular grids.

## hydrophobic_sections
Descriptor generator for sections of molecules likely to be hydrophobic.

## id_chirality
Can identify and enumerate chiral centres. Also useful for identifying
and unmarked and invalid chirality.

## in_database_bdb
An important tool for maintaining a database where the key is a unique smiles.
Databases are built with buildsmidb_bdb.
[in_database](Molecule_Tools/in_database.md)

## in_lilly_database_bdb
## _inventory
Special purpose tools, only relevant inside Lilly.

## isolation_forest
Early implementation of Isolation Forest. The version in sklearn is better today.

## iwbdb_cat
## iwbdb_compare
## iwbdb_delete
## iwbdb_exists
## iwbdb_fetch
## iwbdb_from_tdt
## iwbdb_list
## iwbdb_load
## iwbdb_merge_into_tdt

Command line utilities for BerkeleyDB databases

## iwcut
A version of cut. Understands descriptor files, and writes columns in the order
specified [iwcut](General/iwcut.md)

## iwdemerit
Part of Lilly Medchem Rules.

## iwdescr
Descriptor generator. Generates 250+ mostly interpretable molecular descriptors.
[iwdescr](Molecule_Tools/iwdescr.md)

## iwecfp
Generates EC fingerprints, the -EC to gfp_make.

## iwecfp_database_load
## iwecfp_database_lookup
The synethetic feasibility databases [synthetic feasibility(Molecule_Tools/synthetic_precedent.md).

## iwecfp_intermolecular
An innovative approach to measuring the interaction of a ligand with a protein.

## iwfp
Generates linear fingerprints -FPIW

## iwpathd
Generates T-shaped paths -PATHD. Needs work.

## iwsplit
Split utility. Now Linux split is almost as good, but the naming is better
here and it knows about descriptor files.

## iwstats
Generate statistics for the problem of model prediction - regression only.

## jfilecompare
Compare descriptor files.

## just_columns_with_same_sign
Infrastructure. Some tools need to process features that do not cross zero.

## jwcats
CATS pharmacaphore fetures, also fingerprint generator -CATS

## jwdip
Descriptor generator.

## jwdist
Descriptor generator.

## jwestate
Descriptor generator.

## jwmedv
Descriptor generator.

## jwmolconn
Descriptor generator. Molecular connectivity indices.

## jwmorse
Descriptor generator.

## jwsadb
Descriptor generator. Ideas from publications from Peter Jurs publications.

## kstat_correlated
Identify and remove correlated features

## linear_fingerprint
New version linear fingerprint generator. Still needs work.

## linker_replacement.
Given linkers, two connected fragments, generated by dicer or other tools,
replace interiour parts of a molecule [linker_replacement](Molecule_Tools/linker_replacement.md).

## long_molecules
Identify molecules that have a long max-bond_separator compared to the other
dimesions.

## maccskeys
Descriptor and fingerprint generator.

## make_descriptors.sh
A descriptor generator [make_descriptors](Molecule_Tools/make_descriptors.md).

## make_these_molecules
Given a combinatorial library scheme, make specific molecules from among those possible.

## marvin2gfp
Convert Marvin fingerprint output to gfp form.
## marvin_pka
Post process output from Marvin pKa calculations.

## medchemwizard
Implementation of Abbot paper. The shell wrapper in contrib/bin is set
up to read the reactions from the paper.

## mispredicted
Across a set of train/test splits, examine the predicted files and identify those
molecules that are always predicted poorly. 

## mkfrag
Convert a multi-fragment molecule to separate molecules.

## model_average
Given multiple files containing model predictions, produce an averge prediction.

## mol2qry
Convert a molecule to a substructure query file. A great deal of control is
available control how atoms are translated [mol2qry](Molecule_Tools/mol2qry.md)

## molecular_abstraction
Generate various molecular abstractions [molecular abstraction](/docs/Molecule_Tools/molecular_abstraction.md)

## molecular_grid
Generate a grid around one or more molecules [molecular_grid](Molecule_Tools/molecular_grid.md).

## molecular_merge
Concatentate a file of molecules into one molecule - mostly useful for testing.

## molecular_scaffold
Convert the molecule to scaffold form - molecular_abstraction also does this.

## molecular_transformations
Apply multiple reactions to a molecule.

## molecular_variants
Useful for identifying smaller pharmacaphore equivalents. Driven by external reaction
files that specify transformations.

## molecules_from_reagents

## molecule_subset
Extracts subsets of the atoms in a molecule based on a query.

## msort
## msort_parallel
Sorts files of molecules. Can be very efficient. Also useful for splitting a
set of molecules into non-overlapping sets [msort](Molecule_Tools/msort.md).

## nextdir
Enables cycling through directories in an interactive session.
cd $(nextdir)

## nn2csv
Convert a nearest neighbour file to csv form. Intended as an easier to use nplotnn.

## nn_leader_and_jp
Jarvis Patric clustering and other capabilities.

## nn_merge
Merge multiple nearest neighbour files.

## nn_merge_from_smiles
Mergest nearest neighbour files.

## nn_training
## nn_predictions
Infrastructure. Predictions from nearest neighbours, knn variant.

## nn_single_linkage
Simple clustering

## normalise
Normalise columns in a tabular file.

## notenoughvariance
Identify and remove constant or nearly constnt columns.

## nplotnn
Converts nearest neighbour files to smiles form. Very useful, but can be
complex.

## numbered_smiles
Write a smiles where the isotopic label is one of several molecular properties.
[numbered_smiles](Molecule_Tools/numbered_smiles.md).

## numeric_differences
A difference utility where a threshold for numeric difference can be specified.

## overlapping_fragment_model
A fragment based model that allows atoms to be matched by multiple fragments.

## parallel_nn_search_to_gfp_spread
Multiple individual nearest neighbour searches have been run, and those files

need to be made into a form suitable for the -A option of gfp_spread.
## pharmacophore_2d
Can generate queries for pharmacaphore type searches in 2D.
[pharmacophore_2d](Molecule_Tools/pharmacophore_2d.md).

## ppv
Generates statistics for predictions for a classification model.

## preferred_smiles
Several smiles variants.

## psafp
Toplogical polar surface area

## pubchem_fingerprints
Pubchem Fingerprints

## QED.sh
Command line version for computing QED [QED](Molecule_Tools/QED.md).

## QupKake.sh
Consumes output from atom contribution methods like QupKake and Epok
for pKa computations. [QupKake](Molecule_Tools/QupKake.md).

## r1r2etc
Given a scaffold and substituent points, identify substituents.

## random_fingerprint
Rand Fingerprint generator - useful for testing

## random_geometric_changes
Useful for testing

## random_molecular_permutations
Permutes molecules in chemically sensible ways [random_molecular_permutations](Molecule_Tools/random_molecular_permutations.md).

## random_records
Randomly select records from a file. Preserves the order of the records in
the file, and scalable - unlink `shuf file | head -n 1000`.

## random_smiles
Generate different smiles for a molecule - useful for testing.

## reduced_graph
Generate Reduced Graphs.

## remove_and_label
Remove substituents and label the core atoms.

## remove_matched_atoms
Remove atoms identified by a query.

## retrosynthesis
Infrastructure. Part of reaction informatics.

## rgroup
Yet another R group utilitiy [rgroup](Molecule_Tools/rgroup.md). This is the most supported of the
R group identification tools.

## ring_extraction
## ring_replacement
Infrastructure Ring replacement utility [ring_replacement](Molecule_Tools/ring_replacement.md).

## ring_fingerprint
Generates a fingerprint based on the rings in a molecule.

## ring_size_fingerprint
Generates a fingerprint based on the size of rings.

## ring_substitution
Generates a ring substitution fingerprint.

## ring_trimming
Remove rings and rings that might be part of a fused system.

## rmsigma
Remove outliers from tabular data.

## rotatable_bond_fingerprint
Fingerprint for rotatable bonds

## rotatable_bonds
Compute and filter rotatable bonds [rotatable_bonds](Molecule_Tools/rotatable_bonds.md).

## rule_of_five
Rule of five implementation

## running_average
Running Average

## rxn_fingerprint
## rxn_reverse
## rxn_signature
## rxn_standardize
Infrastructure. Reaction informatics

## rxn_substructure_search
Substructure search reagents or products. Unsophisticated...

## safe_generate
Using the SAFE smiles representation, generate de-novo molecules [SAFE](Molecule_Tools/SAFE.md).

## same_structures

## scaffolds
Generates exhaustive combinations of the rings that comprise a molecular
scaffold [scaffolds](Molecule_Tools/scaffolds.md).

## shuffle_file
Knows how to shuffle some multi-record entities.

## sidechain_switcheroo
Reads a pool of molecules and randomly permutes the sidechains between them.
[sidechain_switcheroo](Molecule_Tools/sidechain_switcheroo.md).

## smi2rings_bdb
Infrastructure: builds a rings database [smi2rings](Molecule_Tools/smi2rings.md).

## smiles_mutation
String mutations on smiles. Mostly generates bad smiles, but sometimes
interesting. More efficient to use random_molecular_permutations.

## sp3_filter
Identify and filter on sp3 character. Can now do the same within iwdescr.

## spearman_rank
Compute Spearman rank for tabular data.

## stratified_samples
Infrastructure. Used when forming train/test splits in svmfp_calibrate.

## substituent_identification
Implements CREM-like functionality [substituent_identification](Molecule_Tools/substituent_identification.md).

## substituent_model

## substitutions

## substructure_match_fraction
Identify the fraction of a molecule implied by matched atoms.

## substructure_mcs

## superimpose_by_matched_atoms
Geometric superimposition based on matched query atoms.

## svmfp_error_vs_distance
Infrastructure. Useful for examining the relationship between model error
and a training set. Now available in gfp_svmfp_score

## svm_lite_to_gfp
Infrastructure. Convert from svm-lite form to gfp.

## tautomer_generation
Not very good tautomer generation. Problem is too hard, do not use. Use
the chemical standardisation functionality instead.

## tcount
Counts lines and columns. The quickest way to ensure that your tabular file
is actually tabular.

## tdt_join
## tdt_sort
## tdt_stats
Utilities for dealing with TDT files - those generated by gfp_make and consumed
by the various gfp_* tools.

## temperature
Infrastructure: computes the MPR fingerprint component

## test_train_split_classification
## test_train_split_cluster
## test_train_split_random
Infrastructure. Used with svmfp_calibrate

## tnass
Implements an idea around dependent fragments - only try one if another
has already matched.

## topotorsion
## topotorsion_fingerprints
Related to Topological Torsion calculations.

## tp1_summarise
## tp_first_pass
Parts of Lilly Medchem Rules

## trxn
Reaction enumeration [trxn](Molecule_Tools/trxn.md)

## tshadow
Descriptor generation - based on 3D structure and projections onto the axes.

## tsmiles
Tester for smiles canonicalisation. Will fail, but failures are interesting to
observe.

## tstandardise
Tester for chemical standardisation. Run across a collection like Chembl, expect
to see a couple of hundred failures. This is an extremely difficult problem.

## tsubstructure
Substructure searching.

## tsubstructure_summarise_hits
Geneate summary reports from tsubstructure.

## tsymmetry
Tester for symmetry, obsolete.

## unbalanced_quotes
Scan a file looking for unbalanced quotes - spans records

## unique_molecules
Discard duplicate molecules [unique_molecules](Molecule_Tools/unique_molecules.md)

## unique_rows
Discard duplicate items in a given column of a tabular file.

## verloop
Infrastructure. Verloop descriptors.

## whatsmissing
In a descriptor file, identify missing values.

## xlogp
Implementation of XLogp logP estimator.

## xray_structure_compare
Compares two structures via their coordintes.
