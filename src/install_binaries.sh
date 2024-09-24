#!/bin/bash

# Copy executables from bazel-bin to another location
# This does all executables in the repo that need to be installed.
# There is no ability to install just a subset.

destination=$1

if [[ -d ${destination} ]] ; then
  true
elif [[ ! -z "${destination}" ]] ; then
  mkdir -l "${destination}"
elif [[ -v LILLYMOL_HOME ]] ; then
  destination="${LILLYMOL_HOME}/bin/${uname}"
else
  echo "Must specify destination directory as argument"
  exit 1
fi
echo "Destination is ${destination}"

# This is called to do the actual copy from bazel-bin/... to a destination.
# Three arguments
#  1. The repository directory name for the source - Molecule_Tools for example
#  2. A bash array of files to copy - fileconv for example
#  3. The destination to where files are copied.

# If the file has not been built, a warning is ussued. No attempt is made
# to rebuild anything.
function copy_files() {
  dir=$1
  local -n to_copy=$2
  destination=$3

  for file in "${to_copy[@]}" ; do
    from="bazel-bin/${dir}/${file}"
    to="${destination}/${file}"
    if [[ ! -s "${from}" ]] ; then
      echo "Source file ${from} not made"
      continue
    fi
    if [[ ! -s ${to} || "${from}" -nt "${to}" ]] ; then
      echo "Copying ${from} to ${to}"
      cp -f ${from} ${to}
    fi
  done
}

# A special version for Go which places the executable in a slightly
# different place
function copy_go_files() {
  dir=$1
  local -n to_copy=$2
  destination=$3

  for file in "${to_copy[@]}" ; do
    from="bazel-bin/${dir}/${file}_/${file}"
    to="${destination}/${file}"
    if [[ ! -s "${from}" ]] ; then
      echo "Source file ${from} not made"
      continue
    fi
    if [[ ! -s ${to} || "${from}" -nt "${to}" ]] ; then
      echo "Copying ${from} to ${to}"
      cp -f ${from} ${to}
    fi
  done
}

# An array of executables to be installed, once for each source directory.
declare -a files=(
        '_fragment_filter'
        'abraham'
        'activity_consistency'
        'align_molecule'
        'alogp'
        'atom_triples'
        'chirality_fingerprint'
        'common_names'
        'dbf'
        'dicer'
        'dicer_to_topological_types'
        'echoqry'
        'echorxn'
        'ec_fingerprint'
        'enumeration'
        'extended_atom_pairs'
        'ez_descriptor'
        'ez_fingerprint'
        'ez_fingerprint_v2'
        'fileconv'
        'fingerprint_substructure'
        'firstatom'
        'get_coordinates'
        'get_substituents'
        'gfp_erg'
        'ghose_crippen'
        'grease'
        'grep_molecule'
        'grid_fingerprint'
        'hydrophobic_sections'
        'id_chirality'
        'iwdemerit'
        'iwdescr'
        'iwecfp'
        'iwecfp_intermolecular'
        'iwfp'
        'iwpathd'
        'jwcats'
        'jwdip'
        'jwdist'
        'jwestate'
        'jwmedv'
        'jwmolconn'
        'jwmorse'
        'jwsadb'
        'linear_fingerprint'
        'long_molecules'
        'maccskeys'
        'make_these_molecules'
        'medchemwizard'
        'minor_changes'
        'mkfrag'
        'mol2qry'
        'mol2SAFE'
        'molecular_abstraction'
        'molecular_grid'
        'molecular_merge'
        'molecular_scaffold'
        'molecular_transformations'
        'molecular_variants'
        'molecule_filter'
        'molecule_subset'
        'molecules_from_reagents'
        'msort'
        'msort_parallel'
        'numbered_smiles'
        'overlapping_fragment_model'
        'parsimonious_set'
        'pharmacophore_2d'
        'preferred_smiles'
        'psafp'
        'pubchem_fingerprints'
        'r1r2etc'
        'random_geometric_changes'
        'random_molecular_permutations'
        'random_smiles'
        'reduced_graph'
        'remove_and_label'
        'remove_matched_atoms'
        'representations'
        'retrosynthesis'
        'rgroup'
        'ring_extraction'
        'ring_fingerprint'
        'ring_replacement'
        'ring_replacement_collate'
        'ring_size_fingerprint'
        'ring_substitution'
        'ring_trimming'
        'rotatable_bond_fingerprint'
        'rotatable_bonds'
        'rule_of_five'
        'rxn_fingerprint'
        'rxn_reverse'
        'rxn_signature'
        'rxn_standardize'
        'rxn_substructure_search'
        'same_structures'
        'smiles_mutation'
        'sp3_filter'
        'substituent_model'
        'substitutions'
        'substructure_match_fraction'
        'substructure_mcs'
        'superimpose_by_matched_atoms'
        'tautomer_generation'
        'temperature'
        'tnass'
        'topotorsion'
        'topotorsion_fingerprints'
        'tp1_summarise'
        'tp_first_pass'
        'trxn'
        'tshadow'
        'tsmiles'
        'tstandardise'
        'tsubstructure'
        'tsubstructure_summarise_hits'
        'tsymmetry'
        'unique_molecules'
        'verloop'
        'xlogp'
        'xray_structure_compare'
)
copy_files "Molecule_Tools" files $destination

files=(
        'buildsmidb_bdb'
        'in_database_bdb'
        'in_lilly_database_bdb'
        'iwecfp_database_load'
        'iwecfp_database_lookup'
        'smi2rings_bdb'
        'substituent_identification'
)
copy_files 'Molecule_Tools_Bdb' files $destination

files=(
        'enough_inventory_bdb'
        'iwbdb_cat'
        'iwbdb_compare'
        'iwbdb_delete'
        'iwbdb_exists'
        'iwbdb_fetch'
        'iwbdb_from_tdt'
        'iwbdb_list'
        'iwbdb_load'
        'iwbdb_merge_into_tdt'
)
copy_files 'BerkeleyDB' files $destination

files=(
        'xgboost_model_evaluate'
)
copy_files 'xgboost' files $destination

files=(
        'grep_sdf'
        'parallel_process_file'
        'no_spaces_in_file_name'
        'regression_to_classification'
        'rxn_reverse'
        'rxnsmiles2smi'
)
copy_go_files 'go' files $destination

files=(
        'complex_chirality'
        'descriptor_file_same_row_order'
        'diff_line_by_line'
        'prediction_bias'
        'rearrange_columns'
)
copy_files 'Obsolete' files $destination

files=(
        'descriptor_file_to_distance_matrix'
)
copy_files 'Obsolete/DistanceMatrix' files $destination

files=(
        'average'
        'bucketise'
        'byte_offset_index'
        'class_label_translation'
        'concat_files'
        'correlate'
        'descriptor_file_cat'
        'descriptor_file_filter'
        'descriptor_file_select_rows'
        'descriptor_file_sort'
        'descriptor_file_to_svm_lite'
        'descriptor_file_transpose'
        'descriptors_to_fingerprint'
        'dfilefilter'
        'dicer_fragments_collate'
        'difference_sort'
        'distribution'
        'exon_correlated'
        'feature_scaling'
        'fetch_sdf'
        'fetch_sdf_quick'
        'fetch_smiles'
        'fetch_smiles_quick'
        'grep_by_column'
        'grid_overlap'
        '_inventory'
        'isolation_forest'
        'iwcut'
        'iwsplit'
        'jfilecompare'
        'just_columns_with_same_sign'
        'kstat_correlated'
        'mispredicted'
        'model_average'
        'nextdir'
        'normalise'
        'notenoughvariance'
        'numeric_differences'
        'nn_single_linkage'
        'ppv'
        'random_column_order'
        'random_records'
        'rmsigma'
        'running_average'
        'shuffle_file'
        'spearman_rank'
        'stratified_samples'
        'svm_lite_to_gfp'
        'svmfp_error_vs_distance'
        'tcount'
        'test_train_split_classification'
        'test_train_split_cluster'
        'test_train_split_random'
        'unbalanced_quotes'
        'unique_rows'
        'whatsmissing'
)
copy_files 'Utilities/General' files $destination

files=(
        'gene_expression_to_proto'
        'gene_expression_nearneighbours'
)
copy_files 'Utilities/GeneExpression' files $destination

files=(
        'descriptor_file_to_01_fingerprints'
        'gfp_add_descriptors'
        'gfp_distance_filter'
        'gfp_distance_matrix'
        'gfp_distance_matrix_iwdm'
        'gfp_flatten_counted'
        'gfp_histogram'
        'gfp_incremental_diversity'
        'gfp_iterative_expansion'
        'gfp_leader'
        'gfp_leader_standard'
        'gfp_standalone'
        'gfp_leader_tbb'
        'gfp_lnearneighbours'
        'gfp_lnearneighbours_standard'
        'gfp_mcs'
        'gfp_naive_bayesian'
        'gfp_nearneighbours'
        'gfp_nearneighbours_single_file'
        'gfp_nearneighbours_single_file_tbb'
        'gfp_pairwise_distances'
        'gfp_profile_activity_by_bits'
        'gfp_single_linkage'
        'gfp_sparse_to_fixed'
        'gfp_spread'
        'gfp_spread_buckets'
        'gfp_spread_omp'
        'gfp_spread_standard'
        'gfp_svmfp_score'
        'gfp_svmfp_score_tbb'
        'gfp_to_descriptors'
        'gfp_to_descriptors_multiple'
        'gfp_to_svm_lite.v2'
        'iwstats'
        'marvin2gfp'
        'nn_leader_and_jp'
        'nn_merge'
        'nn_merge_from_smiles'
        'nn2csv'
        'nplotnn'
        'parallel_nn_search_to_gfp_spread'
        'random_fingerprint'
        'train_test_split_optimise'
)
copy_files 'Utilities/GFP_Tools' files $destination

files=(
        'nn_predictions'
        'nn_training'
        'test_t_test'
)
copy_files 'Utilities/GFP_Knn' files $destination

files=(
        'distance_matrix_activity_difference'
        'distance_matrix_from_distances'
        'distance_matrix_kmedioids'
        'distance_matrix_leader'
        'distance_matrix_nn'
        'distance_matrix_simple_cluster'
        'distance_matrix_spread'
        'distance_matrix_to_distances'
)
copy_files 'Utilities/Distance_Matrix' files $destination

files=(
        'fetch_tdt'
        'fetch_tdt_quick'
        'tdt_join'
        'tdt_sort'
        'tdt_stats'
)
copy_files 'Foundational/iw_tdt' files $destination
