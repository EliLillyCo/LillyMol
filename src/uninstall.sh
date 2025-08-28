#!/bin/bash

# The default build can generate executables that may not be needed.
# This scipt is designed to be run from the Docker container after the
# build in order to remove executables that are not needed for a deployment.

# Some of these are testers that are important. They are not part of
# the GoogleTest based tests, but are nevertheless important, but are not
# run inside the Docker container.

# Others will soon be moved to Obsolete,
# Some have not been used for a while, but may still be useful if a
# need ever recurs.

# This version of the file is tailored to Lilly.

install_dir=$(dirname $0)/../bin/$(uname)
if [[ ! -d ${install_dir} ]] ; then
  echo "Cannot access installation directory ${install_dir}" >&2
  exit 1
fi

declare -a to_remove
to_remove+=('atomic_distance_fingerprint')
to_remove+=('atom_pair_fingerprint')
to_remove+=('descriptor_file_to_mahalanobis_gfp')
to_remove+=('descriptor_file_to_molecular_properties')
to_remove+=('feature_scaling')
to_remove+=('fixed_bit_vector_benchmark')
to_remove+=('geometric_descriptors')
to_remove+=('get_bond_length_distributions')
to_remove+=('get_linkers')
to_remove+=('gfp_difference_fingerprint')
to_remove+=('gfp_group_spread')
to_remove+=('gfp_leader_threaded')
to_remove+=('gfp_leader_opt_v2')
to_remove+=('gfp_leader_standard_threaded')
to_remove+=('gfp_spread_threaded')
to_remove+=('gfp_svmfp_score_pthread')
to_remove+=('gfp_svmfp_score_v2')
to_remove+=('gfp_to_svm_lite_v3')
to_remove+=('gfp_test_train_split')
to_remove+=('gfp_vector_differences')
to_remove+=('ilpd_reduced_graph')
to_remove+=('InformationContent')
to_remove+=('introduction')
to_remove+=('isolation_forest')
to_remove+=('molecule_pca')
to_remove+=('nn_fixed_size_cluster')
to_remove+=('nn_identify_active_outliers')
to_remove+=('nn_id_to_ndx')
to_remove+=('plate_layout')
to_remove+=('ring_closure')
to_remove+=('skeleton')
to_remove+=('smi2linker')
to_remove+=('sparsefp_benchmark')
to_remove+=('structure_database_load')
to_remove+=('structure_database_loookup')
to_remove+=('test_iwdigits')
to_remove+=('test_sparse_bitvector_performance')
to_remove+=('test_train_split_from_leader')
to_remove+=('test_tfdata_record')
to_remove+=('tiwds')
to_remove+=('tsclass')
to_remove+=('tsmiles')
to_remove+=('tspassbyref')
to_remove+=('tstandardise')
to_remove+=('unique_molecules_sorted')

for file in "${to_remove[@]}" ; do
  fname="${install_dir}/${file}"
  if [[ ! -s ${fname} ]] ; then
    continue
  fi
  echo "removing ${fname}"
  rm -f "${fname}"
done


# Not yet decided what to do about these
# multiple_reactions
# grid_proximity
# dicer_fragment_lookup_bdb
# dicer_complementary_fragments_collate
# model_uncertainty
# positional_analogue_scanning
# gfp_evaluate_clustering
# spatial_replacement
# dicer_fragment_overlap
# dicer_to_complemenent_db
# excise_3d_overlap
# gfp_distance_filter_standard
# gfp_nearneighbours_threaded
# ligand_protein_interactions
