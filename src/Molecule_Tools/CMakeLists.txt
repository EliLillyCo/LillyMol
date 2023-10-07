project(molecule_tools)
find_package(Protobuf REQUIRED)
include_directories(SYSTEM ${Protobuf_INCLUDE_DIRS})

set(third_party_libs
  ${IW_EXTRA_LIBS}
)
# Libraries that all executables in this dir will need
set(common_libs
  cmdline
  accumulator
  moleculeio
  iwmolecule
  iwmisc
  iwbits
  tfdatarecord
  protobuf::libprotobuf
  crc32c
  iwstring_data_source
  string_data_source
  iwre2
  iwstring
  xmlParser
)

set(molecule_tools Molecule_Tools)
file(GLOB proto_files ${CMAKE_CURRENT_SOURCE_DIR}/*.proto)
foreach(fname ${proto_files})
  get_filename_component(dir ${fname} DIRECTORY)
  get_filename_component(stem ${fname} NAME_WE)
  set(PROTO_HDRS "${PROTO_HDRS};${dir}/${stem}.pb.h")
  set(PROTO_SRCS "${PROTO_SRCS};${dir}/${stem}.pb.cc")
  set(PROTO_INPUTS "${PROTO_INPUTS};${molecule_tools}/${stem}.proto")
endforeach()

add_custom_command(
  OUTPUT ${PROTO_HDRS} ${PROTO_SRCS}
  COMMAND ${Protobuf_PROTOC_EXECUTABLE}
  ARGS   --cpp_out . -I .
         ${PROTO_INPUTS}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
  COMMENT "Running protoc on ${PROTO_INPUTS}"
)

# Libraries
 
set(maccskeys_fn5_srcs
  maccskeys_fn5.cc
)
add_library(maccskeys_fn5_lib ${maccskeys_fn5_srcs})
target_include_directories(maccskeys_fn5_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(partial_symmetry_srcs
  partial_symmetry.cc
)
add_library(partial_symmetry_lib ${partial_symmetry_srcs})
target_include_directories(partial_symmetry_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(linear_path_fingerprint_srcs
  linear_path_fingerprint.cc
)
add_library(linear_path_fingerprint_lib ${linear_path_fingerprint_srcs})
target_include_directories(linear_path_fingerprint_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(tp_first_pass_lib_srcs
  tp_first_pass_lib.cc
)
add_library(tp_first_pass_lib ${tp_first_pass_lib_srcs})
target_include_directories(tp_first_pass_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# activity_consistency
set(activity_consistency_srcs
  activity_consistency.cc
)
set(activity_consistency_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(activity_consistency ${activity_consistency_srcs})
target_link_libraries(activity_consistency LINK_PUBLIC ${activity_consistency_libs})

# atom_pair_fingerprint
set(atom_pair_fingerprint_srcs
  atom_pair_fingerprint_main.cc
)
set(atom_pair_fingerprint_libs
  sparse_fp_creator
  atom_pair_fingerprint_lib
  ${common_libs}
  ${third_party_libs}
)
add_executable(atom_pair_fingerprint ${atom_pair_fingerprint_srcs})
target_link_libraries(atom_pair_fingerprint LINK_PUBLIC ${atom_pair_fingerprint_libs})

# buildsmidb_leveldb
set(buildsmidb_leveldb_srcs
  buildsmidb_leveldb.cc
  molecule_database_options.cc
  ${PROTO_SRCS}
)
set(buildsmidb_leveldb_libs
  ${common_libs}
  ${third_party_libs}
  leveldb
)
add_executable(buildsmidb_leveldb ${buildsmidb_leveldb_srcs})
target_link_libraries(buildsmidb_leveldb LINK_PUBLIC ${buildsmidb_leveldb_libs})

# common_names
set(common_names_srcs
  common_names.cc
)
set(common_names_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(common_names ${common_names_srcs})
target_link_libraries(common_names LINK_PUBLIC ${common_names_libs})

# dbf
set(dbf_srcs
  dbf.cc
)
set(dbf_libs
  iwhistogram
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(dbf ${dbf_srcs})
target_link_libraries(dbf LINK_PUBLIC ${dbf_libs})

# dicer
set(dicer_srcs
  dicer.cc
)
set(dicer_libs
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(dicer ${dicer_srcs})
target_link_libraries(dicer LINK_PUBLIC ${dicer_libs})

# ec_fingerprint
set(ec_fingerprint_srcs
  ec_fingerprint_main.cc
)
set(ec_fingerprint_libs
  sparse_fp_creator
  ec_fingerprint_lib
  ${common_libs}
  ${third_party_libs}
)
add_executable(ec_fingerprint ${ec_fingerprint_srcs})
target_link_libraries(ec_fingerprint LINK_PUBLIC ${ec_fingerprint_libs})

# echoqry
set(echoqry_srcs
  echoqry.cc
)
set(echoqry_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(echoqry ${echoqry_srcs})
target_link_libraries(echoqry LINK_PUBLIC ${echoqry_libs})

# fileconv
set(fileconv_srcs
  fileconv.cc
  fix_structures.cc
  fix_structures.h
  do_remove_duplicate_fragments.cc
  do_remove_duplicate_fragments.h
  known_fragment_data.cc
  known_fragment_data.h
)

set(fileconv_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(fileconv ${fileconv_srcs})
target_link_libraries(fileconv LINK_PUBLIC ${fileconv_libs})

# geometric_descriptors
set(geometric_descriptors_srcs
  geometric_descriptors.cc
)
set(geometric_descriptors_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(geometric_descriptors ${geometric_descriptors_srcs})
target_link_libraries(geometric_descriptors LINK_PUBLIC ${geometric_descriptors_libs})

# hydrophobic_sections
set(hydrophobic_sections_srcs
  hydrophobic_sections.cc
)
set(hydrophobic_sections_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(hydrophobic_sections ${hydrophobic_sections_srcs})
target_link_libraries(hydrophobic_sections LINK_PUBLIC ${hydrophobic_sections_libs})

# in_database_leveldb
set(in_database_leveldb_srcs
  in_database_leveldb.cc
  molecule_database_options.cc
  ${PROTO_SRCS}
)
set(in_database_leveldb_libs
  ${common_libs}
  ${third_party_libs}
  leveldb
)
add_executable(in_database_leveldb ${in_database_leveldb_srcs})
target_link_libraries(in_database_leveldb LINK_PUBLIC ${in_database_leveldb_libs})

# iwdescr
set(iwdescr_srcs
  iwdescr.cc
  nvrtspsa.cc
  qry_wcharge.cc
)
set(iwdescr_libs
  partial_symmetry_lib
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(iwdescr ${iwdescr_srcs})
target_link_libraries(iwdescr LINK_PUBLIC ${iwdescr_libs})

# iwecfp
set(iwecfp_srcs
  iwecfp.cc
)
set(iwecfp_libs
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(iwecfp ${iwecfp_srcs})
target_link_libraries(iwecfp LINK_PUBLIC ${iwecfp_libs})

# iwfp
set(iwfp_srcs
  iwfp.cc
)
set(iwfp_libs
  sparse_fp_creator
  ${common_libs}
  iwmfingerprint_lib
  ${third_party_libs}
)
add_executable(iwfp ${iwfp_srcs})
target_link_libraries(iwfp LINK_PUBLIC ${iwfp_libs})

# linear_fingerprint
set(linear_fingerprint_srcs
  linear_fingerprint_main.cc
)
set(linear_fingerprint_libs
  sparse_fp_creator
  linear_fingerprint_lib
  ${common_libs}
  ${third_party_libs}
)
add_executable(linear_fingerprint ${linear_fingerprint_srcs})
target_link_libraries(linear_fingerprint LINK_PUBLIC ${linear_fingerprint_libs})

# maccskeys_fn5
set(maccskeys_fn5_srcs
  maccskeys_fn5_main.cc
)
set(maccskeys_fn5_libs
  sparse_fp_creator
  maccskeys_fn5_lib
  ${common_libs}
  ${third_party_libs}
)
add_executable(maccskeys_fn5 ${maccskeys_fn5_srcs})
target_link_libraries(maccskeys_fn5 LINK_PUBLIC ${maccskeys_fn5_libs})

# make_these_molecules
set(make_these_molecules_srcs
  make_these_molecules.cc
)
set(make_these_molecules_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(make_these_molecules ${make_these_molecules_srcs})
target_link_libraries(make_these_molecules LINK_PUBLIC ${make_these_molecules_libs})

# mol2qry
set(mol2qry_srcs
  mol2qry.cc
)
set(mol2qry_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(mol2qry ${mol2qry_srcs})
target_link_libraries(mol2qry LINK_PUBLIC ${mol2qry_libs})

# molecular_abstractions
set(molecular_abstractions_srcs
  mav2.cc
  molecular_abstraction_functions.cc
  molecular_abstraction_specifications.cc
)
set(molecular_abstractions_libs
  iwmfingerprint_lib
  ${common_libs}
  ${third_party_libs}
)
add_executable(molecular_abstractions ${molecular_abstractions_srcs})
target_link_libraries(molecular_abstractions LINK_PUBLIC ${molecular_abstractions_libs})

# molecular_scaffold
set(molecular_scaffold_srcs
  molecular_scaffold.cc
)
set(molecular_scaffold_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(molecular_scaffold ${molecular_scaffold_srcs})
target_link_libraries(molecular_scaffold LINK_PUBLIC ${molecular_scaffold_libs})

# molecular_transformations
set(molecular_transformations_srcs
  molecular_transformations.cc
)
set(molecular_transformations_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(molecular_transformations ${molecular_transformations_srcs})
target_link_libraries(molecular_transformations LINK_PUBLIC ${molecular_transformations_libs})

# molecular_variants
set(molecular_variants_srcs
  molecular_variants.cc
)
set(molecular_variants_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(molecular_variants ${molecular_variants_srcs})
target_link_libraries(molecular_variants LINK_PUBLIC ${molecular_variants_libs})

# molecule_subset
set(molecule_subset_srcs
  molecule_subset.cc
)
set(molecule_subset_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(molecule_subset ${molecule_subset_srcs})
target_link_libraries(molecule_subset LINK_PUBLIC ${molecule_subset_libs})

# msort
set(msort_srcs
  msort.cc
)
set(msort_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(msort ${msort_srcs})
target_link_libraries(msort LINK_PUBLIC ${msort_libs})

# pharmacaphore_2d
set(pharmacaphore_2d_srcs
  pharmacaphore_2d.cc
)
set(pharmacaphore_2d_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(pharmacaphore_2d ${pharmacaphore_2d_srcs})
target_link_libraries(pharmacaphore_2d LINK_PUBLIC ${pharmacaphore_2d_libs})

# preferred_smiles
set(preferred_smiles_srcs
  preferred_smiles.cc
)
set(preferred_smiles_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(preferred_smiles ${preferred_smiles_srcs})
target_link_libraries(preferred_smiles LINK_PUBLIC ${preferred_smiles_libs})

# random_molecular_permutations
set(random_molecular_permutations_srcs
  random_molecular_permutations.cc
  set_of_target_molecules.cc
  set_of_target_molecules.h
  swap_kekule_forms.cc
)
set(random_molecular_permutations_libs
  iwreaction
  gfp
  iw_tdt
  gfp_standard
  linear_path_fingerprint_lib
  maccskeys_fn5_lib
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(random_molecular_permutations ${random_molecular_permutations_srcs})
target_link_libraries(random_molecular_permutations LINK_PUBLIC ${random_molecular_permutations_libs})

# random_smiles
set(random_smiles_srcs
  random_smiles.cc
)
set(random_smiles_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(random_smiles ${random_smiles_srcs})
target_link_libraries(random_smiles LINK_PUBLIC ${random_smiles_libs})

# remove_and_label
set(remove_and_label_srcs
  remove_and_label.cc
)
set(remove_and_label_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(remove_and_label ${remove_and_label_srcs})
target_link_libraries(remove_and_label LINK_PUBLIC ${remove_and_label_libs})

# retrosynthesis
set(retrosynthesis_srcs
  retrosynthesis.cc
)
set(retrosynthesis_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(retrosynthesis ${retrosynthesis_srcs})
target_link_libraries(retrosynthesis LINK_PUBLIC ${retrosynthesis_libs})

# rgroup
set(rgroup_srcs
  rgroup.cc
)
set(rgroup_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(rgroup ${rgroup_srcs})
target_link_libraries(rgroup LINK_PUBLIC ${rgroup_libs})

# ring_substitution
set(ring_substitution_srcs
  ring_substitution.cc
)
set(ring_substitution_libs
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(ring_substitution ${ring_substitution_srcs})
target_link_libraries(ring_substitution LINK_PUBLIC ${ring_substitution_libs})

# ring_extraction
set(ring_extraction_srcs
  ring_extraction.cc
  ring_ext_rep.cc
  reaction_duplicate.cc
)
set(ring_extraction_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(ring_extraction ${ring_extraction_srcs})
target_link_libraries(ring_extraction LINK_PUBLIC ${ring_extraction_libs})

# ring_replacement
set(ring_replacement_srcs
  ring_replacement.cc
  ring_ext_rep.cc
  reaction_duplicate.cc
)
set(ring_replacement_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(ring_replacement ${ring_replacement_srcs})
target_link_libraries(ring_replacement LINK_PUBLIC ${ring_replacement_libs})

# ring_trimming
set(ring_trimming_srcs
  ring_trimming.cc
)
set(ring_trimming_libs
  iwreaction
  ${common_libs}
  ${third_party_libs}
)
add_executable(ring_trimming ${ring_trimming_srcs})
target_link_libraries(ring_trimming LINK_PUBLIC ${ring_trimming_libs})

# rotatable_bonds
set(rotatable_bonds_srcs
  rotatable_bonds.cc
)
set(rotatable_bonds_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(rotatable_bonds ${rotatable_bonds_srcs})
target_link_libraries(rotatable_bonds LINK_PUBLIC ${rotatable_bonds_libs})

# rotatable_bond_fingerprint
set(rotatable_bond_fingerprint_srcs
  rotatable_bond_fingerprint.cc
)
set(rotatable_bond_fingerprint_libs
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(rotatable_bond_fingerprint ${rotatable_bond_fingerprint_srcs})
target_link_libraries(rotatable_bond_fingerprint LINK_PUBLIC ${rotatable_bond_fingerprint_libs})

# rule_of_five
set(rule_of_five_srcs
  rule_of_five.cc
  nvrtspsa.cc
)
set(rule_of_five_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(rule_of_five ${rule_of_five_srcs})
target_link_libraries(rule_of_five LINK_PUBLIC ${rule_of_five_libs})

# rxn_fingerprint
set(rxn_fingerprint_srcs
  rxn_fingerprint.cc
  extended_connectivity_fp.cc
)
set(rxn_fingerprint_libs
  iwreaction
  linear_path_fingerprint_lib
  iwmfingerprint_lib
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(rxn_fingerprint ${rxn_fingerprint_srcs})
target_link_libraries(rxn_fingerprint LINK_PUBLIC ${rxn_fingerprint_libs})

# rxn_reverse
set(rxn_reverse_srcs
  rxn_reverse.cc
)
set(rxn_reverse_libs
  iwreaction
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(rxn_reverse ${rxn_reverse_srcs})
target_link_libraries(rxn_reverse LINK_PUBLIC ${rxn_reverse_libs})

# rxn_signature
set(rxn_signature_srcs
  rxn_signature.cc
)
set(rxn_signature_libs
  iwreaction
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(rxn_signature ${rxn_signature_srcs})
target_link_libraries(rxn_signature LINK_PUBLIC ${rxn_signature_libs})

# rxn_substructure_search
set(rxn_substructure_search_srcs
  rxn_substructure_search.cc
)
set(rxn_substructure_search_libs
  iwreaction
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(rxn_substructure_search ${rxn_substructure_search_srcs})
target_link_libraries(rxn_substructure_search LINK_PUBLIC ${rxn_substructure_search_libs})

# rxn_standardize
set(rxn_standardize_srcs
  rxn_standardize.cc
)
set(rxn_standardize_libs
  iwreaction
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(rxn_standardize ${rxn_standardize_srcs})
target_link_libraries(rxn_standardize LINK_PUBLIC ${rxn_standardize_libs})

# smiles_mutation
set(smiles_mutation_srcs
  smiles_mutation.cc
  random_reactions.cc
)
set(smiles_mutation_libs
  iwreaction
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(smiles_mutation ${smiles_mutation_srcs})
target_link_libraries(smiles_mutation LINK_PUBLIC ${smiles_mutation_libs})

# sp3_filter
set(sp3_filter_srcs
  sp3_filter.cc
)
set(sp3_filter_libs
  iwreaction
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(sp3_filter ${sp3_filter_srcs})
target_link_libraries(sp3_filter LINK_PUBLIC ${sp3_filter_libs})

# tautomer_generation
set(tautomer_generation_srcs
  tautomer_generation.cc
)
set(tautomer_generation_libs
  iwreaction
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(tautomer_generation ${tautomer_generation_srcs})
target_link_libraries(tautomer_generation LINK_PUBLIC ${tautomer_generation_libs})

# temperature
set(temperature_srcs
  temperature.cc
)
set(temperature_libs
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(temperature ${temperature_srcs})
target_link_libraries(temperature LINK_PUBLIC ${temperature_libs})

# tnass
set(tnass_srcs
  tnass.cc
  nass.cc
)
set(tnass_libs
  sparse_fp_creator
 ${common_libs}
  ${third_party_libs}
)
add_executable(tnass ${tnass_srcs})
target_link_libraries(tnass LINK_PUBLIC ${tnass_libs})

# topological_torsion
set(topological_torsion_srcs
  topological_torsion_main.cc
  topological_torsion.cc
)
set(topological_torsion_libs
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(topological_torsion ${topological_torsion_srcs})
target_link_libraries(topological_torsion LINK_PUBLIC ${topological_torsion_libs})

# tp_first_pass
set(tp_first_pass_srcs
  tp_first_pass.cc
)
set(tp_first_pass_libs
  ${common_libs}
  ${third_party_libs}
  tp_first_pass_lib
)
add_executable(tp_first_pass ${tp_first_pass_srcs})
target_link_libraries(tp_first_pass LINK_PUBLIC ${tp_first_pass_libs})

# trxn
set(trxn_srcs
  trxn.cc
)
set(trxn_libs
  iwreaction
  iw_tdt
  ${common_libs}
  ${third_party_libs}
)
add_executable(trxn ${trxn_srcs})
target_link_libraries(trxn LINK_PUBLIC ${trxn_libs})

# tsmiles
set(tsmiles_srcs
  tsmiles.cc
)
set(tsmiles_libs
  ${common_libs}
  ${third_party_libs}
)
add_executable(tsmiles ${tsmiles_srcs})
target_link_libraries(tsmiles LINK_PUBLIC ${tsmiles_libs})

# tsubstructure
set(tsubstructure_srcs
  tsubstructure.cc
  tsubstructure_fp.cc
)
set(tsubstructure_libs
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(tsubstructure ${tsubstructure_srcs})
target_link_libraries(tsubstructure LINK_PUBLIC ${tsubstructure_libs})

# unique_molecules
set(unique_molecules_srcs
  unique_molecules.cc
)
set(unique_molecules_libs
  iw_tdt
  iwreaction
  sparse_fp_creator
  ${common_libs}
  ${third_party_libs}
)
add_executable(unique_molecules ${unique_molecules_srcs})
target_link_libraries(unique_molecules LINK_PUBLIC ${unique_molecules_libs})

set(installable
  fileconv
  tsubstructure
  trxn
  unique_molecules
  common_names
  msort
  in_database_leveldb
  buildsmidb_leveldb
  molecular_abstractions
)

install(
  TARGETS ${installable}
  DESTINATION bin
)
