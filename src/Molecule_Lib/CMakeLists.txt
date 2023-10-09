project(iwmolecule)
find_package(Protobuf REQUIRED)
include_directories(SYSTEM ${Protobuf_INCLUDE_DIRS})

set(Protobuf_IMPORT_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

#set(proto_src
#   ${CMAKE_CURRENT_SOURCE_DIR}/geometric_constraints.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/substructure.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/toggle_kekule_form.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/mol2graph.proto
#)
#protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS ${proto_src})

set(molecule_lib Molecule_Lib)

file(GLOB proto_files ${CMAKE_CURRENT_SOURCE_DIR}/*.proto)
foreach(fname ${proto_files})
  get_filename_component(dir ${fname} DIRECTORY)
  get_filename_component(stem ${fname} NAME_WE)
  set(PROTO_HDRS "${PROTO_HDRS};${dir}/${stem}.pb.h")
  set(PROTO_SRCS "${PROTO_SRCS};${dir}/${stem}.pb.cc")
  set(PROTO_INPUTS "${PROTO_INPUTS};${molecule_lib}/${stem}.proto")
endforeach()

#set(proto_files
#   ${CMAKE_CURRENT_SOURCE_DIR}/geometric_constraints.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/substructure.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/toggle_kekule_form.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/mol2graph.proto
#)

add_custom_command(
  OUTPUT ${PROTO_HDRS} ${PROTO_SRCS}
  COMMAND ${Protobuf_PROTOC_EXECUTABLE}
  ARGS   --cpp_out . -I .
         ${PROTO_INPUTS}
  COMMENT "Running protoc on ${PROTO_INPUTS}"
# COMMENT "Running ${CMAKE_CURRENT_SOURCE_DIR}"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/allowed_elements.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/aromatic.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/atom.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/atom_alias.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/atom_typing.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bond.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bond_list.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/careful_frag.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/charge_assigner.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/charge_calculation.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/charmm.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/chiral_centre.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cif.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cis_trans_bond.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/coordinates.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/coordinate_box.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/count_atoms_in_smiles.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/csubstructure.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dihedral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/donor_acceptor.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/element.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/element_hits_needed.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ematch.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/etrans.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/frag.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/geometric_constraints.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ISIS_Atom_List.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/is_actually_chiral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwrcb.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qry_wstats.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwrnm.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwsubstructure.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mdl_atom_record.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mdl_file_data.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mdl_molecule.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/misc2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mol2graph.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mol2graph_proto.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/molecule.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/molecule_to_query.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/moleculeb.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/moleculed.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/moleculeh.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/moleculer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/molecule_smarts.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mpr.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/numass.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ostream_and_type.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/output.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/parse_smarts_tmp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/path.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/path_around_ring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/path_scoring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/pearlman.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ring_bond_iterator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rmele.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rotbond_common.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rwsubstructure.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rwsubstructure_proto.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/set_of_atoms.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/smiles.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/smiles_support.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_atom_env.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_bonds.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_ring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_ring_base.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_ring_sys.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_chiral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_results.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/symm_class_can_rank.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/standardise.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_a.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_nmab.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_spec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_env.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/symmetry.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/target.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/temp_detach_atoms.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/toggle_kekule_form.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tokenise_atomic_smarts.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/unique.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/u3b.c
  ${CMAKE_CURRENT_SOURCE_DIR}/_substructure.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_qry_wstats.cc
  ${PROTO_SRCS}
)
add_library(iwmolecule ${srcs})

set(moleculeio_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/charmm.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cif.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/csv.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/inchi_dummy.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/marvin.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mdl.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mdl_v30.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mmod.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mrk.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/msi.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/pdb.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rwmolecule.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/smi.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tripos.cc
)
add_library(moleculeio ${moleculeio_srcs})

set(atom_pair_fingerprint_srcs
  atom_pair_fingerprint.cc
)
add_library(atom_pair_fingerprint_lib ${atom_pair_fingerprint_srcs})
target_include_directories(atom_pair_fingerprint_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(ec_fingerprint_srcs
  ec_fingerprint.cc
)
add_library(ec_fingerprint_lib ${ec_fingerprint_srcs})
target_include_directories(ec_fingerprint_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
 
set(linear_fingerprint_srcs
  linear_fingerprint.cc
)
add_library(linear_fingerprint_lib ${linear_fingerprint_srcs})
target_include_directories(linear_fingerprint_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
 
set(iwmfingerprint_srcs
  iwmfingerprint.cc
)
add_library(iwmfingerprint_lib ${iwmfingerprint_srcs})
target_include_directories(iwmfingerprint_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(iwreaction_srcs
  iwreaction.cc
  iwreaction2.cc
  reaction_iterator.cc
  reaction_from_smirks.cc
  reaction_match_conditions.cc
  reaction_proto.cc
  rxnfile.cc
  rxnfile2.cc
  rxnfile3.cc
  rxn_3d_replace.cc
  rxn_bond_angle.cc
  rxn_bond_length.cc
  rxn_dihedral_angle.cc
  rxn_matched_sors.cc
)
add_library(iwreaction ${iwreaction_srcs})
target_include_directories(iwreaction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  
enable_testing()


if(IWTEST)
# tests

set(libs_for_tests
  moleculeio
  iwmolecule
  cmdline
  gmock
  gtest_main
  iwstring_data_source
  string_data_source
  xmlParser
  iwmisc
  iwbits
  iwre2
  iwstring
  protobuf::libprotobuf
  tfdatarecord
  crc32c
  ${IW_EXTRA_LIBS}
)

# aromatic_test
set(aromatic_test_srcs
  aromatic_test.cc
)
set(aromatic_test_libs
  ${libs_for_tests}
)
add_executable(aromatic_test ${aromatic_test_srcs})
target_link_libraries(aromatic_test LINK_PUBLIC ${aromatic_test_libs})

# atom_test
set(atom_test_srcs
  atom_test.cc
)
set(atom_test_libs
  ${libs_for_tests}
)
add_executable(atom_test ${atom_test_srcs})
target_link_libraries(atom_test LINK_PUBLIC ${atom_test_libs})

# coordinate_box_test
set(coordinate_box_test_srcs
  coordinate_box_test.cc
)
set(coordinate_box_test_libs
  ${libs_for_tests}
)
add_executable(coordinate_box_test ${coordinate_box_test_srcs})
target_link_libraries(coordinate_box_test LINK_PUBLIC ${coordinate_box_test_libs})

# geometric_constraints_test
set(geometric_constraints_test_srcs
  geometric_constraints_test.cc
)
set(geometric_constraints_test_libs
  ${libs_for_tests}
)
add_executable(geometric_constraints_test ${geometric_constraints_test_srcs})
target_link_libraries(geometric_constraints_test LINK_PUBLIC ${geometric_constraints_test_libs})

# linear_fingerprint_test
set(linear_fingerprint_test_srcs
  linear_fingerprint_test.cc
)
set(linear_fingerprint_test_libs
  linear_fingerprint_lib
  sparse_fp_creator
  ${libs_for_tests}
)
add_executable(linear_fingerprint_test ${linear_fingerprint_test_srcs})
target_link_libraries(linear_fingerprint_test LINK_PUBLIC ${linear_fingerprint_test_libs})

# misc2_test
set(misc2_test_srcs
  misc2_test.cc
)
set(misc2_test_libs
  ${libs_for_tests}
)
add_executable(misc2_test ${misc2_test_srcs})
target_link_libraries(misc2_test LINK_PUBLIC ${misc2_test_libs})

# ec_fingerprint_test
set(ec_fingerprint_test_srcs
  ec_fingerprint_test.cc
)
set(ec_fingerprint_test_libs
  ec_fingerprint_lib
  sparse_fp_creator
  ${libs_for_tests}
)
add_executable(ec_fingerprint_test ${ec_fingerprint_test_srcs})
target_link_libraries(ec_fingerprint_test LINK_PUBLIC ${ec_fingerprint_test_libs})

# iwreaction_test
set(iwreaction_test_srcs
  iwreaction_test.cc
)
set(iwreaction_test_libs
  iwreaction
  ${libs_for_tests}
)
add_executable(iwreaction_test ${iwreaction_test_srcs})
target_link_libraries(iwreaction_test LINK_PUBLIC ${iwreaction_test_libs})

# molecule_each_index_test
set(molecule_each_index_test_srcs
  molecule_each_index_test.cc
)
set(molecule_each_index_test_libs
  iwreaction
  ${libs_for_tests}
)
add_executable(molecule_each_index_test ${molecule_each_index_test_srcs})
target_link_libraries(molecule_each_index_test LINK_PUBLIC ${molecule_each_index_test_libs})

# moleculeh_test
set(moleculeh_test_srcs
  moleculeh_test.cc
)
set(moleculeh_test_libs
  ${libs_for_tests}
)
add_executable(moleculeh_test ${moleculeh_test_srcs})
target_link_libraries(moleculeh_test LINK_PUBLIC ${moleculeh_test_libs})

# ring_bond_count_test
set(ring_bond_count_test_srcs
  ring_bond_count_test.cc
)
set(ring_bond_count_test_libs
  iwreaction
  ${libs_for_tests}
)
add_executable(ring_bond_count_test ${ring_bond_count_test_srcs})
target_link_libraries(ring_bond_count_test LINK_PUBLIC ${ring_bond_count_test_libs})

# rxnfile_test
set(rxnfile_test_srcs
  rxnfile_test.cc
)
set(rxnfile_test_libs
  iwreaction
  ${libs_for_tests}
)
add_executable(rxnfile_test ${rxnfile_test_srcs})
target_link_libraries(rxnfile_test LINK_PUBLIC ${rxnfile_test_libs})

# rxnfile2_test
set(rxnfile2_test_srcs
  rxnfile2_test.cc
)
set(rxnfile2_test_libs
  iwreaction
  ${libs_for_tests}
)
add_executable(rxnfile2_test ${rxnfile2_test_srcs})
target_link_libraries(rxnfile2_test LINK_PUBLIC ${rxnfile2_test_libs})

# smiles_test
set(smiles_test_srcs
  smiles_test.cc
)
set(smiles_test_libs
  iwreaction
  ${libs_for_tests}
)
add_executable(smiles_test ${smiles_test_srcs})
target_link_libraries(smiles_test LINK_PUBLIC ${smiles_test_libs})

# space_vector_test
set(space_vector_test_srcs
  space_vector_test.cc
)
set(space_vector_test_libs
  iwreaction
  ${libs_for_tests}
)
add_executable(space_vector_test ${space_vector_test_srcs})
target_link_libraries(space_vector_test LINK_PUBLIC ${space_vector_test_libs})

# substructure_env_test
set(substructure_env_test_srcs
  substructure_env_test.cc
)
set(substructure_env_test_libs
  ${libs_for_tests}
)
add_executable(substructure_env_test ${substructure_env_test_srcs})
target_link_libraries(substructure_env_test LINK_PUBLIC ${substructure_env_test_libs})

# substructure_nmab_test
set(substructure_nmab_test_srcs
  substructure_nmab_test.cc
)
set(substructure_nmab_test_libs
  ${libs_for_tests}
)
add_executable(substructure_nmab_test ${substructure_nmab_test_srcs})
target_link_libraries(substructure_nmab_test LINK_PUBLIC ${substructure_nmab_test_libs})

# substructure_ring_test
set(substructure_ring_test_srcs
  substructure_ring_test.cc
)
set(substructure_ring_test_libs
  ${libs_for_tests}
)
add_executable(substructure_ring_test ${substructure_ring_test_srcs})
target_link_libraries(substructure_ring_test LINK_PUBLIC ${substructure_ring_test_libs})

# substructure_ring_system_test
set(substructure_ring_system_test_srcs
  substructure_ring_system_test.cc
)
set(substructure_ring_system_test_libs
  ${libs_for_tests}
)
add_executable(substructure_ring_system_test ${substructure_ring_system_test_srcs})
target_link_libraries(substructure_ring_system_test LINK_PUBLIC ${substructure_ring_system_test_libs})

# substructure_sepma_test
set(substructure_sepma_test_srcs
  substructure_sepma_test.cc
)
set(substructure_sepma_test_libs
  ${libs_for_tests}
)
add_executable(substructure_sepma_test ${substructure_sepma_test_srcs})
target_link_libraries(substructure_sepma_test LINK_PUBLIC ${substructure_sepma_test_libs})

# substructure_single_query_test
set(substructure_single_query_test_srcs
  substructure_single_query_test.cc
)
set(substructure_single_query_test_libs
  ${libs_for_tests}
)
add_executable(substructure_single_query_test ${substructure_single_query_test_srcs})
target_link_libraries(substructure_single_query_test LINK_PUBLIC ${substructure_single_query_test_libs})

# substructure_smarts_test
set(substructure_smarts_test_srcs
  substructure_smarts_test.cc
)
set(substructure_smarts_test_libs
  ${libs_for_tests}
)
add_executable(substructure_smarts_test ${substructure_smarts_test_srcs})
target_link_libraries(substructure_smarts_test LINK_PUBLIC ${substructure_smarts_test_libs})

# substructure_spec_test
set(substructure_spec_test_srcs
  substructure_spec_test.cc
)
set(substructure_spec_test_libs
  ${libs_for_tests}
)
add_executable(substructure_spec_test ${substructure_spec_test_srcs})
target_link_libraries(substructure_spec_test LINK_PUBLIC ${substructure_spec_test_libs})

# substructure_test
set(substructure_test_srcs
  substructure_test.cc
)
set(substructure_test_libs
  ${libs_for_tests}
)
add_executable(substructure_test ${substructure_test_srcs})
target_link_libraries(substructure_test LINK_PUBLIC ${substructure_test_libs})

# symmetry_test
set(symmetry_test_srcs
  symmetry_test.cc
)
set(symmetry_test_libs
  ${libs_for_tests}
)
add_executable(symmetry_test ${symmetry_test_srcs})
target_link_libraries(symmetry_test LINK_PUBLIC ${symmetry_test_libs})

# unique_test
set(unique_test_srcs
  unique_test.cc
)
set(unique_test_libs
  ${libs_for_tests}
)
add_executable(unique_test ${unique_test_srcs})
target_link_libraries(unique_test LINK_PUBLIC ${unique_test_libs})
endif(IWTEST)
