load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")

package(default_visibility = ["//visibility:public"])

#filegroup(
#    name = "charge_assigner_data",
#    srcs = glob([":charge_assigner/**"]),
#)

#genrule(
#  name = 'gfortran_iw',
#  srcs = [
#    'u3b.f',
#  ],
#  outs = ['u3b.o'],
#  cmd = 'gfortran -O3 -c -o $(OUTS) u3b.f'
#)

proto_library(
    name = "atom_type_ext_proto",
    srcs = [
        "atom_type_ext.proto",
    ],
)

cc_proto_library(
    name = "atom_type_ext_cc_proto",
    deps = [
        ":atom_type_ext_proto",
    ],
)

py_proto_library(
    name = "atom_type_ext_py_proto",
    srcs = [
        "atom_type_ext.proto",
    ],
)

proto_library(
    name = "donor_acceptor_proto",
    srcs = [
        "donor_acceptor.proto",
    ],
)

cc_proto_library(
    name = "donor_acceptor_cc_proto",
    deps = [
        ":donor_acceptor_proto",
    ],
)

proto_library(
    name = "mol2graph_proto",
    srcs = [
        "mol2graph.proto",
    ],
)

cc_proto_library(
    name = "mol2graph_cc_proto",
    deps = [
        ":mol2graph_proto",
    ],
)

py_proto_library(
    name = "mol2graph_py_proto",
    srcs = [
        "mol2graph.proto",
    ],
)

proto_library(
    name = "geometric_constraints_proto",
    srcs = [
        "geometric_constraints.proto",
    ],
)

cc_proto_library(
    name = "geometric_constraints_cc_proto",
    deps = [
        ":geometric_constraints_proto",
    ],
)

py_proto_library(
    name = "geometric_constraints_py_proto",
    srcs = [
        "geometric_constraints.proto",
    ],
    visibility = ["//visibility:public"],
)

proto_library(
    name = "molecule_to_query_proto",
    srcs = [
        "molecule_to_query.proto",
    ],
    deps = [
        ":substructure_proto",
    ],
)

cc_proto_library(
    name = "molecule_to_query_cc_proto",
    deps = [
        ":molecule_to_query_proto",
    ],
)

proto_library(
    name = "pharmacophore_proto",
    srcs = [
        "pharmacophore.proto",
    ],
)

cc_proto_library(
    name = "pharmacophore_cc_proto",
    deps = [
        ":pharmacophore_proto",
    ]
)

py_proto_library(
    name = "pharmacophore_py_proto",
    srcs = [
        "pharmacophore.proto",
    ],
    visibility = ["//visibility:public"],
)

proto_library(
    name = "substructure_proto",
    srcs = [
        "substructure.proto",
    ],
    deps = [
        "geometric_constraints_proto",
    ],
)

proto_library(
    name = "standardise_proto",
    srcs = [
        "standardise.proto",
    ],
)

cc_proto_library(
    name = "standardise_cc_proto",
    deps = [
        "standardise_proto",
    ],
)

cc_proto_library(
    name = "substructure_cc_proto",
    deps = [
        ":substructure_proto",
    ],
)
py_proto_library(
    name = "substructure_py_proto",
    srcs = ["substructure.proto"],
    deps = [
        "geometric_constraints_py_proto",
    ],
    visibility = ["//visibility:public"],
)


proto_library(
    name = "reaction_proto",
    srcs = [
        "reaction.proto",
    ],
    deps = [
        ":substructure_proto",
        ":toggle_kekule_form",
    ],
)

cc_proto_library(
    name = "reaction_cc_proto",
    deps = [
        "reaction_proto",
    ],
)

py_proto_library(
    name = "reaction_py_proto",
    srcs = [
        "reaction.proto",
    ],
    deps = [
        "substructure_py_proto",
        "toggle_kekule_form_py_proto",
    ],
    visibility = ["//visibility:public"],
)

proto_library(
    name = "toggle_kekule_form",
    srcs = [
        "toggle_kekule_form.proto",
    ],
    deps = [
        ":substructure_proto",
    ],
)

cc_proto_library(
    name = "toggle_kekule_form_cc_proto",
    deps = [
        ":toggle_kekule_form",
    ],
)

py_proto_library(
    name = "toggle_kekule_form_py_proto",
    srcs = [
        "toggle_kekule_form.proto",
    ],
    deps = [
        ":substructure_py_proto",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "coordinate_box",
    srcs = [
        "coordinate_box.cc",
    ],
    hdrs = [
        "coordinate_box.h",
    ],
    deps = [
        ":space_vector",
        "//Foundational/iwstring",
    ],
)

cc_library(
    name = "iwmolecule",
    srcs = [
        "ISIS_Atom_List.cc",
        "_qry_wstats.cc",
        "_substructure.cc",
        "allowed_elements.cc",
        "aromatic.cc",
        "atom.cc",
        "atom_alias.cc",
        "atom_typing.cc",
        "bond.cc",
        "bond_list.cc",
        "careful_frag.cc",
        "charge_assigner.cc",
        "charge_calculation.cc",
        "chiral_centre.cc",
        "cis_trans_bond.cc",
        "coordinate_box.cc",
        "coordinates.cc",
        "count_atoms_in_smiles.cc",
        "csubstructure.cc",
        "dihedral.cc",
        "donor_acceptor.cc",
        "down_the_bond.cc",
        "element.cc",
        "element_hits_needed.cc",
        "ematch.cc",
        "etrans.cc",
        "frag.cc",
        "geometric_constraints.cc",
        "is_actually_chiral.cc",
        "iwrcb.cc",
        "iwrnm.cc",
        "iwsubstructure.cc",
        "mdl_atom_record.cc",
        "mdl_file_data.cc",
        "mdl_molecule.cc",
        "misc2.cc",
        "moleculeio.cc",
        "mol2graph.cc",
        "mol2graph_proto.cc",
        "molecule.cc",
        "molecule_smarts.cc",
        "molecule_to_query.cc",
        "moleculeb.cc",
        "moleculed.cc",
        "moleculeh.cc",
        "moleculer.cc",
        "numass.cc",
        "ostream_and_type.cc",
        "output.cc",
        "parse_smarts_tmp.cc",
        "path.cc",
        "path_around_ring.cc",
        "path_scoring.cc",
        "pearlman.cc",
        "qry_wstats.cc",
        "ring_bond_iterator.cc",
        "rmele.cc",
        "rotbond_common.cc",
        "rwsubstructure.cc",
        "rwsubstructure_proto.cc",
        "set_of_atoms.cc",
        "smiles.cc",
        "smiles_support.cc",
        "sort_atoms.cc",
        "ss_atom_env.cc",
        "ss_bonds.cc",
        "ss_ring.cc",
        "ss_ring_base.cc",
        "ss_ring_sys.cc",
        "standardise.cc",
        "substructure_a.cc",
        "substructure_chiral.cc",
        "substructure_env.cc",
        "substructure_nmab.cc",
        "substructure_results.cc",
        "substructure_spec.cc",
        "symm_class_can_rank.cc",
        "symmetry.cc",
        "target.cc",
        "temp_detach_atoms.cc",
        "toggle_kekule_form.cc",
        "tokenise_atomic_smarts.cc",
        "u3b.c",
        "unique.cc",
    ],
    hdrs = [
        "allowed_elements.h",
        "aromatic.h",
        "atom.h",
        "atom_alias.h",
        "atom_typing.h",
        "bond.h",
        "bond_list.h",
        "cahn_ingold_prelog.h",
        "charge_assigner.h",
        "charge_calculation.h",
        "chiral_centre.h",
        "cis_trans_bond.h",
        "collection_template.h",
        "coordinate_box.h",
        "coordinates.h",
        "dihedral.h",
        "donor_acceptor.h",
        "element.h",
        "ematch.h",
        "etrans.h",
        "f2c.h",
        "geometric_constraints.h",
        "is_actually_chiral.h",
        "istream_and_type.h",
        "iwmtypes.h",
        "iwrcb.h",
        "iwrnm.h",
        "marvin.h",
        "mdl.h",
        "mdl_atom_record.h",
        "mdl_file_data.h",
        "mdl_molecule.h",
        "misc2.h",
        "moleculeio.h",
        "mol2graph.h",
        "mol2graph_proto.h",
        "molecule.h",
        "molecule_arom.h",
        "molecule_mdl.h",
        "molecule_smi.h",
        "molecule_to_query.h",
        "moleculej.h",
        "numass.h",
        "ostream_and_type.h",
        "output.h",
        "parse_smarts_tmp.h",
        "path.h",
        "path_around_ring.h",
        "path_scoring.h",
        "pearlman.h",
        "qry_wstats.h",
        "rmele.h",
        "rotbond_common.h",
        "rwsubstructure.h",
        "set_of_atoms.h",
        "smiles.h",
        "space_vector.h",
        "standardise.h",
        "substructure.h",
        "symmetry.h",
        "target.h",
        "temp_detach_atoms.h",
        "tmpsssr.h",
        "toggle_kekule_form.h",
        "tokenise_atomic_smarts.h",
        "u3b.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":atom_type_ext_cc_proto",
        ":donor_acceptor_cc_proto",
        ":mol2graph_cc_proto",
        ":molecule_to_query_cc_proto",
        ":pharmacophore_cc_proto",
        ":space_vector",
        ":standardise_cc_proto",
        ":substructure_cc_proto",
        ":toggle_kekule_form_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iw_tdt",
        "//Foundational/iwbits",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwqsort",
        "//Foundational/iwstring",
        "//Foundational/xmlParser:xmlparser",
        "@com_google_protobuf//:protobuf",
    ],
)


cc_library (
    name = "charge_assigner",
    srcs = [
        "charge_assigner.cc",
    ],
    hdrs = [
        "charge_assigner.h",
    ],
    deps = [
        ":iwmolecule",
    ],
)

cc_library (
    name = "donor_acceptor",
    srcs = [
        "donor_acceptor.cc",
    ],
    hdrs = [
        "donor_acceptor.h",
    ],
    deps = [
        ":iwmolecule",
    ],
)

cc_library(
    name = "moleculeio",
    srcs = [
        "charmm.cc",
        "cif.cc",
        "csv.cc",
        "inchi_dummy.cc",
        "marvin.cc",
        "mdl.cc",
        "mdl_v30.cc",
        "mmod.cc",
        "moe.cc",
        "mrk.cc",
        "msi.cc",
        "pdb.cc",
        "rwmolecule.cc",
        "smi.cc",
        "tripos.cc",
        "xyz.cc",
    ],
    hdrs = [
        "readmdl.h",
        "rwmolecule.h",
        "tripos.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":iwmolecule",
    ],
    alwayslink = True,
)

cc_library(
    name = "iwreaction",
    srcs = [
        "iwreaction.cc",
        "iwreaction2.cc",
        "reaction_from_smirks.cc",
        "reaction_iterator.cc",
        "reaction_match_conditions.cc",
        "reaction_proto.cc",
        "rxn_3d_replace.cc",
        "rxn_bond_angle.cc",
        "rxn_bond_length.cc",
        "rxn_dihedral_angle.cc",
        "rxn_matched_sors.cc",
        "rxnfile.cc",
        "rxnfile2.cc",
        "rxnfile3.cc",
    ],
    hdrs = [
        "iwreaction.h",
        "reaction_match_conditions.h",
        "rxn_file.h",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":reaction_cc_proto",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "atom_pair_fingerprint",
    srcs = [
        "atom_pair_fingerprint.cc",
    ],
    hdrs = [
        "atom_pair_fingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "geometric_constraints",
    srcs = [
        "geometric_constraints.cc",
    ],
    hdrs = [
        "geometric_constraints.h",
    ],
    deps = [
        ":geometric_constraints_cc_proto",
        ":iwmolecule",
    ],
)

cc_library(
    name = "linear_fingerprint",
    srcs = [
        "linear_fingerprint.cc",
    ],
    hdrs = [
        "linear_fingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "molecular_formula_lib",
    srcs = [
        "molecular_formula.cc",
    ],
    hdrs = [
        "molecular_formula.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/iwmisc",
        ":iwmolecule",
        ":moleculeio",
    ],
)

cc_library(
    name = "space_vector",
    srcs = [
    ],
    hdrs = [
        "space_vector.h",
    ],
    deps = [
        "//Foundational/iwmisc",
    ],
)

cc_test(
    name = "aromatic_test",
    timeout = "short",
    srcs = ["aromatic_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "atom_test",
    timeout = "short",
    srcs = ["atom_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "cahn_ingold_prelog_test",
    srcs = [
        "cahn_ingold_prelog_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    timeout="short",
)

cc_test(
    name = "dihedral_test",
    srcs = [
        "dihedral_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    timeout="short",
)

cc_test(
    name = "donor_acceptor_test",
    srcs = [
        "donor_acceptor_test.cc",
    ],
    data = [
        "@donor_acceptor//:donor_acceptor_data",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    timeout="short",
)

cc_test(
    name = "element_hits_needed_test",
    srcs = [
        "element_hits_needed_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    timeout="short",
)

cc_test(
    name = "linear_fingerprint_test",
    timeout = "short",
    srcs = ["linear_fingerprint_test.cc"],
    deps = [
        ":iwmolecule",
        ":linear_fingerprint",
        ":moleculeio",
        "//Foundational/iwmisc:sparse_fp_creator",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "frag_test",
    timeout = "short",
    srcs = ["frag_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "geometric_constraints_test",
    timeout = "short",
    srcs = ["geometric_constraints_test.cc"],
    deps = [
        ":geometric_constraints",
        ":geometric_constraints_cc_proto",
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "atom_pair_fingerprint_test",
    timeout = "short",
    srcs = ["atom_pair_fingerprint_test.cc"],
    deps = [
        ":atom_pair_fingerprint",
        ":iwmolecule",
        ":moleculeio",
        "//Foundational/iwmisc:sparse_fp_creator",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "coordinate_box_test",
    timeout = "short",
    srcs = [
        "coordinate_box_test.cc",
    ],
    deps = [
        ":coordinate_box",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "geometry_test",
    timeout = "short",
    srcs = [
        "geometry_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "molecular_formula_test",
    srcs = [
        "molecular_formula_test.cc",
    ],
    deps = [
        "molecular_formula_lib",
        "@googletest//:gtest_main",
    ],
    timeout="short",
)

cc_test(
    name = "molecule_test",
    timeout = "short",
    srcs = [
        "molecule_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "molecule_each_index_test",
    timeout = "short",
    srcs = [
        "molecule_each_index_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "ring_bond_count_test",
    timeout = "short",
    srcs = [
        "ring_bond_count_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "ring_test",
    srcs = [
        "ring_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    timeout="short",
)

cc_test(
    name = "smi_test",
    timeout = "short",
    srcs = [
        "smi_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "smiles_test",
    timeout = "short",
    srcs = [
        "smiles_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "space_vector_test",
    timeout = "short",
    srcs = [
        "space_vector_test.cc",
    ],
    deps = [
        ":space_vector",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_a_test",
    srcs = ["substructure_a_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_test(
    name = "ss_bonds_test",
    srcs = [
        "ss_bonds_test.cc",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
        "@googletest//:gtest_main",
    ],
    timeout="short",
)

cc_test(
    name = "substructure_mam_test",
    timeout = "short",
    srcs = ["substructure_mam_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_sepma_test",
    timeout = "short",
    srcs = ["substructure_sepma_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_test",
    timeout = "short",
    srcs = ["substructure_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_spec_test",
    timeout = "short",
    srcs = ["substructure_spec_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_single_query_test",
    timeout = "short",
    srcs = ["substructure_single_query_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_env_test",
    timeout = "short",
    srcs = ["substructure_env_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_ring_test",
    timeout = "short",
    srcs = ["substructure_ring_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_ring_system_test",
    timeout = "short",
    srcs = ["substructure_ring_system_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_smarts_test",
    timeout = "short",
    srcs = ["substructure_smarts_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "standardise_test",
    timeout = "short",
    srcs = ["standardise_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "reaction_test",
    timeout = "short",
    srcs = ["reaction_test.cc"],
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "rxnfile_test",
    timeout = "short",
    srcs = ["rxnfile_test.cc"],
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "rxnfile2_test",
    timeout = "short",
    srcs = ["rxnfile2_test.cc"],
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "ss_ring_base_test",
    timeout = "short",
    srcs = ["ss_ring_base_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "ss_ring_test",
    timeout = "short",
    srcs = ["ss_ring_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_dtb_test",
    timeout = "short",
    srcs = ["substructure_dtb_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "substructure_nmab_test",
    timeout = "short",
    srcs = ["substructure_nmab_test.cc"],
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "symmetry_test",
    timeout = "short",
    srcs = ["symmetry_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "misc2_test",
    timeout = "short",
    srcs = ["misc2_test.cc"],
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "toggle_kekule_form_test",
    timeout = "short",
    srcs = ["toggle_kekule_form_test.cc"],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "unique_test",
    timeout = "short",
    srcs = ["unique_test.cc"],
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "iwmfingerprint",
    srcs = [
        "iwmfingerprint.cc",
    ],
    hdrs = [
        "iwmfingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "iw_vdw",
    srcs = [
        "van_der_waals.cc",
    ],
    hdrs = [
        "iw_vdw.h",
    ],
    deps = [
        ":iwmolecule",
    ],
)

cc_library(
    name = "ring_data",
    srcs = [
        "ring_data.cc",
    ],
    hdrs = [
        "ring_data.h",
    ],
    deps = [
        ":iwmolecule",
    ],
)
