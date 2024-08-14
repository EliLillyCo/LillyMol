load("@rules_cc//cc:defs.bzl", "cc_binary")
load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":_fragment_filter",
        ":abraham",
        ":activity_consistency",
        ":align_molecule",
        ":alogp",
        ":atom_triples",
        ":chirality_fingerprint",
        ":common_names",
        ":dbf",
        ":dicer",
        ":dicer_to_topological_types",
        ":echoqry",
        ":echorxn",
        ":ec_fingerprint",
        ":enumeration",
        ":extended_atom_pairs",
        ":ez_descriptor",
        ":ez_fingerprint",
        ":ez_fingerprint_v2",
        ":fileconv",
        ":fingerprint_substructure",
        ":firstatom",
        ":get_coordinates",
        ":get_substituents",
        ":gfp_erg",
        ":ghose_crippen",
        ":grease",
        ":grep_molecule",
        ":grid_fingerprint",
        ":hydrophobic_sections",
        ":id_chirality",
        ":iwdemerit",
        ":iwdescr",
        ":iwecfp",
        ":iwecfp_intermolecular",
        ":iwfp",
        ":iwpathd",
        ":jwcats",
        ":jwdip",
        ":jwdist",
        ":jwestate",
        ":jwmedv",
        ":jwmolconn",
        ":jwmorse",
        ":jwsadb",
        ":linear_fingerprint",
        ":long_molecules",
        ":maccskeys",
        ":make_these_molecules",
        ":medchemwizard",
        ":minor_changes",
        ":mkfrag",
        ":mol2qry",
        ":mol2SAFE",
        ":molecular_abstraction",
        ":molecular_grid",
        ":molecular_merge",
        ":molecular_scaffold",
        ":molecular_transformations",
        ":molecular_variants",
        ":molecule_filter",
        ":molecule_subset",
        ":molecules_from_reagents",
        ":msort",
        ":msort_parallel",
        ":numbered_smiles",
        ":overlapping_fragment_model",
        ":parsimonious_set",
        ":pharmacophore_2d",
        ":preferred_smiles",
        ":psafp",
        ":pubchem_fingerprints",
        ":r1r2etc",
        ":random_geometric_changes",
        ":random_molecular_permutations",
        ":random_smiles",
        ":reduced_graph",
        ":remove_and_label",
        ":remove_matched_atoms",
        ":representations",
        ":retrosynthesis",
        ":rgroup",
        ":ring_extraction",
        ":ring_fingerprint",
        ":ring_replacement",
        ":ring_replacement_collate",
        ":ring_size_fingerprint",
        ":ring_substitution",
        ":ring_trimming",
        ":rotatable_bond_fingerprint",
        ":rotatable_bonds",
        ":rule_of_five",
        ":rxn_fingerprint",
        ":rxn_reverse",
        ":rxn_signature",
        ":rxn_standardize",
        ":rxn_substructure_search",
        ":same_structures",
        ":smiles_mutation",
        ":sp3_filter",
        ":substituent_model",
        ":substitutions",
        ":substructure_match_fraction",
        ":substructure_mcs",
        ":superimpose_by_matched_atoms",
        ":tautomer_generation",
        ":temperature",
        ":tnass",
        ":topotorsion",
        ":topotorsion_fingerprints",
        ":tp1_summarise",
        ":tp_first_pass",
        ":trxn",
        ":tshadow",
        ":tsmiles",
        ":tstandardise",
        ":tsubstructure",
        ":tsubstructure_summarise_hits",
        ":tsymmetry",
        ":unique_molecules",
        ":verloop",
        ":xlogp",
        ":xray_structure_compare",
    ],
)

proto_library(
    name = "atomic_distance_fingerprint_proto",
    srcs = [
        "atomic_distance_fingerprint.proto",
    ],
)

cc_proto_library(
    name = "atomic_distance_fingerprint_cc_proto",
    deps = [
        ":atomic_distance_fingerprint_proto",
    ],
)

proto_library(
    name = "common_names_proto",
    srcs = [
        "common_names.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "common_names_cc_proto",
    deps = [
        ":common_names_proto",
    ],
)

proto_library(
    name = "demerit_proto",
    srcs = [
        "demerit.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "demerit_cc_proto",
    deps = [
        ":demerit_proto",
    ],
)

proto_library(
    name = "dicer_fragments_proto",
    srcs = [
        "dicer_fragments.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "dicer_fragments_cc_proto",
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":dicer_fragments_proto",
    ],
)

py_proto_library(
    name = "dicer_fragments_py_proto",
    srcs = [
        "dicer_fragments.proto",
    ],
)

proto_library(
    name = "enumeration_proto",
    srcs = [
        "enumeration.proto",
    ],
)

cc_proto_library(
    name = "enumeration_cc_proto",
    deps = [
        "enumeration_proto",
    ],
)

proto_library(
    name = "iwdescr_proto",
    srcs = [
        "iwdescr.proto",
    ],
)

cc_proto_library(
    name = "iwdescr_cc_proto",
    deps = [
        ":iwdescr_proto",
    ],
)

py_proto_library(
    name = "iwdescr_py_proto",
    srcs = [
        ":iwdescr.proto",
    ],
)

proto_library(
    name = "make_these_molecules_proto",
    srcs = [
        "make_these_molecules.proto",
    ],
)

cc_proto_library(
    name = "make_these_molecules_cc_proto",
    deps = [
        "make_these_molecules_proto",
    ],
)

proto_library(
    name = "minor_changes_proto",
    srcs = [
        "minor_changes.proto",
    ],
)

cc_proto_library(
    name = "minor_changes_cc_proto",
    deps = [
        "minor_changes_proto",
    ],
)

proto_library(
    name = "molecule_database_options_proto",
    srcs = [
        "molecule_database_options.proto",
    ],
    deps = [
        "//Molecule_Lib:mol2graph_proto",
    ],
)

cc_proto_library(
    name = "molecule_database_options_cc_proto",
    deps = [
        ":molecule_database_options_proto",
    ],
)

proto_library(
    name = "molecule_filter_proto",
    srcs = [
        "molecule_filter.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "molecule_filter_cc_proto",
    deps = [
        ":molecule_filter_proto",
    ],
)

proto_library(
    name = "ligand_protein_interactions_proto",
    srcs = [
        "ligand_protein_interactions.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "ligand_protein_interactions_cc_proto",
    deps = [
        ":ligand_protein_interactions_proto",
    ],
)

proto_library(
    name = "multiple_reactions_proto",
    srcs = [
        "multiple_reactions.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "multiple_reactions_cc_proto",
    deps = [
        ":multiple_reactions_proto",
    ],
)

proto_library(
    name = "pharmacophore_2d_proto",
    srcs = [
        "pharmacophore_2d.proto",
    ],
)

cc_proto_library(
    name = "pharmacophore_2d_cc_proto",
    deps = [
        "pharmacophore_2d_proto",
    ],
)

proto_library(
    name = "random_molecular_permutations_proto",
    srcs = [
        "random_molecular_permutations.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "random_molecular_permutations_cc_proto",
    deps = [
        ":random_molecular_permutations_proto",
    ],
)

proto_library(
    name = "reduced_graph_proto",
    srcs = [
        "reduced_graph.proto",
    ],
    deps = [
        "//Molecule_Lib:substructure_proto",
    ],
)

cc_proto_library(
    name = "reduced_graph_cc_proto",
    deps = [
        "reduced_graph_proto",
    ],
)

proto_library(
    name = "retrosynthesis_proto",
    srcs = [
        "retrosynthesis.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "retrosynthesis_cc_proto",
    deps = [
        "retrosynthesis_proto",
    ],
)

proto_library(
    name = "ring_closure_proto",
    srcs = [
        "ring_closure.proto",
    ],
)

cc_proto_library(
    name = "ring_closure_cc_proto",
    deps = [
        ":ring_closure_proto",
    ],
)

proto_library(
    name = "replacement_ring_proto",
    srcs = [
        "replacement_ring.proto",
    ],
)

cc_proto_library(
    name = "replacement_ring_cc_proto",
    deps = [
        ":replacement_ring_proto",
    ],
)

proto_library(
    name = "scaffolds_proto",
    srcs = [
        "scaffolds.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "scaffolds_cc_proto",
    deps = [
        ":scaffolds_proto",
    ],
)

# proto_library(
#     name = "tsubstructure_proto",
#     srcs = [
#         "tsubstructure.proto",
#     ],
#     deps = [
#         "//Molecule_Lib:substructure_proto",
#     ],
# )
#
# cc_proto_library(
#     name = "tsubstructure_cc_proto",
#     deps = [
#         "tsubstructure_proto",
#     ],
# )

proto_library(
    name = "xlogp_proto",
    srcs = [
        "xlogp.proto",
    ],
)

cc_proto_library(
    name = "xlogp_cc_proto",
    deps = [
        "xlogp_proto",
    ],
)

py_proto_library(
    name = "xlogp_py_proto",
    srcs = [
        "xlogp.proto",
    ],
)
cc_binary(
    name = "abraham",
    srcs = [
        "abraham.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "activity_consistency",
    srcs = [
        "activity_consistency.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "aggregate_replacement_rings",
    srcs = [
        "aggregate_replacement_rings.cc",
    ],
    deps = [
        ":replacement_ring_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "align_molecule",
    srcs = [
        "align_molecule.cc",
    ],
    deps = [
        ":spatially_common_matched_atoms",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "alogp",
    srcs = [
        "alogp_main.cc",
    ],
    deps = [
        ":alogp_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "atom_pair_fingerprint",
    srcs = [
        "atom_pair_fingerprint_main.cc",
    ],
    deps = [
        ":fingerprint_writer",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:atom_pair_fingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "atom_triples",
    srcs = [
        "atom_triples.cc",
    ],
    deps = [
        ":fingerprint_writer",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "atomic_distance_fingerprint",
    srcs = [
        "atomic_distance_fingerprint.cc",
    ],
    deps = [
        "atomic_distance_fingerprint_cc_proto",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_absl//absl/algorithm:container",
    ],
)

cc_binary(
    name = "bump_check",
    srcs = [
        "bump_check.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "charge_fingerprint",
    srcs = [
        "charge_fingerprint.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Tools:extended_connectivity_fp",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "chirality_fingerprint",
    srcs = [
        "chirality_fingerprint.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Tools:extended_connectivity_fp",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "common_names",
    srcs = [
        "common_names.cc",
    ],
    deps = [
        "common_names_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dbf",
    srcs = [
        "dbf.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer",
    srcs = [
        "dicer.cc",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        ":dicer_lib",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer_fragment_activity",
    srcs = [
        "dicer_fragment_activity.cc",
    ],
    deps = [
        "dicer_fragments_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer_fragment_overlap",
    srcs = [
        "dicer_fragment_overlap.cc",
    ],
    deps = [
        "dicer_fragments_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer_fragment_replace",
    srcs = [
        "dicer_fragment_replace.cc",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        ":molecule_filter_cc_proto",
        ":molecule_filter_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer_to_topological_types",
    srcs = [
        "dicer_to_topological_types.cc",
    ],
    deps = [
        "dicer_fragments_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "ec_fingerprint",
    srcs = [
        "ec_fingerprint_main.cc",
    ],
    deps = [
        ":ec_fingerprint_lib",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "echoqry",
    srcs = [
        "echoqry.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "echorxn",
    srcs = [
        "echorxn.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "element_density",
    srcs = [
        "element_density.cc",
    ],
    deps = [
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)


cc_binary(
    name = "enumeration",
    srcs = [
        "enumeration.cc",
    ],
    deps = [
        "dicer_fragments_cc_proto",
        "enumeration_cc_proto",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "excise_3d_overlap",
    srcs = [
        "excise_3d_overlap.cc",
    ],
    deps = [
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "extended_atom_pairs",
    srcs = [
        "extended_atom_pairs.cc",
    ],
    deps = [
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ez_descriptor",
    srcs = [
        "ez_descriptor.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ez_fingerprint",
    srcs = [
        "ez_fingerprint.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ez_fingerprint_v2",
    srcs = [
        "ez_fingerprint_v2.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "fileconv",
    srcs = [
        "fileconv.cc",
        "fix_structures.cc",
        "fix_structures.h",
        "known_fragment_data.cc",
        "known_fragment_data.h",
    ],
    deps = [
        ":fileconv_lib",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "fingerprint_substructure",
    srcs = [
        "circular_fingerprint.cc",
        "circular_fingerprint.h",
        "fingerprint_substructure.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "firstatom",
    srcs = [
        "firstatom.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_absl//absl/algorithm:container",
    ],
)

cc_binary(
    name = "geometric_descriptors",
    srcs = [
        "geometric_descriptors.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "get_coordinates",
    srcs = [
        "get_coordinates.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "get_linkers",
    srcs = [
        "get_linkers.cc",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "get_substituents",
    srcs = [
        "get_substituents.cc",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        "//Foundational/data_source:iwtfdata_record",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "gfp_erg",
    srcs = [
        "gfp_erg.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ghose_crippen",
    srcs = [
        "ghose_crippen.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/histogram:iwhistogram",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "grid_fingerprint",
    srcs = [
        "grid_fingerprint.cc",
    ],
    deps = [
        ":insight_grid",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "grid_proximity",
    srcs = [
        "grid_proximity.cc",
    ],
    deps = [
        ":insight_grid",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "hydrophobic_sections",
    srcs = [
        "hydrophobic_sections.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "id_chirality",
    srcs = [
        "id_chirality.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "introduction",
    srcs = [
        "introduction.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "iwdemerit",
    srcs = [
        "iwdemerit.cc",
        "qry_and_demerit.h",
        "substructure_demerits.cc",
        "substructure_demerits.h",
    ],
    deps = [
        ":demerit_cc_proto",
        ":iwdemerit_lib",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "iwdescr",
    srcs = [
        "iwdescr.cc",
        "qry_wcharge.cc",
        "qry_wcharge.h",
    ],
    deps = [
        ":alogp_lib",
        ":iwdescr_cc_proto",
        ":nvrtspsa",
        ":partial_symmetry",
        ":xlogp_lib",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "iwecfp",
    srcs = [
        "iwecfp.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "iwecfp_intermolecular",
    srcs = [
        "iwecfp_intermolecular.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "iwfp",
    srcs = [
        "iwfp.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "iwpathd",
    srcs = [
        "iwpathd.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "jwcats",
    srcs = [
        "jwcats.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "jwsadb",
    srcs = [
        "jw_sa_db_descriptor.cc",
    ],
    tags = ["fortran"],
    deps = [
        ":comma",
        ":read_control_file",
        ":surface_area_molvol",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ligand_protein_interactions",
    srcs = [
        "ligand_protein_interactions.cc",
    ],
    deps = [
        "ligand_protein_interactions_cc_proto",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "linear_fingerprint",
    srcs = [
        "linear_fingerprint_main.cc",
    ],
    deps = [
        ":fingerprint_writer",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:linear_fingerprint",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "maccskeys",
    srcs = [
        "maccskeys_fn5_main.cc",
    ],
    deps = [
        ":maccskeys_fn5_lib",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "make_these_molecules",
    srcs = [
        "make_these_molecules.cc",
    ],
    deps = [
        ":make_these_molecules_cc_proto",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "minor_changes",
    srcs = [
        "minor_changes_main.cc",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        ":minor_changes_cc_proto",
        ":minor_changes_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "mol2qry",
    srcs = [
        "mol2qry.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "mol2SAFE",
    srcs = [
        "mol2SAFE.cc",
    ],
    deps = [
        ":highest_ring_number",
        ":fragment_molecule_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwstring:absl_hash",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_absl//absl/container:node_hash_map",
    ],
)

cc_binary(
    name = "molecular_abstraction",
    srcs = [
        "mav2.cc",
        "molecular_abstraction_functions.cc",
        "molecular_abstraction_functions.h",
        "molecular_abstraction_specifications.cc",
        "molecular_abstraction_specifications.h",
    ],
    deps = [
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "molecular_scaffold",
    srcs = [
        "molecular_scaffold.cc",
    ],
    deps = [
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "molecular_transformations",
    srcs = [
        "molecular_transformations.cc",
    ],
    deps = [
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "molecular_variants",
    srcs = [
        "molecular_variants.cc",
    ],
    deps = [
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "molecule_filter",
    srcs = [
        "molecule_filter.cc",
    ],
    deps = [
        ":molecule_filter_cc_proto",
        ":molecule_filter_lib",
        ":alogp_lib",
        ":xlogp_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:nvrtspsa",
    ]
)

cc_binary(
    name = "molecule_pca",
    srcs = [
        "molecule_pca.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/jama:jama",
        "//Foundational/tnt:tnt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ]
)

cc_binary(
    name = "molecule_subset",
    srcs = [
        "molecule_subset.cc",
    ],
    deps = [
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "molecules_from_reagents",
    srcs = [
        "molecules_from_reagents.cc",
    ],
    deps = [
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "numbered_smiles",
    srcs = [
        "numbered_smiles.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "msort",
    srcs = [
        "msort.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "msort_parallel",
    srcs = [
        "msort_parallel.cc",
    ],
    linkopts = [
        "-ldl",
    ],
    tags = [
        "tbb",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@onetbb//:tbb",
    ],
)

cc_binary(
    name = "multiple_reactions",
    srcs = [
        "multiple_reactions.cc",
    ],
    deps = [
        "multiple_reactions_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "parsimonious_set",
    srcs = [
        "parsimonious_set.cc",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "pharmacophore_2d",
    srcs = [
        "pharmacophore_2d.cc",
    ],
    deps = [
        ":pharmacophore_2d_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "preferred_smiles",
    srcs = [
        "preferred_smiles.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "psafp",
    srcs = [
        "psafp.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:nvrtspsa",
    ],
)

cc_binary(
    name = "pubchem_fingerprints",
    srcs = [
        "pubchem_fingerprints.cc",
    ],
    deps = [
        ":nvrtspsa",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "overlapping_fragment_model",
    srcs = [
        "overlapping_fragment_model.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ]
)

cc_binary(
    name = "random_geometric_changes",
    srcs = [
        "random_geometric_changes.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "random_molecular_permutations",
    srcs = [
        "random_molecular_permutations.cc",
        "set_of_target_molecules.cc",
        "set_of_target_molecules.h",
        "swap_kekule_forms.cc",
    ],
    deps = [
        ":linear_path_fingerprint",
        ":maccskeys_fn5_lib",
        ":mpr",
        ":random_molecular_permutations_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Utilities/GFP_Tools:gfp_standard",
        "//Utilities/GFP_Tools:sparse_collection",
    ],
)

cc_binary(
    name = "random_smiles",
    srcs = [
        "random_smiles.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "reduced_graph",
    srcs = [
        "reduced_graph.cc",
    ],
    deps = [
        "reduced_graph_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "remove_and_label",
    srcs = [
        "remove_and_label.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "remove_matched_atoms",
    srcs = [
        "remove_matched_atoms.cc",
    ],
    deps = [
        ":spatially_common_matched_atoms",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "representations",
    srcs = [
        "representations.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "retrosynthesis",
    srcs = [
        "retrosynthesis.cc",
    ],
    deps = [
        ":retrosynthesis_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "rgroup",
    srcs = [
        "rgroup.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ring_closure",
    srcs = [
        "ring_closure_main.cc",
    ],
    deps = [
        ":ring_closure_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ring_size_fingerprint",
    srcs = [
        "ring_size_fingerprint.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ring_substitution",
    srcs = [
        "ring_substitution.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "ring_extraction",
    srcs = [
        "ring_extraction.cc"
    ],
    deps = [
        ":replacement_ring_cc_proto",
        "//Molecule_Tools:reaction_duplicate",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "ring_fingerprint",
    srcs = [
        "ring_fingerprint.cc",
    ],
    deps = [
        ":reaction_duplicate",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "ring_replacement",
    srcs = [
        "ring_replacement.cc",
    ],
    deps = [
        ":reaction_duplicate",
        ":replacement_ring_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "ring_replacement_collate",
    srcs = [
        "ring_replacement_collate.cc",
    ],
    deps = [
        ":replacement_ring_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "ring_trimming",
    srcs = [
        "ring_trimming.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "rotatable_bonds",
    srcs = [
        "rotatable_bonds.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "rotatable_bond_fingerprint",
    srcs = [
        "rotatable_bond_fingerprint.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "rule_of_five",
    srcs = [
        "rule_of_five.cc",
    ],
    deps = [
        ":nvrtspsa",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "rxn_fingerprint",
    srcs = [
        "ct_cache.h",
        "linear_path_fingerprint.cc",
        "linear_path_fingerprint.h",
        "rxn_fingerprint.cc",
    ],
    deps = [
        ":extended_connectivity_fp",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "rxn_reverse",
    srcs = [
        "rxn_reverse.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "rxn_signature",
    srcs = [
        "rxn_signature.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "rxn_substructure_search",
    srcs = [
        "rxn_substructure_search.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "rxn_standardize",
    srcs = [
        "rxn_standardize.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_library(
    name = "alogp_lib",
    srcs = [
        "alogp.cc",
    ],
    hdrs = [
        "alogp.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    visibility = [
        "//visibility:public",
    ],
)

#cc_library(
#    name = "atomic_distance_fingerprint_lib",
#    srcs = [
#        "atomic_distance_fingerprint.cc",
#    ],
#    hdrs = [
#        "atomic_distance_fingerprint.h",
#    ],
#    deps = [
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:iwreaction",
#        "//Molecule_Lib:moleculeio",
#    ],
#)

cc_library(
    name = "comma",
    srcs = [
        "comma_related_procedures.cc",
        "jw_vector_array_operation.cc",
    ],
    hdrs = [
        "comma_related_procedures.h",
        "jw_vector_array_operation.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "dicer_api",
    hdrs = [
        "dicer_api.h",
    ],
    srcs = [
        "dicer_api.cc",
    ],
    deps = [
        ":dicer_lib",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_library(
    name = "dicer_lib",
    hdrs = [
        "dicer_lib.h",
    ],
    srcs = [
        "dicer_lib.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "ec_fingerprint_lib",
    srcs = [
        "ec_fingerprint.cc",
    ],
    hdrs = [
        "ec_fingerprint.h",
    ],
    deps = [
        ":fingerprint_writer",
        "//Molecule_Lib:iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_library(
    name = "molecule_filter_lib",
    srcs = [
        "molecule_filter_lib.cc",
    ],
    hdrs = [
        "molecule_filter_lib.h",
    ],
    deps = [
        ":molecule_filter_cc_proto",
        ":alogp_lib",
        ":xlogp_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:nvrtspsa",
    ],
)

#cc_library(
#    name = "rxn_to_openrxn_lib",
#    srcs = [
#        "rxn_to_openrxn.cc",
#        "rxn_to_openrxn.h",
#    ],
#    deps = [
#        "//Foundational/cmdline:iwcmdline",
#        "//Foundational/iwaray:iwaray",
#        "//Foundational/iwmisc:iwmisc",
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:iwreaction",
#        "//Molecule_Lib:moleculeio",
#        "//Molecule_Lib:substructure_cc_proto",
#        "@ord_schema//:dataset_cc_proto",
#        "@ord_schema//:reaction_cc_proto",
#    ],
#)

#cc_binary(
#    name = "rxn_to_openrxn",
#    srcs = [
#        "rxn_to_openrxn_main.cc"
#    ],
#    deps = [
#        ":rxn_to_openrxn_lib",
#        "@ord_schema//:dataset_cc_proto",
#        "@ord_schema//:reaction_cc_proto",
#    ],
#)

cc_binary(
    name = "same_structures",
    srcs = [
        "same_structures.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "scaffolds",
    srcs = [
        "scaffolds_main.cc",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        ":scaffolds_lib",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_absl//absl/container:flat_hash_map",
    ],
)

cc_binary(
    name = "skeleton",
    srcs = [
        "skeleton.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "smiles_mutation",
    srcs = [
        "random_reactions.cc",
        "random_reactions.h",
        "smiles_mutation.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "spatial_replacement",
    srcs = [
        "spatial_replacement.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "sp3_filter",
    srcs = [
        "sp3_filter.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)
cc_binary(
    name = "substituent_model",
    srcs = [
        "substituent_model.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "substructure_mcs",
    srcs = [
        "substructure_mcs.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_library(
    name = "spatially_common_matched_atoms",
    srcs = [
        "spatially_common_matched_atoms.cc",
    ],
    hdrs = [
        "spatially_common_matched_atoms.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

#cc_library(
#    name = "surface_area_molvol",
#    srcs = [
#        "surface_area_molvol.cc",
#    ],
#    hdrs = [
#        "surface_area_molvol.h"
#    ],
#    deps = [
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:moleculeio",
#        "//Molecule_Lib:iw_vdw",
#    ],
#)

cc_binary(
    name = "tautomer_generation",
    srcs = [
        "tautomer_generation.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "temperature",
    srcs = [
        "temperature.cc",
    ],
    deps = [
        ":mpr",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:linear_fingerprint",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "tnass",
    srcs = [
        "nass.cc",
        "nass.h",
        "tnass.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:linear_fingerprint",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "topotorsion",
    srcs = [
        "topotorsion.cc",
    ],
    deps = [
        ":topotorsion_support",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwstring:absl_hash",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_absl//absl/container:node_hash_map",
    ],
)

cc_binary(
    name = "topotorsion_fingerprints",
    srcs = [
        "topotorsion_fingerprints.cc",
    ],
    deps = [
        ":topotorsion_support",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "tp_first_pass",
    srcs = [
        "tp_first_pass.cc",
    ],
    deps = [
        ":demerit_cc_proto",
        ":tp_first_pass_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "trxn",
    srcs = [
        "trxn.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "tshadow",
    srcs = [
        "shadow.cc",
        "shadow.h",
        "tshadow.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/jama",
        "//Foundational/tnt",
        "//Molecule_Lib:iw_vdw",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_gitlab_libeigen_eigen//:eigen",
    ],
)

cc_binary(
    name = "tsmiles",
    srcs = [
        "tsmiles.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "tstandardise",
    srcs = [
        "tstandardise.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "tsubstructure",
    srcs = [
        "tsubstructure.cc",
        "tsubstructure_fp.cc",
        "tsubstructure_fp.h",
    ],
    deps = [
        ":demerit_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "tsubstructure_summarise_hits",
    srcs = [
        "tsubstructure_summarise_hits.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_binary(
    name = "tsymmetry",
    srcs = [
        "tsymmetry.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "unique_molecules",
    srcs = [
        "unique_molecules.cc",
    ],
    deps = [
        "//Foundational/iw_tdt",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "verloop",
    srcs = [
        "sterimol.cc",
        "sterimol.h",
        "topological_sterimol.cc",
        "topological_sterimol.h",
        "verloop.cc",
    ],
    deps = [
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iw_vdw",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "grease",
    srcs = [
        "grease.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "grep_molecule",
    srcs = [
        "grep_molecule.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:molecular_formula_lib",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "jwmolconn",
    srcs = [
        "e_state_computation_procedure.cc",
        "e_state_computation_procedure.h",
        "jw_molconn.cc",
        "jw_path_with_values.cc",
        "jw_path_with_values.h",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "mkfrag",
    srcs = [
        "mkfrag.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "_fragment_filter",
    srcs = [
        "fragment_filter.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "jwmorse",
    srcs = [
        "jwmorse.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "long_molecules",
    srcs = [
        "long_molecules.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "substitutions",
    srcs = [
        "substitutions.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "superimpose_by_matched_atoms",
    srcs = [
        "superimpose_by_matched_atoms.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "r1r2etc",
    srcs = [
        "r1r2etc.cc",
    ],
    deps = [
        ":nvrtspsa",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "jwdip",
    srcs = [
        "jw_DiP.cc",
    ],
    deps = [
        ":surface_area_molvol",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "jwdist",
    srcs = [
        "jw_common_alignment.cc",
        "jw_common_alignment.h",
        "jw_distribution_along_longest_axis.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "jwestate",
    srcs = [
        "jwestate.cc",
    ],
    deps = [
        ":estate",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "jwmedv",
    srcs = [
        "jw_MEDV13.cc",
    ],
    deps = [
        ":estate",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "molecular_grid",
    srcs = [
        "molecular_grid.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "medchemwizard",
    srcs = [
        "medchemwizard.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "molecular_merge",
    srcs = [
        "molecular_merge.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "smi2linker",
    srcs = [
        "smi2linker.cc",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "substructure_match_fraction",
    srcs = [
        "substructure_match_fraction.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "tp1_summarise",
    srcs = [
        "tp1_summarise.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "xlogp",
    srcs = [
        "xlogp_main.cc",
    ],
    deps = [
        ":xlogp_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "xray_structure_compare",
    srcs = [
        "xray_structure_compare.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_library(
    name = "demerit",
    srcs = [
        "demerit.cc",
    ],
    hdrs = [
        "demerit.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "estate",
    srcs = [
        "e_state_computation_procedure_ori.cc",
    ],
    hdrs = [
        "e_state_computation_procedure_ori.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "extended_connectivity_fp",
    srcs = [
        "extended_connectivity_fp.cc",
    ],
    hdrs = [
        "extended_connectivity_fp.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "extended_pi",
    srcs = [
        "extended_pi.cc",
    ],
    hdrs = [
        "extended_pi.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "fileconv_lib",
    srcs = [
        "fileconv_opts.cc",
        "fix_structures.cc",
        "known_fragment_data.cc",
        "remove_duplicate_fragments.cc",
    ],
    hdrs = [
        "fileconv_opts.h",
        "fix_structures.h",
        "known_fragment_data.h",
        "remove_duplicate_fragments.h",
    ],
    deps = [
        "//Foundational/iwmisc:compile_time",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "fingerprint_writer",
    srcs = [
        "fingerprint_writer.cc",
    ],
    hdrs = [
        "fingerprint_writer.h",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwbits",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwstring",
    ],
)

cc_library(
    name = "fragment_molecule_lib",
    srcs = [
        "fragment_molecule.cc",
    ],
    hdrs = [
        "fragment_molecule.h",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "highest_ring_number",
    srcs = [
        "highest_ring_number.cc",
    ],
    hdrs = [
        "highest_ring_number.h",
    ],
    deps = [
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_library(
    name = "insight_grid",
    srcs = [
        "insight_grid.cc",
    ],
    hdrs = [
        "insight_grid.h",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "iwdemerit_lib",
    srcs = [
        "iwdemerit_lib.cc",
    ],
    hdrs = [
        "iwdemerit_lib.h",
    ],
    deps = [
        ":demerit",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "linear_path_fingerprint",
    srcs = [
        "linear_path_fingerprint.cc",
    ],
    hdrs = [
        "ct_cache.h",
        "linear_path_fingerprint.h",
    ],
    deps = [
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "maccskeys_fn5_lib",
    srcs = [
        "maccskeys_fn5.cc",
    ],
    hdrs = [
        "maccskeys_fn5.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "minor_changes_lib",
    srcs = [
        "minor_changes.cc",
    ],
    hdrs = [
        "minor_changes.h",
    ],
    deps = [
        ":dicer_fragments_cc_proto",
        ":minor_changes_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "mpr",
    srcs = [
        "mpr.cc",
    ],
    hdrs = [
        "mpr.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_library(
    name = "nvrtspsa",
    srcs = [
        "nvrtspsa.cc",
    ],
    hdrs = [
        "nvrtspsa.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "partial_symmetry",
    srcs = [
        "partial_symmetry.cc",
    ],
    hdrs = [
        "partial_symmetry.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "reaction_duplicate",
    srcs = [
        "reaction_duplicate.cc",
    ],
    hdrs = [
        "reaction_duplicate.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_library(
    name = "read_control_file",
    srcs = [
        "read_control_file.cc",
    ],
    hdrs = [
        "read_control_file.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "ring_closure_lib",
    srcs = [
        "ring_closure.cc",
    ],
    hdrs = [
        "ring_closure.h",
    ],
    deps = [
        ":ring_closure_cc_proto",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "ring_ext_rep",
    srcs = [
        "ring_ext_rep.cc",
    ],
    hdrs = [
        "ring_ext_rep.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_library(
    name = "scaffolds_lib",
    srcs = [
        "scaffolds.cc",
    ],
    hdrs = [
        "scaffolds.h",
    ],
    deps = [
        ":scaffolds_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwstring:iwstring",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_library(
    name = "surface_area_molvol",
    srcs = [
        "jwrandom_preset_array.h",
        "smolvol.c",
        "surface_area_molvol.cc",
    ],
    hdrs = [
        "surface_area_molvol.h",
    ],
    tags = ["fortran"],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Molecule_Lib:iw_vdw",
        "//Molecule_Lib:iwmolecule",
        "@local_f2c//:f2c",
    ],
)

cc_library(
    name = "topotorsion_support",
    srcs = [
        "topotorsion_support.cc",
        "torsion_hash.cc",
    ],
    hdrs = [
        "topotorsion_support.h",
        "torsion_hash.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "tp_first_pass_lib",
    srcs = [
        "tp_first_pass_lib.cc",
    ],
    hdrs = [
        "tp_first_pass_lib.h",
    ],
    deps = [
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
        "//Molecule_Lib:iwmolecule",
    ],
)

#cc_library(
#    name = "surface_area_molvol",
#    srcs = [
#        "surface_area_molvol.cc",
#    ],
#    hdrs = [
#        "surface_area_molvol.h",
#    ],
#    deps = [
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:iw_vdw",
#    ],
#)

cc_library(
    name = "unique_molecules_api",
    srcs = [
        "unique_molecules_api.cc",
    ],
    hdrs = [
        "unique_molecules_api.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "@com_google_absl//absl/container:btree",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_library(
    name = "xlogp_lib",
    srcs = [
        "xlogp.cc",
    ],
    hdrs = [
        "xlogp.h"
    ],
    deps = [
        ":xlogp_cc_proto",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_test(
    name = "alogp_test",
    timeout = "short",
    srcs = ["alogp_test.cc"],
    data = [
        "@charge_assigner//:charge_assigner_data",
    ],
    deps = [
        ":alogp_lib",
        "//Molecule_Lib:charge_assigner",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "dicer_api_test",
    timeout = "short",
    srcs = [
        "dicer_api_test.cc",
    ],
    deps = [
        ":dicer_api",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest_main",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_test(
    name = "ec_fingerprint_test",
    timeout = "short",
    srcs = ["ec_fingerprint_test.cc"],
    deps = [
        ":ec_fingerprint_lib",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/iwmisc:sparse_fp_creator",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "extended_pi_test",
    srcs = [
        "extended_pi_test.cc",
    ],
    deps = [
        ":extended_pi",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "fileconv_test",
    srcs = [
        "fileconv_test.cc",
    ],
    deps = [
        ":fileconv_lib",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "fragment_molecule_test",
    srcs = [
        "fragment_molecule_test.cc",
    ],
    deps = [
        ":fragment_molecule_lib",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "highest_ring_number_test",
    srcs = [
        "highest_ring_number_test.cc",
    ],
    deps = [
        ":highest_ring_number",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "iwdemerit_lib_test",
    srcs = [
        "iwdemerit_lib_test.cc",
    ],
    deps = [
        ":demerit",
        ":iwdemerit_lib",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "minor_changes_test",
    srcs = [
        "minor_changes_test.cc",
    ],
    deps = [
        ":minor_changes_cc_proto",
        ":minor_changes_lib",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "mpr_test",
    srcs = [
        "mpr_test.cc",
    ],
    deps = [
        ":mpr",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "partial_symmetry_test",
    timeout = "short",
    srcs = [
        "partial_symmetry_test.cc",
    ],
    deps = [
        ":partial_symmetry",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "ring_ext_rep_test",
    srcs = [
        "ring_ext_rep_test.cc",
    ],
    deps = [
        ":ring_ext_rep",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_test(
    name = "scaffolds_test",
    srcs = [
        "scaffolds_test.cc",
    ],
    deps = [
        ":scaffolds_lib",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "@com_google_protobuf//:protobuf",
    ],
    timeout="short",
)

cc_test(
    name = "tp_first_pass_test",
    timeout = "short",
    srcs = [
        "tp_first_pass_test.cc",
    ],
    deps = [
        ":tp_first_pass_lib",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "xlogp_test",
    srcs = [
        "xlogp_test.cc",
    ],
    deps = [
        ":xlogp_lib",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ]
)
