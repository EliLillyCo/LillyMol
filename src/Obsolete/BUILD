load("//build_deps:install.bzl", "local_install")

# Some executables are still used in certain systems.
local_install(
    name = "install",
    srcs = [
        ":complex_chirality",
        ":concat_files_large_memory",
        ":descriptor_file_same_row_order",
        ":diff_line_by_line",
        ":gfp_spread_threaded",
        ":prediction_bias",
        ":rearrange_columns",
    ],
)

cc_binary(
    name = "concat_files_large_memory",
    srcs = [
        "concat_files_large_memory.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)     

cc_binary(
    name = "complex_chirality",
    srcs = [
        "complex_chirality.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ]
)   

cc_binary(
    name = "confusion_matrix",
    srcs = [
        "confusion_matrix.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwqsort:iwqsort",
        "//Utilities/GFP_Tools:metric",
    ],
)

# Does not build, interesting idea, maybe return to it someday.
# cc_binary(
#     name = "descriptor_evaluation",
#     srcs = [
#         "descriptor_evaluation.cc",
#     ],
#     tags = [
#         "tbb",
#     ],
#     deps = [
#         "//Foundational/accumulator:accumulator",
#         "//Foundational/cmdline:iwcmdline",
#         "//Foundational/data_source:iwstring_data_source",
#         "//Foundational/iwbits:iwbits",
#         "//Foundational/iwmisc:iwmisc",
#         "//Foundational/iw_tdt:iw_tdt",
#         "//Utilities/GFP_Tools:gfp",
#         "@onetbb//:tbb",
#     ],
# )

cc_binary(
    name = "descriptor_file_to_molecular_properties",
    srcs = [
        "descriptor_file_to_molecular_properties.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_profile_activity",
    srcs = [
        "descriptor_file_profile_activity.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_same_row_order",
    srcs = [
        "descriptor_file_same_row_order.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(

    name = "descriptor_file_to_mahalanobis_gfp",
    srcs = [
        "descriptor_file_to_mahalanobis_gfp.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
)

cc_binary(
    name = "diff_line_by_line",
    srcs = [
        "diff_line_by_line.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)  


cc_binary(
    name = "fingerprint_weights_to_atom_colours",
    srcs = [
        "fingerprint_weights_to_atom_colours.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ]
)

cc_binary(
    name = "folded",
    srcs = [
        "folded.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_bit_filter",
    srcs = [
        "gfp_bit_filter.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_difference_fingerprint",
    srcs = [
        "gfp_difference_fingerprint.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_evaluate_clustering",
    srcs = [
        "gfp_evaluate_clustering.cc",
        "cluster_eval.cc",
        "cluster_eval.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_group_spread",
    srcs = [
        "gfp_group_spread.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
        "//Utilities/GFP_Tools:smiles_id_dist",
    ],
)

cc_binary(
    name = "gfp_leader_opt_v2",
    srcs = [
        "gfp_leader_opt_v2.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_leader_standard_threaded",
    srcs = [
        "gfp_leader_standard_threaded.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp_standard",
        "//Utilities/GFP_Tools:smiles_id_dist",
        "@onetbb//:tbb",
    ],
)

cc_binary(
    name = "gfp_leader_threaded",
    srcs = [
        "gfp_leader_threaded.cc",
        "leader_parallel.cc",
        "leader_parallel.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
        "@onetbb//:tbb",
    ],
)

cc_binary(
    name = "gfp_naive_bayesian_multiclass",
    srcs = [
        "gfp_naive_bayesian_multiclass.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_nearneighbours_threaded",
    srcs = [
        "gfp_nearneighbours_threaded.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_proximity_descriptors",
    srcs = [
        "gfp_proximity_descriptors.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

#  cc_binary(
#      name = "gfp_spread_mpi",
#      srcs = [
#          "gfp_spread_mpi.cc",
#      ],
#      deps = [
#          "//Foundational/accumulator:accumulator",
#          "//Foundational/cmdline:iwcmdline",
#          "//Foundational/data_source:iwstring_data_source",
#          "//Foundational/iwmisc:iwmisc",
#          "//Foundational/iw_tdt:iw_tdt",
#          "//Utilities/GFP_Tools:gfp",
#      ],
#      tags = [
#          "mpi",
#      ],
#  )
#  cc_binary(
#      name = "gfp_spread_mpi_v2",
#      srcs = [
#          "gfp_spread_mpi_v2.cc",
#      ],
#      deps = [
#          "//Foundational/accumulator:accumulator",
#          "//Foundational/cmdline:iwcmdline",
#          "//Foundational/data_source:iwstring_data_source",
#          "//Foundational/iwmisc:iwmisc",
#          "//Foundational/iw_tdt:iw_tdt",
#          "//Utilities/GFP_Tools:gfp",
#      ],
#      tags = [
#          "mpi",
#      ],
#  )

cc_binary(
    name = "gfp_spread_threaded",
    srcs = [
        "gfp_spread_threaded.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
        "//Utilities/GFP_Tools:spread_v2",
    ],
)

cc_binary(
    name = "gfp_test_train_split",
    srcs = [
        "gfp_test_train_split.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
)

cc_binary(
    name = "gfp_vector_differences",
    srcs = [
        "gfp_vector_differences.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "nn_identify_active_outliers",
    srcs = [
        "nn_identify_active_outliers.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
)

cc_binary(
    name = "prediction_bias",
    srcs = [
        "prediction_bias.cc",
        "sgefa.c",
        "sgesl.c",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "@local_f2c//:f2c",
    ],
    linkopts = [
        "-lblas",
    ],
)

cc_binary(
    name = "rearrange_columns",
    srcs = [
        "rearrange_columns.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "ring_combinations",
    srcs = [
        "ring_combinations.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ] 
)

cc_binary(
    name = "rocs_rpt2nn",
    srcs = [
        "rocs_rpt2nn.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Utilities/GFP_Tools:gfp",
        "//Utilities/GFP_Tools:smiles_id_dist",
    ],
)

cc_binary(
    name = "selimsteg_build",
    srcs = [
        "selimsteg_build.cc",  
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ]
)

cc_binary(
    name = "test_train_split_from_leader",
    srcs = [
        "test_train_split_from_leader.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "zero_pad",
    srcs = [
        "zero_pad.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)
