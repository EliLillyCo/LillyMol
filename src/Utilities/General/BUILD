load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":average",
        ":bucketise",
        ":byte_offset_index",
        ":class_label_translation",
        ":concat_files",
        ":concat_files_large_memory",
        ":correlate",
        ":descriptor_file_cat",
        ":descriptor_file_filter",
        ":descriptor_file_select_rows",
        ":descriptor_file_sort",
        ":descriptor_file_to_svm_lite",
        ":descriptor_file_transpose",
        ":descriptors_to_fingerprint",
        ":dfilefilter",
        ":difference_sort",
        ":diff_line_by_line",
        ":distribution",
        ":enough_inventory_bdb",
        ":feature_scaling",
        ":fetch_sdf",
        ":fetch_sdf_quick",
        ":fetch_smiles",
        ":fetch_smiles_quick",
        ":grep_by_column",
        ":grid_overlap",
        ":_inventory",
        ":iwcut",
        ":iwsplit",
        ":jfilecompare",
        ":kstat_correlated",
        ":mispredicted",
        ":model_average",
        ":nextdir",
        ":normalise",
        ":notenoughvariance",
        ":numeric_differences",
        ":random_records",
        ":rearrange_columns",
        ":rmsigma",
        ":running_average",
        ":shuffle_file",
        ":spearman_rank",
        ":stratified_samples",
        ":svm_lite_to_gfp",
        ":tcount",
        ":test_train_split_classification",
        ":test_train_split_cluster",
        ":test_train_split_random",
        ":unbalanced_quotes",
        ":unique_rows",
        ":whatsmissing",
    ],
)

proto_library(
    name = "class_label_translation_proto",
    srcs = [
        "class_label_translation.proto",
    ],
)

cc_proto_library(
    name = "class_label_translation_cc_proto",
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "class_label_translation_proto",
    ],
)

cc_binary(
    name = "average",
    srcs = [
        "average.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "bucketise",
    srcs = [
        "bucketise.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "byte_offset_index",
    srcs = [
        "byte_offset_index.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_cat",
    srcs = [
        "descriptor_file_cat.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_filter",
    srcs = [
        "descriptor_file_filter.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_to_svm_lite",
    srcs = [
        "descriptor_file_to_svm_lite.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_transpose",
    srcs = [
        "descriptor_file_transpose.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "difference_sort",
    srcs = [
        "difference_sort.cc",
    ],
    deps = [
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwqsort:iwqsort",
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
        "//Foundational/iwmisc",
    ],
)

proto_library(
    name = "feature_scaling_proto",
    srcs = [
        "feature_scaling.proto",
    ],
)

cc_proto_library(
    name = "feature_scaling_cc_proto",
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "feature_scaling_proto",
    ],
)

py_proto_library(
    name = "feature_scaling_py_proto",
    srcs = [
        "feature_scaling.proto",
    ],
)

cc_library(
    name = "scaler",
    srcs = [
    ],
    hdrs = [
        "scaler.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "feature_scaling_cc_proto",
    ],
)

cc_library(
    name = "ttitem",
    srcs = [
        "ttitem.cc",
    ],
    hdrs = [
        "ttitem.h",
    ],
    deps = [
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "class_label_translation",
    srcs = [
        "class_label_translation.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "concat_files",
    srcs = [
        "concat_files.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
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
    name = "correlate",
    srcs = [
        "correlate.cc",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:iw_tabular_data",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "descriptor_file_select_rows",
    srcs = [
        "descriptor_file_select_rows.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_sort",
    srcs = [
        "descriptor_file_sort.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "descriptors_to_fingerprint",
    srcs = [
        "descriptors_to_fingerprint.cc",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwstring",
    ],
)

# This is a duplicate, two different versions were merged,
# but both names are in scripts.
cc_binary(
    name = "dfilefilter",
    srcs = [
        "descriptor_file_filter.cc",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:iw_tabular_data",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "dicer_fragments_collate",
    srcs = [
        "dicer_fragments_collate.cc",
    ],
    deps = [
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwstring",
        "//Molecule_Tools:dicer_fragments_cc_proto",
        "@com_google_absl//absl/container:flat_hash_map",
    ],
)

cc_binary(
    name = "distribution",
    srcs = [
        "distribution.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "enough_inventory_bdb",
    srcs = [
        "enough_inventory_bdb.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "feature_scaling",
    srcs = [
        "feature_scaling.cc",
    ],
    deps = [
        ":feature_scaling_cc_proto",
        ":scaler",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwstring",
    ],
)

py_binary(
    name = "feature_scaling_create",
    srcs = [
        "feature_scaling_create.py",
    ],
    deps = [
        "feature_scaling_py_proto",
    ],
)

cc_binary(
    name = "fetch_sdf",
    srcs = [
        "fetch_sdf.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "fetch_sdf_quick",
    srcs = [
        "fetch_sdf_quick.cc",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "fetch_smiles",
    srcs = [
        "fetch_smiles.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "fetch_smiles_quick",
    srcs = [
        "fetch_smiles_quick.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "grep_by_column",
    srcs = [
        "grep_by_column.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "grid_overlap",
    srcs = [
        "grid_overlap.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:compile_time",
        "//Foundational/iwstring",
        "//Utilities/GFP_Tools:nearneighbours_cc_proto",
    ],
)

cc_binary(
    name = "_inventory",
    srcs = [
        "inventory.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "iwsplit",
    srcs = [
        "iwsplit.cc",
    ],
    deps = [
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "iwcut",
    srcs = [
        "iwcut.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "jfilecompare",
    srcs = [
        "jfilecompare.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "kstat_correlated",
    srcs = [
        "kstat_correlated.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "mispredicted",
    srcs = [
        "mispredicted.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "model_average",
    srcs = [
        "model_average.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "nextdir",
    srcs = [
        "nextdir.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwqsort:iwqsort",
    ],
)

cc_binary(
    name = "normalise",
    srcs = [
        "iwtokeniser.cc",
        "iwtokeniser.h",
        "normalise.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "notenoughvariance",
    srcs = [
        "notenoughvariance.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "numeric_differences",
    srcs = [
        "numeric_differences.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)


cc_binary(
    name = "random_records",
    srcs = [
        "random_records.cc",
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
    name = "rmsigma",
    srcs = [
        "rmsigma.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwqsort:iwqsort",
    ],
)

cc_binary(
    name = "running_average",
    srcs = [
        "running_average.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "shuffle_file",
    srcs = [
        "shuffle_file.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "spearman_rank",
    srcs = [
        "spearman_rank.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwqsort:iwqsort",
    ],
)

cc_binary(
    name = "stratified_samples",
    srcs = [
        "stratified_samples.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "svm_lite_to_gfp",
    srcs = [
        "svm_lite_to_gfp.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_binary(
    name = "tcount",
    srcs = [
        "tcount.cc",
        "tokens_in_quoted_string.cc",
        "tokens_in_quoted_string.h",
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
    name = "test_train_split_classification",
    srcs = [
        "test_train_split_classification.cc",
    ],
    deps = [
        ":ttitem",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "test_train_split_cluster",
    srcs = [
        "test_train_split_cluster.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwbits",
        "//Foundational/iwmisc",
    ],
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
    name = "test_train_split_random",
    srcs = [
        "test_train_split_random.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/iwbits:iwbits",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwqsort:iwqsort",
    ],
)

cc_binary(
    name = "unbalanced_quotes",
    srcs = [
        "unbalanced_quotes.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "unique_rows",
    srcs = [
        "iwtokeniser.cc",
        "iwtokeniser.h",
        "unique_rows.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "whatsmissing",
    srcs = [
        "whatsmissing.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
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

cc_test(
    name = "scaler_test",
    srcs = [
        "scaler_test.cc",
    ],
    deps = [
        ":feature_scaling_cc_proto",
        ":scaler",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)
