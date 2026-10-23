load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")

# Note that we have dependencies on Distance_Matrix which is now
# Obsolete. This directory needs to be re-implemented.
cc_binary(
    name = "nn_predictions",
    srcs = [
        "nn_predictions.cc",
    ],
    deps = [
        ":nn_specification",
        ":nn_support",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "nn_training",
    srcs = [
        "nn_training.cc",
    ],
    deps = [
        ":extract_from_tdt_form",
        ":iwpvalue",
        ":nn_support",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Obsolete/Distance_Matrix:distance_matrix_base",
        "//Utilities/GFP_Tools:bsquared",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_library(
    name = "extract_from_tdt_form",
    srcs = [
        "extract_from_tdt_form.cc",
    ],
    hdrs = [
        "extract_from_tdt_form.h",
    ],
    deps = [
        "//Foundational/iwstring",
    ],
)

cc_library(
    name = "iwpvalue",
    srcs = [
        "iwpvalue.cc",
    ],
    hdrs = [
        "iwpvalue.h",
    ],
    visibility = ["//visibility:public"],
)

cc_test(
    name = "iwpvalue_test",
    srcs = [
        "iwpvalue_test.cc",
    ],
    deps = [
        ":iwpvalue",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "nn_specification",
    srcs = [
        "nn_specification.cc",
    ],
    hdrs = [
        "nn_specification.h",
    ],
    deps = [
        ":nn_support",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
    ],
)

cc_library(
    name = "nn_support",
    srcs = [
        "nn_specification.cc",
        "nn_stuff.cc",
        "prediction_weight_function.cc",
    ],
    hdrs = [
        "nn_results.h",
        "nn_results_.h",
        "nn_specification.h",
        "nn_stuff.h",
        "prediction_weight_function.h",
    ],
    deps = [
        ":extract_from_tdt_form",
        "//Foundational/iwbits",
        "//Foundational/iwmisc",
        "//Obsolete/Distance_Matrix:distance_matrix_base",
        "//Utilities/GFP_Tools:gfp",
    ],
)
