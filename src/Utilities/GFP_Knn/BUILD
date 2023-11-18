load("@rules_proto//proto:defs.bzl", "proto_library")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":nn_predictions",
        ":nn_training",
        ":test_t_test",
    ],
)

cc_binary(
    name = "nn_predictions",
    srcs = [
        "nn_predictions.cc",
    ],
    deps = [
        ":nn_specification",
        ":nn_support",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
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
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/Distance_Matrix:distance_matrix_base",
        "//Utilities/GFP_Tools:bsquared",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "test_t_test",
    srcs = [
        "test_t_test.cc",
    ],
    deps = [
        ":iwpvalue",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
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
        "//Foundational/iwstring:iwstring",
    ],
)

cc_library(
    name = "iwpvalue",
    srcs = [
        "iwpvalue.cc",
        "dq1da.c",
        "idamax.c",
        "machcon.c",
        "uni.c",
        "xerror.c",
    ],
    hdrs = [
        "iwpvalue.h",
    ],
    deps = [
        "@local_f2c//:f2c",
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
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
)

cc_library(
    name = "nn_support",
    srcs = [
        "nn_stuff.cc",
        "nn_specification.cc",
        "prediction_weight_function.cc",
    ],
    hdrs = [
        "nn_results_.h",
        "nn_results.h",
        "nn_specification.h",
        "nn_stuff.h",
        "prediction_weight_function.h",
    ],
    deps = [
        ":extract_from_tdt_form",
        "//Foundational/iwbits:iwbits",
        "//Foundational/iwmisc:iwmisc",
        "//Utilities/Distance_Matrix:distance_matrix_base",
        "//Utilities/GFP_Tools:gfp",
    ],
)
