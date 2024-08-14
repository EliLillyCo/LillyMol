load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":gene_expression_to_proto",
        ":gene_expression_nearneighbours",
    ],
)

proto_library(
    name = "gene_expression_proto",
    srcs = [
        "gene_expression.proto",
    ],
)

cc_proto_library(
    name = "gene_expression_cc_proto",
    deps = [
        ":gene_expression_proto",
    ]
)

py_proto_library(
    name = "gene_expression_py_proto",
    srcs = [
        "gene_expression.proto",
    ]
)


cc_library(
    name = "gene_expression_lib",
    srcs = [
        "gene_expression.cc",
    ],
    hdrs = [
        "gene_expression.h",
    ],
    deps = [
        ":gene_expression_cc_proto",
    ]
)

cc_library(
    name = "needle_lib",
    srcs = [
        "needle.cc",
    ],
    hdrs = [
        "needle.h",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/iwstring:iwstring",
        ":gene_expression_cc_proto",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_protobuf//:protobuf",
    ]
)

#cc_binary(
#    name = "gene_expression_benchmark",
#    srcs = [
#        "gene_expression_benchmark.cc",
#    ],
#    deps = [
#        ":gene_expression_cc_proto",
#        ":gene_expression_lib",
#        "@benchmark",
#    ]
#)

cc_binary(
    name = "gene_expression_nearneighbours",
    srcs = [
        "gene_expression_nearneighbours.cc",
    ],
    deps = [
        ":gene_expression_cc_proto",
        ":gene_expression_lib",
        ":needle_lib",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwtfdata_record",
        "@com_google_protobuf//:protobuf",
    ],
    tags = [
        "gene_expression",
    ],
)

cc_binary(
    name = "gene_expression_to_proto",
    srcs = [
        "gene_expression_to_proto.cc",
    ],
    deps = [
        ":gene_expression_cc_proto",
        ":gene_expression_lib",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "@com_google_protobuf//:protobuf",
    ]
)

cc_test (
    name = "needle_test",
    srcs = [
        "needle_test.cc",
    ],
    deps = [
        "needle_lib",
        ":gene_expression_cc_proto",
        "@googletest//:gtest_main",
    ],
)
