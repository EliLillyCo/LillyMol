load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//bazel:py_proto_library.bzl", "py_proto_library")

proto_library(
    name = "ilpd_reduced_graph_proto",
    srcs = [
        "ilpd_reduced_graph.proto",
    ],
)

cc_proto_library(
    name = "ilpd_reduced_graph_cc_proto",
    deps = [
        "ilpd_reduced_graph_proto",
    ],
)

py_proto_library(
    name = "ilpd_reduced_graph_py_proto",
    deps = [
        "ilpd_reduced_graph_proto",
    ],
)

cc_binary(
    name = "ilpd_reduced_graph",
    srcs = [
        "ilpd_reduced_graph.cc",
    ],
    deps = [
        ":ilpd_reduced_graph_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)
