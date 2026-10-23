load("@rules_proto//proto:defs.bzl", "proto_library")
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")

proto_library(
    name = "retrosynthesis_proto",
    srcs = [
        "retrosynthesis.proto",
    ],
)

cc_proto_library( 
    name = "retrosynthesis_cc_proto",
    deps = [
        ":retrosynthesis_proto",
    ],
)

cc_library(
    name = "retrosynthesis_lib",
    srcs = [
        "retrosynthesis.cc",
        "retrosynthesis.h",
    ],
    deps = [
        ":retrosynthesis_cc_proto",
        "//Foundational/iwstring:absl_hash",
        "//Utilities/GFP_Tools:gfp",
        "//Utilities/GFP_Tools:gfp_standard",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary (
    name = "retrosynthesis",
    srcs = [
        "retrosynthesis_main.cc",
    ],
    deps = [
        ":retrosynthesis_lib",
    ],
)
