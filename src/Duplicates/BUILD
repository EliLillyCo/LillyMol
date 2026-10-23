load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
cc_binary(
    name = "topological_torsion",
    srcs = [
        "topological_torsion_main.cc",
    ],
    deps = [
        ":topological_torsion_lib",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
)

cc_library(
    name = "topological_torsion_lib",
    srcs = [
        "topological_torsion.cc",
    ],
    hdrs = [
        "topological_torsion.h",
    ],
    deps = [
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwstring",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_test(
    name = "topological_torsion_test",
    timeout = "short",
    srcs = [
        "topological_torsion_test.cc",
    ],
    deps = [
        ":topological_torsion_lib",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)
