load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":topological_torsion",
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

cc_binary(
    name = "ring_extraction",
    srcs = [
        "ring_extraction.cc"
    ],
    deps = [
        ":replacement_ring_cc_proto",
        "//Molecule_Tools:reaction_duplicate",
        "//Molecule_Tools:ring_ext_rep",
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
    name = "ring_replacement",
    srcs = [
        "ring_replacement.cc"
    ],
    deps = [
        ":replacement_ring_cc_proto",
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
    name = "topological_torsion",
    srcs = [
        "topological_torsion_main.cc",
    ],
    deps = [
        ":topological_torsion_lib",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
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
        "//Foundational/iwstring:iwstring",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule"
    ],
)

cc_test(
    name = "topological_torsion_test",
    srcs = [
        "topological_torsion_test.cc",
    ],
    deps = [
        ":topological_torsion_lib",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "//Molecule_Lib:moleculeio",
    ],
    timeout="short",
)
