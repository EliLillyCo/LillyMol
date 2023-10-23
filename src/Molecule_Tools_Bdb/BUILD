load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":buildsmidb_bdb",
        ":in_database_bdb",
        ":in_lilly_database_bdb",
        ":iwecfp_database_load",
        ":iwecfp_database_lookup",
        ":smi2rings_bdb",
        ":substituent_identification",
    ],
)


proto_library(
    name = "smi2rings_proto",
    srcs = [
        "smi2rings.proto",
    ],
)

cc_proto_library(
    name = "smi2rings_cc_proto",
    deps = [
        "smi2rings_proto",
    ],
)

cc_binary(
    name = "buildsmidb_bdb",
    srcs = [
        "buildsmidb_bdb.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":storage_conditions",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer_fragment_lookup_bdb",
    srcs = [
        "dicer_fragment_lookup_bdb.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Tools:dicer_fragments_cc_proto",
        "@berkeleydb",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "in_database_bdb",
    srcs = [
        "in_database_bdb.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":storage_conditions",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ],
)
cc_binary(  
    name = "in_lilly_database_bdb",   
    srcs = [
        "in_lilly_database_bdb.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":storage_conditions",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "iwecfp_database_load",
    srcs = [
        "iwecfp_database_load.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":iwecfp_database",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "iwecfp_database_lookup",
    srcs = [
        "iwecfp_database_lookup.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":iwecfp_database",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "smi2rings_bdb",
    srcs = [
        "smi2rings_bdb.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":smi2rings_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "substituent_identification",
    srcs = [
        "substituent_identification.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ],
)

cc_library(
    name = "iwecfp_database",
    srcs = [
        "iwecfp_database.cc",
    ],
    hdrs = [
        "iwecfp_database.h",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "@berkeleydb",
    ],
)

cc_library(
    name = "storage_conditions",
    srcs = [
        "storage_conditions.cc",
    ],
    hdrs = [
        "storage_conditions.h",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "@berkeleydb",
    ],
)
