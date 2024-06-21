load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("//build_deps:install.bzl", "local_install")

package(default_visibility = ["//visibility:public"])

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

proto_library(
    name = "substituent_identification_proto",
    srcs = [
      "substituent_identification.proto",
    ],
)

cc_proto_library(
    name = "substituent_identification_cc_proto",
    deps = [
        "substituent_identification_proto",
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
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:dicer_fragments_cc_proto",
        "@berkeleydb",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer_fragment_replace",
    srcs = [
        "dicer_fragment_replace.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:molecular_formula_lib",
        "//Molecule_Tools:dicer_fragments_cc_proto",
        "@berkeleydb",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer2bdb",
    srcs = [
      "dicer2bdb.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwmisc",
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
        ":iwecfp_database_lookup_lib",
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
        "//Foundational/iwstring:absl_hash",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ],
)
cc_binary(
    name = "structure_database_load",
    srcs = [
        "structure_database_load.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":structure_database",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "structure_database_lookup",
    srcs = [
        "structure_database_lookup.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":structure_database",
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
        ":substituent_identification_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@berkeleydb",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_protobuf//:protobuf",
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
    name = "iwecfp_database_shared",
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
        "@berkeleydbshared//:berkeleydb",
    ],
)

cc_library(
    name = "iwecfp_database_lookup_lib",
    srcs = [
        "iwecfp_database_lookup_lib.cc",
    ],
    hdrs = [
        "iwecfp_database_lookup_lib.h",
    ],
    tags = ["berkeleydb"],
    deps = [
        "iwecfp_database",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "@berkeleydb",
    ],
)

cc_library(
    name = "iwecfp_database_lookup_lib_shared",
    srcs = [
        "iwecfp_database_lookup_lib.cc",
    ],
    hdrs = [
        "iwecfp_database_lookup_lib.h",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":iwecfp_database_shared",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "@berkeleydbshared//:berkeleydb",
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

cc_library(
    name = "selimsteg",
    srcs = [
        "selimsteg.cc",
    ],
    hdrs = [
        "selimsteg.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "@berkeleydbshared//:berkeleydb",
    ],
)

cc_library(
    name = "structure_database",
    srcs = [
        "structure_database.cc",
    ],
    hdrs = [
        "structure_database.h",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "@berkeleydb",
    ],
)

cc_library(
    name = "structure_database_shared",
    srcs = [
        "structure_database.cc",
    ],
    hdrs = [
        "structure_database.h",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "@berkeleydbshared//:berkeleydb",
    ],
)
