load("@rules_proto//proto:defs.bzl", "proto_library")
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")

package(default_visibility = ["//visibility:public"])

proto_library(
    name = "conformer_database_proto",
    srcs = [
        "conformer_database.proto",
    ],
)

cc_proto_library(
    name = "conformer_database_cc_proto",
    deps = [
        ":conformer_database_proto",
    ],
)

proto_library(
    name = "dicer_fragment_complement_proto",
    srcs = [
        "dicer_fragment_complement.proto",
    ],
)

cc_proto_library(
    name = "dicer_fragment_complement_cc_proto",
    deps = [
        ":dicer_fragment_complement_proto",
    ],
)

proto_library(
    name = "qupkake_transfer_proto",
    srcs = [
        "qupkake_transfer.proto",
    ],
)

cc_proto_library(
    name = "qupkake_transfer_cc_proto",
    deps = [
        ":qupkake_transfer_proto",
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
        "//Foundational/iwstring:string_with_commas",
        ":storage_conditions",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@third_party//:berkeley_static",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "conformer_database_lookup",
    srcs = [
        "conformer_database_lookup.cc",
    ],
    deps = [
        "conformer_database_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@third_party//:berkeley_static",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "conformer_database_store",
    srcs = [
        "conformer_database_store.cc",
    ],
    deps = [
        "//Foundational/iwstring:string_with_commas",
        "conformer_database_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_static",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer_fragment_complement_replace",
    srcs = [
        "dicer_fragment_complement_replace.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "dicer_fragment_complement_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwstring:absl_hash",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:dicer_fragments_cc_proto",
        "//Molecule_Tools:mformula",
        "@third_party//:berkeley_static",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "dicer_to_complemenent_db",
    srcs = [
        "dicer_to_complemenent_db.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/iwstring:string_with_commas",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwtfdata_record",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Tools:dicer_fragments_cc_proto",
        "@com_google_absl//absl/container:flat_hash_map",
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_static",
    ],
)

cc_binary(
    name = "iwecfp_database_load",
    srcs = [
        "iwecfp_database_load.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/iwstring:string_with_commas",
        ":iwecfp_database",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@third_party//:berkeley_static",
    ],
)

cc_binary(
    name = "iwecfp_database_lookup",
    srcs = [
        "iwecfp_database_lookup.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/iwstring:string_with_commas",
        ":iwecfp_database",
        ":iwecfp_database_lookup_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@third_party//:berkeley_static",
    ],
)

cc_binary(
    name = "qupkake_transfer",
    srcs = [
        "qupkake_transfer.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":iwecfp_database_lookup_lib",
        ":qupkake_transfer_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Foundational/iwstring:absl_hash",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@third_party//:berkeley_static",
        "@com_google_absl//absl/container:flat_hash_set",
    ],
)

cc_binary(
    name = "smi2rings_bdb",
    srcs = [
        "smi2rings_bdb.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/iwstring:string_with_commas",
        ":smi2rings_cc_proto",
        ":ring_database_lookup",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Foundational/iwstring:absl_hash",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_static",
    ],
)

cc_library(
    name = "substituent_identification_lib",
    srcs = [
        "substituent_identification_lib.cc",
    ],
    hdrs = [
        "substituent_identification_lib.h",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/iwstring:string_with_commas",
        ":substituent_identification_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:mformula",
        "//Molecule_Tools:smarts_for_atom_subset",
        "@third_party//:berkeley_static",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_library(
    name = "substituent_identification_lib_shared",
    srcs = [
        "substituent_identification_lib.cc",
    ],
    hdrs = [
        "substituent_identification_lib.h",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/iwstring:string_with_commas",
        ":substituent_identification_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:mformula",
        "//Molecule_Tools:smarts_for_atom_subset",
        "@third_party//:berkeley_shared",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "substituent_identification",
    srcs = [
        "substituent_identification.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":substituent_identification_lib",
    ],
)

cc_test(
    name = "substituent_identification_lib_test",
    srcs = [
        "substituent_identification_lib_test.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":substituent_identification_lib",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_binary(
    name = "substituent_identification_extract",
    srcs = [
      "substituent_identification_extract.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":substituent_identification_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "@com_google_absl//absl/strings:string_view",
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_shared",
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
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_shared",
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
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_shared",
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
        "@third_party//:berkeley_static",
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
        "@third_party//:berkeley_shared",
    ],
)

cc_library(
    name = "ring_database_lookup",
    srcs = ["ring_database_lookup.cc"],
    hdrs = ["ring_database_lookup.h"],
    tags = ["berkeleydb"],
    visibility = ["//visibility:public"],
    deps = [
        ":smi2rings_cc_proto",
        "//Foundational/iwaray",
        "//Foundational/iwstring",
        "@third_party//:berkeley_static",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_library(
    name = "ring_database_lookup_shared",
    srcs = ["ring_database_lookup.cc"],
    hdrs = ["ring_database_lookup.h"],
    tags = ["berkeleydb"],
    visibility = ["//visibility:public"],
    deps = [
        ":smi2rings_cc_proto",
        "//Foundational/iwaray",
        "//Foundational/iwstring",
        "@third_party//:berkeley_shared",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_test(
    name = "ring_database_lookup_test",
    srcs = ["ring_database_lookup_test.cc"],
    tags = ["berkeleydb"],
    deps = [
        ":ring_database_lookup",
        ":smi2rings_cc_proto",
        "@third_party//:berkeley_static",
        "@com_google_protobuf//:protobuf",
        "@googletest//:gtest_main",
    ],
)
