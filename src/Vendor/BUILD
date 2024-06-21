load("@rules_proto//proto:defs.bzl", "proto_library")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":clogp2descriptors_biobyte",
        ":rcorina",
    ],
    tags = [
        "vendor",
    ],
)

cc_binary(
    name = "clogp2descriptors_biobyte",
    srcs = [
        "clogp2descriptors_biobyte.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)


cc_binary(
    name = "iwclogp",
    srcs = [
        "iwclogp.cc",
        "biobyte_stuff.c",
    ],
    tags = [
        "vendor",
    ],
    deps = [
        "@biobyte//:biobyte",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(  
    name = "rcorina",
    srcs = [
        "rcorina.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/string_data_source:string_data_source",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@corina//:corina",
    ],
    tags = [
        "vendor",
    ],
)

