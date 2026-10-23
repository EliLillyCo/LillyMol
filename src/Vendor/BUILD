load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")

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
        "biobyte_stuff.c",
        "iwclogp.cc",
    ],
    tags = [
        "vendor",
        "biobyte",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:alogp_lib",
        "@biobyte",
    ],
)

cc_binary(
    name = "marvin2gfp",
    srcs = [
        "marvin2gfp.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_binary(
    name = "marvin_pka",
    srcs = [
        "marvin_pka.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "rcorina",
    srcs = [
        "rcorina.cc",
    ],
    tags = [
        "vendor",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/string_data_source",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@corina",
    ],
)
#
#cc_library(
#    name = "corina",
#    srcs = [
#       "/lrlhps/apps/gc3tk/c3tk-core/lib/corina/lib/libcorina_linux_rhel6_x86_64_libc2.12_gcc4.4.7.so",
#    ],
#    hdrs = ["/lrlhps/apps/gc3tk/c3tk-core/lib/corina/include/corina_lib.h"],
#    strip_include_prefix = "/lrlhps/apps/gc3tk/c3tk-core/lib/corina",
#    visibility = ["//visibility:public"],
#)
