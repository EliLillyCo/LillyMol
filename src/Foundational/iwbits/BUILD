load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
cc_library(
    name = "iwbits",
    # Linked whole into //Molecule_Lib:molecule_lib_shared, which claims to export
    # it. Without this the linker drops the members Molecule_Lib does not itself
    # reference - fixed_bit_vector.o among them - and extensions that resolve
    # against the shared library cannot find them.
    alwayslink = 1,
    srcs = [
        "bits_in_common.c",
        "bits_in_common_avx.cc",
        "du_bin2ascii.cc",
        "fastprint.cc",
        "fixed_bit_vector.cc",
        "fixed_size_counted_fingerprint.cc",
        "hexbits.cc",
        "iwbits.cc",
        "swap_byte_order.cc",
    ],
    hdrs = [
        "dy_fingerprint.h",
        "fixed_bit_vector.h",
        "iwbits.h",
        "iwbits_support.h",
        "popcount.h"
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwstring",
    ],
)

cc_binary(
    name = "fixed_bit_vector_benchmark",
    srcs = [
        "fixed_bit_vector_benchmark.cc",
    ],
    deps = [
        ":iwbits",
        "@benchmark",
    ],
)

# cc_binary(
#     name = "fppop",
#     srcs = [
#         "fppop.cc",
#     ],
#     deps = [
#         ":iwbits",
#         "//Foundational/cmdline:iwcmdline",
#         "//Foundational/data_source:iwstring_data_source",
#         "//Foundational/iwmisc:iwconfig",
#         "//Foundational/iwmisc:iwmisc",
#         "//Foundational/iwstring:iwstring",
#         "//Molecule_Lib:iwmolecule",
#     ]
# )

cc_test(
    name = "fixed_bit_vector_test",
    srcs = [
        "fixed_bit_vector_test.cc",
    ],
    deps = [
        ":iwbits",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)
