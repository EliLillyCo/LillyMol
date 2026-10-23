load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")

cc_binary(
    name = "fetch_tdt",
    srcs = [
        "fetch_tdt.cc",
    ],
    deps = [
        ":iw_tdt",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
    ],
)

cc_binary(
    name = "fetch_tdt_quick",
    srcs = [
        "fetch_tdt_quick.cc",
    ],
    deps = [
        ":iw_tdt",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
    ],
)

cc_binary(
    name = "tdt_join",
    srcs = [
        "tdt_join.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "tdt_stats",
    srcs = [
        "tdt_stats.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "tdt_sort",
    srcs = [
        "tdt_sort.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwqsort",
    ],
)

cc_library(
    name = "iw_tdt",
    # Linked whole into //Molecule_Lib:molecule_lib_shared, which claims to export
    # it. Without this the linker drops the members Molecule_Lib does not itself
    # reference, and extensions resolving against the shared library cannot find
    # them at load time.
    alwayslink = 1,
    srcs = [
        "_iw_tdt_double.cc",
        "_iw_tdt_float.cc",
        "_iw_tdt_int.cc",
        "_iw_tdt_iwstring.cc",
        "iw_tdt.cc",
    ],
    hdrs = [
        "iw_tdt.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwstring",
    ],
)
