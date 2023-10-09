load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":fetch_tdt",
        ":fetch_tdt_quick",
        ":tdt_join",
        ":tdt_sort",
        ":tdt_stats",
    ],
)

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
