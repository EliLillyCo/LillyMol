load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":iwbdb_cat",
        ":iwbdb_compare",
        ":iwbdb_delete",
        ":iwbdb_exists",
        ":iwbdb_fetch",
        ":iwbdb_from_tdt",
        ":iwbdb_list",
        ":iwbdb_load",
        ":iwbdb_merge_into_tdt",
    ],
)

cc_binary(
    name = "iwbdb_cat",
    srcs = [
        "iwbdb_cat.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "iwbdb_compare",
    srcs = [
        "iwbdb_compare.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":unpack_data",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "iwbdb_delete",
    srcs = [
        "iwbdb_delete.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "iwbdb_exists",
    srcs = [
        "iwbdb_exists.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "iwbdb_fetch",
    srcs = [
        "iwbdb_fetch.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwbits",
        "//Foundational/iwmisc",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "iwbdb_from_tdt",
    srcs = [
        "iwbdb_from_tdt.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "@berkeleydb",
        "@re2",
        "@zlib",
    ],
)

cc_binary(
    name = "iwbdb_list",
    srcs = [
        "iwbdb_list.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        ":unpack_data",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "@berkeleydb",
    ],
)

cc_binary(
    name = "iwbdb_load",
    srcs = [
        "iwbdb_load.cc",
    ],
    tags = ["berkeleydb"],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "@berkeleydb",
    ],
)
cc_binary(
    name = "iwbdb_merge_into_tdt",
    srcs = [
        "iwbdb_merge_into_tdt.cc",
    ],
    tags = ["berkeleydb"],   
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "@berkeleydb",
        "@re2",
    ],
)

cc_library(
    name = "unpack_data",
    srcs = [
        "unpack_data.cc",
    ],
    hdrs = [
        "unpack_data.h",
    ],
    deps = [
        "//Foundational/iwstring",
    ],
)
