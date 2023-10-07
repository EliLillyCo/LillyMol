load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "iwstring_data_source",
    srcs = [
        "iwmmap.cc",
        "iwstring_data_source.cc",
        "string_data_source.cc",
    ],
    hdrs = [
        "iwmmap.h",
        "iwstring_data_source.h",
        "string_data_source.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwaray",
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwstring",
        "@re2",
        "@zlib",
    ],
)

cc_library(
    name = "iwtfdata_record",
    srcs = [
        "tfdatarecord.cc",
    ],
    hdrs = [
        "tfdatarecord.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
        "@com_google_absl//absl/crc:crc32c",
    ],
)

proto_library(
    name = "proto_for_testing",
    srcs = [
        "proto_for_testing.proto",
    ],
)

cc_proto_library(
    name = "proto_for_testing_cc_proto",
    deps = [
        ":proto_for_testing",
    ],
)

cc_test(
    name = "iwstring_data_source_test",
    srcs = [
        "iwstring_data_source_test.cc",
    ],
    linkopts = [
        "-lz",
    ],
    deps = [
        ":iwstring_data_source",
        "@googletest//:gtest_main",
    ],
)

cc_binary(
    name = "test_tfdata_record",
    srcs = [
        "test_tfdata_record.cc",
    ],
    deps = [
        ":iwtfdata_record",
        ":proto_for_testing_cc_proto",
        "//Foundational/cmdline:iwcmdline",
    ],
)

cc_binary(
    name = "tiwds",
    srcs = [
        "tiwds.cc",
    ],
    linkopts = [
        "-lz",
    ],
    deps = [
        ":iwstring_data_source",
        "//Foundational/cmdline:iwcmdline",
    ],
)
