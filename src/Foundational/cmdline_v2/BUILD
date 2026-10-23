load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")

cc_library(
    name = "cmdline_v2",
    srcs = [
        "cmdline.cc",
    ],
    hdrs = [
        "cmdline_v2.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwaray",
        "//Foundational/iwmisc",
        "//Foundational/iwstring",
    ],
)
