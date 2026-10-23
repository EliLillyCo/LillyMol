load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
cc_library(
    name = "iwhistogram",
    srcs = [
        "histogram_fp.cc",
        "iwhistogram.cc",
    ],
    hdrs = [
        "iwhistogram.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Foundational/iwbits",
        "//Foundational/iwstring",
    ],
)
