load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")
cc_library(
    name = "iwcmdline",
    srcs = [
        "cmdline.cc",
        "cmdline_stdstring.cc",
    ],
    hdrs = [
        "cmdline.h",
        "argv_permutation.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwstring",
    ],
)

cc_test(
    name = "argv_permutation_test",
    srcs = ["argv_permutation_test.cc"],
    deps = [
        ":iwcmdline",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "cmdline_test",
    size = "small",
    srcs = ["cmdline_test.cc"],
    deps = [
        ":iwcmdline",
        "@googletest//:gtest_main",
    ],
)
