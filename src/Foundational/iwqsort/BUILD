load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")

cc_library(
    name = "iwqsort",
    hdrs = [
        "iwqsort.h",
        "iwqsort_ptr.h",
        "iwqsort_ptr_.h",
        "iwqsort_tbb.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwaray",
    ],
)

cc_library(
    name = "iwqsort_parallel",
    hdrs = [
        "iwqsort_tbb.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwaray",
        "@onetbb//:tbb",
    ],
)

cc_binary(
    name = "test_iwqsort_fo",
    srcs = [
       "test_iwqsort_fo.cc",
    ],
    deps = [
        ":iwqsort",
        "//Foundational/cmdline:iwcmdline",
    ],
)
