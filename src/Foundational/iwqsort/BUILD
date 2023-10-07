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
