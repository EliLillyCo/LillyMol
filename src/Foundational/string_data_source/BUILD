cc_library(
    name = "string_data_source",
    srcs = [
        "iwstring_string_data_source.cc",
    ],
    hdrs = [
        "iwstring_string_data_source.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwaray",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwstring",
        "@com_google_protobuf//:protobuf",
        "@re2",
    ],
)

cc_test(
    name = "iwstring_string_data_source_test",
    srcs = [
        "iwstring_string_data_source_test.cc",
    ],
    deps = [
        ":string_data_source",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwstring",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
        "@re2",
    ],
)
