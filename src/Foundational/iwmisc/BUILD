load("@rules_proto//proto:defs.bzl", "proto_library")

proto_library(
    name = "for_testing_proto",
    srcs = [
        "for_testing.proto",
    ],
)

cc_proto_library(
    name = "for_testing_cc_proto",
    deps = [
        "for_testing_proto",
    ]
)


cc_library(
    name = "activity_from_file",
    srcs = [
    ],
    hdrs = [
        "activity_data_from_file.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwstring",
    ],
)

cc_library(
    name = "compile_time",
    srcs = [
        "compile_time.cc",
    ],
    hdrs = [
        "compile_time.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
    ],
)

cc_library(
    name = "iwmisc",
    srcs = [
        "_iw_write_array_int.cc",
        "_iwmaxid_double_int.cc",
        "_iwmaxid_float_iwstring.cc",
        "_iwmaxid_int_int.cc",
        "_iwminid_double_int.cc",
        "_iwminmax_float.cc",
        "_iwqsort.cc",
        "_numeric_data_from_file_float.cc",
        "_numeric_data_from_file_int.cc",
        "_set_or_unset_float.cc",
        "append_first_token_of_name.cc",
        "dash_d.cc",
        "dash_f.cc",
        "dash_s.cc",
        "dash_x.cc",
        "endian.cc",
        "fetch_via_regexp.cc",
        "fname.cc",
        "fraction_as_string.cc",
        "hton_stuff.cc",
        "int_comparator.cc",
        "iwdigits.cc",
        "iwdirname.cc",
        "logical_expression.cc",
        "msi_object.cc",
        "new_double.cc",
        "new_float.cc",
        "new_int.cc",
        "new_unsigned_int.cc",
        "proto_support.cc",
        "report_progress.cc",
        "write_first_token_of_string.cc",
        "write_space_suppressed_string.cc",
    ],
    hdrs = [
        "combinations.h",
        "fetch_via_regexp.h",
        "iw_time.h",
        "iwarchive.h",
        "iwdigits.h",
        "iwminmax.h",
        "logical_expression.h",
        "matcher.h",
        "minmaxspc.h",
        "misc.h",
        "msi_object.h",
        "numeric_data_from_file.h",
        "primes.h",
        "proto_support.h",
        "report_progress.h",
        "set_or_unset.h",
        "timsort.hpp",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":internal",
        ":iwconfig",
        ":ltoml",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwqsort",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_library(
    name = "ltoml",
    srcs = [
        "toml_support.h",
    ],
    deps = [
        "@tomlplusplus//:toml++",
        "@com_google_absl//absl/status:status",
    ],
)

cc_library(
    name = "internal",
    hdrs = [
        "new_array_.h",
    ],
)

cc_library(
    name = "iwconfig",
    hdrs = [
        "iwconfig.h",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "sparse_fp_creator",
    srcs = [
        "md5c.c",
        "sparse_fp_creator.cc",
        "sparse_fp_creator_support.cc",
    ],
    hdrs = [
        "md5.h",
        "md5_global.h",
        "sparse_fp_creator.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":iwmisc",
        "//Foundational/iwaray",
        "//Foundational/iwbits",
        "//Foundational/iwqsort",
        "//Foundational/iwstring",
    ],
)

cc_library(
    name = "iw_tabular_data",
    srcs = [
        "iw_tabular_data.cc",
    ],
    hdrs = [
        "iw_tabular_data.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwbits",
    ],
)

cc_library(
    name = "iwre2",
    srcs = [
        "iwre2.cc",
    ],
    hdrs = [
        "iwre2.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwstring",
        "@com_google_absl//absl/strings:strings",
        "@re2",
    ],
)

cc_test(
    name = "combinations_test",
    srcs = [
        "combinations_test.cc",
    ],
    deps = [
        ":iwmisc",
        "//Foundational/iwstring",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "matcher_test",
    srcs = [
        "matcher_test.cc",
    ],
    deps = [
        ":iwmisc",
        "//Foundational/iwstring",
        "@googletest//:gtest_main",
    ],
)

cc_library(
    name = "memoized_floats",
    srcs = [
        "memoized_floats.cc",
    ],
    hdrs = [
        "memoized_floats.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":iwmisc",
    ],
)

cc_library(
    name = "normalisation",
    srcs = [
        "normalisation.cc",
    ],
    hdrs = [
        "normalisation.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
    ],
)

cc_binary(
    name = "test_iwdigits",
    srcs = [
        "test_iwdigits.cc",
    ],
    deps = [
        ":iwmisc",
        "//Foundational/cmdline_v2",
    ],
)

cc_test(
    name = "iwdigits_test",
    srcs = [
        "iwdigits_test.cc",
    ],
    deps = [
        ":iwmisc",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "iwre2_test",
    srcs = [
        "iwre2_test.cc",
    ],
    deps = [
        ":iwre2",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "logical_expression_test",
    srcs = [
        "logical_expression_test.cc",
    ],
    deps = [
        ":iwmisc",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "memoized_floats_test",
    srcs = [
        "memoized_floats_test.cc",
    ],
    deps = [
        ":memoized_floats",
        "@benchmark",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "memoized_floats_benchmark",
    srcs = [
        "memoized_floats_benchmark.cc",
    ],
    deps = [
        ":memoized_floats",
        "@benchmark//:benchmark_main",
    ],
)

cc_test(
    name = "toml_support_test",
    srcs = [
        "toml_support_test.cc",
    ],
    deps = [
        "for_testing_cc_proto",
        ":iwmisc",
        ":ltoml",
        "@googletest//:gtest_main",
        "@tomlplusplus//:toml++",
    ]
)
