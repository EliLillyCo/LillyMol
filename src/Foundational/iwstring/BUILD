load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")

cc_library(
    name = "absl_hash",
    hdrs = [
        "absl_hash.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":iwstring",
        "@com_google_absl//absl/hash",
    ],
)

cc_library(
    name = "iwstring",
    srcs = [
        "IWString_class.cc",
        "append_number.cc",
        "hexprint.cc",
        "iw_stl_hash_map.cc",
        "iwgetline.cc",
        "iwstring.cc",
        "iwstring_and_file_descriptor.cc",
        "iwstring_strstream.cc",
        "iwstrncasecmp.cc",
        "iwwrite.cc",
        "iwwrite_block.cc",
        "iwzlib.cc",
        "remove_suffix.cc",
        "string_change.cc",
        "string_relationals.cc",
        "string_to_float.cc",
        "tokenise_with_quotes.cc",
        "unhtml.cc",
    ],
    hdrs = [
        "iw_stl_hash_map.h",
        "iw_stl_hash_multimap.h",
        "iw_stl_hash_set.h",
        "iwhash.h",
        "iwstring.h",
        "iwzlib.h",
    ],
    linkopts = [
        "-lm",
    ],
    visibility = ["//visibility:public"],
    deps = [
        # not  implemented yet "@fast_float//:fast_float",
        "//Foundational/iwaray",
        "//Foundational/iwmisc:iwconfig",
        # "@third_party//:dragonbox",
        "@zlib",
    ],
)

cc_library(
    name = "string_with_commas",
    hdrs = [
        "string_with_commas.h",
    ],
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "tsclass",
    srcs = [
        "should_match.cc",
        "should_match.h",
        "tsclass.cc",
    ],
    linkopts = [
        "-lm",
    ],
    deps = [
        ":iwstring",
        "//Foundational/cmdline:iwcmdline",
    ],
)

cc_binary(
    name = "tspassbyref",
    srcs = [
        "tspassbyref.cc",
    ],
    deps = [
        ":iwstring",
    ],
)

cc_test(
    name = "iwstring_test",
    timeout = "short",
    srcs = [
        "iwstring_test.cc",
    ],
    deps = [
        ":iwstring",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "string_with_commas_test",
    srcs = [
        "string_with_commas_test.cc",
    ],
    deps = [
        ":string_with_commas",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ]
)

cc_test(
    name = "tokenise_with_quotes_test",
    srcs = [
        "tokenise_with_quotes_test.cc",
    ],
    deps = [
        ":iwstring",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)
