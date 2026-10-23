load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

cc_binary(
    name = "descriptor_similarity",
    srcs = [
        "_iwdescriptors_float.cc",
        "descriptor_similarity.cc",
        "determine_distance_type.cc",
        "iwdescriptor.cc",
        "iwdescriptor.h",
        "set_of_descriptors.h",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "descriptor_similarity_erg",
    srcs = [
        "_iwdescriptors_float.cc",
        "descriptor_similarity_erg.cc",
        "determine_distance_type.cc",
        "iwdescriptor.cc",
        "iwdescriptor.h",
        "set_of_descriptors.h",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_to_distance_matrix",
    srcs = [
        "_iwdescriptors_float.cc",
        "descriptor_file_to_distance_matrix.cc",
        "determine_distance_type.cc",
        "iwdescriptor.cc",
        "iwdescriptor.h",
        "set_of_descriptors.h",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Obsolete/Distance_Matrix:distance_matrix_base",
    ],
)
