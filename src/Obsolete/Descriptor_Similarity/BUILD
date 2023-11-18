load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":descriptor_file_to_distance_matrix",
    ]
)

cc_binary(
    name = "descriptor_similarity",
    srcs = [
        "descriptor_similarity.cc",
        "determine_distance_type.cc",
        "iwdescriptor.cc",
        "iwdescriptor.h",
        "set_of_descriptors.h",
        "_iwdescriptors_float.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "descriptor_similarity_erg",
    srcs = [
        "descriptor_similarity_erg.cc",
        "determine_distance_type.cc",
        "iwdescriptor.cc",
        "iwdescriptor.h",
        "set_of_descriptors.h",
        "_iwdescriptors_float.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "descriptor_file_to_distance_matrix",
    srcs = [
        "descriptor_file_to_distance_matrix.cc",
        "determine_distance_type.cc",
        "iwdescriptor.cc",
        "iwdescriptor.h",
        "set_of_descriptors.h",
        "_iwdescriptors_float.cc",
    ],
    deps = [
        "//Utilities/Distance_Matrix:distance_matrix_base",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)
