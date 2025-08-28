load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":descriptor_file_to_distance_matrix",
    ],
)

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
