load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":distance_matrix_activity_difference",
        ":distance_matrix_from_distances",
        ":distance_matrix_kmedioids",
        ":distance_matrix_leader",
        ":distance_matrix_nn",
        ":distance_matrix_simple_cluster",
        ":distance_matrix_spread",
        ":distance_matrix_to_distances",
    ],
)

cc_binary(
    name = "distance_matrix_activity_difference",
    srcs = [
        "distance_matrix_activity_difference.cc",
    ],
    deps = [
        ":distance_matrix_base",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "distance_matrix_from_distances",
    srcs = [
        "distance_matrix_from_distances.cc",
    ],
    deps = [
        ":distance_matrix_base",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "distance_matrix_kmedioids",
    srcs = [
        "distance_matrix_kmedioids.cc",
    ],
    deps = [
        ":distance_matrix_base",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "distance_matrix_leader",
    srcs = [
        "distance_matrix_leader.cc",
    ],
    deps = [
        ":distance_matrix_base",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "distance_matrix_nn",
    srcs = [
        "distance_matrix_nn.cc",
    ],
    deps = [
        ":distance_matrix_base",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "distance_matrix_spread",
    srcs = [
        "distance_matrix_spread.cc",
    ],
    deps = [
        ":distance_matrix_base",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
    ],
)
cc_binary(
    name = "distance_matrix_simple_cluster",
    srcs = [
        "distance_matrix_simple_cluster.cc",
        "iwhcl.cc",
        "iwhcl.h",
    ],
    deps = [
        ":distance_matrix_base",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "distance_matrix_to_distances",
    srcs = [
        "distance_matrix_to_distances.cc",
    ],
    deps = [
        ":distance_matrix_base",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_library(
    name = "distance_matrix_base",
    srcs = [
        "_dm_to_nn_conditions_float.cc",
        "_iwdistance_matrix_byte.cc",
        "_iwdistance_matrix_double.cc",
        "_iwdistance_matrix_float.cc",
        "_iwdistance_matrix_float_byte.cc",
        "_iwdistance_matrix_int.cc",
        "_masquerading_as_byte_float.cc",
        "iwdmsupport.cc",
    ],
    hdrs = [
        "DM_to_NN_Conditions_.h",
        "IWDistanceMatrixBase.h",
        "Masquerading_as_Byte.h",
        "iwdmsupport.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
    ],
)
