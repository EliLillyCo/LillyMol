load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":xgboost_model_evaluate",
    ]
)

proto_library(
    name = "xgboost_model_proto",
    srcs = [
        "xgboost_model.proto",
    ],
    tags = [
        "xgboost",
    ],
)

cc_proto_library(
    name = "xgboost_model_cc_proto",
    deps = [
        ":xgboost_model_proto",
    ],
    tags = [
        "xgboost",
    ],
    visibility = [
        "//visibility:public",
    ],
)

py_proto_library(
    name = "xgboost_model_py_proto",
    srcs = [
        "xgboost_model.proto",
    ],
    tags = [
        "xgboost",
    ],
)

proto_library(
    name = "random_forest_model_proto",
    srcs = [
        "random_forest_model.proto",
    ],
)

cc_proto_library(
    name = "random_forest_model_cc_proto",
    deps = [
        "random_forest_model_proto",
    ]
)

py_proto_library(
    name = "random_forest_model_py_proto",
    srcs = [
        "random_forest_model.proto",
    ],
)

cc_binary(
    name = "xgboost_model_evaluate",
    srcs = [
        "xgboost_model_evaluate.cc",
    ],
    deps = [
        ":xgboost_model_cc_proto",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "@xgboost//:xgboost",
        "@dlmc//:dlmc",
    ],
    tags = [
        "xgboost",
    ],
)
