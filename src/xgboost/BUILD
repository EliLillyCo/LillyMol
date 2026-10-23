load("@com_google_protobuf//bazel:py_proto_library.bzl", "py_proto_library")
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")

load("@rules_proto//proto:defs.bzl", "proto_library")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")
# load("@rules_proto_grpc_go//:defs.bzl", "go_proto_compile", "go_proto_library")

proto_library(
    name = "xgboost_model_proto",
    srcs = [
        "xgboost_model.proto",
    ],
    tags = [
        "xgboost",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_proto_library(
    name = "xgboost_model_cc_proto",
    tags = [
        "xgboost",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":xgboost_model_proto",
    ],
)

py_proto_library(
    name = "xgboost_model_py_proto",
    deps = [
        ":xgboost_model_proto",
    ],
    tags = [
        "xgboost",
    ],
)

#go_proto_compile(
#    name = "xgboost_model_go_proto",
#    protos = [
#        ":xgboost_model_proto",
#    ],
#)

#go_proto_library(
#    name = "xgboost_model_go_proto",
#    importpath = "xgboost_lib",
#    protos = [
#        ":xgboost_model_proto",
#    ],
#    visibility = ["//visibility:public"],
#)

proto_library(
    name = "random_forest_model_proto",
    srcs = [
        "random_forest_model.proto",
    ],
)

cc_proto_library(
    name = "random_forest_model_cc_proto",
    deps = [
        ":random_forest_model_proto",
    ],
)

py_proto_library(
    name = "random_forest_model_py_proto",
    deps = [
        ":random_forest_model_proto",
    ],
)

cc_binary(
    name = "xgboost_model_evaluate",
    srcs = [
        "xgboost_model_evaluate.cc",
    ],
    tags = [
        "xgboost",
    ],
    deps = [
        ":xgboost_model_cc_proto",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Utilities/General:class_label_translation_cc_proto",
        "//Utilities/General:linear_scaling",
        "@third_party//:dlmc",
        "@third_party//:xgboost",
    ],
)
