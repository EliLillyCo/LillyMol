load("@gazelle//:def.bzl", "gazelle")
load("@rules_go//go:def.bzl", "go_binary", "go_library")
load("//build_deps:install.bzl", "local_install")

gazelle(name = "gazelle")

local_install(
    name = "install",
    srcs = [
        ":grep_sdf",
        ":no_spaces_in_file_name",
        ":parallel_process_file",
        ":regression_to_classification",
        ":rxn_reverse",
        ":rxnsmiles2smi",
    ],
)

# gazelle:prefix xgboost

go_binary(
    name = "parallel_process_file",
    srcs = ["parallel_process_file.go"],
)

go_binary(
    name = "no_spaces_in_file_name",
    srcs = ["no_spaces_in_file_name.go"],
)

go_binary(
    name = "regression_to_classification",
    srcs = ["regression_to_classification.go"],
)

go_binary(
    name = "rxn_reverse",
    srcs = ["rxn_reverse.go"],
)

go_binary(
    name = "rxnsmiles2smi",
    srcs = ["rxnsmiles2smi.go"],
)

go_binary(
    name = "grep_sdf",
    srcs = ["grep_sdf.go"],
)

# Could never get this to work.
# Problem came down to trying to make the import
# google.golang.org/protobuf/proto
# avaialble, and I was never able to figure that out.

# go_binary(
#     name = "xgbd_evaluate",
#     srcs = ["xgbd_evaluate.go"],
#     deps = [
#         "//xgboost:xgboost_model_go_proto",
#       we need something here for the missing dependency.
#     ],
# )
