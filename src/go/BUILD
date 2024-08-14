load("@rules_go//go:def.bzl", "go_binary", "go_test")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":grep_sdf",
        ":parallel_process_file",
        ":no_spaces_in_file_name",
        ":regression_to_classification",
        ":rxn_reverse",
        ":rxnsmiles2smi",
    ],
)

go_binary(
    name = "parallel_process_file",
    srcs = [
        "parallel_process_file.go",
    ],
)

go_binary(
    name = "no_spaces_in_file_name",
    srcs = [
        "no_spaces_in_file_name.go",
    ],
)

go_binary(
    name = "regression_to_classification",
    srcs = [
        "regression_to_classification.go",
    ],
)

go_binary(
    name = "rxn_reverse",
    srcs = [
        "rxn_reverse.go",
    ],
)

go_binary(
    name = "rxnsmiles2smi",
    srcs = [
        "rxnsmiles2smi.go",
    ],
)

go_binary(
    name = "grep_sdf",
    srcs = [
        "grep_sdf.go",
    ],
)
