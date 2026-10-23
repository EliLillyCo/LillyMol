load("@bazel_skylib//rules:common_settings.bzl", "bool_flag")

bool_flag(
    name = "build_inchi",
    build_setting_default = False,
    visibility = ["//visibility:public"],
)

config_setting(
    name = "build_inchi_enabled",
    flag_values = {
        ":build_inchi": "true",
    },
    visibility = ["//visibility:public"],
)
