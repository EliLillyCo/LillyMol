load("@bazel_tools//tools/cpp:cc_toolchain_config_lib.bzl", "tool_path")
load(":cc_toolchain_config.bzl", "cc_toolchain_config")

package(default_visibility = ["//visibility:public"])

cc_toolchain_suite(
    name = "clang_suite",
    toolchains = {
        "skylake": ":skylake_toolchain",
    },
)

filegroup(name = "empty")

cc_toolchain(
    name = "skylake_toolchain",
    all_files = ":empty",
    compiler_files = ":empty",
    dwp_files = ":empty",
    linker_files = ":empty",
    objcopy_files = ":empty",
    strip_files = ":empty",
    supports_param_files = 0,
    toolchain_config = ":skylake_toolchain_config",
    toolchain_identifier = "skylake-toolchain",
)

cc_toolchain_config(name = "skylake_toolchain_config")
