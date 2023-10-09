cc_library(
    name = "iwaray",
    srcs = [
        "_extending_resizable_array_double.cc",
        "_extending_resizable_array_float.cc",
        "_extending_resizable_array_int.cc",
        "_extending_resizable_array_long.cc",
        "_extending_resizable_array_long_long.cc",
        "_extending_resizable_array_uchar.cc",
        "_extending_resizable_array_ushort.cc",
        "_iwaray_float.cc",
        "_iwaray_int.cc",
        "_resizable_array_char.cc",
        "_resizable_array_char2.cc",
        "_resizable_array_double.cc",
        "_resizable_array_float.cc",
        "_resizable_array_int.cc",
        "_resizable_array_int_star.cc",
        "_resizable_array_long.cc",
        "_resizable_array_long_long.cc",
        "_resizable_array_p_int.cc",
        "_resizable_array_uchar.cc",
        "_resizable_array_uint.cc",
        "_resizable_array_unsigned_long.cc",
        "_resizable_array_ushort.cc",
    ],
    hdrs = [
        "iwaray.h",
        "iwaray_np_.h",
        "iwaray_op.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwmisc:iwconfig",
    ],
)

cc_test(
    name = "iwaray_test",
    srcs = [
        "iwaray_test.cc",
    ],
    deps = [
        ":iwaray",
        "@com_google_protobuf//:protobuf",
        "@googletest//:gtest_main",
    ],
)
