cc_library(
    name = "accumulator",
    srcs = [
        "KahanSum.cc",
        "_accumulator_double.cc",
        "_accumulator_double_v2.cc",
        "_accumulator_float.cc",
        "_accumulator_float_v2.cc",
        "_accumulator_int.cc",
        "_accumulator_int_v2.cc",
        "_accumulator_long.cc",
        "_accumulator_long_v2.cc",
        "_accumulator_uint.cc",
        "_accumulator_uint_v2.cc",
        "_accumulator_w_missing_float.cc",
        "_resizable_array_p_accumulator_float.cc",
    ],
    hdrs = [
        "accumulator.h",
        "accumulator_v2.h",
        "kahan_sum.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Foundational/iwaray",
        "//Foundational/iwstring",
    ],
)
