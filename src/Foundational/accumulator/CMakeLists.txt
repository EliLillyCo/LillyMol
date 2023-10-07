project(accumulator)
set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/KahanSum.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_float_v2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_double.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_double_v2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_int_v2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_long.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_long_v2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_uint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_uint_v2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_accumulator_w_missing_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_p_accumulator_float.cc
)

add_library(accumulator ${srcs})
