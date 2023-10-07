project(iw_tdt)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iw_tdt.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_iw_tdt_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_iw_tdt_double.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_iw_tdt_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_iw_tdt_iwstring.cc
)

add_library(iw_tdt ${srcs})


