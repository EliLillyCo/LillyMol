project(cmdline)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cmdline_stdstring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_iwaray_ov.cc
)

add_library(cmdline ${srcs})
target_include_directories(cmdline PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
