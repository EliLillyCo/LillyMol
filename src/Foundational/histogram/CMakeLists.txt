project(iwhistogram)

set(histogram_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwhistogram.cc
)
add_library(iwhistogram ${histogram_srcs})
