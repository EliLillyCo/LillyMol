project(iwqsort)

set(iwqsort_hdr
  ${CMAKE_CURRENT_SOURCE_DIR}/iwqsort.h
)
add_library(iwqsort INTERFACE)
target_sources(iwqsort INTERFACE "$<BUILD_INTERFACE:${iwqsort_header}>")
target_include_directories(iwqsort INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

