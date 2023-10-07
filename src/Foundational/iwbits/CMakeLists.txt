project(iwbits)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/bits_in_common.c
  ${CMAKE_CURRENT_SOURCE_DIR}/du_bin2ascii.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fastprint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fixed_size_counted_fingerprint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/hexbits.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwbits.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/swap_byte_order.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fixed_bit_vector.cc
)

add_library(iwbits ${srcs})


IF(IWTEST)
enable_testing()

set(fixed_bit_vector_test_srcs
  fixed_bit_vector_test.cc
)
set(fixed_bit_vector_test_libs
  iwstring
  iwbits
  gtest_main
)
add_executable(fixed_bit_vector_test ${fixed_bit_vector_test_srcs})
target_link_libraries(fixed_bit_vector_test LINK_PUBLIC ${fixed_bit_vector_test_libs})

add_test(
  NAME fixed_bit_vector_test
  COMMAND fixed_bit_vector_test
)

include(GoogleTest)
gtest_discover_tests(fixed_bit_vector_test)
endif(IWTEST)
