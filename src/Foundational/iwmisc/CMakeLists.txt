project(iwmisc)
set(source_files
  ${CMAKE_CURRENT_SOURCE_DIR}/append_first_token_of_name.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dash_d.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dash_f.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dash_s.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dash_x.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/endian.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fetch_via_regexp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fraction_as_string.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/int_comparator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwdigits.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwdirname.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/logical_expression.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/msi_object.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/new_double.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/new_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/new_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/proto_support.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/report_progress.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/write_space_suppressed_string.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_numeric_data_from_file_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_numeric_data_from_file_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_set_or_unset_float.cc
)

set(header_files
  ${CMAKE_CURRENT_SOURCE_DIR}/fetch_via_regexp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iwarchive.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iwdigits.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iwminmax.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iw_time.h
  ${CMAKE_CURRENT_SOURCE_DIR}/logical_expression.h
  ${CMAKE_CURRENT_SOURCE_DIR}/minmaxspc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/misc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/msi_object.h
  ${CMAKE_CURRENT_SOURCE_DIR}/numeric_data_from_file.h
  ${CMAKE_CURRENT_SOURCE_DIR}/primes.h
  ${CMAKE_CURRENT_SOURCE_DIR}/proto_support.h
  ${CMAKE_CURRENT_SOURCE_DIR}/report_progress.h
  ${CMAKE_CURRENT_SOURCE_DIR}/set_or_unset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/timsort.hpp
)

set(iwconfig_header ${CMAKE_CURRENT_SOURCE_DIR}/iwconfig.h)

add_library(iwconfig INTERFACE)
target_sources(iwconfig INTERFACE "$<BUILD_INTERFACE:${iwconfig_header}>")
target_include_directories(iwconfig INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_library(iwmisc ${source_files} )

set(sparse_fp_creator_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/sparse_fp_creator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sparse_fp_creator_support.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/md5c.c
)

add_library(sparse_fp_creator ${sparse_fp_creator_srcs})

set(iwre2_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwre2.cc
)
add_library(iwre2 ${iwre2_srcs})
target_include_directories(iwre2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} re2)


set(tabular_data_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iw_tabular_data.cc
)
add_library(tabular_data ${tabular_data_srcs})

set(memoized_floats_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/memoized_floats.cc
)
add_library(memoized_floats ${memoized_floats_srcs})
target_include_directories(memoized_floats PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


if(IWTEST)
enable_testing()
include(GoogleTest)

set(iwdigits_test_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwdigits_test.cc
)

set(iwdigits_test_libs
  iwmisc
  iwstring
  gtest
  gmock
  gtest_main
)

add_executable(iwdigits_test ${iwdigits_test_srcs})
target_link_libraries(iwdigits_test LINK_PUBLIC ${iwdigits_test_libs})

add_test(
    NAME iwdigits_test
    COMMAND iwdigits_test
)
gtest_discover_tests(iwdigits_test)

set(memoized_floats_test_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/memoized_floats_test.cc
)

set(memoized_floats_test_libs
  memoized_floats
  iwmisc
  iwstring
  gtest
  gmock
  gtest_main
)
add_executable(memoized_floats_test ${memoized_floats_test_srcs})
target_link_libraries(memoized_floats_test LINK_PUBLIC ${memoized_floats_test_libs})

add_test(
    NAME memoized_floats_test
    COMMAND memoized_floats_test
)

gtest_discover_tests(memoized_floats_test)

set(iwre2_test_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwre2_test.cc
  ${re2}
)

set(iwre2_test_libs
  iwre2
  iwmisc
  iwstring
  gtest
  gmock
  gtest_main
  re2
)

add_executable(iwre2_test ${iwre2_test_srcs})
target_link_libraries(iwre2_test LINK_PUBLIC ${iwre2_test_libs})

add_test(
    NAME iwre2_test
    COMMAND iwre2_test
)

gtest_discover_tests(iwre2_test)
endif(IWTEST)
