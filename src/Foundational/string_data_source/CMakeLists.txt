project(string_data_source)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_string_data_source.cc
)

add_library(string_data_source ${srcs})
target_include_directories(string_data_source PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(IWTEST)
enable_testing()

set(iwstring_string_data_source_test_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_string_data_source_test.cc
)

set(iwstring_string_data_source_test_libs
  iwstring
  string_data_source
  iwre2
  gtest
  gmock
  gtest_main
  re2
)

# Google Test not working yet... I can test with bazel.
add_executable(iwstring_string_data_source_test ${iwstring_string_data_source_test_srcs})
target_link_libraries(iwstring_string_data_source_test LINK_PUBLIC ${iwstring_string_data_source_test_libs})

add_test(
    NAME iwstring_string_data_source_test
    COMMAND iwstring_string_data_source_test
)

include(GoogleTest)
gtest_discover_tests(iwstring_string_data_source_test)
endif(IWTEST)
