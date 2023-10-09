project(iwstring)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/append_number.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/hexprint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwgetline.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/IWString_class.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_and_file_descriptor.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwwrite.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwwrite_block.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstrncasecmp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_strstream.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwzlib.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iw_stl_hash_map.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/remove_suffix.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/string_change.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/string_relationals.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/unhtml.cc
)

add_library(iwstring ${srcs})
target_include_directories(iwstring PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(tsclass_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/tsclass.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/should_match.cc
)
set(tsclass_libs
  iwstring
  cmdline
)

add_executable(tsclass ${tsclass_srcs})
target_link_libraries(tsclass PUBLIC ${tsclass_libs})

if(IWTEST)
enable_testing()

set(iwstring_test_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_test.cc
)
set(iwstring_test_libs
  iwstring
  cmdline
  gtest
  gmock
  gtest_main
)

# Google Test not working yet... I can test with bazel.
add_executable(iwstring_test ${iwstring_test_srcs})
target_link_libraries(iwstring_test LINK_PUBLIC ${iwstring_test_libs})

add_test(
    NAME iwstring_test
    COMMAND iwstring_test
)

include(GoogleTest)
gtest_discover_tests(iwstring_test)
endif(IWTEST)
