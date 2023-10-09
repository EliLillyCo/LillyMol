project(data_source)

set(iwstring_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_data_source.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwmmap.cc
)
add_library(iwstring_data_source ${iwstring_srcs})

set(tfdatarecord_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/tfdatarecord.cc
)
add_library(tfdatarecord ${tfdatarecord_srcs})
target_include_directories(iwstring_data_source PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

IF(IWTEST)
include(GoogleTest)
enable_testing()
find_package(Protobuf REQUIRED)

protobuf_generate_cpp(for_testing_proto_srcs for_testing_proto_hdrs proto_for_testing.proto)

message("hdrs " ${for_testing_proto_srcs})
set(tfdatarecord_test_srcs
  tfdatarecord_test.cc
  ${for_testing_proto_hdrs} 
  ${for_testing_proto_srcs}
)
set(tfdatarecord_test_libs
  tfdatarecord
  iwstring
  iwstring_data_source
  iwbits
  gtest
  gmock
  gtest_main
  snappy
  protobuf -lz
)
# Will not build inside Lilly due to problems in protocol buffers
# add_executable(tfdatarecord_test ${tfdatarecord_test_srcs})
#target_include_directories(tfdatarecord_test
#   PRIVATE
#   ${CMAKE_CURRENT_BINARY_DIR}
#)
#target_link_libraries(tfdatarecord_test LINK_PUBLIC ${tfdatarecord_test_libs})

#add_test(
#  NAME tfdatarecord_test
#  COMMAND tfdatarecord_test
#)
#gtest_discover_tests(tfdatarecord_test)

set(iwstring_data_source_test_srcs
  iwstring_data_source_test.cc
)
set(iwstring_data_source_test_libs
  tfdatarecord
  iwstring_data_source
  iwstring
  iwre2
  gtest
  gmock
  gtest_main
  re2 -lz
)
add_executable(iwstring_data_source_test ${iwstring_data_source_test_srcs})
target_link_libraries(iwstring_data_source_test LINK_PUBLIC ${iwstring_data_source_test_libs})

add_test(
  NAME iwstring_data_source_test
  COMMAND iwstring_data_source_test
)

gtest_discover_tests(iwstring_data_source_test)

endif(IWTEST)
