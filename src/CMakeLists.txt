cmake_minimum_required(VERSION 3.16.1)
cmake_policy(SET CMP0077 NEW)
project(LillyMol)
# TODO:ianwatson get cmake working with this release
message(FATAL_ERROR "cmake infrastructure currently not functional")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# Wrong, causes an extra level
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build)
set(IW_EXTRA_LIBS protobuf snappy -lz re2)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

option(IWTEST "configure for GoogleTest" ON)

if(IWTEST)
include(CTest)
enable_testing()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
#find_package(GTest REQUIRED CONFIG)
endif(IWTEST)

FetchContent_Declare(
  re2
  URL https://github.com/google/re2/archive/refs/heads/main.zip
)
FetchContent_MakeAvailable(re2)
SET(RE2_BUILD_TESTING OFF)
include_directories(${re2_SOURCE_DIR})

find_package(Git QUIET)

FetchContent_Declare(
  snappy
  GIT_REPOSITORY https://github.com/google/snappy
  GIT_SUBMODULES third_party/benchmark third_party/googletest
  GIT_TAG main
)
set(SNAPPY_BUILD_TESTS OFF)
set(SNAPPY_BUILD_BENCHMARKS OFF)
FetchContent_MakeAvailable(snappy)
include_directories(${snappy_SOURCE_DIR})

FetchContent_Declare(
  protobuf
  URL https://github.com/protocolbuffers/protobuf/archive/refs/heads/master.zip
  SOURCE_SUBDIR cmake
  GIT_TAG main
)
set(protobuf_BUILD_TESTS OFF)
#FetchContent_MakeAvailable(protobuf)
include_directories(${protobuf_SOURCE_DIR})
find_package(Protobuf REQUIRED)
include_directories(${protobuf_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

FetchContent_Declare(
  leveldb
  GIT_REPOSITORY https://github.com/google/leveldb
# URL https://github.com/cschuet/leveldb/archive/refs/heads/master.zip
# URL https://github.com/google/leveldb/archive/refs/heads/main.zip
  GIT_TAG master
  GIT_SUBMODULES ""
)
set(LEVELDB_BUILD_TESTS OFF)
set(LEVELDB_BUILD_BENCHMARKS OFF)
FetchContent_MakeAvailable(leveldb)
include_directories(${leveldb_SOURCE_DIR}/include)

# For some reason adding protobuf as a library to targets
# does not work. So, we have a hack.
set(protobuf_build_dir ${CMAKE_CURRENT_BINARY_DIR}/_deps/protobuf-build)
set(protobuf_hack -L${protobuf_build_dir} -lprotobuf)
set(IW_EXTRA_LIBS ${IW_EXTRA_LIBS} ${protobuf_hack})
include_directories(BEFORE SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/_deps/protobuf-src/src)

#FetchContent_Declare(
#  clif
#  GIT_REPOSITORY https://github.com/google/clif
#  GIT_TAG 25a438fad1534f70c758c3d3c30ec361aeca63e2
#)
#FetchContent_MakeAvailable(clif)

# crc32c uses BUILD_TESTING to enable tests
set(CRC32C_BUILD_TESTS OFF)
set(CRC32C_BUILD_BENCHMARKS OFF)
set(CRC32C_USE_GLOG OFF)
FetchContent_Declare(
  crc32c
  URL https://github.com/google/crc32c/archive/refs/heads/main.zip
  GIT_TAG main
)
FetchContent_MakeAvailable(crc32c)
include_directories(${crc32c_SOURCE_DIR}/include)
#execute_process(COMMAND git submodule update --init -- ${CMAKE_CURRENT_BINARY_DIR}/_deps/crc32c-src/third_party
#                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_compile_options(-march=native)
add_link_options(-lz)
SET(CMAKE_EXE_LINKER_FLAGS "-lz")

add_subdirectory(Foundational/iwaray)
add_subdirectory(Foundational/accumulator)
add_subdirectory(Foundational/iwstring)
add_subdirectory(Foundational/cmdline)
add_subdirectory(Foundational/cmdline_v2)
add_subdirectory(Foundational/data_source)
add_subdirectory(Foundational/histogram)
add_subdirectory(Foundational/iwbits)
add_subdirectory(Foundational/iwmisc)
add_subdirectory(Foundational/iwqsort)
add_subdirectory(Foundational/iw_tdt)
add_subdirectory(Foundational/jama)
add_subdirectory(Foundational/string_data_source)
add_subdirectory(Foundational/tnt)
add_subdirectory(Foundational/xmlParser)
add_subdirectory(Utilities/General)
add_subdirectory(Utilities/GFP_Tools)
add_subdirectory(Molecule_Lib)
add_subdirectory(Molecule_Tools)
