project(jama)

file(GLOB hdrs ${PROJECT_SOURCE_DIR}/*.h)

add_library(jama INTERFACE)
target_sources(jama INTERFACE "$<BUILD_INTERFACE:${hdrs}>")
target_include_directories(jama INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
