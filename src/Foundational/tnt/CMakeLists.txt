project(tnt)

file(GLOB hdrs ${PROJECT_SOURCE_DIR}/*.h)

add_library(tnt INTERFACE)
target_sources(tnt INTERFACE "$<BUILD_INTERFACE:${hdrs}>")
target_include_directories(tnt INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
