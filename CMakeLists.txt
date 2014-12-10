cmake_minimum_required(VERSION 2.8.4)
project(sdft)
enable_testing()

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")

set(SOURCE_FILES src/sdft.c)
set(TEST_FILES test/main.c)
include_directories(${CMAKE_SOURCE_DIR}/include)
add_library(sdft ${SOURCE_FILES})
add_executable(test_sdft ${TEST_FILES})
target_link_libraries(test_sdft sdft)
add_test(Tests test_sdft)