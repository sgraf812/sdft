cmake_minimum_required(VERSION 2.8.4)
project(sdft)
enable_testing()

if (CMAKE_CXX_COMPILER_ID STREQUAL GNU OR
    CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    set(warnings "-Wall -Wextra -Werror")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL MSVC)
    set(warnings "/W4 /WX /EHsc")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings}")
set(SOURCE_FILES src/sdft.cpp)
set(TEST_FILES test/main.c)
include_directories(${CMAKE_SOURCE_DIR}/include)
add_library(sdft ${SOURCE_FILES})
add_executable(test_sdft ${TEST_FILES})
target_link_libraries(test_sdft sdft)
add_test(Tests test_sdft)