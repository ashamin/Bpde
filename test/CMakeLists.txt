cmake_minimum_required(VERSION 2.6)
project(basic_test)

add_subdirectory(gtest-1.7.0)
enable_testing()

include_directories(
    ${gtest_SOURCE_DIR}/include
    ${gtest_SOURCE_DIR}
)

add_executable( runUnitTests tests.cpp irregulargrid.cpp )
target_link_libraries(runUnitTests gtest gtest_main bpde)
add_test( runUnitTests runUnitTests )
