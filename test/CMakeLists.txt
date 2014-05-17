cmake_minimum_required(VERSION 2.6)
project(QualityTesting)

add_subdirectory(gtest-1.7.0)
enable_testing()

include_directories(
    ${gtest_SOURCE_DIR}/include
    ${gtest_SOURCE_DIR}
)

add_executable(runQualityTests
    tests.cpp
    irregulargrid.cpp
)

target_link_libraries(runQualityTests
    gtest
    gtest_main
    bpde
)

add_test(runQualityTests
    runQualityTests
)
