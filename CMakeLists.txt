cmake_minimum_required(VERSION 2.8)
project(Bpde)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

option(ENABLE_QUALITY_TESTING "Enabling quality testing" OFF)

add_subdirectory(${PROJECT_SOURCE_DIR}/area)
add_subdirectory(${PROJECT_SOURCE_DIR}/bsolvers)
add_subdirectory(${PROJECT_SOURCE_DIR}/bbuilder)

include_directories(
    ${PROJECT_SOURCE_DIR}/area
    ${PROJECT_SOURCE_DIR}/bsolvers
    ${PROJECT_SOURCE_DIR}/bbuilder
)

set(BPDE_H
    ${AREA_H}
    ${BSOLVERS_H}
    ${BBUILDER_H}
)

set(BPDE_SRC
    ${AREA_SRC}
    ${BSOLVERS_SRC}
    ${BBUILDER_SRC}
)

add_executable(${PROJECT_NAME}
    main.cpp
    ${BPDE_H}
    ${BPDE_SRC}
)

if(ENABLE_QUALITY_TESTING)
    message("Enabling tests")
    add_library(bpde
        SHARED
        ${BPDE_H}
        ${BPDE_SRC}
    )

    add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif()
