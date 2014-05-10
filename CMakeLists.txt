cmake_minimum_required(VERSION 2.8)
project(Bpde)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/area)
add_subdirectory(${PROJECT_SOURCE_DIR}/bsolvers)

include_directories(${PROJECT_SOURCE_DIR}/area)
include_directories(${PROJECT_SOURCE_DIR}/bsolvers)


set(BPDE_H
    ${AREA_H}
    ${BSOLVERS_H}
)

set(BPDE_SRC
    ${AREA_SRC}
    ${BSOLVERS_SRC}
)

add_executable(${PROJECT_NAME}
    main.cpp
    ${BPDE_H}
    ${BPDE_SRC}
)

