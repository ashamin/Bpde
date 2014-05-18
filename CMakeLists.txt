cmake_minimum_required(VERSION 2.8)
project(Bpde)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake
)

find_package(OpenCL REQUIRED)
if(OPENCL_FOUND)
    message("OPENCL FOUND")
    include_directories(${OPENCL_INCLUDE_DIRS})
endif()

option(ENABLE_QUALITY_TESTING "Enabling quality testing" OFF)

add_subdirectory(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/lib)

include_directories(
    ${PROJECT_SOURCE_DIR}/include/BSolvers
    ${PROJECT_SOURCE_DIR}/include
)

add_library(bpde
    SHARED
    ${BPDE_H}
    ${BPDE_SRC}
)

if(OPENCL_FOUND)
    target_link_libraries(bpde
        ${OPENCL_LIBRARIES}
    )
endif()

if(ENABLE_QUALITY_TESTING)
    message("Enabling tests")

    add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif()
