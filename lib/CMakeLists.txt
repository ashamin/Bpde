add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/BSolvers)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/BLogger)

set(BPDE_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/BArea.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BBuilder.cpp
    ${BSOLVERS_SRC}
    ${BLOGGER_SRC}
PARENT_SCOPE
)

