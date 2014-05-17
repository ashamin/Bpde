add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/BSolvers)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/BSolvers
)

set(BPDE_H
    ${CMAKE_CURRENT_SOURCE_DIR}/BArea.h
    ${CMAKE_CURRENT_SOURCE_DIR}/BBuilder.h
    ${BSOLVERS_SRC}
PARENT_SCOPE
)


