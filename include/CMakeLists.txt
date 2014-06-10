add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/BSolvers)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/BLogger)


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/BSolvers
    ${CMAKE_CURRENT_SOURCE_DIR}/BLogger
)

set(BPDE_H
    ${CMAKE_CURRENT_SOURCE_DIR}/BArea.h
    ${CMAKE_CURRENT_SOURCE_DIR}/BBuilder.h
    ${BSOLVERS_SRC}
    ${BLOGGER_H}
PARENT_SCOPE
)


