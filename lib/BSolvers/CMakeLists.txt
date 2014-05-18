if(OPENCL_FOUND)
    set(BSOLVERS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/BSolverOcl.cpp
    )
endif()

set(BSOLVERS_SRC
    ${BSOLVERS_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/BSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BSolverOmp.cpp
PARENT_SCOPE
)

