if(OPENCL_FOUND)
    set(BSOLVERS_H
        ${CMAKE_CURRENT_SOURCE_DIR}/BSolverOcl.h
    )
endif()

set(BSOLVERS_H
    ${BSOLVERS_H}
    ${CMAKE_CURRENT_SOURCE_DIR}/BSolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/BSolverOmp.h
PARENT_SCOPE
)


