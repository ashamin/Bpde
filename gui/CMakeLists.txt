cmake_minimum_required(VERSION 2.8.4)
  project(BpdeGui)

# here u could set your path to cmake directory in qt. to locate FindQt5Windets.cmake module
#set(CMAKE_PREFIX_PATH "/home/ashamin/Qt5.2.1/5.2.1/gcc_64/lib/cmake/")
#set(LD_LIBRARY_PATH /home/ashamin/R/x86_64-pc-linux-gnu-library/3.0/)

  # Tell CMake to run moc when necessary:
  set(CMAKE_AUTOMOC ON)
  # As moc files are generated in the binary dir, tell CMake
  # to always look for includes there:
  set(CMAKE_INCLUDE_CURRENT_DIR ON)

  # Widgets finds its own dependencies.
#  find_package(Qt5Widgets REQUIRED)
find_package(Qt5 COMPONENTS Widgets Svg REQUIRED)



#execute_process(COMMAND R RHOME
#                OUTPUT_VARIABLE R_HOME)

#message(${R_HOME})


set(NUM_TRUNC_CHARS 2)

execute_process(COMMAND R CMD config --cppflags
                OUTPUT_VARIABLE RCPP_FLAGS)
message(${RCPP_FLAGS})
string(SUBSTRING ${RCPP_FLAGS} ${NUM_TRUNC_CHARS} -1 RCPP_FLAGS)
message(${RCPP_FLAGS})

execute_process(COMMAND R CMD config --ldflags
                OUTPUT_VARIABLE RLD_FLAGS)
message(${RLD_FLAGS})

string(LENGTH ${RLD_FLAGS} RLD_FLAGS_LEN)

message(${RLD_FLAGS_LEN})

if (${RLD_FLAGS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLD_FLAGS_L)
    string(STRIP ${RLD_FLAGS_L} RLD_FLAGS_L )
    #link_directories(${RLD_FLAGS_L})
endif()

if (${RLD_FLAGS} MATCHES "[-][l]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLD_FLAGS_l)
    string(STRIP ${RLD_FLAGS_l} RLD_FLAGS_l )
endif()

message(${RLD_FLAGS_L})
message(${RLD_FLAGS_l})

execute_process(COMMAND R CMD config BLAS_LIBS
                OUTPUT_VARIABLE R_BLAS)
string(SUBSTRING ${R_BLAS} ${NUM_TRUNC_CHARS} -1 R_BLAS)

execute_process(COMMAND R CMD config LAPACK_LIBS
                OUTPUT_VARIABLE R_LAPACK)
string(SUBSTRING ${R_LAPACK} ${NUM_TRUNC_CHARS} -1 R_LAPACK)

message(${R_BLAS})
message(${R_LAPACK})

execute_process(COMMAND Rscript -e "Rcpp:::CxxFlags()"
                OUTPUT_VARIABLE RCPP_INCL)
string(SUBSTRING ${RCPP_INCL} ${NUM_TRUNC_CHARS} -1 RCPP_INCL)

message(${RCPP_INCL})


#### NOT REQUIRED. IS EMPTY
#execute_process(COMMAND Rscript -e "Rcpp:::LdFlags()"
#                OUTPUT_VARIABLE RCPP_LIBS)

#message(${RCPP_LIBS})

#if (${RCPP_LIBS} MATCHES "[-][L]([^ ;])+")
#    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPP_LIBS_L)
    #link_directories(${RCPP_LIBS_L} )
#endif()

#if (${RCPP_LIBS} MATCHES "[-][l][R]([^;])+")
#    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPP_LIBS_l)
#endif()

#message(${RCPP_LIBS_L})
#message(${RCPP_LIBS_l})

execute_process(COMMAND Rscript -e "RInside:::CxxFlags()"
                OUTPUT_VARIABLE RINSIDE_INCL)
string(SUBSTRING ${RINSIDE_INCL} ${NUM_TRUNC_CHARS} -1 RINSIDE_INCL)

message(${RINSIDE_INCL})


execute_process(COMMAND Rscript -e "RInside:::LdFlags()"
                OUTPUT_VARIABLE RINSIDE_LIBS)

message(${RINSIDE_LIBS})

if (${RINSIDE_LIBS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDE_LIBS_L)
    #link_directories(${RINSIDELIBS_L})
endif()

## remove part after whitespace if necessary
if (${RINSIDE_LIBS} MATCHES "[-][l][R]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDE_LIBS_l)
endif()
#set(RINSIDE_LIBS_l RInside)

message(${RINSIDE_LIBS_L})
message(${RINSIDE_LIBS_l})


include_directories(${RCPP_FLAGS} ${RCPP_INCL} ${RINSIDE_INCL})

link_directories(${RLD_FLAGS_L} ${RINSIDE_LIBS_L})






  add_executable(BpdeGui main.cpp bpdemainwindow.cpp bpdemainwindow.h)

## remove spaces from R_BLAS and R_LAPACK variables!!
target_link_libraries(BpdeGui
	${RLD_FLAGS_l}
	#${R_BLAS}
	#${R_LAPACK}
	blas
	lapack
	${RINSIDE_LIBS_l}
        bpde
)

  qt5_use_modules(BpdeGui Widgets Svg)
