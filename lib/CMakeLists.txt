set(
        LIB_FILES
        src/common.cpp
        src/primal_dual.cpp
        src/bnf.cpp
)

add_library(nucleolus SHARED ${LIB_FILES})

target_link_libraries(nucleolus glpk)

set(NUCLEOLUS_INCLUDE_DIR include)

target_include_directories(nucleolus PUBLIC ${NUCLEOLUS_INCLUDE_DIR})