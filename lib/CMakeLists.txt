set(
        LIB_FILES
        src/common.cpp
        src/bnf.cpp
)

add_library(nucleolus STATIC ${LIB_FILES})

target_link_libraries(nucleolus glpk)

set(NUCLEOLUS_INCLUDE_DIR include)

target_include_directories(nucleolus PUBLIC ${NUCLEOLUS_INCLUDE_DIR})

set_target_properties(
  nucleolus
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
)