set(
        LIB_FILES
        src/common.cpp
        src/nucleolus.cpp
)

add_library(cgt STATIC ${LIB_FILES})

target_link_libraries(cgt glpk)

set(NUCLEOLUS_INCLUDE_DIR include)

target_include_directories(cgt PUBLIC ${NUCLEOLUS_INCLUDE_DIR})

set_target_properties(
  cgt
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
)