# Add doctest as an interface library
set(DOCTEST_INCLUDE_DIR include)
add_library(doctest INTERFACE)
target_include_directories(doctest INTERFACE ${DOCTEST_INCLUDE_DIR})

set(
        TEST_FILES
        src/nucleolus_test.cpp
)

add_executable (tests ${TEST_FILES})
target_link_libraries(tests doctest nucleolus)

set_target_properties(
  tests
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
)