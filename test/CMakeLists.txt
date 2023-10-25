# Add doctest as an interface library
set(DOCTEST_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/vendor/doctest)
add_library(doctest INTERFACE)

set(JSONCPP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/vendor/jsoncpp/include)
set(JSONCPP_SRC ${CMAKE_SOURCE_DIR}/vendor/jsoncpp/src/jsoncpp.cpp)

set(TEST_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/test/include)
target_include_directories(doctest INTERFACE ${TEST_INCLUDE_DIR} ${DOCTEST_INCLUDE_DIR} ${JSONCPP_INCLUDE_DIR})

add_definitions(-D_TEST_FILES_DIR=\"${CMAKE_SOURCE_DIR}/test/files\")
set(
        TEST_FILES
        ${JSONCPP_SRC}
        src/main_test.cpp
        src/iea_nucleolus_test.cpp
        src/nucleolus_test.cpp
        src/ordinal_shapley_test.cpp
)

add_executable (tests ${TEST_FILES})
target_link_libraries(tests doctest cgt)

set_target_properties(
  tests
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/public"
)

include(${DOCTEST_INCLUDE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(tests)