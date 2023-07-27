# Add doctest as an interface library
set(DOCTEST_INCLUDE_DIR include)
add_library(doctest INTERFACE)
target_include_directories(doctest INTERFACE ${DOCTEST_INCLUDE_DIR})

set(
        TEST_FILES
        src/nucleolus_test.cpp
        src/gen_game.cpp
)

add_executable (tests ${TEST_FILES})
target_link_libraries(tests doctest nucleolus)


add_test(NAME tests COMMAND tests)
