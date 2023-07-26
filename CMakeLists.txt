cmake_minimum_required(VERSION 3.22)
project(nucleolus)

add_definitions(-std=c++17)
set(CMAKE_CXX_STANDARD 17)

if (CMAKE_CXX_COMPILER MATCHES "/em\\+\\+(-[a-zA-Z0-9.])?$")
    message(" * C++ compiler: Emscripten")
    set(CMAKE_CXX_COMPILER_ID "Emscripten")
else ()
    message(" * C++ compiler: ${CMAKE_CXX_COMPILER_ID}")
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Emscripten")
    set(
            INCLUDE_DIRECTORIES
            /home/secci/projects/glpk-wasm/include
    )
    set(
            LIB_DIRECTORIES
            /home/secci/projects/glpk-wasm/lib
    )
else ()
    set(
            INCLUDE_DIRECTORIES
            /home/secci/projects/glpk/include
    )
    set(
            LIB_DIRECTORIES
            /home/secci/projects/glpk/lib
    )
endif ()

include_directories(${INCLUDE_DIRECTORIES})

link_directories(${LIB_DIRECTORIES})

enable_testing()

add_subdirectory(lib)

add_subdirectory(test)
