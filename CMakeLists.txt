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

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "/workspaces.*$")
    message(" * Github workspace detected")
    set(CMAKE_WORKSPACE "github")
elseif (CMAKE_CURRENT_SOURCE_DIR MATCHES "/home/runner/work.*$")
    message(" * Github runner detected")
    set(CMAKE_WORKSPACE "github")
else ()
    set(CMAKE_WORKSPACE "local")
    message(" * Local workspace detected")
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
    if (CMAKE_WORKSPACE STREQUAL "github")
        set(
                INCLUDE_DIRECTORIES
                ${CMAKE_CURRENT_SOURCE_DIR}/vendor/glpk/include
        )
        set(
                LIB_DIRECTORIES
                ${CMAKE_CURRENT_SOURCE_DIR}/vendor/glpk/lib
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
endif ()

include_directories(${INCLUDE_DIRECTORIES})

link_directories(${LIB_DIRECTORIES})

include(CTest)

enable_testing()

add_subdirectory(lib)

add_subdirectory(test)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Emscripten")
    add_subdirectory(module)
endif ()
