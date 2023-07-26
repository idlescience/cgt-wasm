cmake_minimum_required(VERSION 3.22)
project(nucleolus)

add_definitions(-std=c++17)
set(CMAKE_CXX_STANDARD 17)

set(
        INCLUDE_DIRECTORIES
        /home/secci/projects/glpk-5.0/src
)

include_directories(${INCLUDE_DIRECTORIES})

set(
        LIB_DIRECTORIES
        /home/secci/projects/glpk-5.0/src/.libs
)

link_directories(${LIB_DIRECTORIES})

enable_testing()

add_subdirectory(lib)

add_subdirectory(test)
