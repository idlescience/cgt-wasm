cmake_minimum_required(VERSION 3.22)
project(nucleolus)

add_definitions(-std=c++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-s EXIT_RUNTIME=1")

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

set(
        SOURCE_FILES
        main.cpp
        lib/common.h
        lib/common.cpp
        lib/nucleolus.h
        lib/nucleolus.cpp
)

add_executable(nucleolus ${SOURCE_FILES})

target_link_libraries(nucleolus glpk)