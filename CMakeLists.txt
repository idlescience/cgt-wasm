cmake_minimum_required(VERSION 3.22)
project(nucleolus)

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

set(
        SOURCE_FILES
        main.cpp
        common.h
        common.cpp
        PD.h
        PD.cpp
)

add_executable(nucleolus ${SOURCE_FILES})

target_link_libraries(nucleolus glpk)