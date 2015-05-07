cmake_minimum_required(VERSION 2.8.4)
project(SWE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    f_calc.c
    init.c
    io.c
    jacobian.c
    shallow_water.c
    shallow_water.h
    )

add_executable(SWE ${SOURCE_FILES})