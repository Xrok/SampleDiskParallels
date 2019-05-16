cmake_minimum_required(VERSION 3.14)
project(SampleDiskParallels)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

include_directories(.)

add_executable(SampleDiskParallels
        main.cpp
        ParallelSamplingDensePoints.h
        sample.h)
