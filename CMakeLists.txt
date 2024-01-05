cmake_minimum_required(VERSION 3.10)

project(GPSCO
        VERSION "0.0.1")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

add_subdirectory(algorithm)
add_subdirectory(examples)