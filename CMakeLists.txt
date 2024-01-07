cmake_minimum_required(VERSION 3.16)

project(GPSCO
        VERSION "0.0.1")

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

add_subdirectory(algorithm)
add_subdirectory(examples)