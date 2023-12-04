cmake_minimum_required(VERSION 3.20)

find_package(PCL CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

FILE(GLOB INCLUDE "./include/*.h")
FILE(GLOB SRC "./src/*.cpp")

#STATIC library
set(LIB_NAME "gpsco")
add_library(${LIB_NAME} STATIC)
target_include_directories(${LIB_NAME} PUBLIC ./include)
target_sources(${LIB_NAME}
        PUBLIC
        ./include/common.h
        ./include/plane_extraction.h
        PRIVATE
        ./src/common.cpp
        ./src/plane_extraction.cpp)

# link PCL
target_link_libraries(${LIB_NAME} PRIVATE ${PCL_LIBRARIES})

# link spdlog
target_link_libraries(${LIB_NAME} PRIVATE spdlog::spdlog)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${LIB_NAME} PRIVATE /bigobj)
endif ()