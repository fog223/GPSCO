find_package(Eigen3 CONFIG REQUIRED)
find_package(PCL CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(Ceres CONFIG REQUIRED)

FILE(GLOB INCLUDE "./include/*.h")
FILE(GLOB SRC "./src/*.cpp")

#STATIC library
set(LIB_NAME "gpsco")
add_library(${LIB_NAME} STATIC)
target_sources(${LIB_NAME}
        PUBLIC
        ${INCLUDE}
        PRIVATE
        ${SRC})

target_include_directories(${LIB_NAME} PUBLIC ./include)
# link PCL and spdlog
target_link_libraries(${LIB_NAME} PRIVATE ${PCL_LIBRARIES} spdlog::spdlog Ceres::ceres)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${LIB_NAME} PRIVATE /bigobj)
endif ()