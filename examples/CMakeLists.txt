find_package(PCL CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

add_executable(exam_regis exam_regis.cpp)
target_link_libraries(exam_regis PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)
