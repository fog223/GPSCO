find_package(PCL CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

add_executable(exam_plane_extraction exam_plane_extraction.cpp)
target_link_libraries(exam_plane_extraction PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_regis exam_regis.cpp)
target_link_libraries(exam_regis PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_accuracy exam_accuracy.cpp)
target_link_libraries(exam_accuracy PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_overlap exam_overlap.cpp)
target_link_libraries(exam_overlap PRIVATE ${PCL_LIBRARIES} spdlog::spdlog)

