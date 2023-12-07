find_package(PCL CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

add_executable(exam_plane_extraction exam_plane_extraction.cpp)
target_link_libraries(exam_plane_extraction PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_regis exam_regis.cpp)
target_link_libraries(exam_regis PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

#add_executable(exam_registration exam_registration.cpp)
#target_link_libraries(exam_registration PRIVATE gpsco)
#
#add_executable(exam_transmatrix exam_transmatrix.cpp)
#target_link_libraries(exam_transmatrix PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)
#
#add_executable(exam_align exam_align.cpp)
#target_link_libraries(exam_align PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)