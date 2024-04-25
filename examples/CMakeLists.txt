find_package(PCL CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

add_executable(exam_plane_extraction exam_plane_extraction.cpp)
target_link_libraries(exam_plane_extraction PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_regis exam_regis.cpp)
target_link_libraries(exam_regis PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_overlap exam_overlap.cpp)
target_link_libraries(exam_overlap PRIVATE ${PCL_LIBRARIES} spdlog::spdlog)

add_executable(exam_accuracy_hs1 exam_accuracy_hs1.cpp)
target_link_libraries(exam_accuracy_hs1 PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_accuracy_hs2 exam_accuracy_hs2.cpp)
target_link_libraries(exam_accuracy_hs2 PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_accuracy_park exam_accuracy_park.cpp)
target_link_libraries(exam_accuracy_park PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_accuracy_campus exam_accuracy_campus.cpp)
target_link_libraries(exam_accuracy_campus PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_accuracy_apartment exam_accuracy_apartment.cpp)
target_link_libraries(exam_accuracy_apartment PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_accuracy_stairs exam_accuracy_stairs.cpp)
target_link_libraries(exam_accuracy_stairs PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_accuracy_hauptgebaude exam_accuracy_hauptgebaude.cpp)
target_link_libraries(exam_accuracy_hauptgebaude PRIVATE ${PCL_LIBRARIES} spdlog::spdlog gpsco)

add_executable(exam_display exam_display.cpp)
target_link_libraries(exam_display PRIVATE ${PCL_LIBRARIES} spdlog::spdlog)
