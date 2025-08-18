add_executable(test_load_vmf3gr_grid test_load_vmf3_gr_grid.cpp)
add_executable(testAskne zwdAskne.cpp)
add_executable(msiseDataHunter msiseDataHunter.cpp)

target_link_libraries(test_load_vmf3gr_grid PRIVATE rwatmo ${PROJECT_DEPENDENCIES})
target_link_libraries(testAskne PRIVATE rwatmo ${PROJECT_DEPENDENCIES})
target_link_libraries(msiseDataHunter PRIVATE rwatmo ${PROJECT_DEPENDENCIES})

add_test(NAME test_load_vmf3gr_grid COMMAND test_load_vmf3gr_grid)
add_test(NAME testAskne COMMAND testAskne)
add_test(NAME msiseDataHunter COMMAND msiseDataHunter)
