add_executable(test1 vmf3SiteStreamBasic.cpp)
add_executable(testAskne zwdAskne.cpp)
add_executable(msiseDataHunter msiseDataHunter.cpp)

target_link_libraries(test1 PRIVATE rwatmo ${PROJECT_DEPENDENCIES})
target_link_libraries(testAskne PRIVATE rwatmo ${PROJECT_DEPENDENCIES})
target_link_libraries(msiseDataHunter PRIVATE rwatmo ${PROJECT_DEPENDENCIES})

add_test(NAME test1 COMMAND test1)
add_test(NAME testAskne COMMAND testAskne)
add_test(NAME msiseDataHunter COMMAND msiseDataHunter)
