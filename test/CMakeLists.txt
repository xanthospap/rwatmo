add_executable(test1 vmf3SiteStreamBasic.cpp)
add_executable(testAskne zwdAskne.cpp)

target_link_libraries(test1 PRIVATE rwatmo)
target_link_libraries(testAskne PRIVATE rwatmo)

add_test(NAME test1 COMMAND test1)
add_test(NAME testAskne COMMAND testAskne)
