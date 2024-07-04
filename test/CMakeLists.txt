add_executable(test1 vmf3SiteStreamBasic.cpp)

target_link_libraries(test1 PRIVATE rwatmo)

add_test(NAME test1 COMMAND test1)
